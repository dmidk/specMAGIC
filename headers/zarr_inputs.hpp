#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <netcdf.h>

#include "navigation.hpp"
#include "tables.hpp"
#include "types.hpp"

// Fields extracted from the zarr archives by py_utils/extract_zarr_inputs.py, already on
// the output grid. They replace the climatologies (atmosphere) and the land-use albedo pixel
// by pixel, wherever they have a value.
//
// A requested input that is not available at all (missing file or variable, file on another
// grid) is a problem: by default it is reported and the climatology is used, with
// strict = true the run stops. Pixels that a field does not cover only ever fall back.

namespace ZarrInputs {

    // Stands in for NaN once the fields are read: -ffast-math lets the compiler assume
    // that NaN never occurs, so the per-pixel tests must not rely on it
    constexpr float MISSING = -9999.f;

    // Share of pixels falling back to the climatology above which the run warns about it
    constexpr double FALLBACK_WARNING_SHARE = 0.01;

    // LSA SAF gives its visible albedo (AL-VI-DH) below this wavelength and its
    // near-infrared albedo (AL-NI-DH) above it
    constexpr double LSASAF_VIS_NIR_SPLIT_NM = 700.;

    struct Request {
        bool atmosphere = false;   // aerosol and ozone from CAMS, water vapour from DINI and CAMS
        bool albedo = false;       // albedo from LSA SAF
        bool strict = false;       // stop instead of falling back when a requested input is missing
    };

    // One field on the output grid, stored row-major (lat, lon) in a plain array, as in
    // Matrix. It is allocated once while the files are read and never resized, so that the
    // parallel loop only ever reads fixed memory.
    struct Field {

        std::string file;
        std::string provenance;

        Field() noexcept = default;

        Field(int rows, int cols) : rows_(rows), cols_(cols) {
            if (rows_ <= 0 || cols_ <= 0) throw std::runtime_error("Field: dimensions must be positive");
            data_ = new float[static_cast<size_t>(rows_) * cols_];
        }

        ~Field() { delete[] data_; }

        // Owns its data, so it is moved, never copied
        Field(const Field&) = delete;
        Field& operator=(const Field&) = delete;

        Field(Field&& other) noexcept
            : file(std::move(other.file)), provenance(std::move(other.provenance)),
              rows_(other.rows_), cols_(other.cols_), data_(other.data_) {
            other.rows_ = other.cols_ = 0;
            other.data_ = nullptr;
        }

        Field& operator=(Field&& other) noexcept {
            if (this != &other) {
                delete[] data_;
                file = std::move(other.file);
                provenance = std::move(other.provenance);
                rows_ = other.rows_;
                cols_ = other.cols_;
                data_ = other.data_;

                other.rows_ = other.cols_ = 0;
                other.data_ = nullptr;
            }
            return *this;
        }

        bool present() const noexcept { return data_ != nullptr; }

        // Unchecked access. Assumes 0 <= lat < rows and 0 <= lon < cols, which holds for
        // every pixel, since readFile only accepts files on the output grid.
        float operator()(int lat, int lon) const noexcept {
            return data_[static_cast<size_t>(lat) * cols_ + lon];
        }

        float* data() noexcept { return data_; }
        size_t size() const noexcept { return static_cast<size_t>(rows_) * cols_; }

    private:
        int rows_{0}, cols_{0};
        float* data_{nullptr};
    };

    // How the inputs of the day pixels were sourced, summed over threads after the loop
    struct Tally {
        uint64_t primary = 0;
        uint64_t backup = 0;
        uint64_t climatology = 0;
        uint64_t clipped = 0;
        uint64_t water = 0;        // albedo only: no LSA SAF albedo over water, which is normal

        Tally& operator+=(const Tally& o) {
            primary += o.primary; backup += o.backup; climatology += o.climatology;
            clipped += o.clipped; water += o.water;
            return *this;
        }
    };

    struct Usage {
        Tally aod, o3, h2o, albedo;

        Usage& operator+=(const Usage& o) {
            aod += o.aod; o3 += o.o3; h2o += o.h2o; albedo += o.albedo;
            return *this;
        }
    };

    struct Inputs {

        Field aod, o3, h2o, h2o_cams, albedo_vis, albedo_nir, albedo_land;
        Request request;

        // Requested inputs that are not available, reported again after the run
        std::vector<std::string> problems;

        // Concentration axes of the LUTs. Values from the archives are clipped to these,
        // because unlike the smooth monthly climatologies they do leave them
        MAGIC_REAL aod_min = 0, aod_max = 0;
        MAGIC_REAL h2o_min = 0, h2o_max = 0;
        MAGIC_REAL o3_min = 0, o3_max = 0;

        bool requested() const { return request.atmosphere || request.albedo; }

        // Reads every file named in the list (one path per line), then checks that each
        // requested field was found. Throws if any was not and the request is strict.
        void load(const std::string& list_path, const Geography& geo, Request r) {

            request = r;

            if (list_path.empty()) {
                problem("ZARR_INPUTS is not set, so no zarr inputs were given");
            } else if (std::ifstream in{list_path}) {
                std::string path;
                while (std::getline(in, path)) {
                    if (!path.empty()) readFile(path, geo);
                }
            } else {
                problem("could not open the zarr input list " + list_path);
            }

            if (request.atmosphere) {
                expect(aod, "aod (CAMS aod550)", "the aerosol climatology");
                expect(o3, "o3 (CAMS gtco3)", "the ozone climatology");
                expect(h2o, "h2o (DINI wvint)", "CAMS tcwv, then the water vapour climatology");
                expect(h2o_cams, "h2o_cams (CAMS tcwv)", "the water vapour climatology outside DINI");
            }
            if (request.albedo) {
                expect(albedo_vis, "albedo_vis (LSA SAF AL-VI-DH)", "the land-use albedo");
                expect(albedo_nir, "albedo_nir (LSA SAF AL-NI-DH)", "the land-use albedo");
                expect(albedo_land, "albedo_land (LSA SAF QFLAGS)", "counting missing albedo over water too");
            }

            if (request.strict && !problems.empty()) {
                std::string message = "Requested zarr inputs are missing, and strict inputs were asked for:";
                for (const auto& p : problems) message += "\n    " + p;
                throw std::runtime_error(message);
            }
        }

        void setLimits(const Tables::RTM& rtm, const Tables::Correction::Absorber& water,
            const Tables::Correction::Absorber& ozone) {

            aod_min = rtm.aod[0];       aod_max = rtm.aod[rtm.aoddim - 1];
            h2o_min = water.concen[0];  h2o_max = water.concen[water.nrows - 1];
            o3_min = ozone.concen[0];   o3_max = ozone.concen[ozone.nrows - 1];
        }

        // Overrides the climatology values of this pixel wherever an archive has one,
        // and counts where the albedo of this pixel will come from
        void apply(PixelClimate& clim, const Area& a, Usage& u) const {

            if (request.atmosphere) {
                useArchiveValue(clim.aod, aod, nullptr, aod_min, aod_max, u.aod, a);
                useArchiveValue(clim.o3, o3, nullptr, o3_min, o3_max, u.o3, a);

                // DINI first, CAMS outside the DINI domain
                useArchiveValue(clim.h2o, h2o, &h2o_cams, h2o_min, h2o_max, u.h2o, a);
            }

            if (request.albedo) {
                if (hasAlbedo(a)) u.albedo.primary++;
                else if (at(albedo_land, a) == 0.f) u.albedo.water++;
                else u.albedo.climatology++;
            }
        }

        // LSA SAF visible covers 0.4-0.7 um and near-infrared 0.7-4.0 um. Both must be
        // there, so that a pixel never mixes LSA SAF and land-use albedo across bands.
        // Returns -1 where there is no value, which sends the caller to the land-use albedo.
        MAGIC_REAL getKatoAlbedo(int zero_based_kato_band, const Area& a) const {

            if (!hasAlbedo(a)) return (MAGIC_REAL)-1;

            // Kato band 16 (684.2-704.4 nm, zero-based 15) lies mostly below
            // LSASAF_VIS_NIR_SPLIT_NM, so bands 1-16 take the visible albedo
            return at(zero_based_kato_band <= 15 ? albedo_vis : albedo_nir, a);
        }

        MAGIC_REAL getSatelliteAlbedo(double wavelength_nm, const Area& a) const {

            if (!hasAlbedo(a)) return (MAGIC_REAL)-1;
            return at(wavelength_nm < LSASAF_VIS_NIR_SPLIT_NM ? albedo_vis : albedo_nir, a);
        }

        // For the global attributes of the output files: where every requested input came from
        std::vector<std::pair<std::string, std::string>> describe(const Usage& u) const {

            std::vector<std::pair<std::string, std::string>> attrs;

            if (request.atmosphere) {
                attrs.emplace_back("aod_source", source(aod, "aerosol climatology"));
                attrs.emplace_back("o3_source", source(o3, "ozone climatology"));
                attrs.emplace_back("h2o_source", source(h2o, "water vapour climatology"));
                attrs.emplace_back("h2o_backup_source", source(h2o_cams, "water vapour climatology"));
                attrs.emplace_back("aod_usage", summary(u.aod, false));
                attrs.emplace_back("o3_usage", summary(u.o3, false));
                attrs.emplace_back("h2o_usage", summary(u.h2o, false));
            }
            if (request.albedo) {
                attrs.emplace_back("albedo_vis_source", source(albedo_vis, "land-use albedo"));
                attrs.emplace_back("albedo_nir_source", source(albedo_nir, "land-use albedo"));
                attrs.emplace_back("albedo_usage", summary(u.albedo, true));
            }

            std::string all;
            for (const auto& p : problems) all += (all.empty() ? "" : "; ") + p;
            if (!all.empty()) attrs.emplace_back("input_problems", all);

            return attrs;
        }

        // After the run: where the inputs came from, and again every problem, so that a
        // run that fell back is hard to miss
        void report(const Usage& u) const {

            if (request.atmosphere) {
                report("aod", u.aod, false);
                report("o3", u.o3, false);
                report("h2o", u.h2o, false);
            }
            if (request.albedo) report("albedo", u.albedo, true);

            if (problems.empty()) return;

            std::printf("\n[WARNING]: ================ MISSING ZARR INPUTS ================\n");
            for (const auto& p : problems) std::printf("[WARNING]:   %s\n", p.c_str());
            std::printf("[WARNING]: These requested inputs fell back to the climatologies.\n");
            std::printf("[WARNING]: Run with --strict-inputs to stop instead.\n");
            std::printf("[WARNING]: ======================================================\n\n");
        }

    private:

        void problem(const std::string& message) {
            problems.push_back(message);
            std::printf("[WARNING]: %s \n", message.c_str());
        }

        void expect(const Field& f, const char* what, const char* instead) {
            if (!f.present()) problem(std::string(what) + " was requested but not found; using " + instead);
        }

        float at(const Field& f, const Area& a) const {
            if (!f.present()) return MISSING;
            return f(a.nlat, a.nlon);
        }

        bool hasAlbedo(const Area& a) const {
            float vis = at(albedo_vis, a);
            float nir = at(albedo_nir, a);
            return vis >= 0.f && vis <= 1.f && nir >= 0.f && nir <= 1.f;
        }

        // Takes the primary field's value, else the backup's, clipped to the LUT axis.
        // Leaves the climatology value in place when neither has one.
        void useArchiveValue(MAGIC_REAL& target, const Field& primary, const Field* backup,
            MAGIC_REAL lo, MAGIC_REAL hi, Tally& t, const Area& a) const {

            float v = at(primary, a);
            if (v != MISSING) {
                t.primary++;
            } else if (backup && (v = at(*backup, a)) != MISSING) {
                t.backup++;
            } else {
                t.climatology++;
                return;
            }

            MAGIC_REAL x = static_cast<MAGIC_REAL>(v);
            if (x < lo || x > hi) t.clipped++;
            target = std::min(std::max(x, lo), hi);
        }

        static std::string source(const Field& f, const char* fallback) {
            return f.present() ? f.provenance : std::string(fallback) + " (zarr input missing)";
        }

        // Same result as std::isnan, but std::isnan cannot be used: under -ffast-math the
        // compiler assumes NaN never occurs and folds std::isnan(v) to false (GCC 9 does).
        // Testing the IEEE 754 bits directly (exponent all ones, fraction non-zero) survives
        static bool isNan(float v) {
            uint32_t bits;
            std::memcpy(&bits, &v, sizeof bits);
            return (bits & 0x7f800000u) == 0x7f800000u && (bits & 0x007fffffu) != 0;
        }

        // Shares of the day pixels; for albedo, of the day pixels on land
        static double fallbackShare(const Tally& t) {
            uint64_t total = t.primary + t.backup + t.climatology;
            return total == 0 ? 0. : static_cast<double>(t.climatology) / total;
        }

        static std::string summary(const Tally& t, bool albedo) {
            uint64_t total = t.primary + t.backup + t.climatology;
            if (total == 0) return "no day pixels";

            char buf[200];
            if (albedo) {
                std::snprintf(buf, sizeof buf,
                    "%.1f%% of land day pixels from LSA SAF, %.1f%% land-use albedo; "
                    "%.1f%% of all day pixels are water and use the land-use albedo",
                    100. * t.primary / total, 100. * t.climatology / total,
                    100. * t.water / (total + t.water));
            } else {
                std::snprintf(buf, sizeof buf,
                    "%.1f%% zarr, %.1f%% zarr backup, %.1f%% climatology, %.1f%% clipped to LUT",
                    100. * t.primary / total, 100. * t.backup / total,
                    100. * t.climatology / total, 100. * t.clipped / total);
            }
            return buf;
        }

        static void report(const char* name, const Tally& t, bool albedo) {
            if (t.primary + t.backup + t.climatology + t.water == 0) return;

            const char* level = fallbackShare(t) > FALLBACK_WARNING_SHARE ? "[WARNING]: " : "";
            std::printf("%sZarr input %s: %s \n", level, name, summary(t, albedo).c_str());
        }

        static bool sameAxis(const double* file_axis, const MAGIC_EXACT* geo_axis, size_t n) {
            for (size_t i = 0; i < n; i++) {
                if (std::fabs(file_axis[i] - geo_axis[i]) > 1e-3) return false;
            }
            return true;
        }

        static std::string textAttribute(int ncid, int varid, const char* name) {
            size_t len = 0;
            if (nc_inq_attlen(ncid, varid, name, &len) != NC_NOERR) return "";

            std::string value(len, '\0');
            if (nc_get_att_text(ncid, varid, name, value.data()) != NC_NOERR) return "";
            return value;
        }

        void readFile(const std::string& path, const Geography& geo) {

            int ncid;
            if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR) {
                problem("could not open zarr input " + path);
                return;
            }

            // The file must be on exactly the grid we compute on
            int dimid;
            size_t file_nlat = 0, file_nlon = 0;
            bool ok = nc_inq_dimid(ncid, "lat", &dimid) == NC_NOERR
                && nc_inq_dimlen(ncid, dimid, &file_nlat) == NC_NOERR
                && nc_inq_dimid(ncid, "lon", &dimid) == NC_NOERR
                && nc_inq_dimlen(ncid, dimid, &file_nlon) == NC_NOERR
                && file_nlat == geo.nlat && file_nlon == geo.nlon;

            if (ok) {
                double* file_lat = new double[file_nlat];
                double* file_lon = new double[file_nlon];
                int varid;
                ok = nc_inq_varid(ncid, "lat", &varid) == NC_NOERR
                    && nc_get_var_double(ncid, varid, file_lat) == NC_NOERR
                    && nc_inq_varid(ncid, "lon", &varid) == NC_NOERR
                    && nc_get_var_double(ncid, varid, file_lon) == NC_NOERR
                    && sameAxis(file_lat, geo.lat, file_nlat) && sameAxis(file_lon, geo.lon, file_nlon);
                delete[] file_lat;
                delete[] file_lon;
            }

            if (!ok) {
                problem("zarr input " + path + " is not on the output grid of magic-config.asc");
                nc_close(ncid);
                return;
            }

            const std::pair<const char*, Field*> fields[] = {
                {"aod", &aod}, {"o3", &o3}, {"h2o", &h2o}, {"h2o_cams", &h2o_cams},
                {"albedo_vis", &albedo_vis}, {"albedo_nir", &albedo_nir}, {"albedo_land", &albedo_land}};

            for (const auto& f : fields) {

                int varid;
                if (nc_inq_varid(ncid, f.first, &varid) != NC_NOERR) continue;

                if (f.second->present()) {
                    std::printf("[WARNING]: Zarr input %s found again in %s. Keeping the one from %s. \n",
                        f.first, path.c_str(), f.second->file.c_str());
                    continue;
                }

                Field field(static_cast<int>(file_nlat), static_cast<int>(file_nlon));
                if (nc_get_var_float(ncid, varid, field.data()) != NC_NOERR) {
                    problem(std::string("could not read ") + f.first + " from " + path);
                    continue;
                }

                float* values = field.data();
                for (size_t i = 0; i < field.size(); i++) {
                    if (isNan(values[i])) values[i] = MISSING;
                }

                field.file = path;
                field.provenance = textAttribute(ncid, varid, "provenance");
                if (field.provenance.empty()) field.provenance = path;
                *f.second = std::move(field);
            }

            nc_close(ncid);
        }
    };

}
