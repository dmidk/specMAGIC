#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <netcdf.h>

#include "constants.hpp"
#include "navigation.hpp"
#include "tables.hpp"
#include "types.hpp"

// Optional input fields for the time of the satellite image, already on the output grid,
// one per table: aerosol, ozone, water vapour and albedo. They replace the climatologies
// (and the land-use albedo) pixel by pixel, wherever they have a value. The CLI reads them
// from NetCDF files (load), the Python API passes them as arrays (set); either way they end
// up in the same Fields.
//
// A requested input that is not usable at all (missing file or variable, not on the output
// grid) is a problem: by default it is reported and the climatology is used, with
// strict = true the run stops. Pixels that a field does not cover only ever fall back.

namespace InputFields {

    // Stands in for NaN once the fields are read: -ffast-math lets the compiler assume
    // that NaN never occurs, so the per-pixel tests must not rely on it
    constexpr float MISSING = -9999.f;

    // Share of pixels falling back to the climatology above which the run warns about it
    constexpr double FALLBACK_WARNING_SHARE = 0.01;

    // Which tables were given
    struct Request {
        bool aerosol = false;
        bool ozone = false;
        bool water_vapour = false;
        bool albedo = false;
        bool strict = false;       // stop instead of falling back when a given input is unusable

        bool any() const { return aerosol || ozone || water_vapour || albedo; }
    };

    // The CLI's NetCDF file of each table; empty if not given
    struct Paths {
        std::string aerosol, ozone, water_vapour, albedo;
    };

    // One field on the output grid, stored row-major (lat, lon) in a plain array, as in
    // Matrix. It is allocated once while the inputs are read and never resized, so that the
    // parallel loop only ever reads fixed memory.
    struct Field {

        std::string provenance;

        Field() noexcept = default;

        Field(int rows, int cols) : rows_(rows), cols_(cols) {
            if (rows_ <= 0 || cols_ <= 0) throw std::runtime_error("Field: dimensions must be positive");
            data_ = new float[size()];
        }

        ~Field() { delete[] data_; }

        // Owns its data, so it is moved, never copied
        Field(const Field&) = delete;
        Field& operator=(const Field&) = delete;

        Field(Field&& other) noexcept
            : provenance(std::move(other.provenance)),
              rows_(other.rows_), cols_(other.cols_), data_(other.data_) {
            other.rows_ = other.cols_ = 0;
            other.data_ = nullptr;
        }

        Field& operator=(Field&& other) noexcept {
            if (this != &other) {
                delete[] data_;
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
        // every pixel, since only fields on the output grid are accepted.
        float operator()(int lat, int lon) const noexcept {
            return data_[static_cast<long>(lat) * cols_ + lon];
        }

        float* data() noexcept { return data_; }
        long size() const noexcept { return static_cast<long>(rows_) * cols_; }

    private:
        int rows_{0}, cols_{0};
        float* data_{nullptr};
    };

    // How the inputs of the day pixels were sourced, summed over threads after the loop
    struct Tally {
        uint64_t input = 0;
        uint64_t climatology = 0;
        uint64_t clipped = 0;
        uint64_t water = 0;        // albedo only: no albedo over water, which is normal

        Tally& operator+=(const Tally& o) {
            input += o.input; climatology += o.climatology;
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

        Field aod, o3, h2o, albedo_vis, albedo_nir, albedo_land;
        Request request;

        // Given inputs that are not usable, reported again after the run
        std::vector<std::string> problems;

        // Concentration axes of the LUTs. Input values are clipped to these, because unlike
        // the smooth monthly climatologies they do leave them
        MAGIC_REAL aod_min = 0, aod_max = 0;
        MAGIC_REAL h2o_min = 0, h2o_max = 0;
        MAGIC_REAL o3_min = 0, o3_max = 0;

        Inputs() = default;

        // Holds the fields, so it is moved, never copied
        Inputs(const Inputs&) = delete;
        Inputs& operator=(const Inputs&) = delete;
        Inputs(Inputs&&) = default;
        Inputs& operator=(Inputs&&) = default;

        // CLI: reads the file of every table that has a path
        void load(const Paths& paths, const Geography& geo) {

            if (!paths.aerosol.empty()) {
                request.aerosol = true;
                readFile(paths.aerosol, geo, {{"aod", &aod}});
            }
            if (!paths.ozone.empty()) {
                request.ozone = true;
                readFile(paths.ozone, geo, {{"o3", &o3}});
            }
            if (!paths.water_vapour.empty()) {
                request.water_vapour = true;
                readFile(paths.water_vapour, geo, {{"h2o", &h2o}});
            }
            if (!paths.albedo.empty()) {
                request.albedo = true;
                readFile(paths.albedo, geo, {{"albedo_vis", &albedo_vis}, {"albedo_nir", &albedo_nir},
                    {"albedo_land", &albedo_land}});
            }
        }

        // Python: copies one field that the caller has already put on the output grid. The
        // name is that of the variable in the CLI's files, e.g. "aod" or "albedo_vis".
        void set(const std::string& name, const float* values, int rows, int cols,
            const std::string& provenance, const Geography& geo) {

            Field* target = requested(name);

            if (rows != static_cast<int>(geo.nlat) || cols != static_cast<int>(geo.nlon)) {
                problem("input " + name + " (" + provenance + ") is not on the output grid of magic-config.asc");
                return;
            }

            Field f(rows, cols);
            std::copy(values, values + f.size(), f.data());
            f.provenance = provenance;
            *target = cleaned(std::move(f));
        }

        // Python: a given field that could not be put on the output grid, recorded as the
        // CLI records a file it cannot use
        void refuse(const std::string& name, const std::string& message) {
            requested(name);
            problem(message);
        }

        // After load or set: stops if any given input was unusable and the request is strict
        void check() const {
            if (request.strict && !problems.empty()) {
                std::string message = "Given input fields are not usable, and strict inputs were asked for:";
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

        // Overrides the climatology values of this pixel wherever an input has one,
        // and counts where the albedo of this pixel will come from
        void apply(PixelClimate& clim, const Area& a, Usage& u) const {

            if (request.aerosol) useInputValue(clim.aod, aod, aod_min, aod_max, u.aod, a);
            if (request.ozone) useInputValue(clim.o3, o3, o3_min, o3_max, u.o3, a);
            if (request.water_vapour) useInputValue(clim.h2o, h2o, h2o_min, h2o_max, u.h2o, a);

            if (request.albedo) {
                if (hasAlbedo(a)) u.albedo.input++;
                else if (at(albedo_land, a) == 0.f) u.albedo.water++;
                else u.albedo.climatology++;
            }
        }

        // The visible albedo covers 0.4-0.7 um and the near-infrared one 0.7-4.0 um. Both
        // must be there, so that a pixel never mixes input and land-use albedo across bands.
        // Returns -1 where there is no value, which sends the caller to the land-use albedo.
        MAGIC_REAL getKatoAlbedo(int zero_based_kato_band, const Area& a) const {

            if (!hasAlbedo(a)) return (MAGIC_REAL)-1;

            // Kato band 16 (684.2-704.4 nm, zero-based 15) lies mostly below
            // VIS_NIR_SPLIT_NM, so bands 1-16 take the visible albedo
            return at(zero_based_kato_band <= 15 ? albedo_vis : albedo_nir, a);
        }

        MAGIC_REAL getSatelliteAlbedo(double wavelength_nm, const Area& a) const {

            if (!hasAlbedo(a)) return (MAGIC_REAL)-1;
            return at(wavelength_nm < VIS_NIR_SPLIT_NM ? albedo_vis : albedo_nir, a);
        }

        // For the global attributes of the output files: where every given input came from
        std::vector<std::pair<std::string, std::string>> describe(const Usage& u) const {

            std::vector<std::pair<std::string, std::string>> attrs;

            if (request.aerosol) {
                attrs.emplace_back("aod_source", source(aod, "aerosol climatology"));
                attrs.emplace_back("aod_usage", summary(u.aod, false));
            }
            if (request.ozone) {
                attrs.emplace_back("o3_source", source(o3, "ozone climatology"));
                attrs.emplace_back("o3_usage", summary(u.o3, false));
            }
            if (request.water_vapour) {
                attrs.emplace_back("h2o_source", source(h2o, "water vapour climatology"));
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

            if (request.aerosol) report("aod", u.aod, false);
            if (request.ozone) report("o3", u.o3, false);
            if (request.water_vapour) report("h2o", u.h2o, false);
            if (request.albedo) report("albedo", u.albedo, true);

            if (problems.empty()) return;

            std::printf("\n[WARNING]: ================ UNUSABLE INPUT FIELDS ================\n");
            for (const auto& p : problems) std::printf("[WARNING]:   %s\n", p.c_str());
            std::printf("[WARNING]: These given inputs fell back to the climatologies.\n");
            std::printf("[WARNING]: Run with --strict-inputs to stop instead.\n");
            std::printf("[WARNING]: ========================================================\n\n");
        }

    private:

        Field* field(const std::string& name) {
            if (name == "aod") return &aod;
            if (name == "o3") return &o3;
            if (name == "h2o") return &h2o;
            if (name == "albedo_vis") return &albedo_vis;
            if (name == "albedo_nir") return &albedo_nir;
            if (name == "albedo_land") return &albedo_land;
            return nullptr;
        }

        // Marks the table of the named field as given, and returns the field
        Field* requested(const std::string& name) {
            Field* target = field(name);
            if (target == nullptr) throw std::runtime_error("unknown input field " + name);

            if (name == "aod") request.aerosol = true;
            else if (name == "o3") request.ozone = true;
            else if (name == "h2o") request.water_vapour = true;
            else request.albedo = true;
            return target;
        }

        void problem(const std::string& message) {
            problems.push_back(message);
            std::printf("[WARNING]: %s \n", message.c_str());
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

        // Takes the input's value, clipped to the LUT axis. Leaves the climatology value in
        // place where the input has none.
        void useInputValue(MAGIC_REAL& target, const Field& input,
            MAGIC_REAL lo, MAGIC_REAL hi, Tally& t, const Area& a) const {

            float v = at(input, a);
            if (v == MISSING) {
                t.climatology++;
                return;
            }
            t.input++;

            MAGIC_REAL x = static_cast<MAGIC_REAL>(v);
            if (x < lo || x > hi) t.clipped++;
            target = std::min(std::max(x, lo), hi);
        }

        static std::string source(const Field& f, const char* fallback) {
            return f.present() ? f.provenance : std::string(fallback) + " (input missing)";
        }

        // Same result as std::isnan, but std::isnan cannot be used: under -ffast-math the
        // compiler assumes NaN never occurs and folds std::isnan(v) to false (GCC 9 does).
        // Testing the IEEE 754 bits directly (exponent all ones, fraction non-zero) survives
        static bool isNan(float v) {
            uint32_t bits;
            std::memcpy(&bits, &v, sizeof bits);
            return (bits & 0x7f800000u) == 0x7f800000u && (bits & 0x007fffffu) != 0;
        }

        static Field cleaned(Field f) {
            float* values = f.data();
            for (long i = 0; i < f.size(); i++) {
                if (isNan(values[i])) values[i] = MISSING;
            }
            return f;
        }

        // Shares of the day pixels; for albedo, of the day pixels on land
        static double fallbackShare(const Tally& t) {
            uint64_t total = t.input + t.climatology;
            return total == 0 ? 0. : static_cast<double>(t.climatology) / total;
        }

        static std::string summary(const Tally& t, bool albedo) {
            uint64_t total = t.input + t.climatology;
            if (total == 0) return "no day pixels";

            char buf[200];
            if (albedo) {
                std::snprintf(buf, sizeof buf,
                    "%.1f%% of land day pixels from the input, %.1f%% land-use albedo; "
                    "%.1f%% of all day pixels are water and use the land-use albedo",
                    100. * t.input / total, 100. * t.climatology / total,
                    100. * t.water / (total + t.water));
            } else {
                std::snprintf(buf, sizeof buf,
                    "%.1f%% input, %.1f%% climatology, %.1f%% clipped to LUT",
                    100. * t.input / total, 100. * t.climatology / total, 100. * t.clipped / total);
            }
            return buf;
        }

        static void report(const char* name, const Tally& t, bool albedo) {
            if (t.input + t.climatology + t.water == 0) return;

            const char* level = fallbackShare(t) > FALLBACK_WARNING_SHARE ? "[WARNING]: " : "";
            std::printf("%sInput %s: %s \n", level, name, summary(t, albedo).c_str());
        }

        static bool sameAxis(const double* file_axis, const MAGIC_EXACT* geo_axis, unsigned int n) {
            for (unsigned int i = 0; i < n; i++) {
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

        // Values that a variable declares as missing, which become MISSING like NaN does
        static void maskFillValues(int ncid, int varid, Field& field) {
            for (const char* name : {"_FillValue", "missing_value"}) {
                size_t len = 0;
                float fill;
                if (nc_inq_attlen(ncid, varid, name, &len) != NC_NOERR || len != 1
                    || nc_get_att_float(ncid, varid, name, &fill) != NC_NOERR || isNan(fill)) continue;

                float* values = field.data();
                for (long i = 0; i < field.size(); i++) {
                    if (values[i] == fill) values[i] = MISSING;
                }
            }
        }

        // Reads the named variables of one table from its file, which must be on the output
        // grid, each with the dimensions (lat, lon) and stored unpacked
        void readFile(const std::string& path, const Geography& geo,
            std::initializer_list<std::pair<const char*, Field*>> fields) {

            int ncid;
            if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR) {
                problem("could not open input " + path);
                return;
            }

            int lat_dim, lon_dim;
            size_t file_nlat = 0, file_nlon = 0;
            bool ok = nc_inq_dimid(ncid, "lat", &lat_dim) == NC_NOERR
                && nc_inq_dimlen(ncid, lat_dim, &file_nlat) == NC_NOERR
                && nc_inq_dimid(ncid, "lon", &lon_dim) == NC_NOERR
                && nc_inq_dimlen(ncid, lon_dim, &file_nlon) == NC_NOERR
                && file_nlat == geo.nlat && file_nlon == geo.nlon;

            if (ok) {
                double* file_lat = new double[geo.nlat];
                double* file_lon = new double[geo.nlon];
                int varid;
                ok = nc_inq_varid(ncid, "lat", &varid) == NC_NOERR
                    && nc_get_var_double(ncid, varid, file_lat) == NC_NOERR
                    && nc_inq_varid(ncid, "lon", &varid) == NC_NOERR
                    && nc_get_var_double(ncid, varid, file_lon) == NC_NOERR
                    && sameAxis(file_lat, geo.lat, geo.nlat) && sameAxis(file_lon, geo.lon, geo.nlon);
                delete[] file_lat;
                delete[] file_lon;
            }

            if (!ok) {
                problem("input " + path + " is not on the output grid of magic-config.asc");
                nc_close(ncid);
                return;
            }

            for (const auto& f : fields) {

                const std::string name = f.first;
                int varid;
                if (nc_inq_varid(ncid, f.first, &varid) != NC_NOERR) {
                    problem(name + " is not in input " + path);
                    continue;
                }

                // Exactly (lat, lon): the whole variable is read into a field of that size
                int ndims = 0;
                int dimids[NC_MAX_VAR_DIMS];
                if (nc_inq_varndims(ncid, varid, &ndims) != NC_NOERR || ndims != 2
                    || nc_inq_vardimid(ncid, varid, dimids) != NC_NOERR
                    || dimids[0] != lat_dim || dimids[1] != lon_dim) {
                    problem(name + " in input " + path + " must have exactly the dimensions (lat, lon)");
                    continue;
                }

                // nc_get_var_float does not unpack, so packed values would be read raw
                int attid;
                if (nc_inq_attid(ncid, varid, "scale_factor", &attid) == NC_NOERR
                    || nc_inq_attid(ncid, varid, "add_offset", &attid) == NC_NOERR) {
                    problem(name + " in input " + path + " is packed (scale_factor/add_offset); "
                        "write it unpacked, as float");
                    continue;
                }

                Field field(static_cast<int>(geo.nlat), static_cast<int>(geo.nlon));
                if (nc_get_var_float(ncid, varid, field.data()) != NC_NOERR) {
                    problem("could not read " + name + " from " + path);
                    continue;
                }
                maskFillValues(ncid, varid, field);

                field.provenance = textAttribute(ncid, varid, "provenance");
                if (field.provenance.empty()) field.provenance = path;
                *f.second = cleaned(std::move(field));
            }

            nc_close(ncid);
        }
    };

}
