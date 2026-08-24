#pragma once

#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <netcdf.h>

#include "helpers.hpp"
#include "navigation.hpp"
#include "types.hpp"

namespace ModisBrdf {

    enum class Source {
        Band1 = 0,
        Band2,
        Band3,
        Band4,
        Band5,
        Band6,
        Band7,
        Vis,
        Nir,
        Shortwave,
        Count
    };

    struct BrdfTriplet {
        MAGIC_REAL fiso = 0;
        MAGIC_REAL fvol = 0;
        MAGIC_REAL fgeo = 0;
        bool valid = false;
    };

    inline const char* sourceName(Source source) {
        switch (source) {
            case Source::Band1: return "Band1";
            case Source::Band2: return "Band2";
            case Source::Band3: return "Band3";
            case Source::Band4: return "Band4";
            case Source::Band5: return "Band5";
            case Source::Band6: return "Band6";
            case Source::Band7: return "Band7";
            case Source::Vis: return "vis";
            case Source::Nir: return "nir";
            case Source::Shortwave: return "shortwave";
            default: return "";
        }
    }

    inline std::string filenameFor(Source source) {
        return std::string("Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_")
            + sourceName(source) + ".nc";
    }

    inline Source sourceForSatelliteWavelength(double wavelength_nm) {
        if (wavelength_nm >= 620 && wavelength_nm <= 670) return Source::Band1;
        if (wavelength_nm >= 841 && wavelength_nm <= 876) return Source::Band2;
        if (wavelength_nm >= 459 && wavelength_nm <= 479) return Source::Band3;
        if (wavelength_nm >= 545 && wavelength_nm <= 565) return Source::Band4;
        if (wavelength_nm >= 1230 && wavelength_nm <= 1250) return Source::Band5;
        if (wavelength_nm >= 1628 && wavelength_nm <= 1652) return Source::Band6;
        if (wavelength_nm >= 2105 && wavelength_nm <= 2155) return Source::Band7;
        if (wavelength_nm < 700) return Source::Vis;
        if (wavelength_nm < 2100) return Source::Nir;
        return Source::Shortwave;
    }

    inline Source sourceForKatoBand(int band) {
        // band is expected as zero-based Kato index, 0..31.
        // This mapping is deliberately conservative:
        // narrow MODIS bands where suitable, broad VIS/NIR/SW otherwise.
        if (band <= 7) return Source::Vis;          // UV/blue/visible
        if (band <= 11) return Source::Vis;         // green/yellow visible
        if (band <= 14) return Source::Band1;       // red, around 620-670 nm
        if (band <= 19) return Source::Band2;       // near 865 nm
        if (band <= 22) return Source::Nir;
        if (band == 23) return Source::Band5;
        if (band == 24 || band == 25) return Source::Band6;
        if (band == 26 || band == 27) return Source::Band7;
        return Source::Shortwave;
    }

    inline MAGIC_REAL blackSkyAlbedo(
        MAGIC_REAL fiso,
        MAGIC_REAL fvol,
        MAGIC_REAL fgeo,
        MAGIC_REAL solar_zenith_angle
    ) {
        MAGIC_REAL theta = solar_zenith_angle;
        MAGIC_REAL theta2 = theta * theta;
        MAGIC_REAL theta3 = theta2 * theta;

        MAGIC_REAL vol_kernel =
            (MAGIC_REAL)-0.007574
            - (MAGIC_REAL)0.070987 * theta2
            + (MAGIC_REAL)0.307588 * theta3;

        MAGIC_REAL geo_kernel =
            (MAGIC_REAL)-1.284909
            - (MAGIC_REAL)0.166314 * theta2
            + (MAGIC_REAL)0.041840 * theta3;

        return fiso + fvol * vol_kernel + fgeo * geo_kernel;
    }

    inline MAGIC_REAL whiteSkyAlbedo(
        MAGIC_REAL fiso,
        MAGIC_REAL fvol,
        MAGIC_REAL fgeo
    ) {
        return fiso
            + (MAGIC_REAL)0.189184 * fvol
            - (MAGIC_REAL)1.377622 * fgeo;
    }

    inline MAGIC_REAL diffuseFraction(MAGIC_REAL cos_sza) {
        if (cos_sza <= (MAGIC_REAL)0) return (MAGIC_REAL)1;

        MAGIC_REAL direct =
            (MAGIC_REAL)1098.0
            * std::exp((MAGIC_REAL)-0.057 / cos_sza)
            * cos_sza;

        MAGIC_REAL diffuse =
            (MAGIC_REAL)94.23 * std::sqrt(cos_sza);

        MAGIC_REAL denom = direct + diffuse;
        if (denom <= (MAGIC_REAL)0) return (MAGIC_REAL)1;

        MAGIC_REAL fraction = (MAGIC_REAL)2.0 * diffuse / denom;
        return std::min((MAGIC_REAL)1.0, std::max((MAGIC_REAL)0.0, fraction));
    }

    inline MAGIC_REAL surfaceAlbedo(
        MAGIC_REAL fiso,
        MAGIC_REAL fvol,
        MAGIC_REAL fgeo,
        MAGIC_REAL cos_sza
    ) {
        MAGIC_REAL theta = std::acos(std::min((MAGIC_REAL)1, std::max((MAGIC_REAL)0, cos_sza)));

        MAGIC_REAL black = blackSkyAlbedo(fiso, fvol, fgeo, theta);
        MAGIC_REAL white = whiteSkyAlbedo(fiso, fvol, fgeo);
        MAGIC_REAL diffuse = diffuseFraction(cos_sza);

        MAGIC_REAL albedo = ((MAGIC_REAL)1 - diffuse) * black + diffuse * white;

        if (!std::isfinite(albedo)) return (MAGIC_REAL)-1;
        if (albedo < (MAGIC_REAL)0 || albedo > (MAGIC_REAL)1) return (MAGIC_REAL)-1;

        return albedo;
    }

    class ModisBrdfAlbedo {
    public:
        static constexpr int SourceCount = static_cast<int>(Source::Count);
        static constexpr int FullLon = 7200;
        static constexpr int FullLat = 3600;
        static constexpr MAGIC_REAL Scale = (MAGIC_REAL)0.001;
        static constexpr short FillValue = 32767;

        bool enabled = false;
        bool loaded = false;

        int nlon = 0;
        int nlat = 0;
        int lon_start = 0;
        int lat_start = 0;

        MAGIC_REAL lon_min = 0;
        MAGIC_REAL lat_max = 0;

        std::array<std::vector<MAGIC_REAL>, SourceCount> fiso;
        std::array<std::vector<MAGIC_REAL>, SourceCount> fvol;
        std::array<std::vector<MAGIC_REAL>, SourceCount> fgeo;
        std::array<bool, SourceCount> source_loaded{};

        ModisBrdfAlbedo() = default;

        void load(Config c, int month) {
            enabled = c.use_modis_brdf_albedo == 1;

            if (!enabled) {
                std::printf("MODIS BRDF albedo disabled; using land-use albedo.\n");
                return;
            }

            if (month < 1 || month > 12) {
                std::printf("Invalid month for MODIS BRDF albedo: %d\n", month);
                enabled = false;
                return;
            }

            defineRegion(c);

            for (int s = 0; s < SourceCount; ++s) {
                Source source = static_cast<Source>(s);
                std::string filepath = c.path + c.modis_brdf_dir + filenameFor(source);

                bool ok = readSource(filepath, source, month);
                source_loaded[s] = ok;

                if (!ok) {
                    std::printf("WARNING: Could not read MODIS BRDF source %s from %s\n",
                        sourceName(source), filepath.c_str());
                }
            }

            loaded = true;
        }

        MAGIC_REAL getAlbedo(Source source, Area a, MAGIC_REAL cos_sza) const {
            int s = static_cast<int>(source);
            if (!enabled || !loaded || !source_loaded[s]) return (MAGIC_REAL)-1;

            BrdfTriplet brdf = interpolate(source, a.degrees_lat, a.degrees_lon);
            if (!brdf.valid) return (MAGIC_REAL)-1;

            return surfaceAlbedo(brdf.fiso, brdf.fvol, brdf.fgeo, cos_sza);
        }

        MAGIC_REAL getKatoAlbedo(int zero_based_kato_band, Area a, MAGIC_REAL cos_sza) const {
            return getAlbedo(sourceForKatoBand(zero_based_kato_band), a, cos_sza);
        }

        MAGIC_REAL getSatelliteAlbedo(double wavelength_nm, Area a, MAGIC_REAL cos_sza) const {
            return getAlbedo(sourceForSatelliteWavelength(wavelength_nm), a, cos_sza);
        }

    private:
        void defineRegion(Config c) {
            MAGIC_REAL region_lon_min = c.lonbegin + c.deltalon;
            MAGIC_REAL region_lon_max = region_lon_min + c.dxy * (c.londim - 1);
            MAGIC_REAL region_lat_min = c.latbegin;
            MAGIC_REAL region_lat_max = c.latbegin + c.dxy * (c.latdim - 1);

            lon_start = lonIndex(region_lon_min) - 2;
            int lon_end = lonIndex(region_lon_max) + 2;

            lat_start = latIndex(region_lat_max) - 2;
            int lat_end = latIndex(region_lat_min) + 2;

            lon_start = std::max(0, lon_start);
            lat_start = std::max(0, lat_start);
            lon_end = std::min(FullLon - 1, lon_end);
            lat_end = std::min(FullLat - 1, lat_end);

            nlon = lon_end - lon_start + 1;
            nlat = lat_end - lat_start + 1;

            lon_min = lonAt(lon_start);
            lat_max = latAt(lat_start);

            std::printf("MODIS BRDF region: lon_start=%d nlon=%d lat_start=%d nlat=%d\n",
                lon_start, nlon, lat_start, nlat);
        }

        static int lonIndex(MAGIC_REAL lon) {
            MAGIC_REAL wrapped = lon;
            while (wrapped < (MAGIC_REAL)-180) wrapped += (MAGIC_REAL)360;
            while (wrapped > (MAGIC_REAL)180) wrapped -= (MAGIC_REAL)360;

            int index = static_cast<int>(std::lround((wrapped - (MAGIC_REAL)-179.975) / (MAGIC_REAL)0.05));
            return std::max(0, std::min(FullLon - 1, index));
        }

        static int latIndex(MAGIC_REAL lat) {
            int index = static_cast<int>(std::lround(((MAGIC_REAL)89.975 - lat) / (MAGIC_REAL)0.05));
            return std::max(0, std::min(FullLat - 1, index));
        }

        static MAGIC_REAL lonAt(int index) {
            return (MAGIC_REAL)-179.975 + (MAGIC_REAL)0.05 * index;
        }

        static MAGIC_REAL latAt(int index) {
            return (MAGIC_REAL)89.975 - (MAGIC_REAL)0.05 * index;
        }

        int localIndex(int local_lon, int local_lat) const {
            return local_lon * nlat + local_lat;
        }

        bool readSource(const std::string& filepath, Source source, int month) {
            int ncid = 0;
            int retval = nc_open(filepath.c_str(), NC_NOWRITE, &ncid);
            if (retval != NC_NOERR) return false;

            bool ok = readVariable(ncid, "FISO_MEAN", fiso[static_cast<int>(source)], month)
                && readVariable(ncid, "FVOL_MEAN", fvol[static_cast<int>(source)], month)
                && readVariable(ncid, "FGEO_MEAN", fgeo[static_cast<int>(source)], month);

            nc_close(ncid);
            return ok;
        }

        bool readVariable(int ncid, const char* varname, std::vector<MAGIC_REAL>& out, int month) {
            int varid = 0;
            int retval = nc_inq_varid(ncid, varname, &varid);
            if (retval != NC_NOERR) return false;

            std::vector<short> raw(nlon * nlat);

            size_t start[3] = {
                static_cast<size_t>(month - 1),
                static_cast<size_t>(lon_start),
                static_cast<size_t>(lat_start)
            };

            size_t count[3] = {
                1,
                static_cast<size_t>(nlon),
                static_cast<size_t>(nlat)
            };

            retval = nc_get_vara_short(ncid, varid, start, count, raw.data());
            if (retval != NC_NOERR) return false;

            out.resize(raw.size());

            for (size_t i = 0; i < raw.size(); ++i) {
                if (raw[i] == FillValue) out[i] = (MAGIC_REAL)-999;
                else out[i] = static_cast<MAGIC_REAL>(raw[i]) * Scale;
            }

            return true;
        }

        BrdfTriplet interpolate(Source source, MAGIC_REAL lat, MAGIC_REAL lon) const {
            int s = static_cast<int>(source);

            int global_lon = lonIndex(lon);
            int global_lat = latIndex(lat);

            int local_lon = global_lon - lon_start;
            int local_lat = global_lat - lat_start;

            if (local_lon < 0 || local_lon >= nlon || local_lat < 0 || local_lat >= nlat) {
                return {};
            }

            int idx = localIndex(local_lon, local_lat);

            BrdfTriplet result;
            result.fiso = fiso[s][idx];
            result.fvol = fvol[s][idx];
            result.fgeo = fgeo[s][idx];

            result.valid =
                result.fiso >= 0
                && result.fvol >= 0
                && result.fgeo >= 0
                && std::isfinite(result.fiso)
                && std::isfinite(result.fvol)
                && std::isfinite(result.fgeo);

            return result;
        }
    };

}