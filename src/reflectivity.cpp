#include "../headers/reflectivity.hpp"


namespace Reflectivity {

    int wavelengthToKatoBand(double lambda) {

        // Input: wavelength IN NANOMETRES

        // Idiot check
        if (lambda < 240.1 || lambda > 4605.7)
            return -1;

        if (lambda < 272.5) return 1;
        if (lambda < 283.4) return 2;
        if (lambda < 306.8) return 3;
        if (lambda < 327.8) return 4;
        if (lambda < 362.5) return 5;
        if (lambda < 407.5) return 6;
        if (lambda < 452.0) return 7;
        if (lambda < 517.7) return 8;
        if (lambda < 540.0) return 9;
        if (lambda < 549.5) return 10;
        if (lambda < 566.6) return 11;
        if (lambda < 605.0) return 12;
        if (lambda < 625.0) return 13;
        if (lambda < 666.7) return 14;
        if (lambda < 684.2) return 15;
        if (lambda < 704.4) return 16;
        if (lambda < 742.6) return 17;
        if (lambda < 791.5) return 18;
        if (lambda < 844.5) return 19;
        if (lambda < 889.0) return 20;
        if (lambda < 974.9) return 21;
        if (lambda < 1045.7) return 22;
        if (lambda < 1194.2) return 23;
        if (lambda < 1515.9) return 24;
        if (lambda < 1613.5) return 25;
        if (lambda < 1964.8) return 26;
        if (lambda < 2153.5) return 27;
        if (lambda < 2275.2) return 28;
        if (lambda < 3001.9) return 29;
        if (lambda < 3635.4) return 30;
        if (lambda < 3991.0) return 31;

        return 32;
    }

    MAGIC_REAL cloudIndex(MAGIC_REAL radi, MAGIC_REAL rmin, MAGIC_REAL rmax) {

        if (rmax == rmin) return 0.0;
        return (radi - rmin) / (rmax - rmin);

    }

    MAGIC_REAL clearSkyIndex(MAGIC_REAL CAL) {

        if (CAL <= (MAGIC_REAL) CAL_MIN_CLEAR)
            return (MAGIC_REAL) K_CLEAR_EXTREME;

        if (CAL <= (MAGIC_REAL) CAL_LINEAR_MAX) {
            if (CAL <= (MAGIC_REAL) CAL_NOISE_THRESHOLD)
                return (MAGIC_REAL) 1.0;

            return (MAGIC_REAL) 1.0 - CAL;
        }

        if (CAL <= (MAGIC_REAL) CAL_TRANSITION_MAX) {
            return (MAGIC_REAL) K_A0
                + (MAGIC_REAL) K_A1 * CAL
                + (MAGIC_REAL) K_A2 * CAL * CAL;
        }

        return (MAGIC_REAL) K_OVERCAST;
    }


    // Some kind of Lambert corrections
    // Reflection can seem different at low angles near horizon

    // Used over ocean
    MAGIC_REAL zenithAdjustment(int land, MAGIC_REAL cos_zen) {

        if (cos_zen < (MAGIC_REAL)0.0 || cos_zen > (MAGIC_REAL)1.0)
            throw std::runtime_error("cos_zen out of range");

        static constexpr MAGIC_REAL d[20] = {
            0.40, 0.44, 0.32, 0.39, 0.22,
            0.28, 0.40, 0.47, 0.53, 0.53,
            0.35, 0.41, 0.10, 0.40, 0.10,
            0.40, 0.41, 0.58, 0.10, 0.10
        };

        if (land < 0 || land >= 20)
            throw std::runtime_error("invalid land class");

        return ((MAGIC_REAL)1.0 + d[land]) / ((MAGIC_REAL) 1.0 + (MAGIC_REAL) 2.0 * d[land] * cos_zen);
    }

    // Used over anything non ocean
    MAGIC_REAL zenithAdjustmentVis(int land, MAGIC_REAL sza) {
        if (sza < (MAGIC_REAL) 0.0 || sza > (MAGIC_REAL) 1.571)
            throw std::runtime_error("sza out of range");

        static constexpr MAGIC_REAL B1[20] = {
            0.092,0.092,0.052,0.068,0.067,
            0.106,0.128,0.095,0.1,0.136,
            0.35,0.135,0.10,0.135,0.10,
            0.135,0.10,0.10,0.10,0.10
        };

        static constexpr MAGIC_REAL B2[20] = {
            0.077,0.077,0.1,0.1,0.086,
            0.062,0.056,0.090,0.058,0.018,
            0.1,0.024,0.10,0.024,0.10,
            0.005,0.10,0.10,0.10,0.10
        };

        if (land < 0 || land >= 20)
            throw std::runtime_error("invalid land class");

        auto g1 = (MAGIC_REAL) -0.007574 - (MAGIC_REAL) 0.070987 * sza * sza + (MAGIC_REAL) 0.307588 * sza * sza * sza;
        auto g2 = (MAGIC_REAL) -1.284909 - (MAGIC_REAL) 0.166314 * sza * sza + (MAGIC_REAL) 0.04184 * sza * sza * sza;

        constexpr MAGIC_REAL g1n = (MAGIC_REAL) 0.26781;
        constexpr MAGIC_REAL g2n = (MAGIC_REAL) -1.4192;

        constexpr MAGIC_REAL g1m = (MAGIC_REAL) 0.06244;
        constexpr MAGIC_REAL g2m = (MAGIC_REAL) -1.3517;

        MAGIC_REAL renorm =
            (MAGIC_REAL) 1.0 / ((MAGIC_REAL) 1.0
                + B1[land]*(g1m - g1n)
                + B2[land]*(g2m - g2n));

        return renorm *
            ((MAGIC_REAL) 1.0
                + B1[land]*(g1 - g1n)
                + B2[land]*(g2 - g2n));
    }

    MAGIC_REAL zenithCorrection(int land_class, MAGIC_REAL cos_zenith_angle) {

        MAGIC_REAL z = 0; 

        // TODO: Check the indexing here. 
        // For consistency, it should be land_class -1, but this 
        // differs from specmagic now. 
        // There may be a bug in specmagic now.

        if (land_class == 17) z = zenithAdjustment(land_class-1, cos_zenith_angle);
        else z = zenithAdjustmentVis(land_class-1, acos(cos_zenith_angle));

        assert(z > 0);

        return z;
    }

    void makeAlbedoCorrection(Climate& climatologies, PixelClimate& c, Area a, MAGIC_REAL cos_sza) {

        int idx_lu = climatologies.lu.index(a.degrees_lat, a.degrees_lon);
        int land_class = climatologies.lu.grid[idx_lu];   // 1..20

        c.albedo_correction = zenithCorrection(land_class, cos_sza);
    }

    MAGIC_REAL getSurfaceAlbedo(Climate& climatologies, Area a, GroundAlbedo& alb, 
        int band) {
        
        // Index in land use table
        int idx_lu = climatologies.lu.index(a.degrees_lat, a.degrees_lon);
        int land_class = climatologies.lu.grid[idx_lu];   // 1..20
        int land_index = land_class - 1;                  // 0..19
        float albedo = alb.getAlbedo(land_index, band);

        return albedo;


    }

    MAGIC_REAL getBestKatoSurfaceAlbedo(
       Climate& climatologies,
       Area a,
       GroundAlbedo& alb,
       const ModisBrdf::ModisBrdfAlbedo& modis,
       int band,
       MAGIC_REAL cos_sza,
       MAGIC_REAL fallback_correction
   ) {
        MAGIC_REAL modis_albedo = modis.getKatoAlbedo(band, a, cos_sza);

        if (modis_albedo >= (MAGIC_REAL)0 && modis_albedo <= (MAGIC_REAL)1) {
            return modis_albedo;
        }

        return getSurfaceAlbedo(climatologies, a, alb, band) * fallback_correction;
    }

    MAGIC_REAL getBestSatelliteSurfaceAlbedo(
        Climate& climatologies,
        Area a,
        GroundAlbedo& alb,
        const ModisBrdf::ModisBrdfAlbedo& modis,
        double wavelength_nm,
        MAGIC_REAL cos_sza,
        MAGIC_REAL fallback_correction
    ) {
        MAGIC_REAL modis_albedo = modis.getSatelliteAlbedo(wavelength_nm, a, cos_sza);

        if (modis_albedo >= (MAGIC_REAL)0 && modis_albedo <= (MAGIC_REAL)1) {
            return modis_albedo;
        }

        int band = Reflectivity::wavelengthToKatoBand(wavelength_nm);
        band--;

        return getSurfaceAlbedo(climatologies, a, alb, band) * fallback_correction;
    }

  
}

// Big monster func to calculate CAL for this pixel
MAGIC_REAL effectiveCloudAlbedo(
    Image& img,
    SolarParameters sun,
    Climate& climatologies,
    GroundAlbedo& alb,
    const ModisBrdf::ModisBrdfAlbedo& modis,
    Area a,
    PixelClimate& clim,
    int line,
    int col
) {

    MAGIC_REAL Rmax = 3500;

    // Radiance, scaled by angle
    MAGIC_REAL radi = (static_cast<MAGIC_REAL>(img.im.at(line, col)) - (MAGIC_REAL) DARK_OFFSET) / sun.corrected_cos_sza;

    MAGIC_REAL albedo = Reflectivity::getBestSatelliteSurfaceAlbedo(
            climatologies,
            a,
            alb,
            modis,
            img.info.wavelength,
            sun.cos_sza,
            clim.albedo_correction
        );

    MAGIC_REAL Rmin = 0;
    Rmin = (MAGIC_REAL) 1.1 * Rmax * (albedo + (MAGIC_REAL) 0.05);

    MAGIC_REAL CAL_continuous = Reflectivity::cloudIndex(radi, Rmin, Rmax);

    return CAL_continuous;
}

