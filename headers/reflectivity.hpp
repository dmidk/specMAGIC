#pragma once

#include <stdexcept>
#include "constants.hpp"
#include "tables.hpp"
#include "image.hpp"
#include "sun_geometry.hpp"
#include "modis_brdf.hpp"

using namespace Tables;

namespace Reflectivity {

    int wavelengthToKatoBand(double lambda);

    MAGIC_REAL cloudIndex(MAGIC_REAL radi, MAGIC_REAL rmin, MAGIC_REAL rmax);

    MAGIC_REAL clearSkyIndex(MAGIC_REAL CAL);

    MAGIC_REAL zenithAdjustment(int land, MAGIC_REAL cos_zen);

    MAGIC_REAL zenithAdjustmentVis(int land, MAGIC_REAL sza);

    MAGIC_REAL zenithCorrection(int land_class, MAGIC_REAL cos_zenith_angle);

    void makeAlbedoCorrection(Climate& climatologies, PixelClimate& c, Area a, MAGIC_REAL cos_sza);

    MAGIC_REAL getSurfaceAlbedo(Climate& climatologies, Area a, GroundAlbedo& alb, int band);

    MAGIC_REAL getBestKatoSurfaceAlbedo(
        Climate& climatologies,
        Area a,
        GroundAlbedo& alb,
        const ModisBrdf::ModisBrdfAlbedo& modis,
        int band,
        MAGIC_REAL cos_sza,
        MAGIC_REAL fallback_correction
    );

    MAGIC_REAL getBestSatelliteSurfaceAlbedo(
        Climate& climatologies,
        Area a,
        GroundAlbedo& alb,
        const ModisBrdf::ModisBrdfAlbedo& modis,
        double wavelength_nm,
        MAGIC_REAL cos_sza,
        MAGIC_REAL fallback_correction
    );

}

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
    );