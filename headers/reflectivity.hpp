#pragma once

#include <stdexcept>
#include "constants.hpp"
#include "tables.hpp"
#include "image.hpp"
#include "sun_geometry.hpp"
#include "modis_brdf.hpp"
#include "input_fields.hpp"

using namespace Tables;

namespace Reflectivity {

    int wavelengthToKatoBand(double lambda);

    MAGIC_REAL cloudIndex(MAGIC_REAL radi, MAGIC_REAL rmin, MAGIC_REAL rmax);

    MAGIC_REAL clearSkyIndex(MAGIC_REAL CAL);

    MAGIC_REAL zenithAdjustment(int land, MAGIC_REAL cos_zen);

    MAGIC_REAL zenithAdjustmentVis(int land, MAGIC_REAL sza);

    MAGIC_REAL zenithCorrection(int land_class, MAGIC_REAL cos_zenith_angle);

    void makeAlbedoCorrection(Climate& climatologies, PixelClimate& c, Area a, MAGIC_REAL cos_sza, MAGIC_REAL declination,
        bool input_albedo);

    MAGIC_REAL getSurfaceAlbedo(Climate& climatologies, Area a, GroundAlbedo& alb, int band);

    MAGIC_REAL getBestKatoSurfaceAlbedo(Climate& climatologies, Area a, AlbedoType alb_source,
        GroundAlbedo& alb, const ModisBrdf::ModisBrdfAlbedo& modis, const InputFields::Inputs& inputs, int band, MAGIC_REAL cos_sza, MAGIC_REAL fallback_correction,
        MAGIC_REAL noon_correction);

    MAGIC_REAL getBestSatelliteSurfaceAlbedo(Climate& climatologies, Area a, AlbedoType alb_source, GroundAlbedo& alb,
        const ModisBrdf::ModisBrdfAlbedo& modis, const InputFields::Inputs& inputs, double wavelength_nm, MAGIC_REAL cos_sza, MAGIC_REAL fallback_correction,
        MAGIC_REAL noon_correction);

}

MAGIC_REAL effectiveCloudAlbedo(Image& img, SolarParameters sun, Climate& climatologies, AlbedoType alb_source, GroundAlbedo& alb,
    const ModisBrdf::ModisBrdfAlbedo& modis, const InputFields::Inputs& inputs, Area a, PixelClimate& clim, int line, int col);