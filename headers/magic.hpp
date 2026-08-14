#pragma once
#include <stdexcept>
#include <cmath>
#include <array>
#include <cstdio>
#include "tables.hpp"
#include "radiation.hpp"
#include "navigation.hpp"


using namespace Tables;
using namespace Correction;


namespace magicHelpers {

    int borders(MAGIC_REAL* xa, int n, MAGIC_REAL x);

    MAGIC_REAL upperWeight(MAGIC_REAL* axis, int hi, MAGIC_REAL x);

    int nearestGgIndex(RTM& rtm, MAGIC_REAL gg_value);

    MAGIC_REAL interpolate(Absorber& lut, MAGIC_REAL* axis,
        MAGIC_REAL* const* table, MAGIC_REAL x, int band);

    MAGIC_REAL absorberCorrection(Absorber& lut, IrradianceMode mode,
        int band, MAGIC_REAL concentration, MAGIC_REAL cos_sza);
}

using namespace magicHelpers;

MAGIC_REAL magic(RTM& rtm, Absorber& waterLut, Absorber& ozoneLut,
    IrradianceMode mode, int band, MAGIC_REAL cos_sza, MAGIC_REAL f_ext,
    PixelClimate climate);

void wavelengthCorrectionCloudySky(Correction::Spectral& spectral,
    MAGIC_REAL k, MAGIC_REAL* GHI);