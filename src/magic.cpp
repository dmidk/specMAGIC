#include "../headers/magic.hpp"


namespace magicHelpers {

    // Identical to borders member function of Climatology, maybe we want to combine somewhere
    inline int borders(MAGIC_REAL* xa, int n, MAGIC_REAL x) {
        
        int klo = 0;
        int khi = n - 1;
        MAGIC_REAL sign = (xa[n - 1] < xa[0]) ? -1.0f : 1.0f;

        while (khi - klo > 1) {
            int k = (khi + klo) / 2;
            if (sign * xa[k] > sign * x) {
                khi = k;
            } else {
                klo = k;
            }
        }
        return khi;
    }

    inline MAGIC_REAL upperWeight(MAGIC_REAL* axis, int hi, MAGIC_REAL x) {
        
        MAGIC_REAL denom = axis[hi] - axis[hi - 1];
        if (denom == 0.0f) {
            throw std::runtime_error("upperWeight: repeated adjacent axis values");
        }

        return 1.0f - ((axis[hi] - x) / denom);
    }

    inline int nearestGgIndex(RTM& rtm, MAGIC_REAL gg_value) {
        
        int gg_hi = borders(rtm.gg, rtm.ggdim, gg_value);

        MAGIC_REAL denom = rtm.gg[gg_hi] - rtm.gg[gg_hi - 1];
        if (denom == 0.0f) {
            throw std::runtime_error("nearestGgIndex: repeated adjacent gg values");
        }

        MAGIC_REAL w = (rtm.gg[gg_hi] - gg_value) / denom;

        // Specmagic now uses the rule:
        // if closer to upper node -> gg_hi
        // else -> gg_hi - 1
        return (w < 0.5f) ? gg_hi : (gg_hi - 1);
    }

    // Not completely identical to the member function of Climatology
    // potentially worth checking...
    inline MAGIC_REAL interpolate(Absorber& lut, MAGIC_REAL* axis,
        MAGIC_REAL* const* table, MAGIC_REAL x, int band) {

        int klo = 0;
        int khi = lut.nrows - 1;

        while (khi - klo > 1) {
            int k = (khi + klo) / 2;
            if (axis[k] > x) {
                khi = k;
            } else {
                klo = k;
            }
        }

        MAGIC_REAL denom = axis[khi] - axis[klo];
        if (denom == 0.0f) {
            throw std::runtime_error("interpolate: repeated adjacent absorber axis values");
        }

        // This is the weight in specmagic-now:
        // w = 1 at lower node, 0 at upper node
        MAGIC_REAL w = (axis[khi] - x) / denom;

        if (w >= 0.0f && w <= 1.0f) {
            // Exact specmagic formula:
            return (1.0f - w) * table[khi][band] + w * table[klo][band];
        }

        printf("[MAGIC] WARNING: Dodgy interpolation in Absorber calculation. Quality may be degraded.\n");
        return -1000.0f;
    }


    inline MAGIC_REAL absorberCorrection(Absorber& lut, IrradianceMode mode,
        int band, MAGIC_REAL concentration, MAGIC_REAL cos_sza) {

        // DNI calculation
        if (mode == IrradianceMode::Beam) {
            return interpolate(lut, lut.concen, lut.delta_beam,
                                            concentration, band)
                * std::pow(cos_sza, lut.ab[band]);
            
        // GHI Calculation
        } else if (mode == IrradianceMode::Global) {
            return interpolate(lut, lut.concen, lut.delta_global,
                                            concentration, band)
                * std::pow(cos_sza, lut.ag[band]);
        } else { 
            assert(!"Got undefined IrradianceMode!");
        }

    }

}


using namespace magicHelpers;

MAGIC_REAL magic(RTM& rtm, Absorber& waterLut, Absorber& ozoneLut,
    IrradianceMode mode, int band, MAGIC_REAL cos_sza, MAGIC_REAL f_ext,
    PixelClimate climate) {


    // dumb sanity checks
    if (band < 0 || band >= rtm.bandsdim) {
        throw std::out_of_range("magic: band index out of range");
    }
    if (waterLut.ncols != rtm.bandsdim || ozoneLut.ncols != rtm.bandsdim) {
        throw std::runtime_error("magic: incompatible LUT band dimensions");
    }

    MAGIC_REAL Gtmp[2];
    MAGIC_REAL Gmlbmc = 0.0f;

    // Choose nearest neighbour
    int i_gg = nearestGgIndex(rtm, climate.gg);

    // AOD interpolation fraction 
    int aod_kh = borders(rtm.aod, rtm.aoddim, climate.aod);
    MAGIC_REAL w = upperWeight(rtm.aod, aod_kh, climate.aod);

    // SSA interpolation fraction 
    int ssa_kh = borders(rtm.ssa, rtm.ssadim, climate.ssa);
    MAGIC_REAL ws = upperWeight(rtm.ssa, ssa_kh, climate.ssa);

    // Calculate rtm at these two nodes
    // Mapping:
    //   i = aod_kh     -> Gtmp[0] (upper AOD)
    //   i = aod_kh - 1 -> Gtmp[1] (lower AOD)
    for (int i = aod_kh - 1; i <= aod_kh; ++i) {
        const int out_idx = aod_kh - i;

        if (mode == IrradianceMode::Beam) {
            Gtmp[out_idx] =
                ws * rtm.Im[i_gg][ssa_kh][i][band]
                   * std::exp(
                         rtm.btau[i_gg][ssa_kh][i][band] /
                         std::pow(cos_sza, rtm.ab[i_gg][ssa_kh][i][band])
                     )
              + (1.0f - ws) * rtm.Im[i_gg][ssa_kh - 1][i][band]
                   * std::exp(
                         rtm.btau[i_gg][ssa_kh - 1][i][band] /
                         std::pow(cos_sza, rtm.ab[i_gg][ssa_kh - 1][i][band])
                     );
        } else {
            Gtmp[out_idx] =
                ws * rtm.Im[i_gg][ssa_kh][i][band]
                   * std::exp(
                         rtm.gtau[i_gg][ssa_kh][i][band] /
                         std::pow(cos_sza, rtm.ag[i_gg][ssa_kh][i][band])
                     )
              + (1.0f - ws) * rtm.Im[i_gg][ssa_kh - 1][i][band]
                   * std::exp(
                         rtm.gtau[i_gg][ssa_kh - 1][i][band] /
                         std::pow(cos_sza, rtm.ag[i_gg][ssa_kh - 1][i][band])
                     );
        }
    }


    // Interpolate in AOD and apply f_ext * cos_sza
    //   Gtmp[0] = upper AOD
    //   Gtmp[1] = lower AOD
    // and w is weight toward upper AOD,
    Gmlbmc = f_ext * cos_sza * (w * Gtmp[0] + (1.0f - w) * Gtmp[1]);

    // Add on the absorber corrections
    Gmlbmc += absorberCorrection(waterLut, mode, band, climate.h2o, cos_sza);
    Gmlbmc += absorberCorrection(ozoneLut, mode, band, climate.o3, cos_sza);

    // Add the surface albedo correction
    // Not relevant for the Beam case
    if (mode == IrradianceMode::Global) {
        Gmlbmc *= (0.98 + 0.1 * static_cast<MAGIC_REAL>( climate.surface_albedo) );
    }

    // Unpleasant clamps
    // TODO: Do we need these?
    if (cos_sza < (MAGIC_REAL)0.025) {
        Gmlbmc = 0.0;
    }

    if (Gmlbmc < (MAGIC_REAL)0.0) {
        Gmlbmc = 0.0;
    }

    return Gmlbmc;
}



void wavelengthCorrectionCloudySky(Correction::Spectral& spectral,
    MAGIC_REAL k, MAGIC_REAL* GHI) {

    int kh = 0;
    int kl = 0;
    MAGIC_REAL mkd = 0.0;

    // Find the interval in which k lies
    if (k < spectral.k_COD[6]) {
        // COD > 80: assume same wavelength correction as for 80
        kh = 6;
        kl = 6;
        mkd = 0.0;
    }
    else if (k >= spectral.k_COD[2]) {
        // Note: preserves original >= boundary behavior
        kh = 1;
        kl = 2;
    }
    else if (k < spectral.k_COD[5] && k >= spectral.k_COD[6]) {
        kh = 5;
        kl = 6;
    }
    else if (k < spectral.k_COD[4] && k >= spectral.k_COD[5]) {
        kh = 4;
        kl = 5;
    }
    else if (k < spectral.k_COD[3] && k >= spectral.k_COD[4]) {
        kh = 3;
        kl = 4;
    }
    else if (k < spectral.k_COD[2] && k >= spectral.k_COD[3]) {
        kh = 2;
        kl = 3;
    }
    else {
        // Preserve original behavior: print error and continue.
        // In the legacy code this falls through with kh = kl = 0.
        std::fprintf(
            stderr,
            "error in wavelength correction, k (%5.4f) not within if loop\n",
            k
        );
    }

    // Interpolate between the two correction values.
    MAGIC_REAL distanceToUpper =
        std::fabs(spectral.k_COD[kh] - static_cast<MAGIC_REAL>(k));

    MAGIC_REAL intervalWidth =
        spectral.k_COD[kh] - spectral.k_COD[kl];

    if (kh == kl) {
        mkd = (MAGIC_REAL)0.0;
    }
    else {
        mkd = distanceToUpper / intervalWidth;
    }

    // Apply wavelength-dependent correction to each band.
    // NOTE:
    // This preserves the original indexing: correction row = kk + 1.
    for (int kk = KATO_MIN - 1; kk < KATO_MAX; kk++) {
        const MAGIC_REAL factor =
            mkd * spectral.cor[kk + 1][kl] +
            ((MAGIC_REAL)1.0 - mkd) * spectral.cor[kk + 1][kh];

        GHI[kk] = static_cast<MAGIC_REAL>(factor * GHI[kk]);
    }
}