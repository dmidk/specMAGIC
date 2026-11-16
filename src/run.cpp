#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <cassert>
#include <stdexcept>
#include <string>
#include "../headers/run.hpp"
#include "../headers/tables.hpp"
#include "../headers/magic.hpp"
#include "../headers/satellite.hpp"
#include "../headers/sun_geometry.hpp"
#include "../headers/reflectivity.hpp"
#include "../headers/modis_brdf.hpp"


using namespace Tables;

Geography outputGrid(Config& c, MAGIC_EXACT subsatellite_longitude_deg) {

    // Global coords across ENTIRE image
    Geography geo = Geography(c.latdim, c.londim);
    geo.makeCoordinates(c, subsatellite_longitude_deg);
    return geo;
}

MagicRun runMagic(Config& c, Image& img, Geography& geo, AlbedoType albedo_source,
    InputFields::Inputs& inputs) {

    const MAGIC_EXACT satellite_longitude_deg = img.info.subsatellite_longitude_deg;

    // ------ Read lookup tables ---------
    GroundAlbedo alb(c.num_bands);
    alb.fromFile(c);

    ModisBrdf::ModisBrdfAlbedo modis;
    // Empty if MODIS is not enabled
    if (albedo_source == AlbedoType::MODIS) {
        modis.load(c, satellite_longitude_deg, img.timestamp.month, albedo_source);
    }

    // For the cloudy sky correction
    Correction::Spectral lc(KATO_MAX + 1);
    lc.fromFile(c);

    Correction::Absorber water(c.h2odim, c.num_bands);
    water.fromFile(c.path+c.clut_spec_h2o);

    Correction::Absorber ozone(c.o3dim, c.num_bands);
    ozone.fromFile(c.path+c.clut_spec_o3);

    RTM monster(c.ggdim, c.ssadim, c.aoddim, c.num_bands);
    monster.fromFile(c.path+c.clut_spec);
    // ------ Done reading climatologies and lookup tables ---------


    // Read climatologies
    Climate climatologies = Climate(c, img.timestamp.month);

    // Optional input fields, which override the climatologies per pixel
    const bool use_inputs = inputs.request.any();
    if (use_inputs) {
        inputs.setLimits(monster, water, ozone);
    }

    InputFields::Usage input_usage;

    // Allocate the radiation matrices
    MagicRun run(geo.nlat, geo.nlon);
    Output& radiation = run.radiation;

    // Counters
    uint64_t day = 0;
    uint64_t night = 0;
    uint64_t invalid = 0;

    // Number of loop iterations
    int pixels = geo.nlat * geo.nlon;

    int kats = KATO_MAX + 1;

    // The first error in the loop over the pixels, if any
    int failed = 0;
    std::string failure;

    // Begin central loop
    // Geography is not written to here, so should be shared
    // Radiation matrix written to PER ELEMENT, so should be shared
    // im info read-only, should be shared
    double t1 = omp_get_wtime();
    #pragma omp parallel
    {

        // Make the scratch memory
        // Where we do the calculations
        // This MUST be thread local
        Scratch sc(kats);

        // Climate variables for this pixel ONLY
        PixelClimate clim = PixelClimate();

        // How this thread's pixels were sourced from the input fields
        InputFields::Usage local_usage;

        #pragma omp for schedule(static, 4096) reduction(+:invalid,night,day)
        for (int pix = 0; pix < pixels; pix++) {

            // An error cannot leave a parallel region (it would terminate the program, and
            // with it Python), so the first one is kept, the remaining pixels are skipped,
            // and it is thrown again after the loop
            int stop;
            #pragma omp atomic read
            stop = failed;
            if (stop) continue;

            try {

                // Coordinates for this pixel only
                Area a;
                a.makeArea(geo, pix, satellite_longitude_deg);

                // Get the pixel position in line and columns
                unsigned int col, lin;

                MAGIC_EXACT lon_for_projection = a.lon - (satellite_longitude_deg / RADTODEGREE);
                Satellite::geo2Image(a.lat, lon_for_projection, img.info, col, lin);

                int line = lin; int column = col;
                if (img.info.flip_vertical) line = Satellite::flipVertical(lin, img.info.num_lines);
                if (img.info.flip_horizontal) column = Satellite::flipHorizontal(col, img.info.num_columns);

                // Get the time (as experienced by satellite)
                MAGIC_EXACT UTC_time = Satellite::calcObsTime(line, img.info, img.timestamp);

                // Calculate all the other sun geometry things
                SolarParameters sun = SunGeometry::solarParameters(img.timestamp, UTC_time, a);

                // Distance to subsatellite point
                MAGIC_EXACT nadir = sqrt(pow(a.lat, 2) + pow((a.lon - a.deltalon_rad), 2));

                // Check if the pixel is valid
                PixelState s = SunGeometry::calcPixelValidity(nadir, sun.cos_sza);

                switch (s) {

                    case PixelState::Invalid:
                    {
                        invalid++;
                        radiation.GHI(a.nlat, a.nlon) = -1;

                        // Many possible reasons for pixel being invalid
                        // Sun close to horizon, or viewing angle too far away
                        // Expected behaviour ~5% of pixels are invalid

                        #pragma omp simd
                        for (int band = KATO_MIN; band < KATO_MAX; band++) {
                            sc.GHI_spectral[band] = -1;
                            sc.DNI_spectral[band] = -1;
                        }
                    }
                    break;

                    case PixelState::Night:
                    {
                        night++;
                        radiation.GHI(a.nlat, a.nlon) = 0;

                        // We see nothing if nighttime, so skip


                        #pragma omp simd
                        for (int band = KATO_MIN; band < KATO_MAX; band++) {
                            sc.GHI_spectral[band] = 0;
                            sc.DNI_spectral[band] = 0;
                        }
                    }
                    break;

                    case PixelState::Day:
                    {
                        day++;

                        // Put this pixel's climatology information into
                        // the PixelClimate object
                        climatologies.makeLocalClimate(clim, a);
                        if (use_inputs) {
                            inputs.apply(clim, a, local_usage);
                        }

                        // Add the albedo correction
                        Reflectivity::makeAlbedoCorrection(climatologies, clim, a, sun.cos_sza, sun.declination,
                            albedo_source == AlbedoType::INPUT);

                        // Correction for an FCI measuring fault
                        SunGeometry::correctZenithAngle(sun);

                        // Calculate the CAL for this pixel
                        MAGIC_REAL cal = effectiveCloudAlbedo(img, sun, climatologies, albedo_source, alb,
                            modis, inputs, a, clim, line, column);

                        // Put the CAL into the radiation matrix
                        radiation.encodeCAL(cal, a);

                        // Get the clear sky index
                        MAGIC_REAL k = Reflectivity::clearSkyIndex(cal);

                        MAGIC_REAL CSR = 0;
                        for (int band = KATO_MIN; band < KATO_MAX; band++) {

                            // Calculate surface albedo at this point
                            // Prefer MODIS BRDF albedo if available, otherwise use old land-use albedo.
                            clim.surface_albedo = static_cast<MAGIC_EXACT>(Reflectivity::getBestKatoSurfaceAlbedo(
                                climatologies, a, albedo_source, alb, modis, inputs, band, sun.cos_sza, clim.albedo_correction,
                                clim.noon_albedo_correction));

                            // Calculate the clear-sky ghi at this pixel
                            sc.GHI_spectral[band] = magic(monster, water, ozone, IrradianceMode::Global, band,
                                sun.cos_sza, sun.distance, clim);

                            // Calculate the clear-sky DNI at this pixel
                            sc.DNI_spectral[band] = magic(monster, water, ozone, IrradianceMode::Beam, band,
                                sun.cos_sza, sun.distance, clim);

                            // Get the clearsky value BEFORE cloudysky corrections
                            CSR += sc.GHI_spectral[band];

                            // Eq. 13 of Mueller et al 2012
                            sc.GHI_spectral[band] *= k;

                            // Empirical correction for the beam component
                            MAGIC_REAL beam_check = k - (MAGIC_REAL)0.38 * ((MAGIC_REAL)1. - k);
                            if (beam_check < 0) sc.DNI_spectral[band] = 0.;
                            else sc.DNI_spectral[band] *= pow(beam_check, (MAGIC_REAL)2.5);


                        }

                        // Wavelength correction for the cloudy sky
                        wavelengthCorrectionCloudySky(lc, k, sc.GHI_spectral);

                        // Per-pixel sum
                        MAGIC_REAL GHI = 0;
                        MAGIC_REAL DNI = 0;

                        // Sum up the band contributions
                        // Gets the total GHI/DNI for this pixel
                        for (int band = KATO_MIN; band < KATO_MAX; band++) {

                            // To avoid indexing errors
                            int idx = band - 1;

                            GHI += sc.GHI_spectral[idx];
                            DNI += sc.DNI_spectral[idx];

                        }

                        // Put these values in the matrices
                        radiation.encodeGHI(GHI, a);
                        radiation.encodeCSR(CSR, a);
                        radiation.encodeDNI(DNI, sun.corrected_cos_sza, a);

                    }
                    break;

                    default:
                        // This is undefined behaviour
                        assert(!"Got an undefined pixel state in the satellite image!");

                }   // End PixelState switch

            } catch (const std::exception& e) {
                #pragma omp critical(magic_failure)
                {
                    if (!failed) failure = e.what();
                    #pragma omp atomic write
                    failed = 1;
                }
            }

        }   // End for loop over pixels

        #pragma omp critical
        input_usage += local_usage;

    }    // End parallel region
    run.loop_seconds = omp_get_wtime() - t1;

    if (failed) throw std::runtime_error(failure);

    // For safety
    printPixels(day, night, invalid);

    // Record which input fields were used
    if (use_inputs) {
        inputs.report(input_usage);
        run.attributes = inputs.describe(input_usage);
    }

    return run;
}
