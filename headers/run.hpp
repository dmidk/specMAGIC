#pragma once
#include <string>
#include <utility>
#include <vector>
#include "helpers.hpp"
#include "image.hpp"
#include "input_fields.hpp"
#include "modis_brdf.hpp"
#include "navigation.hpp"
#include "radiation.hpp"
#include "tables.hpp"
#include "types.hpp"

// One run of magic on one image: what the CLI (main.cpp) and the Python API
// (python/bindings.cpp) share, so that both give the same results.

struct MagicRun {
    Output radiation;

    // Global attributes for the output: where the input fields came from
    std::vector<std::pair<std::string, std::string>> attributes;

    // How many day pixels used each input, behind the usage attributes
    InputFields::Usage usage;

    // Time spent in the loop over the pixels
    double loop_seconds = 0;

    MagicRun(int nlat, int nlon) : radiation(nlat, nlon) {}

    // Holds the output matrices, so it is moved, never copied
    MagicRun(const MagicRun&) = delete;
    MagicRun& operator=(const MagicRun&) = delete;
    MagicRun(MagicRun&&) = default;
    MagicRun& operator=(MagicRun&&) = default;
};

// The lookup tables and the climatologies of one month, read from disk before a run. A run
// only reads them (its threads share them), so the images of a series in that month can
// share them too. Built in place, as the tables cannot be copied or moved.
struct MagicTables {
    unsigned int month;
    bool modis;            // whether the MODIS albedo of the month is loaded

    Tables::GroundAlbedo alb;
    ModisBrdf::ModisBrdfAlbedo modis_albedo;
    Tables::Correction::Spectral lc;          // for the cloudy sky correction
    Tables::Correction::Absorber water;
    Tables::Correction::Absorber ozone;
    Tables::RTM monster;
    Tables::Climate climatologies;

    MagicTables(Config& c, MAGIC_EXACT satellite_longitude_deg, unsigned int month, bool modis);

    MagicTables(const MagicTables&) = delete;
    MagicTables& operator=(const MagicTables&) = delete;
};

// The output grid of the config, shifted with the satellite's longitude
Geography outputGrid(Config& c, MAGIC_EXACT subsatellite_longitude_deg);

// GHI, DNI, CAL and CSR of every pixel of the output grid. The input fields must already
// have been loaded or set for that grid (InputFields::Inputs::load or set).
// The tables must be those of the image's month, with MODIS loaded if albedo_source is MODIS.
MagicRun runMagic(Config& c, Image& img, Geography& geo, MagicTables& tables, AlbedoType albedo_source,
    InputFields::Inputs& inputs);
