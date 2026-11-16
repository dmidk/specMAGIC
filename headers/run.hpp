#pragma once
#include <string>
#include <utility>
#include <vector>
#include "helpers.hpp"
#include "image.hpp"
#include "input_fields.hpp"
#include "navigation.hpp"
#include "radiation.hpp"
#include "types.hpp"

// One run of magic on one image: what the CLI (main.cpp) and the Python API
// (python/bindings.cpp) share, so that both give the same results.

struct MagicRun {
    Output radiation;

    // Global attributes for the output: where the input fields came from
    std::vector<std::pair<std::string, std::string>> attributes;

    // Time spent in the loop over the pixels
    double loop_seconds = 0;

    MagicRun(int nlat, int nlon) : radiation(nlat, nlon) {}

    // Holds the output matrices, so it is moved, never copied
    MagicRun(const MagicRun&) = delete;
    MagicRun& operator=(const MagicRun&) = delete;
    MagicRun(MagicRun&&) = default;
    MagicRun& operator=(MagicRun&&) = default;
};

// The output grid of the config, shifted with the satellite's longitude
Geography outputGrid(Config& c, MAGIC_EXACT subsatellite_longitude_deg);

// GHI, DNI, CAL and CSR of every pixel of the output grid. The input fields must already
// have been loaded or set for that grid (InputFields::Inputs::load or set).
MagicRun runMagic(Config& c, Image& img, Geography& geo, AlbedoType albedo_source,
    InputFields::Inputs& inputs);
