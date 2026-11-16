// The Python module specmagic._magic: one run of magic with the image and the input fields
// given from memory, through the same runMagic as the CLI. The Python package (specmagic)
// prepares the arrays; see python/specmagic/__init__.py.

#include <cstdio>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "../headers/read.hpp"
#include "../headers/run.hpp"

namespace py = pybind11;

using Counts = py::array_t<MAGIC_INT, py::array::c_style | py::array::forcecast>;
using Values = py::array_t<float, py::array::c_style | py::array::forcecast>;
using Grid = std::tuple<int, int, float, float, float>;

static py::array_t<MAGIC_INT> toNumpy(const Matrix& m) {
    py::array_t<MAGIC_INT> out({m.rows(), m.cols()});
    std::copy(m.data(), m.data() + static_cast<long>(m.rows()) * m.cols(), out.mutable_data());
    return out;
}

static py::array_t<MAGIC_EXACT> toNumpy(const MAGIC_EXACT* values, unsigned int n) {
    py::array_t<MAGIC_EXACT> out(n);
    std::copy(values, values + n, out.mutable_data());
    return out;
}

// The config, with the output grid replaced if one is given
static Config config(const std::string& home, const std::optional<Grid>& grid) {
    Config c = loadConfig("magic-config.asc", home);
    if (grid) {
        std::tie(c.latdim, c.londim, c.latbegin, c.lonbegin, c.dxy) = *grid;
    }
    return c;
}

// The lat/lon axes of the output grid, onto which Python regrids the input fields
static py::tuple outputGridAxes(const std::string& home, const std::string& satellite, std::optional<Grid> grid) {
    Config c = config(home, grid);
    Geography geo = outputGrid(c, loadMetadata("", satellite, home).subsatellite_longitude_deg);
    return py::make_tuple(toNumpy(geo.lat, geo.nlat), toNumpy(geo.lon, geo.nlon));
}

static py::dict run(const std::string& home, const std::string& channel, const std::string& satellite,
    Counts counts, int resolution, std::tuple<int, int, int, int, int> start,
    const std::map<std::string, std::pair<Values, std::string>>& fields,
    const std::map<std::string, std::string>& problems,
    const std::string& albedo, bool strict, std::optional<Grid> grid) {

    if (counts.ndim() != 2) throw std::invalid_argument("the image must be 2-D (y, x)");

    Config c = config(home, grid);

    Image img;
    auto [year, month, day, hour, minute] = start;
    img.fromArray(channel, satellite, home, counts.data(), static_cast<int>(counts.shape(0)),
        static_cast<int>(counts.shape(1)), resolution, DateTime(year, month, day, hour, minute));

    if (albedo != "MODIS" && albedo != "LANDMAP") throw std::invalid_argument("unknown albedo " + albedo);
    AlbedoType albedo_source = albedo == "MODIS" ? AlbedoType::MODIS : AlbedoType::LANDMAP;

    Geography geo = outputGrid(c, img.info.subsatellite_longitude_deg);

    // The given fields, and those that Python could not put on the output grid, recorded
    // as the CLI records its files
    InputFields::Inputs inputs;
    inputs.request.strict = strict;
    for (const auto& [name, message] : problems) {
        inputs.refuse(name, message);
    }
    for (const auto& [name, field] : fields) {
        const Values& values = field.first;
        if (values.ndim() != 2) throw std::invalid_argument("input " + name + " must be 2-D (lat, lon)");
        inputs.set(name, values.data(), static_cast<int>(values.shape(0)), static_cast<int>(values.shape(1)),
            field.second, geo);
    }
    std::fflush(stdout);
    inputs.check();

    if (inputs.request.albedo) {
        if (albedo_source == AlbedoType::MODIS) {
            throw std::invalid_argument("albedo MODIS and an albedo input cannot be used together");
        }
        albedo_source = AlbedoType::INPUT;
    }

    MagicRun result = [&] {
        py::gil_scoped_release release;
        return runMagic(c, img, geo, albedo_source, inputs);
    }();

    // magic reports with printf; flush it now, so that it appears with the run in Python
    // (e.g. in a notebook or under pytest) rather than when the process exits
    std::fflush(stdout);

    py::dict attributes;
    for (const auto& [key, value] : result.attributes) attributes[py::str(key)] = value;

    py::dict out;
    out["GHI"] = toNumpy(result.radiation.GHI);
    out["DNI"] = toNumpy(result.radiation.DNI);
    out["CAL"] = toNumpy(result.radiation.CAL);
    out["CSR"] = toNumpy(result.radiation.CSR);
    out["lat"] = toNumpy(geo.lat, geo.nlat);
    out["lon"] = toNumpy(geo.lon, geo.nlon);
    out["attributes"] = attributes;
    out["loop_seconds"] = result.loop_seconds;
    return out;
}

PYBIND11_MODULE(_magic, m) {
    m.doc() = "specMAGIC's C++ core; use specmagic.run instead of calling this directly";
    m.def("run", &run, py::arg("home"), py::arg("channel"), py::arg("satellite"), py::arg("counts"),
        py::arg("resolution"), py::arg("start"), py::arg("fields"), py::arg("problems"), py::arg("albedo"),
        py::arg("strict"), py::arg("grid"));
    m.def("output_grid", &outputGridAxes, py::arg("home"), py::arg("satellite"), py::arg("grid"));
}
