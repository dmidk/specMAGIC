"""
The CLI and the Python API must give the same results for the same image and inputs: both
run the same C++ (runMagic), one reading files, the other taking arrays.

    cmake -S . -B build && cmake --build build      # the CLI, build/magic
    uv sync                                         # the Python module
    uv run pytest tests/test_cli_python.py
"""
import os
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

import specmagic

REPO = Path(__file__).resolve().parents[1]
MAGIC = REPO / "build" / "magic"
OUTPUTS = ("GHI", "DNI", "CAL", "CSR")

# A small output grid over Europe, so that a run takes seconds
GRID = {27: 31, 28: 41, 29: 50.0, 30: 5.0, 31: 0.2}     # latdim, londim, latbegin, lonbegin, dxy
START = np.datetime64("2026-06-21T12:00")
RESOLUTION = 2000
PIXELS = 5568

pytestmark = pytest.mark.skipif(not MAGIC.exists(), reason=f"the CLI is not built ({MAGIC})")


def make_home(home):
    """A home for magic: the repository's tables, a config with GRID, and output directories."""
    for name in ("climatologies", "luts", "satellites"):
        (home / name).symlink_to(REPO / name)
    for name in OUTPUTS:
        (home / "out" / name).mkdir(parents=True)

    lines, value_index = [], 0
    for line in (REPO / "magic-config.asc").read_text().splitlines(keepends=True):
        value, _, comment = line.partition("#")
        if value.strip():
            if value_index in GRID:
                line = f"{GRID[value_index]}\t\t# {comment.strip()}\n"
            value_index += 1
        lines.append(line)
    (home / "magic-config.asc").write_text("".join(lines))
    return home


def make_image():
    """A full disc of counts with smooth brightness and a band of clouds."""
    y, x = np.mgrid[0:PIXELS, 0:PIXELS]
    counts = 400 + 300 * np.sin(x / 300.0) * np.cos(y / 500.0)
    counts[(y // 150) % 7 == 0] = 2800
    return xr.DataArray(counts.astype(np.int16), dims=("y", "x"), name="vis_06",
                        attrs={"resolution": RESOLUTION, "calibration": "counts",
                               "start_time": str(START).replace("T", " ") + ":00"})


def make_inputs(lat, lon):
    """Input fields on the output grid, with gaps and values outside the LUT axes."""
    rng = np.random.default_rng(1)
    shape = (len(lat), len(lon))

    def field(name, values, provenance):
        return xr.DataArray(values.astype(np.float32), dims=("lat", "lon"), coords={"lat": lat, "lon": lon},
                            name=name, attrs={"provenance": provenance})

    h2o = rng.uniform(5, 40, shape)
    h2o[:10, :10] = np.nan                        # outside a regional model's domain
    vis, nir = rng.uniform(0.02, 0.3, shape), rng.uniform(0.1, 0.5, shape)
    land = (rng.uniform(size=shape) > 0.3).astype(float)
    vis[land == 0] = np.nan
    nir[land == 0] = np.nan
    return {"aerosol": field("aod", rng.uniform(0.02, 3.0, shape), "test aod"),
            "ozone": field("o3", rng.uniform(190, 560, shape), "test o3"),
            "water_vapour": field("h2o", h2o, "test h2o"),
            "albedo_vis": field("albedo_vis", vis, "test albedo_vis"),
            "albedo_nir": field("albedo_nir", nir, "test albedo_nir"),
            "land": field("albedo_land", land, "test land")}


@pytest.fixture(scope="module")
def case(tmp_path_factory):
    """The home with the image for the CLI, the image and inputs for Python, and the input files."""
    home = make_home(tmp_path_factory.mktemp("home"))
    image = make_image()
    image_path = specmagic.write_image(image, home)
    (home / "image_list.txt").write_text(f"{image_path}\n")

    lat, lon = specmagic.output_grid(home, "mtg")
    inputs = make_inputs(lat, lon)
    files = specmagic.write_inputs(home / "inputs", **inputs)
    return SimpleNamespace(home=home, image=image, inputs=inputs, files=files)


@pytest.fixture(scope="module")
def without_inputs(case):
    """The Python run without any input fields, which several tests compare against."""
    return specmagic.run(case.image, home=case.home)


def run_cli(case, files):
    """The CLI's outputs and their global attributes, for the given input files."""
    env = {**os.environ, "SATELLITE": "mtg", "ALBEDO": "LANDMAP", "STRICT_INPUTS": "0", "OMP_NUM_THREADS": "4",
           "AEROSOL_FILE": str(files.get("aerosol", "")), "OZONE_FILE": str(files.get("ozone", "")),
           "WATER_VAPOUR_FILE": str(files.get("water_vapour", "")), "ALBEDO_FILE": str(files.get("albedo", ""))}
    subprocess.run([str(MAGIC), f"{case.home}/", "vis_06", "0"], cwd=case.home, env=env, check=True,
                   stdout=subprocess.DEVNULL)

    stamp = str(START).replace("-", "").replace("T", "").replace(":", "")
    out = {}
    for name in OUTPUTS:
        with xr.open_dataset(case.home / "out" / name / f"{name}{stamp}.nc", mask_and_scale=False) as ds:
            out[name] = ds[name].values[0]
            attrs = dict(ds.attrs)
    return out, attrs


def assert_outputs_equal(a, b):
    for name in OUTPUTS:
        np.testing.assert_array_equal(np.asarray(a[name]), np.asarray(b[name]), err_msg=name)


def assert_same(case, inputs, files):
    """The CLI with the files and Python with the inputs give the same outputs and attributes."""
    cli, cli_attrs = run_cli(case, files)
    py = specmagic.run(case.image, home=case.home, **inputs)
    assert_outputs_equal(py, cli)
    assert py.attrs == cli_attrs
    return py


## ------ The CLI and Python agree ------

def test_without_inputs(case):
    py = assert_same(case, {}, {})
    assert (py.GHI.values > 0).mean() > 0.5


def test_with_every_input(case):
    py = assert_same(case, case.inputs, case.files)
    assert "clipped to LUT" in py.attrs["aod_usage"]
    assert py.attrs["h2o_source"] == "test h2o"


def test_inputs_change_the_result(case, without_inputs):
    with_inputs = specmagic.run(case.image, home=case.home, **case.inputs)
    assert not np.array_equal(without_inputs.GHI.values, with_inputs.GHI.values)


## ------ Input files the CLI must not use ------

def write_aod(case, name, values, dims=("lat", "lon"), encoding=None):
    """An aerosol file on the output grid, written as some other tool might."""
    aod = case.inputs["aerosol"]
    ds = xr.Dataset({"aod": (dims, values, {"provenance": "test aod"})}, coords={"lat": aod.lat, "lon": aod.lon})
    path = case.home / "inputs" / name
    ds.to_netcdf(path, encoding={"aod": encoding or {"dtype": "float32", "_FillValue": None}})
    return path


@pytest.mark.parametrize("layout, problem", [
    ("time", "must have exactly the dimensions (lat, lon)"),
    ("lon_lat", "must have exactly the dimensions (lat, lon)"),
    ("packed", "is packed"),
])
def test_cli_refuses_a_malformed_file(case, without_inputs, layout, problem):
    aod = case.inputs["aerosol"].values
    if layout == "time":
        path = write_aod(case, "aod_time.nc", np.stack([aod] * 3), ("time", "lat", "lon"))
    elif layout == "lon_lat":
        path = write_aod(case, "aod_lon_lat.nc", aod.T, ("lon", "lat"))
    else:
        path = write_aod(case, "aod_packed.nc", aod,
                         encoding={"dtype": "int16", "scale_factor": 0.001, "_FillValue": -1})

    # The CLI does not use the file, and falls back to the climatology for it
    cli, attrs = run_cli(case, {"aerosol": path})
    assert_outputs_equal(cli, without_inputs)
    assert attrs["aod_source"] == "aerosol climatology (input missing)"
    assert problem in attrs["input_problems"]


def test_cli_fill_value_falls_back_like_nan(case):
    aod = case.inputs["aerosol"].values.copy()
    aod[5:15, 5:20] = np.nan
    path = write_aod(case, "aod_fill.nc", aod, encoding={"dtype": "float32", "_FillValue": -999.0})
    with xr.open_dataset(path, mask_and_scale=False) as ds:
        assert (ds.aod.values == -999).any()

    py = assert_same(case, {"aerosol": case.inputs["aerosol"].copy(data=aod)}, {"aerosol": path})
    assert "% climatology" in py.attrs["aod_usage"]
    assert "0.0% climatology" not in py.attrs["aod_usage"]


## ------ Unusable inputs and errors ------

def unusable_aerosol(case):
    """An aerosol field that Python cannot use: off the output grid, in units it does not know."""
    return case.inputs["aerosol"].isel(lat=slice(0, -1)).assign_attrs(units="mm")


def test_unusable_input_is_recorded_alike(case):
    # The CLI cannot read a file without aod, and Python cannot convert units it does not know
    cli, cli_attrs = run_cli(case, {"aerosol": case.files["ozone"]})
    py = specmagic.run(case.image, home=case.home, aerosol=unusable_aerosol(case))

    assert_outputs_equal(py, cli)
    without_problems = lambda attrs: {k: v for k, v in attrs.items() if k != "input_problems"}
    assert without_problems(py.attrs) == without_problems(cli_attrs)
    assert py.attrs["aod_source"] == "aerosol climatology (input missing)"
    assert "units 'mm'" in py.attrs["input_problems"]
    assert "aod is not in input" in cli_attrs["input_problems"]


def test_strict_python_raises_for_an_unusable_input(case):
    with pytest.raises(RuntimeError, match="units 'mm'"):
        specmagic.run(case.image, home=case.home, aerosol=unusable_aerosol(case), strict=True)


def test_error_in_the_loop_is_a_python_exception(case):
    # Antarctic ozone in October lies outside the ozone LUT, which magic does not allow
    october = case.image.assign_attrs(start_time="2026-10-15 12:00:00")
    with pytest.raises(RuntimeError, match="outside"):
        specmagic.run(october, home=case.home, grid=(-68, -64, -20, 20, 0.5))
