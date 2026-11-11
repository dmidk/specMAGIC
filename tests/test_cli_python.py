"""
The CLI and the Python API must give the same results for the same image and inputs: both
run the same C++ (runMagic), one reading files, the other taking arrays.

    cmake -S . -B build && cmake --build build      # the CLI, build/magic
    uv sync                                         # the Python module
    uv run python -m unittest discover -s tests
"""
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import xarray as xr

import specmagic

REPO = Path(__file__).resolve().parents[1]
MAGIC = REPO / "build" / "magic"

# A small output grid over Europe, so that a run takes seconds
GRID = {27: 31, 28: 41, 29: 50.0, 30: 5.0, 31: 0.2}     # latdim, londim, latbegin, lonbegin, dxy
START = np.datetime64("2026-06-21T12:00")
RESOLUTION = 2000
PIXELS = 5568


def make_home(directory):
    """A home for magic: the repository's tables, a config with GRID, and output directories."""
    home = Path(directory)
    for name in ("climatologies", "luts", "satellites"):
        (home / name).symlink_to(REPO / name)
    for name in ("GHI", "DNI", "CAL", "CSR"):
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


@unittest.skipUnless(MAGIC.exists(), f"the CLI is not built ({MAGIC})")
class CliAndPythonAgree(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.home = make_home(cls.tmp.name)
        cls.image = make_image()
        image_path = specmagic.write_image(cls.image, cls.home)
        (cls.home / "image_list.txt").write_text(f"{image_path}\n")

        lat, lon = specmagic.output_grid(cls.home, "mtg")
        cls.inputs = make_inputs(lat, lon)
        cls.files = specmagic.write_inputs(cls.home / "inputs", **cls.inputs)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_cli(self, files):
        env = {**os.environ, "SATELLITE": "mtg", "ALBEDO": "LANDMAP", "STRICT_INPUTS": "0", "OMP_NUM_THREADS": "4",
               "AEROSOL_FILE": str(files.get("aerosol", "")), "OZONE_FILE": str(files.get("ozone", "")),
               "WATER_VAPOUR_FILE": str(files.get("water_vapour", "")), "ALBEDO_FILE": str(files.get("albedo", ""))}
        subprocess.run([str(MAGIC), f"{self.home}/", "vis_06", "0"], cwd=self.home, env=env, check=True,
                       stdout=subprocess.DEVNULL)

        stamp = str(START).replace("-", "").replace("T", "").replace(":", "")
        out = {}
        for name in ("GHI", "DNI", "CAL", "CSR"):
            with xr.open_dataset(self.home / "out" / name / f"{name}{stamp}.nc", mask_and_scale=False) as ds:
                out[name] = ds[name].values[0]
                attrs = dict(ds.attrs)
        return out, attrs

    def assert_same(self, inputs, files):
        cli, cli_attrs = self.run_cli(files)
        py = specmagic.run(self.image, home=self.home, **inputs)
        for name in ("GHI", "DNI", "CAL", "CSR"):
            np.testing.assert_array_equal(py[name].values, cli[name], err_msg=name)
        self.assertEqual(py.attrs, cli_attrs)
        return py

    def test_without_inputs(self):
        py = self.assert_same({}, {})
        self.assertGreater((py.GHI.values > 0).mean(), 0.5)

    def test_with_every_input(self):
        py = self.assert_same(self.inputs, self.files)
        self.assertIn("clipped to LUT", py.attrs["aod_usage"])
        self.assertEqual(py.attrs["h2o_source"], "test h2o")

    def test_inputs_change_the_result(self):
        without = specmagic.run(self.image, home=self.home)
        with_inputs = specmagic.run(self.image, home=self.home, **self.inputs)
        self.assertFalse(np.array_equal(without.GHI.values, with_inputs.GHI.values))

    def write_aod(self, name, values, dims=("lat", "lon"), encoding=None):
        """An aerosol file on the output grid, written as some other tool might."""
        aod = self.inputs["aerosol"]
        ds = xr.Dataset({"aod": (dims, values, {"provenance": "test aod"})}, coords={"lat": aod.lat, "lon": aod.lon})
        path = self.home / "inputs" / name
        ds.to_netcdf(path, encoding={"aod": encoding or {"dtype": "float32", "_FillValue": None}})
        return path

    def assert_refused(self, path, problem):
        """The CLI does not use the file, and falls back to the climatology for it."""
        cli, attrs = self.run_cli({"aerosol": path})
        without = specmagic.run(self.image, home=self.home)
        for name in ("GHI", "DNI", "CAL", "CSR"):
            np.testing.assert_array_equal(cli[name], without[name].values, err_msg=name)
        self.assertEqual(attrs["aod_source"], "aerosol climatology (input missing)")
        self.assertIn(problem, attrs["input_problems"])

    def test_cli_refuses_more_dimensions_than_lat_lon(self):
        aod = self.inputs["aerosol"].values
        self.assert_refused(self.write_aod("aod_time.nc", np.stack([aod] * 3), ("time", "lat", "lon")),
                            "must have exactly the dimensions (lat, lon)")

    def test_cli_refuses_lon_lat(self):
        aod = self.inputs["aerosol"].values
        self.assert_refused(self.write_aod("aod_lonlat.nc", aod.T, ("lon", "lat")),
                            "must have exactly the dimensions (lat, lon)")

    def test_cli_refuses_packed_values(self):
        self.assert_refused(self.write_aod("aod_packed.nc", self.inputs["aerosol"].values,
                                           encoding={"dtype": "int16", "scale_factor": 0.001, "_FillValue": -1}),
                            "is packed")

    def test_cli_fill_value_falls_back_like_nan(self):
        aod = self.inputs["aerosol"].values.copy()
        aod[5:15, 5:20] = np.nan
        path = self.write_aod("aod_fill.nc", aod, encoding={"dtype": "float32", "_FillValue": -999.0})
        with xr.open_dataset(path, mask_and_scale=False) as ds:
            self.assertTrue((ds.aod.values == -999).any())
        py = self.assert_same({"aerosol": self.inputs["aerosol"].copy(data=aod)}, {"aerosol": path})
        self.assertIn("% climatology", py.attrs["aod_usage"])
        self.assertNotIn("0.0% climatology", py.attrs["aod_usage"])

    def test_unusable_input_is_recorded_alike(self):
        # The CLI cannot read a file without aod, and Python cannot convert units it does not know
        cli, cli_attrs = self.run_cli({"aerosol": self.files["ozone"]})
        bad = self.inputs["aerosol"].isel(lat=slice(0, -1)).assign_attrs(units="mm")
        py = specmagic.run(self.image, home=self.home, aerosol=bad)
        for name in ("GHI", "DNI", "CAL", "CSR"):
            np.testing.assert_array_equal(py[name].values, cli[name], err_msg=name)
        self.assertEqual({k: v for k, v in py.attrs.items() if k != "input_problems"},
                         {k: v for k, v in cli_attrs.items() if k != "input_problems"})
        self.assertEqual(py.attrs["aod_source"], "aerosol climatology (input missing)")
        self.assertIn("units 'mm'", py.attrs["input_problems"])
        self.assertIn("aod is not in input", cli_attrs["input_problems"])

    def test_strict_python_raises_for_an_unusable_input(self):
        bad = self.inputs["aerosol"].isel(lat=slice(0, -1)).assign_attrs(units="mm")
        with self.assertRaisesRegex(RuntimeError, "units 'mm'"):
            specmagic.run(self.image, home=self.home, aerosol=bad, strict=True)

    def test_error_in_the_loop_is_a_python_exception(self):
        # Antarctic ozone in October lies outside the ozone LUT, which magic does not allow
        october = self.image.assign_attrs(start_time="2026-10-15 12:00:00")
        with self.assertRaisesRegex(RuntimeError, "outside"):
            specmagic.run(october, home=self.home, grid=(-68, -64, -20, 20, 0.5))


if __name__ == "__main__":
    unittest.main()
