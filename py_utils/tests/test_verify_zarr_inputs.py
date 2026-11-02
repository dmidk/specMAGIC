"""
Tests for the archive layout checks of verify_zarr_inputs.py, on small synthetic datasets.

    uv run python -m unittest discover -s py_utils/tests
"""
import sys
import unittest
from pathlib import Path

import numpy as np
import xarray as xr
from pyproj import CRS

# The module under test lives in py_utils, one level up
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from verify_zarr_inputs import check_cams, check_dini, check_lsasaf

HOURS = np.arange(np.datetime64("2026-08-23T00:00"), np.datetime64("2026-08-23T03:00"),
                  np.timedelta64(1, "h")).astype("datetime64[ns]")


def cams():
    lat, lon = np.linspace(90, -90, 5), np.arange(0.0, 360.0, 90.0)
    dims = ("valid_time", "latitude", "longitude")
    data = np.zeros((len(HOURS), len(lat), len(lon)), dtype=np.float32)
    return xr.Dataset(
        {"aod550": (dims, data, {"units": "~"}),
         "gtco3": (dims, data, {"units": "kg m**-2"}),
         "tcwv": (dims, data, {"units": "kg m**-2"})},
        coords={"valid_time": HOURS, "latitude": lat, "longitude": lon})


def dini():
    y, x = np.arange(3) * 2500.0, np.arange(4) * 2500.0
    latlon = np.zeros((len(y), len(x)))
    data = np.zeros((len(HOURS), 2, len(y), len(x)))
    return xr.Dataset(
        {"wvint": (("forecast_reference_time", "step", "y", "x"), data, {"units": "kg m**-2"}),
         "latitude": (("y", "x"), latlon), "longitude": (("y", "x"), latlon),
         "dini_projection": ((), 0.0, {"crs_wkt": CRS.from_epsg(3035).to_wkt()})},
        coords={"forecast_reference_time": HOURS, "y": y, "x": x,
                "step": np.array([0, 1], dtype="timedelta64[h]").astype("timedelta64[ns]")})


def lsasaf():
    days = np.array(["2026-08-22", "2026-08-23"], dtype="datetime64[ns]")
    angle = np.linspace(-0.15, 0.15, 4)
    dims = ("time", "y", "x")
    shape = (len(days), len(angle), len(angle))
    geostationary = {"grid_mapping_name": "geostationary", "perspective_point_height": 35786400.0,
                     "semi_major_axis": 6378137.0, "semi_minor_axis": 6356752.31424518,
                     "longitude_of_projection_origin": 0.0, "sweep_angle_axis": "y"}
    return xr.Dataset(
        {"AL-VI-DH": (dims, np.zeros(shape)), "AL-NI-DH": (dims, np.zeros(shape)),
         "QFLAGS": (dims, np.zeros(shape, dtype=np.int8)), "Z_Age": (dims, np.zeros(shape, dtype=np.int8)),
         "geostationary": ((), b"", geostationary)},
        coords={"time": days, "y": angle[::-1], "x": angle})


class Cams(unittest.TestCase):

    def test_matching_archive_passes(self):
        self.assertEqual(check_cams(cams()), [])

    def test_a_missing_variable_only_drops_that_variable(self):
        errors = check_cams(cams().drop_vars("gtco3"))
        self.assertEqual([var for var, _ in errors], ["gtco3"])

    def test_wrong_units_only_drop_that_variable(self):
        ds = cams()
        ds.tcwv.attrs["units"] = "mm"
        self.assertEqual([var for var, _ in check_cams(ds)], ["tcwv"])

    def test_longitudes_from_minus_180_are_refused(self):
        ds = cams().assign_coords(longitude=np.arange(-180.0, 180.0, 90.0))
        self.assertIn(None, [var for var, _ in check_cams(ds)])

    def test_a_regional_grid_is_refused(self):
        ds = cams().isel(longitude=slice(0, 3))
        self.assertIn(None, [var for var, _ in check_cams(ds)])


class Dini(unittest.TestCase):

    def test_matching_archive_passes(self):
        self.assertEqual(check_dini(dini()), [])

    def test_steps_in_hours_pass(self):
        self.assertEqual(check_dini(dini().assign_coords(step=[0, 1])), [])

    def test_descending_y_is_refused(self):
        ds = dini().isel(y=slice(None, None, -1))
        self.assertIn(None, [var for var, _ in check_dini(ds)])

    def test_missing_projection_is_refused(self):
        self.assertIn(None, [var for var, _ in check_dini(dini().drop_vars("dini_projection"))])


class LsaSaf(unittest.TestCase):

    def test_matching_archive_passes(self):
        self.assertEqual(check_lsasaf(lsasaf()), [])

    def test_undecoded_albedo_is_refused(self):
        ds = lsasaf()
        ds["AL-VI-DH"] = ds["AL-VI-DH"].astype(np.int16)
        self.assertNotEqual(check_lsasaf(ds), [])

    def test_days_must_start_at_midnight(self):
        ds = lsasaf()
        ds = ds.assign_coords(time=ds.time + np.timedelta64(12, "h"))
        self.assertNotEqual(check_lsasaf(ds), [])

    def test_axes_in_metres_are_refused(self):
        ds = lsasaf()
        ds = ds.assign_coords(x=ds.x * 35786400.0)
        self.assertNotEqual(check_lsasaf(ds), [])

    def test_every_problem_drops_the_whole_archive(self):
        errors = check_lsasaf(lsasaf().drop_vars("Z_Age"))
        self.assertEqual([var for var, _ in errors], [None])


if __name__ == "__main__":
    unittest.main()
