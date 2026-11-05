"""
Tests for the helpers of extract_zarr_inputs.py.

    uv run python -m unittest discover -s py_utils/tests
"""
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

# The module under test lives in py_utils, one level up
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from extract_zarr_inputs import (axis_weights, bilinear_points, bilinear_regular, forecast,
                                 lsasaf_file, open_lsasaf_files, preferred_order)
from test_verify_zarr_inputs import lsasaf
from verify_zarr_inputs import check_lsasaf

HOUR = np.timedelta64(1, "h")


class PreferredOrder(unittest.TestCase):

    times = np.arange(np.datetime64("2026-08-23T00:00"), np.datetime64("2026-08-23T06:00"), HOUR)

    def order(self, target, max_before=3 * HOUR, max_after=np.timedelta64(30, "m")):
        idx = preferred_order(self.times, np.datetime64(target, "s"), max_before, max_after)
        return [str(self.times[i])[11:16] for i in idx]

    def test_later_hour_first_when_nearer(self):
        self.assertEqual(self.order("2026-08-23T00:50"), ["01:00", "00:00"])

    def test_earlier_wins_a_tie(self):
        self.assertEqual(self.order("2026-08-23T01:30")[:2], ["01:00", "02:00"])

    def test_window(self):
        self.assertEqual(self.order("2026-08-23T09:00"), [])
        self.assertEqual(self.order("2026-08-23T04:10"), ["04:00", "03:00", "02:00"])

    def test_daily_never_later(self):
        days = np.arange(np.datetime64("2026-08-20"), np.datetime64("2026-08-25")).astype("datetime64[s]")
        idx = preferred_order(days, np.datetime64("2026-08-23T00:00", "s"),
                              np.timedelta64(10, "D"), np.timedelta64(0, "s"))
        self.assertEqual([str(days[i])[:10] for i in idx],
                         ["2026-08-23", "2026-08-22", "2026-08-21", "2026-08-20"])


class Forecast(unittest.TestCase):

    # Not "run", which would hide TestCase.run
    reference_time = np.datetime64("2026-08-23T12:00", "ns")

    def test_timedelta_step(self):
        self.assertEqual(forecast(self.reference_time, np.timedelta64(3600 * 10**9, "ns")),
                         "run 2026-08-23T12:00 +1 h")

    def test_step_in_hours(self):
        self.assertEqual(forecast(self.reference_time, 11), "run 2026-08-23T12:00 +11 h")


class Interpolation(unittest.TestCase):

    def test_axis_weights(self):
        i, w, inside = axis_weights(np.array([0.0, 1.0, 2.0]), np.array([0.25, 2.0, 3.0]))
        np.testing.assert_array_equal(i, [0, 1, 1])
        np.testing.assert_allclose(w, [0.25, 1.0, 2.0])
        np.testing.assert_array_equal(inside, [True, True, False])

    def test_descending_axis_is_refused(self):
        with self.assertRaises(ValueError):
            axis_weights(np.array([2.0, 1.0, 0.0]), np.array([0.5]))

    def test_bilinear_points_reproduces_a_plane(self):
        y, x = np.arange(4.0), np.arange(5.0)
        values = 2 * y[:, None] + 3 * x[None, :]
        out = bilinear_points(y, x, values, np.array([0.5, 2.25, 9.0]), np.array([1.5, 3.75, 1.0]))
        np.testing.assert_allclose(out[:2], [2 * 0.5 + 3 * 1.5, 2 * 2.25 + 3 * 3.75])
        self.assertTrue(np.isnan(out[2]))

    def test_bilinear_regular_wraps_and_flips(self):
        # Descending latitude and 0..360 longitude, as in CAMS
        src_lat = np.array([10.0, 0.0, -10.0])
        src_lon = np.arange(0.0, 360.0, 10.0)
        values = np.tile(src_lon, (3, 1)) + src_lat[:, None]
        out = bilinear_regular(src_lat, src_lon, values, np.array([5.0]), np.array([-5.0, 15.0]))

        # -5 lies between 350 and the wrapped 360, whose value is that of 0
        np.testing.assert_allclose(out, [[(350 + 0) / 2 + 5, 15 + 5]], rtol=1e-6)


class LsaSafDailyFiles(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

        # One file per day, as LSA SAF distributes them, from 2026-08-01 to 2026-08-25
        day = lsasaf().isel(time=[0])
        for d in np.arange(np.datetime64("2026-08-01"), np.datetime64("2026-08-26")):
            day.assign_coords(time=[d.astype("datetime64[ns]")]).to_netcdf(lsasaf_file(self.dir, d))

    def tearDown(self):
        self.tmp.cleanup()

    def test_file_name(self):
        self.assertEqual(os.path.basename(lsasaf_file("/d", np.datetime64("2025-01-01T13:40"))),
                         "LSA-105_MTG_MTDAL_MTG-FD_202501010000.nc")

    def test_opens_the_image_day_and_the_ten_days_before(self):
        ds = open_lsasaf_files(self.dir, np.datetime64("2026-08-20T13:40", "s"))
        self.assertEqual([str(t)[:10] for t in ds.time.values],
                         [str(d) for d in np.arange(np.datetime64("2026-08-10"), np.datetime64("2026-08-21"))])
        self.assertEqual(check_lsasaf(ds), [])

    def test_takes_what_there_is_before_a_gap(self):
        ds = open_lsasaf_files(self.dir, np.datetime64("2026-08-30T12:00", "s"))
        self.assertEqual(str(ds.time.values[-1])[:10], "2026-08-25")

    def test_no_file_in_the_window(self):
        with self.assertRaises(FileNotFoundError):
            open_lsasaf_files(self.dir, np.datetime64("2026-09-30T12:00", "s"))


if __name__ == "__main__":
    unittest.main()
