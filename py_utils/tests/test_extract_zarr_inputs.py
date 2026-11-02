"""
Tests for the pure helpers of extract_zarr_inputs.py.

    uv run python -m unittest discover -s py_utils/tests
"""
import sys
import unittest
from pathlib import Path

import numpy as np

# The module under test lives in py_utils, one level up
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from extract_zarr_inputs import axis_weights, bilinear_points, bilinear_regular, preferred_order

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


if __name__ == "__main__":
    unittest.main()
