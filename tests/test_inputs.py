"""
Tests for picking the time slice of an input and regridding it, and for preparing the image.

    uv run pytest tests/test_inputs.py
"""
from datetime import datetime

import numpy as np
import pytest
import xarray as xr
from pyorbital.astronomy import sun_earth_distance_correction
from pyproj import CRS, Transformer

from specmagic import lsasaf_albedo, write_inputs
from specmagic._image import FULL_DISC_PIXELS, FILL, image_counts, to_counts
from specmagic._inputs import candidates, forecast, positions, preferred_order, regrid, regrid_input, sample

HOUR = np.timedelta64(1, "h")
LAT, LON = np.arange(50.0, 56.0, 0.5), np.arange(5.0, 15.0, 0.5)


## ------ Picking the time ------

def order(target, max_before=3 * HOUR, max_after=np.timedelta64(30, "m")):
    """The hours of 2026-08-23 00:00..05:00 that suit target, best first."""
    times = np.arange(np.datetime64("2026-08-23T00:00"), np.datetime64("2026-08-23T06:00"), HOUR)
    idx = preferred_order(times, np.datetime64(target, "s"), max_before, max_after)
    return [str(times[i])[11:16] for i in idx]


def test_later_hour_first_when_nearer():
    assert order("2026-08-23T00:50") == ["01:00", "00:00"]


def test_earlier_wins_a_tie():
    assert order("2026-08-23T01:30")[:2] == ["01:00", "02:00"]


def test_window():
    assert order("2026-08-23T09:00") == []
    assert order("2026-08-23T04:10") == ["04:00", "03:00", "02:00"]


def test_daily_never_later():
    days = np.arange(np.datetime64("2026-08-20"), np.datetime64("2026-08-25")).astype("datetime64[ns]")
    da = xr.DataArray(np.zeros((5, 2, 2)), dims=("time", "y", "x"), coords={"time": days})
    labels = [label for _, label in candidates(da, np.datetime64("2026-08-23T13:40", "s"), daily=True)]
    assert labels == [f"day 2026-08-{d}" for d in (23, 22, 21, 20)]


def test_newest_run_for_the_nearest_valid_time():
    runs = np.array(["2026-08-23T10:00", "2026-08-23T11:00"], dtype="datetime64[ns]")
    steps = np.arange(4).astype("timedelta64[h]").astype("timedelta64[ns]")
    da = xr.DataArray(np.zeros((2, 4, 2, 2)), dims=("forecast_reference_time", "step", "y", "x"),
                      coords={"forecast_reference_time": runs, "step": steps})
    _, label = next(candidates(da, np.datetime64("2026-08-23T12:10", "s"), daily=False))
    assert label == "run 2026-08-23T11:00 +1 h, valid 2026-08-23T12:00"


def test_forecast_along_the_time_axis():
    times = np.array(["2026-08-23T12:00"], dtype="datetime64[ns]")
    da = xr.DataArray(np.zeros((1, 2, 2)), dims=("valid_time", "y", "x"),
                      coords={"valid_time": times, "forecast_reference_time": ("valid_time", times - HOUR),
                              "forecast_period": ("valid_time", [HOUR])})
    _, label = next(candidates(da, np.datetime64("2026-08-23T12:00", "s"), daily=False))
    assert label == "run 2026-08-23T11:00 +1 h, valid 2026-08-23T12:00"


def test_step_in_hours():
    assert forecast(np.datetime64("2026-08-23T12:00"), 11) == "run 2026-08-23T12:00 +11 h"


## ------ Regridding ------

def test_positions_either_direction():
    np.testing.assert_allclose(positions([0.0, 1.0, 2.0], [0.5, 3.0]), [0.5, np.nan])
    np.testing.assert_allclose(positions([2.0, 1.0, 0.0], [0.5]), [1.5])


def test_bilinear_reproduces_a_plane():
    y, x = np.arange(4.0), np.arange(5.0)
    da = xr.DataArray(2 * y[:, None] + 3 * x[None, :], dims=("y", "x"))
    out = sample(da, np.array([0.5, 2.25, np.nan]), np.array([1.5, 3.75, 1.0]), "bilinear")
    np.testing.assert_allclose(out, [2 * 0.5 + 3 * 1.5, 2 * 2.25 + 3 * 3.75, np.nan], rtol=1e-6)


def test_global_grid_in_0_360_wraps():
    # Descending latitude and 0..360 longitude, as in CAMS
    src_lat, src_lon = np.array([60.0, 50.0, 40.0]), np.arange(0.0, 360.0, 10.0)
    da = xr.DataArray(np.tile(src_lon, (3, 1)) + src_lat[:, None], dims=("latitude", "longitude"),
                      coords={"latitude": src_lat, "longitude": src_lon})
    out = regrid(da, np.array([55.0]), np.array([-5.0, 15.0]), "bilinear")
    # -5 lies between 350 and the wrapped 360, whose value is that of 0
    np.testing.assert_allclose(out, [[(350 + 0) / 2 + 55, 15 + 55]], rtol=1e-6)


def test_projected_grid_reproduces_a_plane():
    crs = CRS.from_epsg(3035)
    x, y = np.arange(3.8e6, 4.8e6, 2.5e4), np.arange(2.8e6, 3.8e6, 2.5e4)
    da = xr.DataArray(2e-6 * y[:, None] + 3e-6 * x[None, :], dims=("y", "x"), name="h2o",
                      coords={"y": y, "x": x, "crs": ((), 0, {"crs_wkt": crs.to_wkt()})})
    out = regrid(da, LAT, LON, "bilinear")
    px, py = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(*np.meshgrid(LON, LAT))
    np.testing.assert_allclose(out, 2e-6 * py + 3e-6 * px, rtol=1e-5)


def test_projected_grid_needs_its_grid_mapping():
    da = xr.DataArray(np.zeros((2, 2)), dims=("y", "x"), coords={"y": [0.0, 1.0], "x": [0.0, 1.0]}, name="h2o")
    with pytest.raises(ValueError, match="grid mapping"):
        regrid(da, LAT, LON, "bilinear")


## ------ The input fields ------

def lsasaf_dataset():
    """A small LSA SAF MTDAL dataset on the geostationary grid, in scan angles as LSA SAF gives it."""
    angle = np.linspace(-0.155, 0.155, 400)
    days = np.array(["2026-08-22", "2026-08-23"], dtype="datetime64[ns]")
    dims = ("time", "y", "x")
    shape = (2, 400, 400)
    mapping = {"grid_mapping_name": "geostationary", "perspective_point_height": 35786400.0,
               "semi_major_axis": 6378137.0, "semi_minor_axis": 6356752.31424518,
               "longitude_of_projection_origin": 0.0, "sweep_angle_axis": "y"}
    return xr.Dataset(
        {"AL-VI-DH": (dims, np.full(shape, 0.2)), "AL-NI-DH": (dims, np.full(shape, 0.3)),
         "QFLAGS": (dims, np.ones(shape, dtype=np.int8)), "Z_Age": (dims, np.r_[[np.full((400, 400), 3)],
                                                                             [np.full((400, 400), 20)]]),
         "geostationary": ((), b"", mapping)},
        coords={"time": days, "y": ("y", angle[::-1], {"units": "radian"}), "x": ("x", angle, {"units": "radian"})})


def test_ozone_in_kg_m2_becomes_du():
    da = xr.DataArray(np.full((2, 2), 300 * 2.1415e-5), dims=("latitude", "longitude"), name="gtco3",
                      coords={"latitude": [40.0, 60.0], "longitude": [0.0, 20.0]}, attrs={"units": "kg m**-2"})
    field = regrid_input(da, "ozone", np.datetime64("2026-08-23T12:00", "s"), LAT, LON)
    np.testing.assert_allclose(field.values, 300, rtol=1e-5)
    assert field.name == "o3"


def test_unknown_units_are_refused():
    da = xr.DataArray(np.ones((2, 2)), dims=("latitude", "longitude"), name="tcwv",
                      coords={"latitude": [40.0, 60.0], "longitude": [0.0, 20.0]}, attrs={"units": "mm"})
    with pytest.raises(ValueError, match="units"):
        regrid_input(da, "water_vapour", np.datetime64("2026-08-23T12:00", "s"), LAT, LON)


def test_lsasaf_albedo_drops_old_pixels():
    albedo = lsasaf_albedo(lsasaf_dataset())
    when = np.datetime64("2026-08-23T12:00", "s")

    # The image's day is 20 days old everywhere, so nothing is usable then; the day before is
    vis = regrid_input(albedo["albedo_vis"], "albedo_vis", when, LAT, LON)
    np.testing.assert_allclose(vis.values, 0.2)
    assert "day 2026-08-22" in vis.attrs["provenance"]
    land = regrid_input(albedo["land"], "land", when, LAT, LON)
    np.testing.assert_array_equal(land.values, 1)


def test_albedo_file_gets_a_land_mask(tmp_path):
    vis = xr.DataArray(np.full((2, 3), 0.2, np.float32), dims=("lat", "lon"), name="albedo_vis",
                       coords={"lat": [1.0, 2.0], "lon": [1.0, 2.0, 3.0]}, attrs={"provenance": "v"})
    paths = write_inputs(tmp_path, albedo_vis=vis, albedo_nir=vis.rename("albedo_nir"))
    with xr.open_dataset(paths["albedo"]) as ds:
        assert sorted(ds.data_vars) == ["albedo_land", "albedo_nir", "albedo_vis"]
        np.testing.assert_array_equal(ds.albedo_land.values, 1)


## ------ The image ------

def crop_coordinates(rows, cols, pixels, step):
    """x/y of the given rows and columns of the full disc, whose pixel centres lie at (i - pixels / 2 + 0.5) * step."""
    return {"y": (np.asarray(rows) - pixels / 2 + 0.5) * step, "x": (np.asarray(cols) - pixels / 2 + 0.5) * step}


def test_crop_goes_to_its_place_on_the_disc():
    pixels = FULL_DISC_PIXELS[2000]
    crop = xr.DataArray(np.arange(6.0).reshape(2, 3), dims=("y", "x"), name="vis_06",
                        coords=crop_coordinates(range(4000, 4002), range(100, 103), pixels, 2000.0),
                        attrs={"resolution": 2000, "calibration": "counts", "start_time": "2026-08-23 12:00:00"})
    counts, resolution, when, channel = image_counts(crop, "mtg")
    assert (counts.shape, resolution, when, channel) == ((pixels, pixels), 2000, datetime(2026, 8, 23, 12), "vis_06")
    np.testing.assert_array_equal(counts[4000:4002, 100:103], [[0, 1, 2], [3, 4, 5]])
    assert counts[0, 0] == FILL


def test_channel_from_satpy_attributes():
    # satpy names the DataArray after the L1c variable, and the channel in its attributes
    image = xr.DataArray(np.zeros((FULL_DISC_PIXELS[2000],) * 2, np.int16), dims=("y", "x"),
                         name="effective_radiance",
                         attrs={"name": "vis_06", "resolution": 2000, "calibration": "counts",
                                "start_time": datetime(2026, 8, 23, 12)})
    assert image_counts(image, "mtg")[3] == "vis_06"


def test_unpacked_counts_are_refused(tmp_path):
    # As preprocessMTG.py writes the image: int16 counts packed with the radiance calibration
    counts = np.array([[300, 400, 500], [600, 700, 800]])
    scale, offset = 0.00709761, -1.4479115
    radiance = xr.DataArray(counts * scale + offset, dims=("y", "x"), name="vis_06",
                            coords=crop_coordinates(range(4000, 4002), range(100, 103), FULL_DISC_PIXELS[2000], 2000.0),
                            attrs={"resolution": 2000, "calibration": "counts", "start_time": "2026-08-23 12:00:00"})
    path = tmp_path / "image.nc"
    radiance.to_netcdf(path, encoding={"vis_06": {"dtype": "int16", "scale_factor": scale, "add_offset": offset,
                                                  "_FillValue": FILL}})

    with xr.open_dataset(path) as ds:
        with pytest.raises(ValueError, match="mask_and_scale=False"):
            image_counts(ds.vis_06, "mtg")
    with xr.open_dataset(path, mask_and_scale=False) as ds:
        disc = image_counts(ds.vis_06, "mtg")[0]
        np.testing.assert_array_equal(disc[4000:4002, 100:103], counts)


def test_reflectance_back_to_counts():
    when = datetime(2026, 8, 23, 12)
    counts = np.array([[300.0, 2000.0]])
    scale, offset, irradiance = 0.00709761, -1.4479115, 65.72778320
    d = sun_earth_distance_correction(datetime(2026, 8, 23, 12, 5))
    reflectance = 100 * (counts * scale + offset) * np.pi * d ** 2 / irradiance
    np.testing.assert_allclose(to_counts(reflectance, "vis_06", when), counts, rtol=1e-9)
