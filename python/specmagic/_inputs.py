"""
The input fields: the time slice of a source that suits the image, regridded onto specMAGIC's
output grid. magic itself only takes fields already on that grid, so this is the only place
where sources on other grids are handled.
"""
from pathlib import Path

import numpy as np
import xarray as xr
from pyproj import CRS, Transformer

from . import _magic

# 1 Dobson Unit of ozone in kg m-2
DU_IN_KG_M2 = 2.1415e-5

KG_M2 = ("kg m**-2", "kg m-2")
DIMENSIONLESS = ("~", "1", "")
DOBSON = ("DU", "Dobson", "Dobson units")

# How far a source time may lie from the image time. The nearest available time wins, and a
# later one is only used if it is already in the source. Albedo is daily, so it is the
# image's own day or one before it.
HOURLY_MAX_BEFORE = np.timedelta64(3, "h")
HOURLY_MAX_AFTER = np.timedelta64(30, "m")
DAILY_MAX_BEFORE = np.timedelta64(10, "D")

# What each kind of input is called in magic (and in the CLI's files), whether it is daily,
# and how it is regridded: albedo and land by the nearest pixel, the atmosphere bilinearly
KINDS = {
    "aerosol": {"name": "aod", "daily": False, "method": "bilinear"},
    "ozone": {"name": "o3", "daily": False, "method": "bilinear"},
    "water_vapour": {"name": "h2o", "daily": False, "method": "bilinear"},
    "albedo_vis": {"name": "albedo_vis", "daily": True, "method": "nearest"},
    "albedo_nir": {"name": "albedo_nir", "daily": True, "method": "nearest"},
    "land": {"name": "albedo_land", "daily": True, "method": "nearest"},
}


## ------ The output grid ------

def grid_definition(grid):
    """
    (latdim, londim, latbegin, lonbegin, dxy) of the output grid for magic, from
    grid = (lat_min, lat_max, lon_min, lon_max, resolution) in degrees.
    """
    lat_min, lat_max, lon_min, lon_max, dxy = grid
    latdim = round((lat_max - lat_min) / dxy) + 1
    londim = round((lon_max - lon_min) / dxy) + 1
    return latdim, londim, lat_min, lon_min, dxy


def output_grid(home, satellite, grid=None):
    """
    The lat/lon axes of the output grid, as magic builds them: from magic-config.asc in home,
    or from grid = (lat_min, lat_max, lon_min, lon_max, resolution), shifted by the satellite's
    longitude.
    """
    return _magic.output_grid(f"{home}/", satellite, grid_definition(grid) if grid is not None else None)


def on_grid(da, lat, lon):
    """Whether da is already a field on the output grid, as magic takes it."""
    return (da.dims == ("lat", "lon") and da.sizes["lat"] == len(lat) and da.sizes["lon"] == len(lon)
            and np.allclose(da.lat.values, lat, atol=1e-3) and np.allclose(da.lon.values, lon, atol=1e-3))


## ------ The time slice ------

def preferred_order(times, target, max_before, max_after):
    """Indices of the usable times, best first: nearest, and the earlier of two equally near."""
    times = np.asarray(times).astype("datetime64[s]")
    diff = times - target
    idx = np.flatnonzero((diff >= -max_before) & (diff <= max_after))
    return idx[np.lexsort((times[idx], np.abs(diff[idx])))]


def stamp(when):
    return str(np.datetime_as_string(np.datetime64(when, "m"), unit="m"))


def forecast(run, step):
    """The run and step of a forecast field, for its provenance."""
    if not np.issubdtype(np.asarray(step).dtype, np.timedelta64):
        step = np.timedelta64(int(step), "h")
    return f"run {stamp(run)} +{int(step / np.timedelta64(1, 'h'))} h"


def candidates(da, target, daily):
    """
    The 2-D slices of da that suit an image at target, best first, each with a description of
    its time. Without a time dimension, da itself is the only one.
    """
    if daily:
        target = target.astype("datetime64[D]").astype("datetime64[s]")
        max_before, max_after = DAILY_MAX_BEFORE, np.timedelta64(0, "s")
    else:
        max_before, max_after = HOURLY_MAX_BEFORE, HOURLY_MAX_AFTER

    # Forecasts by run and step, as DINI: nearest valid time, then the most recent run
    if "forecast_reference_time" in da.dims and "step" in da.dims:
        runs = da.forecast_reference_time.values.astype("datetime64[s]")
        steps = da.step.values
        if not np.issubdtype(steps.dtype, np.timedelta64):
            steps = steps.astype("timedelta64[h]")
        run_idx, step_idx = (a.ravel() for a in np.meshgrid(np.arange(len(runs)), np.arange(len(steps)),
                                                           indexing="ij"))
        valid = runs[run_idx] + steps[step_idx]
        diff = valid - target
        ok = np.flatnonzero((diff >= -max_before) & (diff <= max_after))
        for k in ok[np.lexsort((-runs[run_idx[ok]].astype(np.int64), valid[ok], np.abs(diff[ok])))]:
            yield (da.isel(forecast_reference_time=run_idx[k], step=step_idx[k]),
                   f"{forecast(runs[run_idx[k]], steps[step_idx[k]])}, valid {stamp(valid[k])}")
        return

    time_dim = next((d for d in da.dims if d in da.coords and np.issubdtype(da[d].dtype, np.datetime64)), None)
    if time_dim is None:
        yield da, "without time"
        return

    # One time axis; the forecast it comes from, if recorded along it (as in CAMS)
    times = da[time_dim].values
    runs = da.coords.get("forecast_reference_time")
    steps = da.coords.get("forecast_period")
    has_forecast = runs is not None and steps is not None and runs.dims == steps.dims == (time_dim,)
    for i in preferred_order(times, target, max_before, max_after):
        label = f"valid {stamp(times[i])}" if not daily else f"day {str(times[i])[:10]}"
        if has_forecast:
            label = f"{forecast(runs.values[i], steps.values[i])}, {label}"
        yield da.isel({time_dim: i}), label


## ------ Regridding ------

def positions(axis, points):
    """Fractional index of points along a monotonic axis; NaN outside it."""
    axis = np.asarray(axis, dtype=np.float64)
    if axis[-1] > axis[0]:
        return np.interp(points, axis, np.arange(len(axis)), left=np.nan, right=np.nan)
    return len(axis) - 1 - np.interp(points, axis[::-1], np.arange(len(axis)), left=np.nan, right=np.nan)


def sample(da, iy, ix, method):
    """
    Values of the 2-D da (y, x) at the fractional indices iy, ix: bilinearly, or of the
    nearest pixel. NaN where an index is NaN. Only the window that the points need is read.
    """
    out = np.full(np.shape(iy), np.nan, dtype=np.float32)
    ok = np.isfinite(iy) & np.isfinite(ix)
    if not ok.any():
        return out
    ny, nx = da.shape

    if method == "nearest":
        jy, jx = np.rint(iy[ok]).astype(int), np.rint(ix[ok]).astype(int)
        y0, x0 = jy.min(), jx.min()
        window = da.isel({da.dims[0]: slice(y0, jy.max() + 1), da.dims[1]: slice(x0, jx.max() + 1)}).values
        out[ok] = window[jy - y0, jx - x0]
        return out

    jy = np.clip(np.floor(iy[ok]).astype(int), 0, max(ny - 2, 0))
    jx = np.clip(np.floor(ix[ok]).astype(int), 0, max(nx - 2, 0))
    wy, wx = iy[ok] - jy, ix[ok] - jx
    y0, x0 = jy.min(), jx.min()
    window = da.isel({da.dims[0]: slice(y0, min(jy.max() + 2, ny)),
                      da.dims[1]: slice(x0, min(jx.max() + 2, nx))}).values.astype(np.float64)
    jy, jx = jy - y0, jx - x0
    jy1, jx1 = np.minimum(jy + 1, window.shape[0] - 1), np.minimum(jx + 1, window.shape[1] - 1)
    out[ok] = ((1 - wy) * (1 - wx) * window[jy, jx] + (1 - wy) * wx * window[jy, jx1]
               + wy * (1 - wx) * window[jy1, jx] + wy * wx * window[jy1, jx1])
    return out


def grid_mapping(da):
    """The CRS of a projected da, from its grid mapping, which must be one of its coordinates."""
    name = da.attrs.get("grid_mapping")
    if name not in da.coords:
        name = next((c for c in da.coords if {"grid_mapping_name", "crs_wkt"} & set(da[c].attrs)), None)
    if name is None:
        raise ValueError(f"{da.name} is on an x/y grid, but its grid mapping is not one of its coordinates; "
                         "open the source with decode_coords='all', or assign the grid mapping variable "
                         "as a coordinate")
    attrs = da[name].attrs
    return (CRS.from_wkt(attrs["crs_wkt"]) if "crs_wkt" in attrs else CRS.from_cf(attrs)), attrs


def regrid(da, lat, lon, method):
    """A 2-D da on its own grid, as float32 on the output grid (lat, lon); NaN where it has no value."""
    if on_grid(da, lat, lon):
        return da.values.astype(np.float32)

    lat_name = next((n for n in ("latitude", "lat") if n in da.dims), None)
    lon_name = next((n for n in ("longitude", "lon") if n in da.dims), None)
    lat2d, lon2d = np.broadcast_arrays(lat[:, None], lon[None, :])

    # A regular latitude/longitude grid, possibly global and in 0..360
    if lat_name and lon_name:
        da = da.transpose(lat_name, lon_name)
        src_lon = da[lon_name].values
        if src_lon.max() > 180:
            lon2d = np.mod(lon2d, 360.0)
        step = abs(src_lon[1] - src_lon[0])
        if np.isclose(src_lon[-1] - src_lon[0] + step, 360.0, rtol=1e-3):
            # Global: repeat the first column after the last, so the seam interpolates
            first = da.isel({lon_name: [0]}).assign_coords({lon_name: [src_lon[0] + 360.0]})
            da = xr.concat([da, first], dim=lon_name)
        return sample(da, positions(da[lat_name].values, lat2d), positions(da[lon_name].values, lon2d), method)

    # A projected grid, such as DINI's Lambert conformal or the geostationary one of LSA SAF
    if "x" in da.dims and "y" in da.dims:
        crs, attrs = grid_mapping(da)
        da = da.transpose("y", "x")
        x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(lon2d, lat2d)

        # Geostationary coordinates can be scan angles in radians rather than metres
        if attrs.get("grid_mapping_name") == "geostationary" and str(da.x.attrs.get("units", "")).startswith("rad"):
            x, y = x / attrs["perspective_point_height"], y / attrs["perspective_point_height"]
        x = np.where(np.isfinite(x), x, np.nan)
        y = np.where(np.isfinite(y), y, np.nan)
        return sample(da, positions(da.y.values, y), positions(da.x.values, x), method)

    raise ValueError(f"cannot tell the grid of {da.name} with dims {da.dims}: it needs latitude/longitude "
                     "or x/y dimensions")


## ------ The input fields ------

def scale_to_magic(kind, da):
    """The factor that takes da to the units magic uses for this kind of input."""
    units = str(da.attrs.get("units", ""))
    if kind == "aerosol":
        accepted = {u: 1.0 for u in DIMENSIONLESS}
    elif kind == "ozone":
        accepted = {**{u: 1.0 / DU_IN_KG_M2 for u in KG_M2}, **{u: 1.0 for u in DOBSON}}
    elif kind == "water_vapour":
        accepted = {u: 1.0 for u in KG_M2}
    else:
        return 1.0
    if units not in accepted:
        raise ValueError(f"{da.name} has units '{units}', expected one of {sorted(accepted)}")
    return accepted[units]


def source_of(da):
    return da.encoding.get("source") or da.attrs.get("source") or "memory"


def regrid_input(da, kind, image_time, lat, lon):
    """
    One input field for an image at image_time, on the output grid lat/lon, in magic's units,
    with its provenance. Raises ValueError if da cannot provide it.
    """
    if kind not in KINDS:
        raise ValueError(f"unknown kind of input {kind}, expected one of {sorted(KINDS)}")
    if on_grid(da, lat, lon):
        field = da.astype(np.float32)
        field.attrs = {"provenance": da.attrs.get("provenance", f"{da.name} from {source_of(da)}")}
        return field.rename(KINDS[kind]["name"])

    target = np.datetime64(image_time, "s")
    scale = scale_to_magic(kind, da)
    for candidate, when in candidates(da, target, KINDS[kind]["daily"]):
        if candidate.ndim != 2:
            raise ValueError(f"{da.name} has dims {candidate.dims} left after picking the time; expected 2")
        values = regrid(candidate, lat, lon, KINDS[kind]["method"]) * np.float32(scale)
        if np.isfinite(values).any():
            break
    else:
        raise ValueError(f"{da.name} from {source_of(da)} has no value on the output grid at a time "
                         f"within the window of {stamp(target)}")

    if kind in ("albedo_vis", "albedo_nir"):
        values[(values < 0) | (values > 1)] = np.nan
    elif kind == "land":
        values = np.where(np.isfinite(values) & (values > 0.5), 1, 0).astype(np.float32)

    return xr.DataArray(values, dims=("lat", "lon"), coords={"lat": lat, "lon": lon},
                        name=KINDS[kind]["name"],
                        attrs={"provenance": f"{da.name} {when} from {source_of(da)}"})


def lsasaf_albedo(ds, max_age_days=10):
    """
    The albedo inputs from an LSA SAF MTDAL dataset (an archive or the daily files): the
    visible and near-infrared black-sky albedo, dropped where older than max_age_days, and
    the land mask from bit 0 of QFLAGS. Pass them on as run(..., **lsasaf_albedo(ds)).
    """
    mapping = ds["geostationary"]
    fresh = (ds["Z_Age"] >= 0) & (ds["Z_Age"] <= max_age_days)
    out = {}
    for kind, var in (("albedo_vis", "AL-VI-DH"), ("albedo_nir", "AL-NI-DH")):
        out[kind] = ds[var].where(fresh).assign_coords(geostationary=mapping)
        out[kind].attrs = {**ds[var].attrs, "grid_mapping": "geostationary"}
        out[kind].encoding = ds[var].encoding
    land = (ds["QFLAGS"].astype("uint8") & 1).assign_coords(geostationary=mapping)
    land.attrs = {"grid_mapping": "geostationary"}
    land.encoding = ds["QFLAGS"].encoding
    out["land"] = land.rename("QFLAGS")
    return out


def write_inputs(directory, **fields):
    """
    Writes on-grid input fields (from regrid_input) as the NetCDF files the CLI reads, one per
    table: aerosol.nc, ozone.nc, water_vapour.nc and albedo.nc. Returns their paths by table.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    tables = {"aerosol": ["aerosol"], "ozone": ["ozone"], "water_vapour": ["water_vapour"],
              "albedo": ["albedo_vis", "albedo_nir", "land"]}
    paths = {}
    for table, kinds in tables.items():
        given = {KINDS[k]["name"]: fields[k] for k in kinds if fields.get(k) is not None}
        if not given:
            continue
        if table == "albedo":
            given.setdefault("albedo_land", all_land(next(iter(given.values()))))

        path = directory / f"{table}.nc"
        xr.Dataset(given).to_netcdf(path, encoding={n: {"dtype": "float32", "_FillValue": None} for n in given})
        paths[table] = path
    return paths


def all_land(field):
    """The land mask to use when none is given: everything is land."""
    return xr.ones_like(field).rename("albedo_land").assign_attrs(
        provenance="no land mask given, so all pixels count as land")
