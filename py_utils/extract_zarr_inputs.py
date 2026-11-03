"""
Extracts inputs for one satellite image from the zarr archives, and writes them on the
specMAGIC output grid, so that magic can use them instead of the climatologies and albedo maps.

    python extract_zarr_inputs.py <image_list> <config> <satellite_toml> <out_dir> <list_out>
        [--cams PATH] [--dini PATH] [--lsasaf PATH] [--strict]

Every archive is optional. Whatever cannot be extracted is reported, and magic falls back to
the climatology (or the land-use albedo) for it. With --strict, any such problem makes this
script exit with 1 instead.

One NetCDF file is written per archive and source time, and their paths are listed in
<list_out>, which magic reads through the ZARR_INPUTS environment variable. A file that an
earlier run already wrote for the same source time and grid is reused, not extracted again.
"""
import os
import sys
import argparse
import tomllib
from datetime import datetime, timezone

import numpy as np
import xarray as xr
from pyproj import CRS, Transformer

from verify_zarr_inputs import check_cams, check_dini, check_lsasaf

# 1 Dobson Unit of ozone in kg m-2
DU_IN_KG_M2 = 2.1415e-5

# How far a source time may lie from the image time. The nearest available time wins,
# and a later one is only used if it is already in the store.
HOURLY_MAX_BEFORE = np.timedelta64(3, "h")
HOURLY_MAX_AFTER = np.timedelta64(30, "m")
DAILY_MAX_BEFORE = np.timedelta64(10, "D")

# LSA SAF pixels whose albedo was last updated longer ago than this are dropped
MAX_ALBEDO_AGE_DAYS = 10

# Requested inputs that could not be extracted
problems = []


def warn(message):
    print(f"[WARNING] zarr inputs: {message}", file=sys.stderr, flush=True)


def problem(message):
    """A requested input that cannot be used. Fatal with --strict."""
    problems.append(message)
    warn(message)


## ------ What the model grid and image look like ------

def config_values(config_path):
    """The value lines of magic-config.asc, as read positionally by loadConfig in src/read.cpp."""
    values = []
    with open(config_path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                values.append(line)
    return values


def output_grid(config_path, subsatellite_longitude_deg):
    """
    The lat/lon axes that Geography::makeCoordinates builds in headers/navigation.hpp.
    magic checks the axes of every file it reads, so any drift here is caught there.
    """
    v = config_values(config_path)
    latdim, londim = int(v[27]), int(v[28])

    # The config is parsed into floats, then widened, so mirror that rounding
    latbegin, lonbegin, dxy = (float(np.float32(x)) for x in v[29:32])

    lat = latbegin + dxy * np.arange(latdim)
    lon = lonbegin + dxy * np.arange(londim) + subsatellite_longitude_deg
    lon = np.where(lon > 180.0, lon - 360.0, lon)

    attrs = {"latbegin": latbegin, "lonbegin": lonbegin, "dxy": dxy,
             "subsatellite_longitude_deg": subsatellite_longitude_deg}
    return lat, lon, attrs


def image_time(image_list):
    """Start time of the image, parsed from its file name as Image::read_from_list does."""
    with open(image_list) as f:
        name = next(line.strip() for line in f if line.strip())

    stamp = os.path.basename(name).split("-")[2][:12]
    return np.datetime64(datetime.strptime(stamp, "%Y%m%d%H%M"), "s")


def subsatellite_longitude(satellite_toml):
    with open(satellite_toml, "rb") as f:
        return float(tomllib.load(f)["navigation"]["subsatellite_longitude_deg"])


## ------ Choosing the source time ------

def preferred_order(times, target, max_before, max_after):
    """Indices of the usable times, best first: nearest, and the earlier of two equally near."""
    times = np.asarray(times).astype("datetime64[s]")
    diff = times - target
    idx = np.flatnonzero((diff >= -max_before) & (diff <= max_after))
    return idx[np.lexsort((times[idx], np.abs(diff[idx])))]


def first_with_data(order, load, cached, label):
    """
    The first candidate that an earlier run already extracted, or that has any data on the
    output grid. Returns its index and fields, where the fields are None if cached.
    """
    for i in order:
        if cached(i):
            return i, None
        fields = load(i)
        if any(np.isfinite(f).any() for f in fields.values() if f.dtype.kind == "f"):
            return i, fields
    problem(f"{label}: no usable time within the allowed window")
    return None, None


def unusable(source, path, errors):
    """
    Reports where an archive does not match the specification (verify_zarr_inputs.py).
    Returns the variables that cannot be used, or None if the whole archive cannot.
    """
    for _, message in errors:
        problem(f"{source}: {message} in {path}")
    if any(var is None for var, _ in errors):
        return None
    return {var for var, _ in errors}


def stamp(when):
    return str(np.datetime_as_string(when, unit="m"))


def forecast(run, step):
    """The run and step of a forecast field, for its provenance."""
    if not np.issubdtype(np.asarray(step).dtype, np.timedelta64):
        step = np.timedelta64(int(step), "h")
    return f"run {stamp(run)} +{int(step / np.timedelta64(1, 'h'))} h"


## ------ Earlier extractions ------

def source_path(out_dir, source, when):
    return os.path.join(out_dir, f"{source}-{stamp(when).replace('-', '').replace('T', '').replace(':', '')}.nc")


def on_grid(ds, lat, lon, grid_attrs):
    """Whether a file was written for this output grid."""
    return (ds.sizes.get("lat") == len(lat) and ds.sizes.get("lon") == len(lon)
            and all(k in ds.attrs and np.isclose(ds.attrs[k], v) for k, v in grid_attrs.items()))


def has_cached(path, names, lat, lon, grid_attrs, provenance=None):
    """
    Whether the file at path already holds all of names on this grid. If provenance is
    given, it must match too, so that a newer DINI run replaces an older one.
    """
    if not os.path.exists(path):
        return False
    try:
        with xr.open_dataset(path) as ds:
            if not on_grid(ds, lat, lon, grid_attrs):
                return False
            return all(n in ds and (provenance is None or ds[n].attrs.get("provenance") == provenance)
                       for n in names)
    except (OSError, ValueError):
        return False


## ------ Regridding ------

def axis_weights(src, dst):
    """Lower neighbour and upper weight of each dst point on an ascending src axis."""
    if not np.all(np.diff(src) > 0):
        raise ValueError("axis_weights needs a strictly ascending source axis")

    i = np.clip(np.searchsorted(src, dst, side="right") - 1, 0, len(src) - 2)
    w = (dst - src[i]) / (src[i + 1] - src[i])
    inside = (dst >= src[0]) & (dst <= src[-1])
    return i, w, inside


def bilinear_regular(src_lat, src_lon, values, lat, lon):
    """Bilinear interpolation from a global lat/lon grid with 0..360 longitudes."""

    # Ascending latitude, and a wrapped column so the antimeridian interpolates
    order = np.argsort(src_lat)
    src_lat, values = src_lat[order], values[order]
    src_lon = np.append(src_lon, src_lon[0] + 360.0)
    values = np.concatenate([values, values[:, :1]], axis=1)

    iy, wy, oky = axis_weights(src_lat, lat)
    ix, wx, okx = axis_weights(src_lon, np.mod(lon, 360.0))
    wy, wx = wy[:, None], wx[None, :]

    out = ((1 - wy) * (1 - wx) * values[np.ix_(iy, ix)]
           + (1 - wy) * wx * values[np.ix_(iy, ix + 1)]
           + wy * (1 - wx) * values[np.ix_(iy + 1, ix)]
           + wy * wx * values[np.ix_(iy + 1, ix + 1)])

    out[~oky, :] = np.nan
    out[:, ~okx] = np.nan
    return out.astype(np.float32)


def bilinear_points(src_y, src_x, values, y, x):
    """Bilinear interpolation of scattered points on a regular, ascending y/x grid."""
    iy, wy, oky = axis_weights(src_y, y)
    ix, wx, okx = axis_weights(src_x, x)

    out = ((1 - wy) * (1 - wx) * values[iy, ix]
           + (1 - wy) * wx * values[iy, ix + 1]
           + wy * (1 - wx) * values[iy + 1, ix]
           + wy * wx * values[iy + 1, ix + 1])

    out[~(oky & okx)] = np.nan
    return out


## ------ The sources ------
#
# Each returns a list of (name, field, source time, attrs), where field is None if an
# earlier run already wrote it for that source time.

def extract_cams(path, target, lat, lon, grid_attrs, out_dir):
    """Aerosol and ozone, plus water vapour as the backup for DINI."""
    ds = xr.open_zarr(path)
    broken = unusable("cams", path, check_cams(ds))
    if broken is None:
        return []
    fields = [("aod", "aod550", 1.0), ("o3", "gtco3", 1.0 / DU_IN_KG_M2), ("h2o_cams", "tcwv", 1.0)]

    times = ds.valid_time.values
    order = preferred_order(times, target, HOURLY_MAX_BEFORE, HOURLY_MAX_AFTER)

    # The forecast each valid time comes from, if the archive records it
    has_forecast = "forecast_reference_time" in ds.coords and "forecast_period" in ds.coords
    runs = ds.forecast_reference_time.values if has_forecast else None
    steps = ds.forecast_period.values if has_forecast else None

    results = []
    for name, var, scale in fields:
        label = f"cams {var}"
        if var in broken:
            continue

        def provenance(i):
            run = f"{forecast(runs[i], steps[i])}, " if has_forecast else ""
            return f"cams {var} {run}valid {stamp(times[i])} from {path}"

        def load(i):
            da = ds[var].isel(valid_time=i)
            values = da.values.astype(np.float64) * scale
            return {name: bilinear_regular(da.latitude.values, da.longitude.values, values, lat, lon)}

        # Matching the provenance too, so that a valid time whose forecast was replaced
        # in the archive is extracted again
        def cached(i):
            return has_cached(source_path(out_dir, "cams", times[i]), [name], lat, lon, grid_attrs,
                              provenance(i))

        i, loaded = first_with_data(order, load, cached, label)
        if i is None:
            continue

        results.append((name, None if loaded is None else loaded[name], times[i],
                        {"provenance": provenance(i)}))
    return results


def extract_dini(path, target, lat, lon, grid_attrs, out_dir):
    """Water vapour from the most recent DINI run valid at the chosen time."""
    ds = xr.open_zarr(path)
    if unusable("dini", path, check_dini(ds)) != set():
        return []
    label = "dini wvint"

    # Only transform the output points that can lie inside the DINI domain
    dini_lon = np.mod(ds.longitude.values + 180.0, 360.0) - 180.0
    dini_lat = ds.latitude.values
    in_box = ((lat[:, None] >= dini_lat.min()) & (lat[:, None] <= dini_lat.max())
              & (lon[None, :] >= dini_lon.min()) & (lon[None, :] <= dini_lon.max()))
    if not in_box.any():
        problem(f"{label}: the output grid does not overlap the DINI domain")
        return []

    # Every (run, step) pair, and the time it is valid at
    runs = ds.forecast_reference_time.values.astype("datetime64[s]")
    steps = ds.step.values
    if not np.issubdtype(steps.dtype, np.timedelta64):
        steps = steps.astype("timedelta64[h]")
    run_idx, step_idx = np.meshgrid(np.arange(len(runs)), np.arange(len(steps)), indexing="ij")
    run_idx, step_idx = run_idx.ravel(), step_idx.ravel()
    valid = runs[run_idx] + steps[step_idx]

    # Nearest valid time first, and for the same valid time the most recent run
    diff = valid - target
    ok = np.flatnonzero((diff >= -HOURLY_MAX_BEFORE) & (diff <= HOURLY_MAX_AFTER))
    order = ok[np.lexsort((-runs[run_idx[ok]].astype(np.int64), valid[ok], np.abs(diff[ok])))]

    def provenance(k):
        return (f"dini wvint {forecast(runs[run_idx[k]], steps[step_idx[k]])}, "
                f"valid {stamp(valid[k])} from {path}")

    lat2d, lon2d = np.broadcast_arrays(lat[:, None], lon[None, :])
    projected = []

    def load(k):
        if not projected:
            crs = CRS.from_wkt(ds.dini_projection.attrs["crs_wkt"])
            projected.extend(Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(
                lon2d[in_box], lat2d[in_box]))
        x, y = projected

        values = ds.wvint.isel(forecast_reference_time=run_idx[k], step=step_idx[k]).values
        field = np.full(lat2d.shape, np.nan, dtype=np.float32)
        field[in_box] = bilinear_points(ds.y.values, ds.x.values, values, y, x)
        return {"h2o": field}

    def cached(k):
        return has_cached(source_path(out_dir, "dini", valid[k]), ["h2o"], lat, lon, grid_attrs,
                          provenance(k))

    k, loaded = first_with_data(order, load, cached, label)
    if k is None:
        return []

    return [("h2o", None if loaded is None else loaded["h2o"], valid[k],
             {"provenance": provenance(k)})]


def extract_lsasaf(path, target, lat, lon, grid_attrs, out_dir):
    """
    Visible and near-infrared black-sky albedo, used as they are, plus the land/water mask
    from QFLAGS, so that magic can tell normal gaps over water from missing land albedo.
    """
    ds = xr.open_zarr(path)
    if unusable("lsasaf", path, check_lsasaf(ds)) is None:
        return []

    names = ["albedo_vis", "albedo_nir", "albedo_land"]
    times = ds.time.values

    # Daily product: the image's own day, or the most recent day before it
    day = target.astype("datetime64[D]").astype("datetime64[s]")
    order = preferred_order(times, day, DAILY_MAX_BEFORE, np.timedelta64(0, "s"))

    def cached(i):
        return has_cached(source_path(out_dir, "lsasaf", times[i]), names, lat, lon, grid_attrs)

    lat2d, lon2d = np.broadcast_arrays(lat[:, None], lon[None, :])
    disc = {}

    def locate():
        """Nearest LSA SAF pixel of every output point that lies on the disc."""
        crs = CRS.from_cf(ds.geostationary.attrs)
        height = ds.geostationary.attrs["perspective_point_height"]
        x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(lon2d, lat2d)

        # Geostationary scan angles in radians; points off the disc come back non-finite
        sx, sy = ds.x.values, ds.y.values
        ix = np.rint((x / height - sx[0]) / ((sx[-1] - sx[0]) / (len(sx) - 1)))
        iy = np.rint((y / height - sy[0]) / ((sy[-1] - sy[0]) / (len(sy) - 1)))
        on_disc = (np.isfinite(ix) & np.isfinite(iy)
                   & (ix >= 0) & (ix < len(sx)) & (iy >= 0) & (iy < len(sy)))
        disc.update(on_disc=on_disc, ix=ix[on_disc].astype(int), iy=iy[on_disc].astype(int))

    def load(i):
        if not disc:
            locate()
        on_disc, ix, iy = disc["on_disc"], disc["ix"], disc["iy"]
        if not on_disc.any():
            raise ValueError("the output grid is not on the LSA SAF disc")

        # Only read the window that the output grid covers
        window = dict(time=i, y=slice(iy.min(), iy.max() + 1), x=slice(ix.min(), ix.max() + 1))

        def sample(var):
            values = ds[var].isel(**window).values
            return values[iy - iy.min(), ix - ix.min()]

        age = sample("Z_Age")
        fresh = (age >= 0) & (age <= MAX_ALBEDO_AGE_DAYS)

        fields = {}
        for name, var in (("albedo_vis", "AL-VI-DH"), ("albedo_nir", "AL-NI-DH")):
            v = sample(var).astype(np.float32)
            v[~fresh | (v < 0) | (v > 1)] = np.nan
            fields[name] = np.full(lat2d.shape, np.nan, dtype=np.float32)
            fields[name][on_disc] = v

        # Bit 0 of QFLAGS is set on land. Off the disc counts as not land.
        land = np.zeros(lat2d.shape, dtype=np.int8)
        land[on_disc] = sample("QFLAGS").astype(np.uint8) & 1
        fields["albedo_land"] = land
        return fields

    i, loaded = first_with_data(order, load, cached, "lsasaf albedo")
    if i is None:
        return []

    when = times[i]
    if when.astype("datetime64[D]") != target.astype("datetime64[D]"):
        warn(f"lsasaf: no product for the image's day, using {str(when)[:10]}")

    results = []
    for name, var in (("albedo_vis", "AL-VI-DH"), ("albedo_nir", "AL-NI-DH"), ("albedo_land", "QFLAGS")):
        attrs = {"provenance": f"lsasaf {var} day {str(when)[:10]} from {path}"}
        results.append((name, None if loaded is None else loaded[name], when, attrs))
    return results


## ------ Output ------

def write_source_files(source, results, lat, lon, grid_attrs, target, out_dir):
    """
    One file per source time, written atomically, keeping what an earlier run already
    wrote there for this grid. Returns the paths of all files used, written or reused.
    """
    used = []
    for when in sorted({r[2] for r in results}):
        path = source_path(out_dir, source, when)
        new = {name: (field, attrs) for name, field, t, attrs in results if t == when and field is not None}

        if not new:
            print(f" Reusing {path}")
            used.append(path)
            continue

        data_vars = {}
        if os.path.exists(path):
            with xr.open_dataset(path) as old:
                if on_grid(old, lat, lon, grid_attrs):
                    data_vars = {n: old[n].load().variable for n in old.data_vars if n not in new}

        for name, (field, attrs) in new.items():
            data_vars[name] = xr.Variable(("lat", "lon"), field, attrs=attrs)

        ds = xr.Dataset(data_vars, coords={"lat": ("lat", lat), "lon": ("lon", lon)})
        ds.attrs = {**grid_attrs, "source": source, "source_time": stamp(when),
                    "image_time": stamp(target),
                    "created": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}

        encoding = {name: {"zlib": True, "complevel": 1, "_FillValue": None,
                           "dtype": "int8" if name == "albedo_land" else "float32"}
                    for name in data_vars}
        tmp = path + ".tmp"
        ds.to_netcdf(tmp, engine="netcdf4", encoding=encoding)
        os.replace(tmp, path)

        print(f" Wrote {', '.join(new)} valid {stamp(when)} to {path}")
        used.append(path)
    return used


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("image_list")
    parser.add_argument("config")
    parser.add_argument("satellite_toml")
    parser.add_argument("out_dir")
    parser.add_argument("list_out")
    parser.add_argument("--cams", default="")
    parser.add_argument("--dini", default="")
    parser.add_argument("--lsasaf", default="")
    parser.add_argument("--strict", action="store_true",
                        help="exit with 1 if a requested input cannot be extracted")
    args = parser.parse_args()

    target = image_time(args.image_list)
    lat, lon, grid_attrs = output_grid(args.config, subsatellite_longitude(args.satellite_toml))
    os.makedirs(args.out_dir, exist_ok=True)

    print(f" Extracting zarr inputs for image time {target}...")

    used = []
    for source, path, extract in (("cams", args.cams, extract_cams),
                                  ("dini", args.dini, extract_dini),
                                  ("lsasaf", args.lsasaf, extract_lsasaf)):
        if not path:
            continue
        try:
            results = extract(path, target, lat, lon, grid_attrs, args.out_dir)
        except Exception as e:
            problem(f"{source}: could not extract from {path} ({type(e).__name__}: {e})")
            continue
        used += write_source_files(source, results, lat, lon, grid_attrs, target, args.out_dir)

    with open(args.list_out, "w") as f:
        f.writelines(p + "\n" for p in used)

    if problems and args.strict:
        print(f"[ERROR] zarr inputs: {len(problems)} requested input(s) could not be extracted "
              "and --strict was given:", file=sys.stderr)
        for p in problems:
            print(f"    {p}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
