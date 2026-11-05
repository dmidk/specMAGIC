"""
Runs specMAGIC on every MTG image in a time range, read from the MTG zarr archive, and
collects the outputs of all images in one zarr store.

    uv run python py_utils/run_specmagic.py START END OUT_STORE [--channel vis_06]
        [--cams-dini-atmosphere] [--albedo LANDMAP|MODIS|LSASAF] [--lsasaf-nc DIR] [--strict-inputs]

START and END are UTC times such as 2026-07-01T00:00, both included. magic must already be
built (build/magic). If OUT_STORE exists, only images after its last time are run, so an
interrupted run continues with the same command.

Everything magic reads and writes for an image (the image as NetCDF, the extracted zarr
inputs and the output files) lives on the RAM disk /dev/shm, and only as long as it is needed.
The outputs and the image are deleted as soon as the outputs are in the store. The extracted
inputs serve several images (CAMS and DINI an hour, LSA SAF a day), and are deleted once an
image no longer uses them: the images are run in time order, and a later image never goes
back to an earlier input.

magic reads FCI counts on the full disc, while the MTG archive holds reflectance on a crop of
it. The reflectance is turned back into counts with the FCI calibration, and the crop is put
in its place on an otherwise empty full disc. The output grid must lie inside the crop.
"""
import os
import re
import sys
import shutil
import argparse
import tempfile
import subprocess
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr
from pyproj import CRS, Transformer
from pyorbital.astronomy import sun_earth_distance_correction

import extract_zarr_inputs as zarr_inputs

REPO = Path(__file__).resolve().parents[1]
MAGIC = REPO / "build" / "magic"
SATELLITE = "mtg"
RAM_DISK = "/dev/shm"

DATA = "/dmidata/projects/energivejr-data"
MTG_ZARR = {"vis": f"{DATA}/mtg/level1c_normal_res_cropped/consolidated/vis.zarr",
            "nir": f"{DATA}/mtg/level1c_normal_res_cropped/consolidated/nir.zarr"}
CAMS_ZARR = f"{DATA}/cams/cams.zarr"
DINI_ZARR = f"{DATA}/dini/dini.zarr"
LSASAF_ZARR = f"{DATA}/mtg/lsasaf/lsasaf_data.zarr"

# FCI calibration of the VNIR channels: effective radiance = counts * scale + offset, the
# channel effective solar irradiance, and the radiance unit conversion coefficient. Read from
# the L1c files of 2026-08-23 in test_data (data/<channel>/measured). The last one is also in
# the attributes of the MTG archive, and is checked against it, as a sign that the
# calibration of the archive is the same.
FCI_CALIBRATION = {
    # channel: (scale, offset, solar irradiance, radiance unit conversion coefficient)
    "vis_04": (0.00411406, -0.83926755, 38.09846115, 50.26570892),
    "vis_05": (0.00536305, -1.0940619, 49.66482162, 38.74919891),
    "vis_06": (0.00709761, -1.4479115, 65.72778320, 24.66601944),
    "vis_08": (0.0077124, -1.5733292, 71.42111206, 13.51755047),
    "vis_09": (0.00522184, -1.0652555, 72.53574371, 11.98462009),
    "nir_13": (0.00485798, -0.9910284, 67.48144531, 5.26138878),
    "nir_16": (0.00556013, -1.1342671, 61.78792191, 3.84947896),
    "nir_22": (0.00340031, -0.6936624, 37.78647614, 1.96577501),
}

# Pixels across the full disc for each resolution in metres, and the value satpy gives
# pixels without data
FULL_DISC_PIXELS = {500: 22272, 1000: 11136, 2000: 5568}
FILL = -1

# The largest box at whole degrees that lies inside the crop of the MTG archive
DEFAULT_EXTENT = (31.0, 72.0, -19.0, 34.0)
DEFAULT_RESOLUTION = 0.025

OUTPUTS = ("GHI", "DNI", "CAL", "CSR")
ATMOSPHERE_ATTRS = ("aod_source", "o3_source", "h2o_source", "h2o_backup_source",
                    "aod_usage", "o3_usage", "h2o_usage")
ALBEDO_ATTRS = ("albedo_vis_source", "albedo_nir_source", "albedo_usage")


## ------ The MTG image ------

def open_mtg(channel):
    group = channel.split("_")[0]
    if channel not in FCI_CALIBRATION or group not in MTG_ZARR:
        sys.exit(f"[ERROR] Only the VNIR channels can be run from the MTG archive, not {channel}.")
    mtg = xr.open_zarr(MTG_ZARR[group])

    archived = re.search(r"array\(([-0-9.eE+]+)", mtg[channel].attrs.get("radiance_unit_conversion_coefficient", ""))
    if not archived or not np.isclose(float(archived.group(1)), FCI_CALIBRATION[channel][3], rtol=1e-4):
        print(f"[WARNING] The calibration of {channel} in {MTG_ZARR[group]} differs from FCI_CALIBRATION, "
              "so the counts may be off.", file=sys.stderr)
    return mtg


def disc_position(mtg, channel):
    """Size of the full disc, and the row and column of the crop's first pixel in it."""
    resolution = int(mtg[channel].attrs["resolution"])
    pixels = FULL_DISC_PIXELS[resolution]
    step = (mtg.x.values[-1] - mtg.x.values[0]) / (mtg.sizes["x"] - 1)

    # Pixel centres of the full disc lie at (i - pixels / 2 + 0.5) * step
    col, row = (v[0] / step + pixels / 2 - 0.5 for v in (mtg.x.values, mtg.y.values))
    if abs(col - round(col)) > 0.01 or abs(row - round(row)) > 0.01:
        sys.exit("[ERROR] The MTG crop does not lie on the pixels of the full disc.")
    return resolution, pixels, round(row), round(col)


def to_counts(reflectance, channel, when):
    """Inverts satpy's calibration: reflectance = 100 * radiance * pi * d^2 / solar irradiance."""
    scale, offset, irradiance, _ = FCI_CALIBRATION[channel]

    # Sun-Earth distance in AU, near the middle of the 10-minute scan
    d = sun_earth_distance_correction(when + timedelta(minutes=5))
    radiance = reflectance / 100.0 * irradiance / (np.pi * d ** 2)
    counts = np.rint((radiance - offset) / scale)
    return np.where(np.isfinite(counts), counts, FILL).astype(np.int16)


def write_image(mtg, channel, i, when, image_dir):
    """The image as preprocessMTG.py writes it: counts on the full disc. None if it is empty."""
    reflectance = mtg[channel].isel(start_time=i).values
    if not np.isfinite(reflectance).any():
        return None

    resolution, pixels, row, col = disc_position(mtg, channel)
    disc = np.full((pixels, pixels), FILL, dtype=np.int16)
    disc[row:row + reflectance.shape[0], col:col + reflectance.shape[1]] = to_counts(reflectance, channel, when)

    end = when + timedelta(minutes=10)
    cycle = when.hour * 6 + when.minute // 10 + 1
    path = image_dir / f"fill-{channel}-{when:%Y%m%d%H%M%S}-{end:%Y%m%d%H%M%S}-{cycle}.nc"
    attrs = {"resolution": np.int64(resolution), "start_time": f"{when:%Y-%m-%d %H:%M:%S}",
             "end_time": f"{end:%Y-%m-%d %H:%M:%S}", "calibration": "counts",
             "comment": f"Counts recalculated from the reflectance in {MTG_ZARR[channel.split('_')[0]]}"}
    xr.Dataset({channel: (("y", "x"), disc, attrs)}).to_netcdf(
        path, engine="netcdf4", encoding={channel: {"_FillValue": FILL}})
    return path


## ------ The run directory on the RAM disk ------

def make_run_dir(extent, resolution):
    """
    A home directory for magic: its config with the output grid replaced, the tables and
    satellite files linked from the repository, and empty output directories.
    """
    run_dir = Path(tempfile.mkdtemp(prefix="specmagic-", dir=RAM_DISK))
    for name in ("climatologies", "luts", "satellites"):
        (run_dir / name).symlink_to(REPO / name)
    for name in OUTPUTS:
        (run_dir / "out" / name).mkdir(parents=True)
    (run_dir / "mtg").mkdir()

    lat_min, lat_max, lon_min, lon_max = extent
    grid = {27: round((lat_max - lat_min) / resolution) + 1,   # latdim
            28: round((lon_max - lon_min) / resolution) + 1,   # londim
            29: lat_min, 30: lon_min, 31: resolution}          # latbegin, lonbegin, dxy

    # magic reads the config by the position of its value lines, as zarr_inputs.config_values does
    lines, value_index = [], 0
    for line in open(REPO / "magic-config.asc"):
        value, _, comment = line.partition("#")
        if value.strip():
            if value_index in grid:
                line = f"{grid[value_index]}\t\t# {comment.strip()}\n"
            value_index += 1
        lines.append(line)
    (run_dir / "magic-config.asc").write_text("".join(lines))
    return run_dir


def check_grid_on_crop(mtg, channel, lat, lon):
    """Stops if any point of the output grid lies outside the MTG crop."""
    crs = CRS.from_cf(mtg[mtg[channel].attrs["grid_mapping"]].attrs)
    edge_lat = np.concatenate([lat, lat, np.full(len(lon), lat[0]), np.full(len(lon), lat[-1])])
    edge_lon = np.concatenate([np.full(len(lat), lon[0]), np.full(len(lat), lon[-1]), lon, lon])
    x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(edge_lon, edge_lat)
    inside = (np.isfinite(x) & (x >= mtg.x.values[0]) & (x <= mtg.x.values[-1])
              & (y >= mtg.y.values[0]) & (y <= mtg.y.values[-1]))
    if not inside.all():
        sys.exit("[ERROR] The output grid reaches outside the MTG crop, where magic would have no "
                 f"image. --extent {' '.join(map(str, DEFAULT_EXTENT))} fits.")


## ------ The output store ------

def provenance_names(atmosphere, albedo):
    names = (ATMOSPHERE_ATTRS if atmosphere else ()) + (ALBEDO_ATTRS if albedo else ())
    return names + ("input_problems",) if names else ()


def last_stored_time(store, settings):
    """The last time in the store, or None if there is no store. Stops if it was run differently."""
    if not os.path.exists(store):
        return None
    with xr.open_zarr(store) as ds:
        different = {k: (ds.attrs.get(k), v) for k, v in settings.items() if ds.attrs.get(k) != v}
        if different:
            sys.exit(f"[ERROR] {store} was made with other settings (stored, now): {different}")
        return ds.time.values[-1]


def append_outputs(store, run_dir, when, names, settings):
    """Adds the four outputs of one image, and its provenance, to the store."""
    data_vars, lat, lon, attrs = {}, None, None, {}
    for name in OUTPUTS:
        path = run_dir / "out" / name / f"{name}{when:%Y%m%d%H%M}.nc"
        with xr.open_dataset(path, mask_and_scale=False, decode_times=False) as ds:
            data_vars[name] = (("time", "lat", "lon"), ds[name].values, {"units": ds[name].attrs["units"]})
            lat, lon, attrs = ds.lat.values, ds.lon.values, dict(ds.attrs)

    # The global attributes of the output files differ per image, so they become variables
    for name in names:
        data_vars[name] = (("time",), np.array([attrs.get(name, "")], dtype=object))

    ds = xr.Dataset(data_vars, coords={"time": [np.datetime64(when, "ns")], "lat": lat, "lon": lon},
                    attrs=settings)
    if not os.path.exists(store):
        encoding = {name: {"_FillValue": FILL, "chunks": (1, len(lat), len(lon))} for name in OUTPUTS}
        encoding["time"] = {"units": "minutes since 1970-01-01", "dtype": "int64"}
        ds.to_zarr(store, mode="w-", encoding=encoding)
    else:
        ds.drop_vars(["lat", "lon"]).to_zarr(store, append_dim="time")


## ------ Running ------

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("start", type=datetime.fromisoformat)
    parser.add_argument("end", type=datetime.fromisoformat)
    parser.add_argument("store", help="output zarr store, created or appended to")
    parser.add_argument("--channel", default="vis_06")
    parser.add_argument("--cams-dini-atmosphere", action="store_true",
                        help="take aerosol and ozone from CAMS and water vapour from DINI, then CAMS")
    parser.add_argument("--albedo", default="LANDMAP", choices=["LANDMAP", "MODIS", "LSASAF"])
    parser.add_argument("--lsasaf-nc", default="", metavar="DIR",
                        help="read LSA SAF from this directory of daily LSA SAF files instead of "
                             f"{LSASAF_ZARR}")
    parser.add_argument("--strict-inputs", action="store_true",
                        help="stop if a requested zarr input is missing, instead of falling back")
    parser.add_argument("--extent", type=float, nargs=4, default=DEFAULT_EXTENT,
                        metavar=("LAT_MIN", "LAT_MAX", "LON_MIN", "LON_MAX"))
    parser.add_argument("--resolution", type=float, default=DEFAULT_RESOLUTION, help="in degrees")
    parser.add_argument("--threads", type=int, default=8, help="OpenMP threads of magic")
    args = parser.parse_args()

    # Keep this script's lines in order with those of magic, which share the terminal
    sys.stdout.reconfigure(line_buffering=True)

    if not MAGIC.exists():
        sys.exit(f"[ERROR] {MAGIC} does not exist. Build magic first, e.g. with scripts/mtg.sh.")
    if args.albedo == "MODIS" and not (REPO / "climatologies" / "modis-brdf").is_dir():
        print("[WARNING] climatologies/modis-brdf is missing; run scripts/mtg.sh --albedo MODIS "
              "once to download it.", file=sys.stderr)

    settings = {"channel": args.channel, "albedo": args.albedo,
                "cams_dini_atmosphere": int(args.cams_dini_atmosphere),
                "extent": list(args.extent), "resolution": args.resolution}
    names = provenance_names(args.cams_dini_atmosphere, args.albedo == "LSASAF")
    sources = {"cams": CAMS_ZARR if args.cams_dini_atmosphere else "",
               "dini": DINI_ZARR if args.cams_dini_atmosphere else "",
               "lsasaf": LSASAF_ZARR if args.albedo == "LSASAF" else "",
               "lsasaf_nc": args.lsasaf_nc if args.albedo == "LSASAF" else ""}

    mtg = open_mtg(args.channel)
    times = mtg.start_time.values
    selected = np.flatnonzero((times >= np.datetime64(args.start)) & (times <= np.datetime64(args.end)))
    last = last_stored_time(args.store, settings)
    if last is not None:
        print(f"{args.store} already holds images up to {str(last)[:16]}; continuing after it.")
        selected = selected[times[selected] > last]
    print(f"{len(selected)} images to run.")
    if not len(selected):
        return

    run_dir = make_run_dir(args.extent, args.resolution)
    try:
        config = run_dir / "magic-config.asc"
        lat, lon, grid_attrs = zarr_inputs.output_grid(
            config, zarr_inputs.subsatellite_longitude(REPO / "satellites" / f"{SATELLITE}.toml"))
        check_grid_on_crop(mtg, args.channel, lat, lon)

        env = {**os.environ, "SATELLITE": SATELLITE, "ALBEDO": args.albedo,
               "CAMS_DINI_ATMOSPHERE": "1" if args.cams_dini_atmosphere else "0",
               "STRICT_INPUTS": "1" if args.strict_inputs else "0",
               "ZARR_INPUTS": str(run_dir / "zarr_inputs.txt"),
               "OMP_NUM_THREADS": str(args.threads), "HDF5_USE_FILE_LOCKING": "FALSE"}

        inputs_in_use = set()
        for i in selected:
            when = times[i].astype("datetime64[s]").astype(datetime)
            print(f" ------ {when:%Y-%m-%d %H:%M} UTC ------", flush=True)

            image = write_image(mtg, args.channel, i, when, run_dir / "mtg")
            if image is None:
                print(f"[WARNING] The MTG archive has no data for {when}; skipping it.", file=sys.stderr)
                continue
            (run_dir / "image_list.txt").write_text(f"{image}\n")

            zarr_inputs.problems.clear()
            used = []
            if any(sources.values()):
                used = zarr_inputs.extract_inputs(np.datetime64(when, "s"), lat, lon, grid_attrs,
                                                  run_dir / "zarr_inputs", **sources)
                if zarr_inputs.problems and args.strict_inputs:
                    sys.exit("[ERROR] Requested zarr inputs could not be extracted, and "
                             "--strict-inputs was given:\n    " + "\n    ".join(zarr_inputs.problems))
            (run_dir / "zarr_inputs.txt").write_text("".join(f"{p}\n" for p in used))

            # The earlier images are in the store, and inputs this image does not use, no
            # later image will use either
            for path in inputs_in_use - set(used):
                os.remove(path)
            inputs_in_use = set(used)

            result = subprocess.run([str(MAGIC), f"{run_dir}/", args.channel, "0"], cwd=run_dir, env=env)
            if result.returncode != 0:
                sys.exit(f"[ERROR] magic failed for {when} with exit code {result.returncode}. "
                         f"{args.store} holds every image before it.")

            append_outputs(args.store, run_dir, when, names, settings)

            # Done with this image
            image.unlink()
            for name in OUTPUTS:
                (run_dir / "out" / name / f"{name}{when:%Y%m%d%H%M}.nc").unlink()
    finally:
        shutil.rmtree(run_dir)


if __name__ == "__main__":
    main()
