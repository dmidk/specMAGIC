"""
specMAGIC from Python: GHI, DNI, CAL and CSR from a satellite image, with optional input
fields for its time, all given as xarray DataArrays.

    import xarray as xr
    import specmagic

    image = xr.open_zarr("vis.zarr").vis_06.sel(start_time="2026-07-15T12:20")
    cams = xr.open_zarr("cams.zarr")
    out = specmagic.run(image, aerosol=cams.aod550, ozone=cams.gtco3)

The input fields can be on their own grids, as from CAMS, DINI or LSA SAF. For each, the time
slice that suits the image is picked and regridded onto the output grid here, in Python, and
then passed to the same C++ code that the CLI runs.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from . import _magic
from ._image import image_counts, write_image
from ._inputs import KINDS, all_land, grid_definition, lsasaf_albedo, output_grid, write_inputs
from ._inputs import regrid_input as _regrid_input

__all__ = ["run", "regrid_input", "lsasaf_albedo", "write_inputs", "write_image", "output_grid"]

OUTPUT_UNITS = {"GHI": "W m-2", "DNI": "W m-2", "CAL": "%", "CSR": "W m-2"}
FILL = -1


def default_home():
    """The directory with magic-config.asc, luts, climatologies and satellites: the repository."""
    if os.environ.get("SPECMAGIC_HOME"):
        return Path(os.environ["SPECMAGIC_HOME"])
    repo = Path(__file__).resolve().parents[2]
    if (repo / "magic-config.asc").exists():
        return repo
    raise ValueError("cannot find specMAGIC's tables; pass home= or set SPECMAGIC_HOME to the repository")


def regrid_input(da, kind, image_time, grid=None, satellite="mtg", home=None):
    """
    One input field for an image at image_time, picked from da and regridded onto the output
    grid, in the units magic uses. kind is one of aerosol, ozone, water_vapour, albedo_vis,
    albedo_nir or land. The result can be passed to run as it is, combined with another one
    first (e.g. DINI water vapour with CAMS outside the DINI domain, dini.fillna(cams)), or
    written for the CLI with write_inputs.
    """
    home = Path(home) if home else default_home()
    lat, lon = output_grid(home, satellite, grid)
    return _regrid_input(da, kind, np.datetime64(pd.Timestamp(image_time).to_datetime64(), "s"), lat, lon)


def run(image, *, aerosol=None, ozone=None, water_vapour=None, albedo_vis=None, albedo_nir=None, land=None,
        albedo="LANDMAP", strict=False, grid=None, satellite="mtg", home=None):
    """
    Runs specMAGIC on one image, and returns GHI, DNI, CAL and CSR on the output grid.

    image: the satellite channel as satpy gives it, named after the channel, with the
        attributes resolution and calibration ("counts", or "reflectance" for MTG), and a
        start time; the full disc, or a crop of it with its x/y coordinates.
    aerosol, ozone, water_vapour: optional fields that replace the climatologies, e.g. CAMS
        aod550 and gtco3, and DINI wvint. Pixels without a value use the climatology.
    albedo_vis, albedo_nir, land: an optional black-sky albedo for the local noon, as from
        LSA SAF (see lsasaf_albedo), visible (0.4-0.7 um) and near-infrared (0.7-4 um), with
        a land mask (1 on land); without one, all pixels count as land. Replaces the land-use
        albedo, so albedo must then be LANDMAP.
    albedo: LANDMAP or MODIS, the albedo where there is no albedo input.
    strict: raise instead of falling back when a given input is not usable. Otherwise an
        unusable input is reported, and recorded in the attribute input_problems.
    grid: (lat_min, lat_max, lon_min, lon_max, resolution) in degrees instead of the grid of
        magic-config.asc.
    home: the directory with magic-config.asc, luts, climatologies and satellites.

    Each input can be on its own grid with a time dimension: the slice for the image is
    picked and regridded as in regrid_input. Fields already on the output grid (dims lat,
    lon) are used as they are.
    """
    home = Path(home) if home else default_home()
    counts, resolution, when, channel = image_counts(image, satellite)
    lat, lon = output_grid(home, satellite, grid)

    given = {"aerosol": aerosol, "ozone": ozone, "water_vapour": water_vapour,
             "albedo_vis": albedo_vis, "albedo_nir": albedo_nir, "land": land}
    given = {kind: da for kind, da in given.items() if da is not None}
    if ("albedo_vis" in given) != ("albedo_nir" in given):
        raise ValueError("albedo_vis and albedo_nir must be given together")
    if "land" in given and "albedo_vis" not in given:
        raise ValueError("land is the land mask of the albedo input, so it needs albedo_vis and albedo_nir")

    # An input that cannot be put on the output grid is handed to magic as a problem, which
    # falls back to the climatology, or stops the run if strict, as for the CLI's files
    fields, problems = {}, {}
    for kind, da in given.items():
        try:
            fields[kind] = _regrid_input(da, kind, np.datetime64(when, "s"), lat, lon)
        except ValueError as e:
            problems[KINDS[kind]["name"]] = f"{kind}: {e}"
    if "albedo_vis" in fields and "land" not in given:
        fields["land"] = all_land(fields["albedo_vis"])

    arrays = {f.name: (np.ascontiguousarray(f.values, dtype=np.float32), f.attrs["provenance"])
              for f in fields.values()}
    raw = _magic.run(f"{home}/", channel, satellite, counts, resolution,
                     (when.year, when.month, when.day, when.hour, when.minute), arrays, problems, albedo, strict,
                     grid_definition(grid) if grid is not None else None)

    attrs = dict(raw["attributes"])

    out = xr.Dataset({name: (("lat", "lon"), raw[name], {"units": units}) for name, units in OUTPUT_UNITS.items()},
                     coords={"lat": raw["lat"], "lon": raw["lon"], "time": np.datetime64(when, "ns")},
                     attrs=attrs)
    for name in OUTPUT_UNITS:
        out[name].encoding["_FillValue"] = FILL
    return out
