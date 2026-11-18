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
then passed to the same C++ code that the CLI runs. The image can also be a series, such as a
day of images, with one more dimension: specmagic.run(images.sel(start_time="2026-07-15"), ...).
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from . import _magic, _mtg
from ._inputs import KINDS, Regridded, all_land, grid_definition, lsasaf_albedo, output_grid, write_inputs
from ._inputs import regrid_input as _regrid_input

__all__ = ["run", "regrid_input", "lsasaf_albedo", "write_inputs", "write_image", "output_grid"]

OUTPUT_UNITS = {"GHI": "W m-2", "DNI": "W m-2", "CAL": "%", "CSR": "W m-2"}
FILL = -1

# The module that brings an image of each satellite to the counts magic reads. Another
# satellite needs its own module like _mtg.py, and a satellites/<name>.toml.
IMAGES = {"mtg": _mtg}


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
    Pick the time of an input field for an image, and regrid it onto the output grid.

    This is what ``run`` does with each input. Call it yourself to combine two sources
    before the run, e.g. DINI water vapour with CAMS around the DINI domain, or to write
    files for the CLI with ``write_inputs``.

    Parameters
    ----------
    da : xarray.DataArray
        The input field, as described for the inputs of ``run``.
    kind : {"aerosol", "ozone", "water_vapour", "albedo_vis", "albedo_nir", "land"}
        Which input it is.
    image_time : str, datetime or numpy.datetime64
        When the image was taken.
    grid, satellite, home
        As for ``run``.

    Returns
    -------
    xarray.DataArray
        The field on the output grid, dims ``(lat, lon)``, float32, in magic's units, with
        NaN where it has no value. Its ``provenance`` attribute says where it came from.
        ``run`` takes it as it is.

    Raises
    ------
    ValueError
        If the field cannot be used: unknown grid or units, or no value on the output grid
        at a suitable time.

    Examples
    --------
    >>> when = image.start_time.values
    >>> h2o = specmagic.regrid_input(dini.wvint, "water_vapour", when).fillna(
    ...     specmagic.regrid_input(cams.tcwv, "water_vapour", when))
    >>> out = specmagic.run(image, water_vapour=h2o)
    """
    home = Path(home) if home else default_home()
    lat, lon = output_grid(home, satellite, grid)
    return _regrid_input(da, kind, np.datetime64(pd.Timestamp(image_time).to_datetime64(), "s"), lat, lon)


def satellite_image(satellite):
    """The module that reads images of the satellite."""
    if satellite not in IMAGES:
        raise NotImplementedError(f"specmagic cannot read images of {satellite!r} yet, only of {', '.join(IMAGES)}")
    return IMAGES[satellite]


def write_image(image, directory, satellite="mtg"):
    """
    Write the image as the NetCDF file that the CLI reads.

    Parameters
    ----------
    image : xarray.DataArray
        One image, as described for ``run``.
    directory : str or pathlib.Path
        Where to write the file. It is named as ``preprocessMTG.py`` names its files.
    satellite : str, default "mtg"
        The satellite that took the image.

    Returns
    -------
    pathlib.Path
        The path of the file.
    """
    return satellite_image(satellite).write_image(image, directory)


def run(image, *, aerosol=None, ozone=None, water_vapour=None, albedo_vis=None, albedo_nir=None, land=None,
        albedo="LANDMAP", strict=False, grid=None, satellite="mtg", home=None):
    """
    Run specMAGIC on a satellite image, or on a series of images.

    All arrays are xarray DataArrays, from anything xarray can open: a satpy Scene, zarr
    archives or NetCDF files. Lazily opened data is only read where it is needed, and
    nothing is written to disk. The optional input fields replace the climatologies. They
    can be on their own grids and have a time dimension: for each image, the right time is
    picked and regridded onto the output grid.

    Parameters
    ----------
    image : xarray.DataArray
        One channel of the satellite, as satpy gives it.

        - dims: ``(y, x)`` for one image. For a series, add one dimension in front, of any
          name, e.g. ``(start_time, y, x)``.
        - coords:

          - ``x`` and ``y`` in metres from the centre of the disc, as satpy gives them.
            Needed for a crop of the disc; optional for the full disc.
          - ``start_time`` or ``time``: when the image was taken. For a series, a
            coordinate along the series dimension, with a different time for each image.
        - attrs:

          - ``name``: the channel, e.g. ``"vis_06"``. If missing, the DataArray's name is used.
          - ``resolution``: in metres, e.g. 2000. MTG has 500, 1000 and 2000.
          - ``calibration``: ``"counts"`` (the default) or ``"reflectance"``.
          - ``start_time``: only used if there is no time coordinate.
        - values: counts or reflectance, see below. NaN marks pixels without data.

        magic computes the cloud index from the raw FCI counts, so the image must give them,
        in one of two ways:

        - Counts (recommended): load them with ``scene.load([...], calibration="counts")``.
          They are used as they are. A NetCDF file of counts packed with the radiance
          calibration, as ``preprocessMTG.py`` writes them, must be opened with
          ``xr.open_dataset(path, mask_and_scale=False)``. Otherwise xarray unpacks the
          counts to radiance, and the image is refused.
        - Reflectance in %, satpy's default, for the MTG channels vis_04 to nir_22. It is
          turned back into counts with the FCI calibration (``_mtg.FCI_CALIBRATION``) and
          the Sun-Earth distance, then rounded. A warning is given if satpy's
          ``radiance_unit_conversion_coefficient`` attribute does not match that
          calibration.
    aerosol : xarray.DataArray, optional
        Aerosol optical depth at 550 nm, e.g. CAMS ``aod550``.

        - dims: a grid and an optional time; see *Input grids* and *Input times* below.
          E.g. CAMS: ``(valid_time, latitude, longitude)``.
        - attrs: ``units`` must be dimensionless: ``""``, ``"1"`` or ``"~"``.
    ozone : xarray.DataArray, optional
        Total column ozone, e.g. CAMS ``gtco3``.

        - dims: a grid and an optional time, as for ``aerosol``.
        - attrs: ``units`` must be ``"DU"`` or ``"kg m-2"`` (also written ``"kg m**-2"``).
          kg m-2 is converted to DU.
    water_vapour : xarray.DataArray, optional
        Total column water vapour, e.g. DINI ``wvint`` or CAMS ``tcwv``.

        - dims: a grid and an optional time, as for ``aerosol``. E.g. DINI:
          ``(forecast_reference_time, step, y, x)``.
        - attrs: ``units`` must be ``"kg m-2"`` (also written ``"kg m**-2"``).
    albedo_vis, albedo_nir : xarray.DataArray, optional
        Black-sky albedo for the sun at local noon: visible (0.4-0.7 um) and near-infrared
        (0.7-4 um). E.g. LSA SAF ``AL-VI-DH`` and ``AL-NI-DH``; ``lsasaf_albedo`` makes
        both, and ``land``, from an LSA SAF dataset. Give both or neither.

        - dims: a grid and an optional daily time; see *Input grids* and *Input times*.
          E.g. LSA SAF: ``(time, y, x)``.
        - coords: those of the grid and the time, e.g. for LSA SAF ``time``, ``x`` and
          ``y`` (in radians) and ``geostationary``. Optionally ``age_days``, see below.
        - attrs: optionally ``max_age_days``, see below. Units are not checked.
        - values: fractions from 0 to 1. Values outside are dropped.

        A pixel uses the albedo input where both albedos have a value, and the land-use
        albedo everywhere else. Both are taken from the nearest pixel of the source.

        Age (optional): the coordinate ``age_days``, with the same dims as the albedo,
        gives the age of each value in days; a negative age means no value. The attribute
        ``max_age_days`` gives the limit. A value older than the limit on the image's day
        is dropped. If an earlier day is used because the image's own day is missing, the
        days in between are added to its age.
    land : xarray.DataArray, optional
        Land mask for the albedo input: 1 on land and 0 on water, e.g. LSA SAF
        ``QFLAGS & 1``. Only used together with ``albedo_vis`` and ``albedo_nir``.

        - dims and coords: as for ``albedo_vis``.
        - values: above 0.5 counts as land. A mask already on the output grid is used as
          it is, so it must be exactly 0 or 1.

        The mask changes no output value. It only decides how pixels without an albedo
        input are reported: on water this is expected (LSA SAF has no values over the
        sea), on land it is a fallback to the land-use albedo, counted in
        ``albedo_usage`` and warned about above 1%. Without a mask, every pixel counts as
        land, so the sea counts as fallback.
    albedo : {"LANDMAP", "MODIS"}, default "LANDMAP"
        The albedo to use where there is no albedo input. MODIS needs the MODIS maps to be
        downloaded, and cannot be combined with an albedo input.
    strict : bool, default False
        If True, raise when a given input cannot be used at all. If False, that input falls
        back to the climatology, and the problem is reported and recorded in
        ``input_problems``. Single pixels without a value always fall back.
    grid : tuple of float, optional
        The output grid, ``(lat_min, lat_max, lon_min, lon_max, resolution)`` in degrees.
        By default, the grid of ``magic-config.asc``.
    satellite : str, default "mtg"
        The satellite that took the image. Only MTG is supported.
    home : str or pathlib.Path, optional
        The directory with ``magic-config.asc``, ``luts``, ``climatologies`` and
        ``satellites``. By default ``$SPECMAGIC_HOME``, or else the repository.

    Returns
    -------
    xarray.Dataset
        - data_vars: int16, with dims ``(lat, lon)`` for one image and
          ``(time, lat, lon)`` for a series.

          - ``GHI``: global horizontal irradiance in W m-2.
          - ``DNI``: direct normal irradiance in W m-2.
          - ``CAL``: effective cloud albedo in %, from 0 to 100.
          - ``CSR``: clear-sky GHI in W m-2.

          All are 0 at night. Pixels that cannot be computed (sun close to the horizon, or
          too far from the point below the satellite) have GHI -1 and the others 0.
        - coords: ``lat`` and ``lon`` in degrees, and ``time``, when the image was taken.
        - attrs (only when inputs are given): for each input, ``*_source`` says where it came
          from and which time was used, and ``*_usage`` how many of the day pixels used it.
          ``input_problems`` lists any problems. For a series, ``*_usage`` counts the day
          pixels of all images together, and ``*_source`` and ``input_problems`` become
          coordinates along ``time``, one value per image.

    Raises
    ------
    ValueError
        If the image cannot be used (wrong dims, attributes or calibration, or two images
        of a series with the same time), if only one of ``albedo_vis`` and ``albedo_nir``
        is given, if ``albedo`` is unknown, or if ``albedo="MODIS"`` is given together with
        an albedo input.
    NotImplementedError
        If the satellite is not supported.
    RuntimeError
        If ``strict`` is True and an input cannot be used, or if magic fails, e.g. on a
        climatology value outside its lookup table. In a series, this stops the whole series.

    Notes
    -----
    **Input grids.** Each input field is on one of three kinds of grid:

    - A regular latitude/longitude grid: dims ``latitude`` and ``longitude`` (or ``lat``
      and ``lon``), in any order, with 1-D coords in degrees. It may be global and use
      longitudes from 0 to 360. E.g. CAMS.
    - A projected grid: dims ``y`` and ``x``, with 1-D coords in metres, or in radians for
      a geostationary grid whose ``x`` has ``units="radian"``. The grid mapping variable
      must be a coordinate too, which needs the data opened with ``decode_coords="all"``.
      Other coords, such as 2-D latitude and longitude, are ignored. E.g. DINI (coord
      ``dini_projection``) and LSA SAF (coord ``geostationary``).
    - The output grid itself: dims exactly ``(lat, lon)``, with the coords of
      ``output_grid``, as ``regrid_input`` returns it. Such a field has no time, is used
      as it is, and its units are not checked, so it must already be in magic's units
      (dimensionless AOD, ozone in DU, water vapour in kg m-2).

    Aerosol, ozone and water vapour are interpolated bilinearly, so a regional field such
    as DINI only covers its own domain. The albedo and land mask take the nearest pixel.

    **Input times.** The time of an input field is either one ``datetime64`` dimension
    with a coordinate, of any name (e.g. ``time`` or ``valid_time``), or the two forecast
    dimensions ``forecast_reference_time`` and ``step`` (a ``timedelta64``, or a number of
    hours). After the time is picked, the field must be 2-D, so it can have no other
    dimensions. A field without a time is used for every image.

    - Aerosol, ozone and water vapour: the time nearest to the image, from 3 hours before
      it to 30 minutes after it. For forecasts, the latest run for that time.
    - Albedo and land mask (daily): the image's own day, or else the latest day up to 10
      days before it. Never a later day.

    If the picked time has no value on the output grid, the next best time is tried.

    **Sources.** ``*_source`` names the file an input was opened from. xarray records this
    for NetCDF files but not for zarr; there, set a ``source`` attribute to name it.
    Otherwise it says "memory". CAMS-style ``forecast_reference_time`` and
    ``forecast_period`` coords along the time dimension add the run and step.

    **Series.** The images of a series are run one by one, each as if run alone. Work
    they share is done only once: the tables are read once per month, each source's grid
    is matched to the output grid once, and a regridded time slice is kept for as long as
    the next images use it.

    Examples
    --------
    >>> cams = xr.open_zarr("cams.zarr")
    >>> lsasaf = xr.open_zarr("lsasaf.zarr", decode_coords="all")
    >>> images = xr.open_zarr("vis.zarr").vis_06.sel(start_time="2025-03-13")
    >>> out = specmagic.run(images, aerosol=cams.aod550, ozone=cams.gtco3,
    ...                     **specmagic.lsasaf_albedo(lsasaf), grid=(52, 59, 6, 17, 0.025))
    """
    series = [dim for dim in image.dims if dim not in ("y", "x")] if image.ndim > 2 else []
    if len(series) > 1:
        raise ValueError(f"the image has the dimensions {image.dims}; expected y, x and at most one more, "
                         "the time of a series")

    home = Path(home) if home else default_home()
    reader = satellite_image(satellite)
    lat, lon = output_grid(home, satellite, grid)

    given = {"aerosol": aerosol, "ozone": ozone, "water_vapour": water_vapour,
             "albedo_vis": albedo_vis, "albedo_nir": albedo_nir, "land": land}
    given = {kind: da for kind, da in given.items() if da is not None}
    if ("albedo_vis" in given) != ("albedo_nir" in given):
        raise ValueError("albedo_vis and albedo_nir must be given together")
    if "land" in given and "albedo_vis" not in given:
        raise ValueError("land is the land mask of the albedo input, so it needs albedo_vis and albedo_nir")

    # What is the same for several images of a series is done once: the tables are read once
    # per month, and the inputs regridded once per slice of their source
    session = _magic.Session(f"{home}/", satellite, albedo, grid_definition(grid) if grid is not None else None)
    regridded = Regridded()
    done = set()

    def run_one(one):
        counts, resolution, when, channel = reader.image_counts(one)
        if when in done:
            # A slice of a series without its own time coordinate keeps the series' start_time
            # attribute, the time of one image
            raise ValueError(f"two images of the series have the time {when}; each needs its own time, "
                             "as a coordinate start_time or time along the series dimension")
        done.add(when)
        regridded.next_image()
        return run_image(session, counts, resolution, when, channel, given, lat, lon, strict, regridded)

    if not series:
        return run_one(image)[0]
    return join_series([run_one(image.isel({series[0]: i})) for i in range(image.sizes[series[0]])])


def run_image(session, counts, resolution, when, channel, given, lat, lon, strict, regridded):
    """
    One image through magic: the output as a Dataset, and the day pixel counts behind its
    usage attributes.
    """
    # An input that cannot be put on the output grid is handed to magic as a problem, which
    # falls back to the climatology, or stops the run if strict, as for the CLI's files
    fields, problems = {}, {}
    for kind, da in given.items():
        try:
            fields[kind] = _regrid_input(da, kind, np.datetime64(when, "s"), lat, lon, regridded)
        except ValueError as e:
            problems[KINDS[kind]["name"]] = f"{kind}: {e}"
    if "albedo_vis" in fields and "land" not in given:
        fields["land"] = all_land(fields["albedo_vis"])

    arrays = {f.name: (np.ascontiguousarray(f.values, dtype=np.float32), f.attrs["provenance"])
              for f in fields.values()}
    raw = session.run(channel, counts, resolution, (when.year, when.month, when.day, when.hour, when.minute),
                      arrays, problems, strict)

    out = xr.Dataset({name: (("lat", "lon"), raw[name], {"units": units}) for name, units in OUTPUT_UNITS.items()},
                     coords={"lat": raw["lat"], "lon": raw["lon"], "time": np.datetime64(when, "ns")},
                     attrs=dict(raw["attributes"]))
    for name in OUTPUT_UNITS:
        out[name].encoding["_FillValue"] = FILL
    return out, dict(raw["usage"])


def join_series(runs):
    """
    The runs of a series of images joined along time. The usage attributes are over the day
    pixels of all the images, so an image with more day pixels weighs more. What else differs
    from image to image, the source of each input (with the time of the slice used) and any
    input problems, becomes a coordinate along time; an image without problems has "".
    """
    out = xr.concat([ds for ds, _ in runs], dim="time")

    totals = {}
    for _, usage in runs:
        for key, counts in usage.items():
            totals[key] = [a + b for a, b in zip(totals.get(key, (0, 0, 0, 0)), counts)]
    out.attrs = {key: _magic.usage_summary(*counts, albedo=key == "albedo_usage") for key, counts in totals.items()}

    for key in dict.fromkeys(key for ds, _ in runs for key in ds.attrs if key not in totals):
        out.coords[key] = ("time", [ds.attrs.get(key, "") for ds, _ in runs])
    for name in OUTPUT_UNITS:
        out[name].encoding["_FillValue"] = FILL
    return out
