"""
The MTG image as magic reads it: FCI counts on the full disc, laid out as satpy saves them
(rows from south to north). A DataArray of counts or of reflectance, of the full disc or of
a crop of it, is brought to that form here.

Everything here is specific to MTG/FCI: its calibration, the size of its full disc, and the
file names of preprocessMTG.py. Another satellite needs a module like this one, listed in
IMAGES in __init__.py.
"""
import re
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from pyorbital.astronomy import sun_earth_distance_correction

# FCI calibration of the VNIR channels: effective radiance = counts * scale + offset, the
# channel effective solar irradiance, and the radiance unit conversion coefficient. Read from
# the L1c files of 2026-08-23 in test_data (data/<channel>/measured). satpy keeps the last
# one in the attributes of reflectance, so it is checked against it, as a sign that the
# calibration of the data is the same.
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


def start_time(image):
    """The nominal start of the image's repeat cycle, from satpy's attribute or a time coordinate."""
    for value in (image.attrs.get("start_time"), image.coords.get("start_time"), image.coords.get("time")):
        if value is None:
            continue
        if isinstance(value, xr.DataArray):
            if value.size != 1:
                continue
            value = value.values.reshape(-1)[0]
        return pd.Timestamp(value).to_pydatetime().replace(tzinfo=None)
    raise ValueError("the image has no start_time attribute or coordinate")


def to_counts(reflectance, channel, when):
    """Inverts satpy's calibration: reflectance = 100 * radiance * pi * d^2 / solar irradiance."""
    scale, offset, irradiance, _ = FCI_CALIBRATION[channel]

    # Sun-Earth distance in AU, near the middle of the 10-minute scan
    d = sun_earth_distance_correction(when + timedelta(minutes=5))
    radiance = reflectance / 100.0 * irradiance / (np.pi * d ** 2)
    return (radiance - offset) / scale


def image_counts(image):
    """
    The counts of the full disc as int16 (y, x), the resolution in metres, the start time and
    the channel, from an image DataArray as satpy gives it: named after its channel, with
    the attributes resolution and calibration ("counts" or "reflectance").
    """
    # satpy names the channel in the attributes; the DataArray itself may carry the name of
    # the variable in the L1c files (effective_radiance)
    channel = image.attrs.get("name") or image.name
    resolution = int(image.attrs["resolution"])
    when = start_time(image)
    calibration = image.attrs.get("calibration", "counts")

    if "y" in image.dims and "x" in image.dims:
        image = image.transpose("y", "x")
    values = image.values.astype(np.float64)

    if calibration == "reflectance":
        if channel not in FCI_CALIBRATION:
            raise ValueError(f"reflectance can only be turned into counts for the FCI VNIR channels, not "
                             f"{channel}; pass counts instead")
        archived = re.search(r"array\(([-0-9.eE+]+)", str(image.attrs.get("radiance_unit_conversion_coefficient", "")))
        if archived and not np.isclose(float(archived.group(1)), FCI_CALIBRATION[channel][3], rtol=1e-4):
            warnings.warn(f"The calibration of {channel} differs from FCI_CALIBRATION, so the counts may be off.")
        values = to_counts(values, channel, when)
    elif calibration != "counts":
        raise ValueError(f"the image's calibration is {calibration}; expected counts or reflectance")
    else:
        # A packed file of counts, as preprocessMTG.py writes, still says counts after xarray
        # has unpacked it to radiance
        finite = values[np.isfinite(values)]
        if "scale_factor" in image.encoding or "add_offset" in image.encoding or np.any(finite != np.rint(finite)):
            raise ValueError(f"{channel} says it is counts, but xarray has unpacked it with its scale_factor "
                             "and add_offset (or it is not whole numbers); open the file with "
                             "mask_and_scale=False to keep the counts")

    counts = np.where(np.isfinite(values), np.rint(values), FILL).astype(np.int16)

    pixels = FULL_DISC_PIXELS[resolution]
    if counts.shape == (pixels, pixels) and not ("x" in image.coords and "y" in image.coords):
        return counts, resolution, when, channel

    # A crop, or a full disc with coordinates: put every pixel in its place on the full disc,
    # whose pixel centres lie at (i - pixels / 2 + 0.5) * step
    step = abs(float(image.x.values[1] - image.x.values[0]))
    cols, rows = ((v / step + pixels / 2 - 0.5) for v in (image.x.values, image.y.values))
    if np.abs(cols - np.rint(cols)).max() > 0.01 or np.abs(rows - np.rint(rows)).max() > 0.01:
        raise ValueError("the image's x/y coordinates do not lie on the pixels of the full disc")
    disc = np.full((pixels, pixels), FILL, dtype=np.int16)
    disc[np.ix_(np.rint(rows).astype(int), np.rint(cols).astype(int))] = counts
    return disc, resolution, when, channel


def write_image(image, directory):
    """Writes the image as the NetCDF file that the CLI reads, as preprocessMTG.py does, and returns its path."""
    counts, resolution, when, channel = image_counts(image)
    end = when + timedelta(minutes=10)
    cycle = when.hour * 6 + when.minute // 10 + 1

    path = Path(directory) / f"fill-{channel}-{when:%Y%m%d%H%M%S}-{end:%Y%m%d%H%M%S}-{cycle}.nc"
    attrs = {"resolution": np.int64(resolution), "start_time": f"{when:%Y-%m-%d %H:%M:%S}",
             "end_time": f"{end:%Y-%m-%d %H:%M:%S}", "calibration": "counts"}
    xr.Dataset({channel: (("y", "x"), counts, attrs)}).to_netcdf(path, encoding={channel: {"_FillValue": FILL}})
    return path
