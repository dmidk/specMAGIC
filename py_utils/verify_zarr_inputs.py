"""
Checks that zarr archives have the layout that extract_zarr_inputs.py expects, as specified
in docs/inputs.md, without reading any data.

    python verify_zarr_inputs.py [--cams PATH] [--dini PATH] [--lsasaf PATH]

Exits with 1 if any archive does not match. extract_zarr_inputs.py runs the same checks on
the archives it has already opened, so a pipeline run needs no separate verification.
"""
import sys
import argparse

import numpy as np
import xarray as xr
from pyproj import CRS

KG_M2 = ("kg m**-2", "kg m-2")
DIMENSIONLESS = ("~", "1", "")

# Each check returns a list of (variable, message). The variable is None where the problem
# makes the whole archive unusable, and names the variable where only that one is.


def check_variable(ds, var, dims, units=None):
    if var not in ds:
        return [(var, f"{var} is missing")]
    errors = []
    if ds[var].dims != dims:
        errors.append((var, f"{var} has dims {ds[var].dims}, expected {dims}"))
    if units is not None and ds[var].attrs.get("units", "") not in units:
        errors.append((var, f"{var} has units '{ds[var].attrs.get('units', '')}', expected one of {units}"))
    return errors


def check_axis(ds, name, ascending=False, regular=False):
    """A 1-D numeric coordinate, strictly monotonic, optionally ascending and evenly spaced."""
    if name not in ds.coords or ds[name].ndim != 1:
        return [(None, f"{name} is missing or not a 1-D coordinate")]
    values = ds[name].values
    if values.dtype.kind not in "fiu" or len(values) < 2:
        return [(None, f"{name} must be numeric with at least 2 values")]
    step = np.diff(values)
    if not (np.all(step > 0) or (not ascending and np.all(step < 0))):
        return [(None, f"{name} must be strictly {'ascending' if ascending else 'monotonic'}")]
    if regular and not np.allclose(step, step[0], rtol=1e-3):
        return [(None, f"{name} must be evenly spaced")]
    return []


def check_times(ds, name, kind="M"):
    if name not in ds.coords or ds[name].ndim != 1:
        return [(None, f"{name} is missing or not a 1-D coordinate")]
    if ds[name].dtype.kind != kind:
        expected = "datetime64" if kind == "M" else "timedelta64"
        return [(None, f"{name} must decode to {expected}, not {ds[name].dtype}")]
    return []


def check_cams(ds):
    errors = check_times(ds, "valid_time")
    errors += check_axis(ds, "latitude")
    errors += check_axis(ds, "longitude", ascending=True, regular=True)
    if not errors:
        lon = ds.longitude.values
        if lon[0] < 0 or lon[-1] >= 360:
            errors.append((None, "longitude must lie in 0..360"))
        elif not np.isclose(lon[0] + 360 - lon[-1], lon[1] - lon[0], rtol=1e-3):
            errors.append((None, "longitude must cover the globe, so that it wraps around"))

    # Optional: the forecast each valid time comes from, recorded in the provenance
    forecast = ("forecast_reference_time", "forecast_period")
    if any(name in ds.coords for name in forecast):
        if not all(name in ds.coords and ds[name].dims == ("valid_time",) for name in forecast):
            errors.append((None, "forecast_reference_time and forecast_period must both be "
                                 "coordinates along valid_time, or both be absent"))
        else:
            errors += check_times(ds, "forecast_reference_time")
            # The forecast step is a timedelta, or a number of hours
            if ds.forecast_period.dtype.kind not in "fiu":
                errors += check_times(ds, "forecast_period", kind="m")

    dims = ("valid_time", "latitude", "longitude")
    errors += check_variable(ds, "aod550", dims, DIMENSIONLESS)
    errors += check_variable(ds, "gtco3", dims, KG_M2)
    errors += check_variable(ds, "tcwv", dims, KG_M2)
    return errors


def check_dini(ds):
    errors = check_times(ds, "forecast_reference_time")
    # The forecast step is a timedelta, or a number of hours
    if not ("step" in ds.coords and ds.step.dtype.kind in "fiu"):
        errors += check_times(ds, "step", kind="m")
    errors += check_axis(ds, "x", ascending=True)
    errors += check_axis(ds, "y", ascending=True)

    for name in ("latitude", "longitude"):
        if name not in ds or ds[name].dims != ("y", "x"):
            errors.append((None, f"{name} must be a 2-D variable on (y, x)"))

    if "dini_projection" not in ds or "crs_wkt" not in ds.dini_projection.attrs:
        errors.append((None, "dini_projection with a crs_wkt attribute is missing"))
    else:
        try:
            CRS.from_wkt(ds.dini_projection.attrs["crs_wkt"])
        except Exception as e:
            errors.append((None, f"dini_projection crs_wkt is not a valid CRS ({e})"))

    errors += check_variable(ds, "wvint", ("forecast_reference_time", "step", "y", "x"), KG_M2)
    return errors


def check_lsasaf(ds):
    errors = check_times(ds, "time")
    if not errors:
        t = ds.time.values
        if not np.all(t == t.astype("datetime64[D]")):
            errors.append((None, "time must be at 00:00 of each product day"))

    # Scan angles in radians: the full disc spans about +-0.16
    for name in ("x", "y"):
        axis_errors = check_axis(ds, name, regular=True)
        if not axis_errors and np.abs(ds[name].values).max() > 0.2:
            axis_errors = [(None, f"{name} must be the scan angle in radians")]
        errors += axis_errors

    if "geostationary" not in ds:
        errors.append((None, "geostationary grid mapping is missing"))
    else:
        attrs = ds.geostationary.attrs
        if "perspective_point_height" not in attrs:
            errors.append((None, "geostationary has no perspective_point_height"))
        try:
            CRS.from_cf(attrs)
        except Exception as e:
            errors.append((None, f"geostationary is not a valid CF grid mapping ({e})"))

    dims = ("time", "y", "x")
    for var in ("AL-VI-DH", "AL-NI-DH"):
        var_errors = check_variable(ds, var, dims)
        if not var_errors and ds[var].dtype.kind != "f":
            var_errors = [(var, f"{var} must decode to albedo in 0-1 (scale_factor), not {ds[var].dtype}")]
        errors += var_errors
    for var in ("QFLAGS", "Z_Age"):
        var_errors = check_variable(ds, var, dims)
        if not var_errors and ds[var].dtype.kind not in "iu":
            var_errors = [(var, f"{var} must be an integer, not {ds[var].dtype}")]
        errors += var_errors

    # Every one of them is needed for the albedo
    return [(None, message) for _, message in errors]


CHECKS = {"cams": check_cams, "dini": check_dini, "lsasaf": check_lsasaf}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for source in CHECKS:
        parser.add_argument(f"--{source}", default="")
    args = parser.parse_args()

    failed = False
    for source, check in CHECKS.items():
        path = getattr(args, source)
        if not path:
            continue
        try:
            errors = check(xr.open_zarr(path))
        except Exception as e:
            errors = [(None, f"cannot be opened ({type(e).__name__}: {e})")]
        print(f"{source} {path}: {'OK' if not errors else 'does not match the specification'}")
        for _, message in errors:
            print(f"    {message}")
        failed |= bool(errors)

    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
