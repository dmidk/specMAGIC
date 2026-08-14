import os
os.environ["DASK_NUM_WORKERS"] = "1"
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import re
from pyproj import CRS, Transformer
import warnings
from contextlib import contextmanager

@contextmanager
def suppress_known_warnings():

    ## Get rid of known but harmless warnings during preprocessing
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"dtype uint16 not compatible with CF-1\.7\.",
            category=UserWarning,
            module=r"satpy\.writers\.cf_writer",
        )

        warnings.filterwarnings(
            "ignore",
            message=r"invalid value encountered in cast",
            category=RuntimeWarning,
            module=r"numpy\..*",
        )

        yield


## Checking for missing data

## Pulls geographic chunk out of the long filename
def get_chunk_number(filename):

    base = os.path.basename(filename)

    match = re.search(r"_(\d{4})_(\d{4})\.nc$", base)

    if match is None:
        raise ValueError(
            f"[PYTHON] Could not extract MTG chunk number from filename: {base}"
        )

    return int(match.group(2))


def check_required_chunks(filenames):

    # Run the actual checks for missing data 
    # Missing chunk 38 gives immediate failure, this is directly over DK 
    # Chunks 37 and 39 give warnings.

    found_chunks = {get_chunk_number(f) for f in filenames}

    if 38 not in found_chunks:
        print(
            "[PYTHON] FATAL: Required MTG chunk 0038 is missing. "
            "Data over Denmark is not available!",
            file=sys.stderr,
        )
        sys.exit(1)

    # Surrounding chunks are useful but not fatal
    for chunk in [37, 39]:
        if chunk not in found_chunks:
            msg = (
                f"[PYTHON] MTG chunk {chunk:04d} is missing. "
                "We are missing data surrounding Denmark."
            )
            warnings.warn(msg, RuntimeWarning)
            print("[PYTHON] WARNING:", msg)


### -------------- PLOTTERS 

def plot_hemisphere(da, figs_folder):
     
    # Very simple plot that just shows if data is missing
    # Plots entire northern hemisphere

    # fill what may be missing
    v = da.values
    fill = da.attrs.get("_FillValue", None)
    if fill is not None:
        v = np.where(v == fill, np.nan, v) 

    plt.figure()
    plt.imshow(v, cmap='viridis', origin="lower")
    plt.title("Seen by MTG")
    plt.colorbar(label="counts")
    plt.tight_layout()
    plt.savefig(figs_folder+"/hemisphere.png", dpi=300)



## Helper for checking
def find_coord_name(possible_names, ds):
    """"Fetches coordinate names frm a ncdf file"""
    for name in possible_names:
        if name in ds.coords:
            return name
        if name in ds.variables:
            return name
    return None

## Helper to make plot subsets
def slice_by_bounds_1d(coord_1d, vmin, vmax):
    """Slicer with check for which way the data runs"""
    cmin = float(coord_1d.min())
    cmax = float(coord_1d.max())
    ascending = cmax > cmin
    return slice(vmin, vmax) if ascending else slice(vmax, vmin)


def centers_to_edges(c):
    """Convert 1D coordinate centres to edges for pcolormesh with shading='flat'."""
    """For making plots look nice"""
    c = np.asarray(c)
    dc = np.diff(c) / 2.0
    edges = np.empty(c.size + 1, dtype=c.dtype)
    edges[1:-1] = c[:-1] + dc
    edges[0] = c[0] - dc[0]
    edges[-1] = c[-1] + dc[-1]
    return edges


def extent_lonlat_to_geos_bounds(tf: Transformer, lon_min, lon_max, lat_min, lat_max, n=200):
    """For the projection
    Compute a safe x/y box in projected space by sampling the boundary
    of lat/lon extent. `tf` transforms lon/lat to x/y (meters).
    """
    lons_top = np.linspace(lon_min, lon_max, n)
    lats_top = np.full(n, lat_max)

    lons_bot = np.linspace(lon_min, lon_max, n)
    lats_bot = np.full(n, lat_min)

    lats_left = np.linspace(lat_min, lat_max, n)
    lons_left = np.full(n, lon_min)

    lats_right = np.linspace(lat_min, lat_max, n)
    lons_right = np.full(n, lon_max)

    lons = np.concatenate([lons_top, lons_bot, lons_left, lons_right])
    lats = np.concatenate([lats_top, lats_bot, lats_left, lats_right])

    xs, ys = tf.transform(lons, lats)

    # drop any non-finite 
    xs = np.asarray(xs)
    ys = np.asarray(ys)
    ok = np.isfinite(xs) & np.isfinite(ys)
    if not np.any(ok):
        raise ValueError("No finite projected points found for extent boundary.")

    x0, x1 = float(xs[ok].min()), float(xs[ok].max())
    y0, y1 = float(ys[ok].min()), float(ys[ok].max())
    return x0, x1, y0, y1

## Plots the MTG data that has been processed into one file
def plot_handled_MTG(data: xr.Dataset, channel, out_dir: str, extent):
    """
    Plot MTG using native geostationary x/y coordinates 
    extent = [lat_min, lon_min, lat_max, lon_max, stride]
    """

    lat_min, lon_min, lat_max, lon_max, stride = extent

    # In case there are multiple channels, we only care about one
    da = data[channel]

    # Read gridmapp variable directly from the file
    gm_name = da.attrs.get("grid_mapping")
    if gm_name is None or gm_name not in data:
        raise ValueError("No grid_mapping variable found (expected da.attrs['grid_mapping']).")

    # dumb but works
    gm = data[gm_name]
    gm_attrs = gm.attrs

    # Build a pyproj CRS 
    if "crs_wkt" in gm_attrs:
        geos_proj = CRS.from_wkt(gm_attrs["crs_wkt"])
    else:
        # fallback: build CRS from CF attrs 
        geos_proj = CRS.from_cf(gm_attrs)

    # Transformer: lon/lat (EPSG:4326) -> geostationary projected x/y (meters)
    tf = Transformer.from_crs("EPSG:4326", geos_proj, always_xy=True)

    x0, x1, y0, y1 = extent_lonlat_to_geos_bounds(tf, lon_min, lon_max, lat_min, lat_max, n=200)

    # some padding
    dx = float(np.median(np.diff(data["x"].values)))
    dy = float(np.median(np.diff(data["y"].values)))
    pad = 5  # num pixels for padding

    x0 -= pad * dx
    x1 += pad * dx
    y0 -= pad * dy
    y1 += pad * dy

    # Subset in x/y
    # The preprocessing does not keep lats/lons, for speed
    x = data["x"]
    y = data["y"]
    x_slice = slice_by_bounds_1d(x, x0, x1)
    y_slice = slice_by_bounds_1d(y, y0, y1)

    sub = da.sel(x=x_slice, y=y_slice)

    # Optional striding for speed
    if stride and stride > 1:
        sub = sub.isel(x=slice(None, None, stride), y=slice(None, None, stride))

    # Build Cartopy CRS for the data coordinates (x/y)
    # Use parameters from attrs to construct Geostationary CRS.
    # probably overkill?
    geos_crs = ccrs.Geostationary(
        central_longitude=float(gm_attrs.get("longitude_of_projection_origin", 0.0)),
        satellite_height=float(gm_attrs["perspective_point_height"]),
        sweep_axis=str(gm_attrs.get("sweep_angle_axis", "y")),
        false_easting=float(gm_attrs.get("false_easting", 0.0)),
        false_northing=float(gm_attrs.get("false_northing", 0.0)),
        globe=ccrs.Globe(
            semimajor_axis=float(gm_attrs["semi_major_axis"]),
            semiminor_axis=float(gm_attrs["semi_minor_axis"]),
        ),
    )

    # Plot
    ax_proj = ccrs.Mercator()
    data_crs_latlon = ccrs.PlateCarree()

    fig = plt.figure(figsize=(10, 7))
    ax = plt.axes(projection=ax_proj)

    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=data_crs_latlon)

    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor="#3a7d44", alpha=0.25)
    ax.add_feature(cfeature.OCEAN, facecolor="#1f4e79", alpha=0.35)  # deep-ish blue

    ax.gridlines(draw_labels=True, linewidth=0.3, color="gray", alpha=0.6)

    # Just to make it extra pretty
    x_edges = centers_to_edges(sub["x"].values)
    y_edges = centers_to_edges(sub["y"].values)
    vmin = float(sub.quantile(0.02))
    vmax = float(sub.quantile(0.98))

    mesh = ax.pcolormesh(x_edges, y_edges, sub.values, transform=geos_crs,
        shading="flat", cmap="Greys_r", alpha=0.7, vmin=vmin, vmax=vmax)

    cb = fig.colorbar(mesh, ax=ax, shrink=0.8, pad=0.08, fraction=0.05)
    cb.set_label(f"Seen by MTG ({da.attrs.get('units', 'unknown units')})")

    plt.title(f"Raw MTG data (channel: {channel})")

    outname = "MTG_raw.png"
    plt.savefig(f"{out_dir}/{outname}", dpi=300, bbox_inches="tight", pad_inches=0.1)


## Searches through data files and plots the most recent one from timestamp
## For use on the output data produced by specmagic now
def plot_latest_product(product, data_dir, figs_dir, extent):
    """
    product: 'GHI', 'DNI', 'CAL' or 'CSR'
    """
    lat_min, lon_min, lat_max, lon_max, stride = extent

    ## Find newest file 
    candidates = sorted(data_dir.glob(f"{product}*.nc"), key=lambda p: p.stat().st_mtime)
    if not candidates:
        print(f"[WARN] No files matching '{product}*.nc' found in {data_dir}")
        return

    latest_file = candidates[-1]
    print(f"[{product}] Plotting latest file: {latest_file.name}")

    ## Parse timestamp from filename like GHIhrYYYYMMDDHHMM...
    m = re.search(rf"{product}(\d{{12}})", latest_file.name)
    if not m:
        print(f"[WARN] Could not parse timestamp from filename: {latest_file.name}")
        stamp = None
        dt_str = "unknown time"
    else:
        stamp = m.group(1)
        dt = datetime.strptime(stamp, "%Y%m%d%H%M")
        dt_str = dt.strftime("%d/%m/%y, %H:%M")

    ds = xr.open_dataset(latest_file)

    ## Pick variable
    data_vars = list(ds.data_vars)
    if not data_vars:
        raise ValueError(f"[{product}] Dataset contains no data variables.")

    preferred = [v for v in data_vars if product.lower() in v.lower()]
    varname = preferred[0] if preferred else data_vars[0]
    da = ds[varname]


    ## If there is a time dimension, choose the last time step
    ## We do NOT anticipate a time dimension
    ## This will give a printout if it finds one
    for tdim in ("time", "valid_time", "datetime"):
        if tdim in da.dims:
            da = da.isel({tdim: -1})
            break

    ## Get coords
    lat_name = find_coord_name(["lat", "latitude", "y"], ds)
    lon_name = find_coord_name(["lon", "longitude", "x"], ds)
    if lat_name is None or lon_name is None:
        raise ValueError(
            f"[{product}] Could not find latitude/longitude coords. "
            f"Inspect ds.coords and adjust names."
        )

    lat = ds[lat_name]
    lon = ds[lon_name]

    ## Convert longitude to -180...180
    if np.issubdtype(lon.dtype, np.number):
        lon_vals = lon.values
        if np.nanmax(lon_vals) > 180:
            lon_wrapped = (((lon + 180) % 360) - 180)
            ds = ds.assign_coords({lon_name: lon_wrapped}).sortby(lon_name)
            da = ds[varname]
            for tdim in ("time", "valid_time", "datetime"):
                if tdim in da.dims:
                    da = da.isel({tdim: -1})
                    break
            lon = ds[lon_name]
            lat = ds[lat_name]

    lat_slice = slice_by_bounds_1d(lat, lat_min, lat_max)
    lon_slice = slice_by_bounds_1d(lon, lon_min, lon_max)

    subset = da.sel({lat_name: lat_slice, lon_name: lon_slice})
    if stride and stride > 1:
        subset = subset.isel({
            lat_name: slice(None, None, stride),
            lon_name: slice(None, None, stride)
        })

    ## Check the subsetting
    if subset.size == 0:
        print(f"[WARN] [{product}] Subset is empty for extent {extent}.")
        ds.close()
        return

    ## Fixed colour limits, for giffing
    PRODUCT_LIMITS = {
        "GHI": (0, 1000),   
        "DNI": (0, 1000),   
        "CAL": (None, None),
        "CSR": (0, 1000),      
    }
    vmin, vmax = PRODUCT_LIMITS.get(product, (None, None))

    COLOURMAPS = {
        "GHI": "bone_r",
        "DNI": "bone_r", 
        "CAL": "bone",
        "CSR": "hot",
    }
    colours = COLOURMAPS.get(product, (None))

    ## Actual plotting
    proj = ccrs.Mercator()
    data_crs = ccrs.PlateCarree()

    fig = plt.figure(figsize=(10, 7))

    ## Fixed axes positions: [left, bottom, width, height]
    ax = fig.add_axes([0.08, 0.10, 0.72, 0.80], projection=proj)
    cax = fig.add_axes([0.84, 0.15, 0.03, 0.70])
    
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=data_crs)

    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor="#3a7d44", alpha=0.25)
    ax.add_feature(cfeature.OCEAN, facecolor="#1f4e79", alpha=0.35)  # deep-ish blue

    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="gray", alpha=0.6)
    gl.top_labels = False
    gl.right_labels = False

    mesh = ax.pcolormesh(
        subset[lon_name], subset[lat_name], subset,
        transform=data_crs,
        shading="auto",
        cmap=colours,
        vmin=vmin,
        vmax=vmax
    )

    cb = fig.colorbar(mesh, cax=cax)
    cb.set_label(f"{varname} ({subset.attrs.get('units', 'unknown units')})")

    ax.set_title(f"{product} at {dt_str}")

    outname = f"{product.lower()}_{stamp if stamp else 'unknown'}.png"
    plt.savefig(figs_dir + "/" + outname, dpi=400)
    plt.close(fig)

    ds.close()

    ## END function plot latest product
