# specMAGIC

## Attributions

The SpecMAGIC code has been developed by Richard Müller, Uwe Pfeifroth and Elloise Fangel-Lloyd. Input data and climatologies are provided by Richard Müller.

Annette Hammer and Axel Kemper contributed during their time at the University of Oldenburg to the development of the concepts of the methods and code. Tanja Behrendt contributed to the evaluation of SPECMAGIC and supported the project with her knowledge about the physics of spectral resolved irradiance.

The 2012 version of SPECMAGIC is described in
Mueller, R.; Behrendt, T.; Hammer, A.; Kemper, A. [A New Algorithm for the Satellite-Based Retrieval of Solar Surface Irradiance in Spectral Bands](https://doi.org/10.3390/rs4030622).

## Requirements

- A Linux installation 
- CMake >= 3.16
- C++ >= 17
- GCC >= 9.0 
- A uv installation, for dependency management. [Please see the uv docs for more information.](https://docs.astral.sh/uv/)
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```
- A Git Large File System (LFS) installation. Without it, the test data `.zip` will not be downloaded correctly and the demo will return a `BadZipFile` error. Install Git LFS for your platform:

```bash
sudo apt install git-lfs
```

Then initialize it and pull the test data:
```bash
git lfs install
git lfs pull
```

## Installation

```bash
git clone [url]
cd specMAGIC
uv sync
```
sets up the Python environment and builds the Python package `specmagic`, including its C++ code. Note that the cloning step can take a minute or two, as the test data available for the demo are large. The command-line program is built by the driver script `mtg.sh` when it runs.

Currently, only Meteosat Third Generation (MTG) is supported, but support for other satellites is forthcoming. The satellite data are not supplied, apart from one image for the demo.

## Using specMAGIC from Python

The Python package `specmagic` is the easiest way to use specMAGIC. It takes the satellite image, and optionally input fields for its time, as xarray DataArrays from anything xarray can open (a satpy Scene, zarr archives, NetCDF files), and returns GHI, DNI, CAL (cloud index) and CSR (clear-sky radiance) as an xarray Dataset. It runs the same C++ code as the command line, in memory: nothing is written to disk.

### One image

The image can come straight from the MTG level 1c files with satpy:

```python
import glob
import hdf5plugin
import specmagic
from satpy import Scene

scene = Scene(filenames=glob.glob("path/to/*FCI-1C-RRAD-FDHSI-FD--CHK-BODY*_O_0079_*.nc"),
              reader="fci_l1c_nc", reader_kwargs={"engine": "h5netcdf"})
scene.load(["vis_06"], calibration="counts")

out = specmagic.run(scene["vis_06"], grid=(52, 59, 6, 17, 0.025))   # lat_min, lat_max, lon_min, lon_max, resolution
out.GHI.plot()
out.to_netcdf("magic.nc")
```

Without `grid`, the output grid is that of `magic-config.asc` (see [choosing a geographic region](#choosing-a-geographic-region)).

### Input fields for the time of the image

Instead of the monthly climatologies and the land-use albedo, the aerosol, ozone, water vapour and surface albedo can be given for the time of the image, for example from CAMS, DINI and LSA SAF. Each is optional, can be on its own grid and have a time dimension: the slice for the image is picked and regridded for you. Pixels without a value use the climatology. These data are not supplied with specMAGIC; you must provide them yourself.

```python
import xarray as xr

image = xr.open_zarr("vis.zarr").vis_06.sel(start_time="2025-03-13T12:00")
cams = xr.open_zarr("cams.zarr")
dini = xr.open_zarr("dini.zarr", decode_coords="all")       # projected grids need decode_coords="all"
lsasaf = xr.open_zarr("lsasaf.zarr", decode_coords="all")

out = specmagic.run(image, aerosol=cams.aod550, ozone=cams.gtco3, water_vapour=dini.wvint,
                    **specmagic.lsasaf_albedo(lsasaf), grid=(52, 59, 6, 17, 0.025))
```

A given input that cannot be used falls back to the climatology with a warning, or raises with `strict=True`. The attributes of the output record where each input came from and how much of it was used.

### A series of images

An image with one more dimension, such as a day of images from an archive, is a series:

```python
images = xr.open_zarr("vis.zarr").vis_06.sel(start_time="2025-03-13")
out = specmagic.run(images, aerosol=cams.aod550, ozone=cams.gtco3, water_vapour=dini.wvint,
                    **specmagic.lsasaf_albedo(lsasaf), grid=(52, 59, 6, 17, 0.025))
```

Each image is run as if alone, and the output has a `time` dimension. What the images share is done once for the whole series (reading the tables, finding the output grid in each source's grid, regridding a slice of a source that several images use), so a series is several times faster than running its images one by one.

See [the input documentation](docs/inputs.md#from-python) for the image formats, the input fields, how they are picked and regridded, and the output attributes.

## Using specMAGIC from the command line

The driver script `scripts/mtg.sh` pre-processes MTG level 1c files, builds and runs the specMAGIC program on them, and plots the result. It writes GHI, DNI, CAL and CSR as NetCDF files under `out/`.

### To run the code in demo mode

`mtg.sh` can run in demo mode, using the data that comes pre-packaged in the `test_data` directory. These test data originate from EUMETSAT Meteosat products published in 2026. The single `.zip` file corresponds to a single level 1c satellite image taken by MTG/FCI on the 23rd August 2026 at 13:00 UTC. These exist under a CC-BY-4.0 free license; for details, please [the EUMETSAT website.](https://user.eumetsat.int/catalogue/EO:EUM:DAT:0989)

```bash
./scripts/mtg.sh --demo
```

The default MTG channel is `vis_06`, which is 640 nm. Other channels can be specified using `--channel`, e.g. 

```bash
./scripts/mtg.sh --demo --channel "nir_13"
```

To see a full list of input options, 

```bash
./scripts/mtg.sh -h
```
Figures will appear in the `/figs` directory. By default, plots of GHI, DNI, CAL (cloud index) and CSR (clear-sky radiance) are produced. 

### To run the code with other MTG data

For the user to run specMAGIC for dates/times other than that supplied in the demo, the user must supply their own data as well as a date and time of interest. The program will then search for `*.nc` files matching the basic EUMETSAT file name convention using that specified date and time. 

First, specify the path to the level 1c MTG satellite data 

```bash
export MTG_RAW_DIR='path/to/data/archive'
```

Then call the code with a timestamp specified

```bash
./scripts/mtg.sh --ref-time "Yesterday 08:00"
```

The timestamp is in 24-hour hour time (so use "Yesterday 13:00" rather than "Yesterday 1pm"). Explicit calendar dates can all be passed. All of the below are accepted datetime formats. 
```bash
2026-08-31 08:00
2026/08/31 08:00
"Today 08:00"        # assuming today is the 31st August
```

Remember that all times are in UTC, to match EUMETSAT satellite conventions. Passing `08:00` will search for data acquired at `08:00` as clocked by the satellite.

### Choosing a geographic region

The region that specMAGIC produces output for is defined in `magic-config.asc` (from Python, `grid=` overrides it). The user can specify the starting latitude and longitude, and the number of points that the simulated grid should extend from these points. The default is the MTG-view of the Northern Hemisphere, spanning from 0°N, 60°W to 65.025°N, 60.025°E. The spatial resolution of the output is also specified here as 0.025°. At present, altering the spatial resolution may result in unexpected behaviour.

It is possible to specify a smaller region than seen by the satellite, which may be beneficial for performance reasons. Attempting to specify a simulated area outside of the satellite's field of view will raise a runtime error. 

The post-processing routine to plot maps of the generated parameters plots a region over Denmark by default. It is also possible to plot the entire Northern Hemisphere region seen by MTG by passing `--whole-domain`.

### Configuring SpecMAGIC

The SpecMAGIC config is found in `root/magic-config.asc`. Here also is a short overview of the available lookup tables and input data. Alterations here are at the user's own risk.

### Multiple timesteps and animating

From Python, `specmagic.run` takes a series of images directly (see above). On the command line, the scripts `scripts/multi.sh` and `scripts/animate.sh` can be used to run the code for multiple satellite images and create an animation of the result. The user must provide their own MTG 
data for this functionality, as the testing data includes only one image.

### MODIS BRDF albedo maps

By default, specMAGIC uses some older landmaps whose resolution is coarse. It is also possible to use MODIS hi-res land albedo maps for the clearsky calculation. 
These are too large to be stored here, but the user can toggle whether to download them using the `--albedo` keyword.

```
./scripts/mtg.sh --albedo MODIS
```

The total file size is ~ 6 GB. The download may take some minutes.
A destination check is implemented, so that if the files are found to already exist in the correct directory, they will not be overwritten. 

From Python, pass `albedo="MODIS"` to `specmagic.run` once the maps are downloaded (by the command above, or `uv run python py_utils/download_MODIS_maps.py climatologies/modis-brdf`).

### Input fields for the time of the image

From the command line, each input is a NetCDF file already on the output grid; the command line does not regrid them or pick their time (`specmagic.write_inputs` can make such files):

```
./scripts/mtg.sh --ref-time "2026-08-23 13:00" --aerosol-file aerosol.nc --ozone-file ozone.nc \
    --water-vapour-file water_vapour.nc --albedo-file albedo.nc
```

A given input that cannot be used falls back to the climatology with a warning, or stops the run with `--strict-inputs`. See [the input documentation](docs/inputs.md#from-the-cli) for the file format.

## Tests

The tests use pytest, which `uv sync` installs with the dev dependencies. The tests comparing the CLI with Python need the CLI built first, and are skipped otherwise:

```bash
cmake -S . -B build && cmake --build build
uv run pytest
```

## Further documentation

There exist some notes on various aspects of the program, linked below.

:sunny: [Physics](https://github.com/dmidk/specMAGIC/blob/main/docs/physics.md)
:cloud: [Including other satellites](https://github.com/dmidk/specMAGIC/blob/main/docs/other_satellites.md)
:zap: [OpenMP for parallelism](https://github.com/dmidk/specMAGIC/blob/main/docs/omp.md)
:key: [Input data](https://github.com/dmidk/specMAGIC/blob/main/docs/inputs.md)
:artificial_satellite: [Meteosat Third Generation](https://github.com/dmidk/specMAGIC/blob/main/docs/mtg.md)
