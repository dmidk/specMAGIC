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

## Input data

This code is assuming the presence of input satellite data. Currently, only Meteosat Third Generation (MTG) is supported, but support for other satellites is forthcoming. The user must provide a file path to available satellite data. 

## To run the code in demo mode

```bash
git clone [url]
cd specMAGIC
uv sync
```
will set up the python environment and necessary plugins. Note that the cloning step can take a minute or two, as the test data available for the demo are large.

The driver script, which pre-processes the MTG data, calls the SpecMAGIC code, and does some small postprocessing, is called `mtg.sh`. It can run in demo mode, using the data that comes pre-packaged in the `test_data` directory. These test data originate from EUMETSAT Meteosat products published in 2026. The single `.zip` file corresponds to a single level 1c satellite image taken by MTG/FCI on the 23rd August 2026 at 13:00 UTC. These exist under a CC-BY-4.0 free license; for details, please [the EUMETSAT website.](https://user.eumetsat.int/catalogue/EO:EUM:DAT:0989)

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

## Choosing a geographic region 

The region that specMAGIC produces output for is defined in `magic-config.asc`. The user can specify the starting latitude and longitude, and the number of points that the simulated grid should extend from these points. The default is the MTG-view of the Northern Hemisphere, spanning from 0°N, 60°W to 65.025°N, 60.025°E. The spatial resolution of the output is also specified here as 0.025°. At present, altering the spatial resolution may result in unexpected behaviour.

It is possible to specify a smaller region than seen by the satellite, which may be beneficial for performance reasons. Attempting to specify a simulated area outside of the satellite's field of view will raise a runtime error. 

The post-processing routine to plot maps of the generated parameters plots a region over Denmark by default. It is also possible to plot the entire Northern Hemisphere region seen by MTG by passing `--whole-domain`.

## To run the code with other MTG data 

As above, 
```bash
git clone [url]
uv sync
```
will set up the python environment and necessary plugins. 

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

## Configuring SpecMAGIC

The SpecMAGIC config is found in `root/magic-config.asc`. Here also is a short overview of the available lookup tables and input data. Alterations here are at the user's own risk.

## Multiple timesteps and animating

The scripts `scripts/multi.sh` and `scripts/animate.sh` can be used to run the code for multiple satellite images and create an animation of the result. The user must provide their own MTG 
data for this functionality, as the testing data includes only one image.

## MODIS BRDF albedo maps

By default, specMAGIC uses some older landmaps whose resolution is coarse. It is also possible to use MODIS hi-res land albedo maps for the clearsky calculation. 
These are too large to be stored here, but the user can toggle whether to download them using the `--albedo` keyword.

```
./scripts/mtg.sh --albedo MODIS
``` 

The total file size is ~ 6 GB. The download may take some minutes.
A destination check is implemented, so that if the files are found to already exist in the correct directory, they will not be overwritten. 

## Input fields for the time of the image

Instead of the monthly climatologies and the land-use albedo, the aerosol, ozone, water vapour and surface albedo can be given for the time of the satellite image, for example from CAMS, DINI and LSA SAF. Each is optional, and pixels without a value use the climatology. These data are not supplied with specMAGIC; you must provide them yourself.

From the command line, each is a NetCDF file already on the output grid:

```
./scripts/mtg.sh --ref-time "2026-08-23 13:00" --aerosol-file aerosol.nc --ozone-file ozone.nc \
    --water-vapour-file water_vapour.nc --albedo-file albedo.nc
```

From Python, the `specmagic` package (built by `uv sync`) takes the image and the inputs as xarray DataArrays on their own grids, picks the time and regrids them, and runs the same C++ code:

```python
import xarray as xr
import specmagic

cams = xr.open_zarr("cams.zarr")
image = xr.open_zarr("vis.zarr").vis_06.sel(start_time="2025-01-01T12:00")
out = specmagic.run(image, aerosol=cams.aod550, ozone=cams.gtco3, grid=(31, 72, -19, 34, 0.025))
```

A given input that cannot be used falls back to the climatology with a warning, or stops the run with `--strict-inputs` (`strict=True` in Python). See [the input documentation](docs/inputs.md#input-fields-for-the-time-of-the-image) for the file format, the Python API, and how the inputs are used.

## Further documentation

There exist some notes on various aspects of the program, linked below.

:sunny: [Physics](https://github.com/dmidk/specMAGIC/blob/main/docs/physics.md)
:cloud: [Including other satellites](https://github.com/dmidk/specMAGIC/blob/main/docs/other_satellites.md)
:zap: [OpenMP for parallelism](https://github.com/dmidk/specMAGIC/blob/main/docs/omp.md)
:key: [Input data](https://github.com/dmidk/specMAGIC/blob/main/docs/inputs.md)
:artificial_satellite: [Meteosat Third Generation](https://github.com/dmidk/specMAGIC/blob/main/docs/mtg.md)
