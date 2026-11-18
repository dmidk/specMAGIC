# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added 
- Added the option to download and run with MODIS high-resolution land albedo maps, instead of the coarse landmaps available by default. [!1](https://github.com/dmidk/specMAGIC/pull/1) @irenelivia [!12](https://github.com/dmidk/specMAGIC/pull/12) @SimonKamuk
- Added the satellite `.toml` file to avoid hardcoded satellite-specific values. [!18](https://github.com/dmidk/specMAGIC/pull/18) @ElloiseFangelLloyd @KristianHMoller
- Added optional input fields for the time of the satellite image (aerosol, ozone, water vapour and a noon black-sky albedo such as LSA SAF's), replacing the climatologies and the land-use albedo where they have values: as on-grid NetCDF files from the CLI (`--aerosol-file`, `--ozone-file`, `--water-vapour-file`, `--albedo-file`, and `--strict-inputs`), or as xarray DataArrays on any grid from the new Python package `specmagic` (`specmagic.run`), which runs the same C++ through pybind11. [!34](https://github.com/dmidk/specMAGIC/pull/34) @SimonKamuk
- The MODIS albedo can now be the fallback under an albedo input: with `--albedo-fallback MODIS-LANDMAP` and `--albedo-file` (or `albedo_fallback="MODIS-LANDMAP"` and `albedo_vis`/`albedo_nir` in Python), a pixel uses the albedo input, then the MODIS albedo, then the land-use albedo. Before, the two could not be combined. @SimonKamuk
- Added series of images to the Python package: `specmagic.run` takes an image with one more (time) dimension, runs each image as if alone, and returns the output with a `time` dimension. What the images share is done once: the lookup tables and climatologies are read once per month, the output grid is located in each input's grid once, and a slice of an input is regridded once for all images that use it. @SimonKamuk

### Changed
- Throw error message rather than silent failure in absorber interpolation. [!9](https://github.com/dmidk/specMAGIC/pull/9) @SimonKamuk
- Get the climatologies' grid spacing directly from axis rather than assume the cell-centred grid assumption. [!10](https://github.com/dmidk/specMAGIC/pull/10) @SimonKamuk
- Reworked `geo2MTGImage` function to work for other satellites as well as MTG. [!18](https://github.com/dmidk/specMAGIC/pull/18) @ElloiseFangelLloyd. 
- Use longitude of the subsatellite point in the `geo2Image` function to correctly calculate for satellites located at nonzero longitudes. [!20](https://github.com/dmidk/specMAGIC/pull/20) @KristianHMoller
- Removed the automatic vertical flipping applied to image, as was necessary for MTG. Instead, whether or not each image should be flipped can now be specified in the `satellite.toml` file, and the flipping will be applied conditionally. [!22](https://github.com/dmidk/specMAGIC/pull/22) @ElloiseFangelLloyd @KristianHMoller
- Renamed the driver script's `--albedo LANDMAP|MODIS` to `--albedo-fallback LANDMAP|MODIS-LANDMAP` (the environment variable `ALBEDO` of `magic` to `ALBEDO_FALLBACK`), and `specmagic.run`'s `albedo` to `albedo_fallback`, since it is the albedo where there is no albedo input, and MODIS falls back to the land-use albedo. @SimonKamuk
- The LSA SAF albedo's age limit (`max_age_days` of `specmagic.lsasaf_albedo`) now counts from the image's day: when an earlier day is used because the image's own is missing, the days in between are added to its age. @SimonKamuk
- `specmagic` reads only the forecast run and step of the slice it uses, not the whole `forecast_reference_time` and `forecast_period` coordinates along CAMS's `valid_time`, which saves about a second per CAMS field and image. @SimonKamuk
- `specmagic` takes an image's time from its `start_time` or `time` coordinate before its `start_time` attribute, which a slice of an archive keeps from the whole archive. @SimonKamuk
- The README leads with the Python package, and the input documentation describes it before the command line. @SimonKamuk


## [v1.0.0]

Initial release of specMAGIC. The application calculates gloval horizontal irradiance (GHI) from input MTG-FCI satellite data. OpenMP is included for parallelism. A test dataset is provided alongside the requisite pre-processing. Post-processing utilities make simple plots of the available outputs: GHI, Direct Normal Irradiance (the beam component of GHI), Clear Sky Radiance, and Effective Cloud Albedo.