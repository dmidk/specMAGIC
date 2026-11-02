# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added 
- Added the option to download and run with MODIS high-resolution land albedo maps, instead of the coarse landmaps available by default. [!1](https://github.com/dmidk/specMAGIC/pull/1) @irenelivia [!12](https://github.com/dmidk/specMAGIC/pull/12) @SimonKamuk
- Added the satellite `.toml` file to avoid hardcoded satellite-specific values. [!18](https://github.com/dmidk/specMAGIC/pull/18) @ElloiseFangelLloyd @KristianHMoller
- Added the option `--zarr-atmosphere-inputs` to take aerosol and ozone from CAMS and water vapour from DINI (then CAMS), and the albedo option `--albedo LSASAF`, for the time of the satellite image, falling back to the climatologies wherever these are missing, or stopping with `--strict-inputs`. The archives are supplied by the user, and `py_utils/verify_zarr_inputs.py` checks their format. [!34](https://github.com/dmidk/specMAGIC/pull/34) @SimonKamuk

### Changed
- Throw error message rather than silent failure in absorber interpolation. [!9](https://github.com/dmidk/specMAGIC/pull/9) @SimonKamuk
- Get the climatologies' grid spacing directly from axis rather than assume the cell-centred grid assumption. [!10](https://github.com/dmidk/specMAGIC/pull/10) @SimonKamuk
- Reworked `geo2MTGImage` function to work for other satellites as well as MTG. [!18](https://github.com/dmidk/specMAGIC/pull/18) @ElloiseFangelLloyd. 
- Use longitude of the subsatellite point in the `geo2Image` function to correctly calculate for satellites located at nonzero longitudes. [!20](https://github.com/dmidk/specMAGIC/pull/20) @KristianHMoller
- Removed the automatic vertical flipping applied to image, as was necessary for MTG. Instead, whether or not each image should be flipped can now be specified in the `satellite.toml` file, and the flipping will be applied conditionally. [!22](https://github.com/dmidk/specMAGIC/pull/22) @ElloiseFangelLloyd @KristianHMoller


## [v1.0.0]

Initial release of specMAGIC. The application calculates gloval horizontal irradiance (GHI) from input MTG-FCI satellite data. OpenMP is included for parallelism. A test dataset is provided alongside the requisite pre-processing. Post-processing utilities make simple plots of the available outputs: GHI, Direct Normal Irradiance (the beam component of GHI), Clear Sky Radiance, and Effective Cloud Albedo.