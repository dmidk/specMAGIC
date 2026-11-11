# Inputs

SpecMAGIC needs four things: a **satellite image**, the **clear-sky atmospheric state**, a
**surface albedo**, and a **cloud spectral response**. Only the image is dynamic. Everything
else comes from the static tables described below, all of which are named in
[`magic-config.asc`](../magic-config.asc).

That config file is read **positionally** by `loadConfig` in `src/read.cpp`, so the order of
its value lines is significant: inserting or removing one shifts every subsequent field. The
tables are loaded in `src/main.cpp`, in the block marked `Read lookup tables`, and in the
`Tables::Climate` constructor.

## What the model needs

| Physical input | Files | Where it enters the physics |
|---|---|---|
| Aerosol | `macc-ecwmf-clim03-06rich.dat` + `magic-clear_spectral.lut` | Produces the clear-sky irradiance itself, per band, for both global and beam |
| Water vapour | `h2o_ecmwf.dat` + `wvcorr-l.lut` | Additive absorption loss, per band |
| Ozone | `ERAint-o3-clim.asc` + `o3corr-l.lut` | Additive absorption loss, per band |
| Surface albedo | *either* `IGBPa_2006.map` + `reflec.lut` *or* `modis-brdf/*.nc` | Multiplicative on clear-sky global, and sets the dark limit of the cloud index |
| Cloud spectral response | `lambdacor.lut` | Multiplicative on the all-sky band GHI |

The cloud amount itself is **not** a static input — it is derived per pixel from the
satellite radiance as the effective cloud albedo (CAL), then converted to a clear-sky index.
The only static table involved is the spectral correction in the last row.

## Clear-sky atmosphere

The three atmospheric constituents all follow the same two-part pattern:

1. A **climatology** giving the spatial distribution of the constituent, as long-term
   monthly means on a regular lat/lon grid. Read by `Tables::Climatology`, which retains
   only the month of the image being processed and bilinearly interpolates to each pixel.
2. A **lookup table** giving the radiative response as a function of that concentration,
   pre-computed by a radiative transfer model and read at 32 Kato bands.

The consequence of that split is that **the climatology's units must match the LUT's
concentration axis**, which is why both are given for each pair below. This matters when
substituting a climatology: a wrongly-scaled field is still a valid file.

### Aerosol

`climatologies/macc-ecwmf-clim03-06rich.dat` — MACC/ECMWF aerosol climatology for 2003-06,
smoothed, read into `Tables::Aerosol::aod`.

* Grid: 720 x 360, 0.5 degrees, cell-centred (-179.75 to 179.75).
* Format: header line, then `lon lat`, 12 monthly AOD values, then a further 12 columns
  which the code **skips** (`extra_cols_to_skip = 12`).
* Units: AOD is dimensionless, nominally at 550 nm.

The skipped columns are all `0.950` and appear to be a single-scattering albedo. The code
does not use them: single-scattering albedo and asymmetry parameter are hardcoded as
`ssa = 0.95` and `gg = 0.7` in `Tables::AerosolOptics`, on the grounds that the existing
climatology has no spatial variation in either.

`luts/magic-clear_spectral.lut` — the RTM output, read by `Tables::RTM`. This is the core
clear-sky table: unlike water vapour and ozone, which only apply corrections, this one
generates the clear-sky irradiance. See the DWD article for how it is derived.

* Shape: 2 `gg` x 3 `ssa` x 11 `aod` x 32 bands. Six blocks, each opened by a
  `ssa= <value> gg= <value>` line, each block 11 x 32 rows of seven columns:
  `aod lambda Im gtau ag btau ab`.
* Axes: `ssa` in {0.70, 0.85, 1.00}, `gg` in {0.60, 0.78}, `aod` in {0.00, 0.10, 0.20,
  0.30, 0.45, 0.60, 0.80, 1.00, 1.20, 1.50, 2.00}.
* Units: `aod` dimensionless; `lambda` in nm; `Im` is a fitted Lambert-Beer prefactor with
  dimensions of W m^-2 per band; `gtau` and `btau` are dimensionless optical-depth-like
  terms (negative) for the global and beam components; `ag` and `ab` are dimensionless
  exponents applied to the cosine of the solar zenith angle.

Interpolation is linear in `aod` and `ssa` but **nearest-neighbour in `gg`**. Since `gg` is
a hardcoded constant of 0.7, this always selects the 0.78 node.

### Water vapour

`climatologies/h2o_ecmwf.dat` — ERA-Interim total-column water vapour.

* Grid: 1441 x 721, 0.25 degrees, node-centred (-180 to 180 inclusive).
* Units: **kg m^-2**, equivalently mm of precipitable water. Observed range is 0.2 to 73.4,
  which matches the 0-70 axis of the LUT below.
* Tracked with Git LFS, so `git lfs pull` is required. Until then it is a ~130 byte pointer
  file and cannot be read.

`luts/wvcorr-l.lut` — water vapour absorption, read by `Tables::Correction::Absorber`.

* Shape: header line, then 18 x 32 records of six columns,
  `concen lambda delta_Ig ag delta_Ib ab`.
* Axis: 0, 2.5, 5.0 ... 15 then 20, 25 ... 70 kg m^-2. The finer sampling at low
  concentrations reflects the logarithmic saturation of water vapour absorption.
* Units: `lambda` in nm; `delta_Ig` and `delta_Ib` are **negative** irradiance losses in
  W m^-2 per band, global and beam, added to the clear-sky result; `ag` and `ab` are
  dimensionless exponents applied to the cosine of the solar zenith angle.

Note that the unused legacy `h2oclim.dat` is in **cm**, not kg m^-2, so it is not a drop-in
replacement — substituting it would put the LUT lookup out by a factor of ten.

### Ozone

`climatologies/ERAint-o3-clim.asc` — ERA-Interim total-column ozone.

* Grid: 361 x 181, 1 degree, node-centred (-180 to 180 inclusive). The config file records
  this as 360 x 181, but the recorded dimensions are unused: `Tables::Climatology::read`
  infers the grid by scanning the file.
* Units: **Dobson Units**. Observed range is 157 to 457.

`luts/o3corr-l.lut` — ozone absorption. Same format and column meanings as `wvcorr-l.lut`.

* Shape: header line, then 8 x 32 records.
* Axis: 210, 255, 300 ... 525 DU.

Note that the climatology reaches below the 210 DU floor of that axis. Every such cell is
Antarctic (1.7 percent of the file, all between 66 and 90 degrees south); the default output
grid starts at the equator, where the minimum is 234 DU.

## Surface albedo

Surface albedo is consumed in two distinct places, at two different spectral resolutions:

* Per Kato band, in `magic`, as a multiplicative correction `0.98 + 0.1 * albedo` on the
  clear-sky global irradiance.
* At the satellite channel's own wavelength, in `effectiveCloudAlbedo`, to set `Rmin` — the
  radiance a cloud-free pixel would show. This feeds the cloud index directly, so an albedo
  error here propagates into the retrieved cloud amount rather than just scaling the result.

There are two interchangeable sources, selected by `--albedo` on the driver script and
resolved by `Reflectivity::getBestKatoSurfaceAlbedo` and
`Reflectivity::getBestSatelliteSurfaceAlbedo`.

### Default: land use plus a spectral albedo LUT

`climatologies/IGBPa_2006.map` — IGBP land-use classification, read by `Tables::LandUse`.

* Grid: 2160 x 1080, 6 cells per degree, indexed by `Tables::LandUse::index`.
* Format: raw bytes, one `unsigned char` per cell, no header. Must be exactly 2332800 bytes
  or the read is rejected.
* Units: land class index, 1 to 20, dimensionless. The class table is in
  [`physics.md`](physics.md).

`luts/reflec.lut` — spectral albedo per land class, read by `Tables::GroundAlbedo`.

* Shape: header line, then 20 rows (one per land class) of 32 columns (one per band).
* Units: albedo, dimensionless, 0 to 1.

The two are used together: the land-use map gives a class, the class indexes a row of
`reflec.lut`. The result is then multiplied by a zenith-angle reflectance correction from
`Reflectivity::zenithCorrection`, which accounts for the non-Lambertian brightening of
surfaces at grazing angles. That correction's coefficients are **not** in any table — they
are hardcoded per land class in `src/reflectivity.cpp`, with a separate formula for ocean
(class 17) and everything else.

### Alternative: MODIS BRDF

`climatologies/modis-brdf/` replaces **both** files above, deriving albedo from measured BRDF
kernel weights rather than a class lookup. Only read when the driver is given
`--albedo MODIS`; otherwise `ModisBrdf::ModisBrdfAlbedo` stays empty. The ten NetCDF files
are fetched by `py_utils/download_MODIS_maps.py` and named
`Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_<source>.nc`, where `<source>` is `Band1`
to `Band7`, `vis`, `nir` or `shortwave`.

* Grid: 7200 x 3600, 0.05 degrees, by 12 months. Only the region covered by the output grid
  is read into memory.
* Variables: `FISO_MEAN`, `FVOL_MEAN`, `FGEO_MEAN`.
* Units: stored as `short`, scaled by 0.001, giving dimensionless BRDF kernel weights.
  32767 is the fill value.

The three weights are combined into a black-sky and a white-sky albedo and mixed by diffuse
fraction, giving an albedo in 0-1. Because this is an angular model rather than a class
lookup, the hardcoded `zenithCorrection` is **not** applied on this path — the BRDF kernels
already carry the angular dependence.

Which of the ten files is used depends on the consumer: `sourceForKatoBand` maps each Kato
band to a MODIS band or a broad vis/nir/shortwave product, while
`sourceForSatelliteWavelength` maps the satellite channel's wavelength. So the same pixel
can draw on different files for its clear-sky albedo and its cloud-index albedo.

### How the two sources combine

MODIS is a preference, not a commitment. A pixel whose BRDF parameters are fill-valued, out
of the loaded region, or non-finite yields no usable MODIS albedo, and that pixel falls back
to the land-use albedo for that band. Since MCD43C1 is a land product, this is the normal
case over water rather than an exceptional one, so a `--albedo MODIS` run is in practice a
mixture: MODIS over land, land-use classes over sea.

## Cloud spectral response

`luts/lambdacor.lut` — corrects the broadband cloud transmission to a wavelength-dependent
transmission, i.e. from the clear-sky to the all-sky case. Read by
`Tables::Correction::Spectral` and applied by `wavelengthCorrectionCloudySky` after the
per-band clear-sky irradiance has been scaled by the clear-sky index.

* Format: header line; then one row of six clear-sky index `k*` values corresponding to six
  cloud optical depths; then 32 rows of `lambda` plus six correction factors.
* Units: `k*` and the correction factors are dimensionless, the latter multiplicative on the
  per-band GHI; `lambda` in nm.
* Coverage: the file's own header notes the corrections are only calculated over
  317.3-2638.5 nm. Rows outside that range repeat the nearest computed value.
* The file has 34 non-empty lines, all of which are read; the remainder is trailing blanks.

The table is indexed by the clear-sky index `k`, which comes from the satellite-derived cloud
index, so this is the one static table whose lookup depends on the image rather than on
geography.

## Input fields for the time of the image

Four optional input fields replace static tables with values for the time of the satellite
image, for example from CAMS, DINI and LSA SAF:

| Input | Field in magic | Units | Replaces |
|---|---|---|---|
| Aerosol | `aod` | dimensionless, at 550 nm | `macc-ecwmf-clim03-06rich.dat` |
| Ozone | `o3` | DU | `ERAint-o3-clim.asc` |
| Water vapour | `h2o` | kg m^-2 | `h2o_ecmwf.dat` |
| Albedo | `albedo_vis`, `albedo_nir`, `albedo_land` | 0-1, and 1 on land | land use and `reflec.lut` |

They are independent, so any combination can be given, from any source, and each pixel
without a value in a field uses the climatology (or the land-use albedo). magic takes them
on its output grid: either as NetCDF files (the CLI), or as arrays from Python (the
`specmagic` package), which picks the time and regrids them first. Both run the same C++
(`runMagic` in `src/run.cpp`, the fields in `headers/input_fields.hpp`), and
`tests/test_cli_python.py` checks that they give the same results.

The data are not supplied with specMAGIC; you must provide them yourself.

### From the CLI

`scripts/mtg.sh` takes one NetCDF file per input, already on the output grid of
`magic-config.asc`, and passes them to `magic` as environment variables:

| Option | Environment variable | Variables in the file |
|---|---|---|
| `--aerosol-file FILE` | `AEROSOL_FILE` | `aod` |
| `--ozone-file FILE` | `OZONE_FILE` | `o3` |
| `--water-vapour-file FILE` | `WATER_VAPOUR_FILE` | `h2o` |
| `--albedo-file FILE` | `ALBEDO_FILE` | `albedo_vis`, `albedo_nir`, `albedo_land` |

Each variable has exactly the dimensions `(lat, lon)`, in that order (no time dimension),
with the `lat` and `lon` coordinates of the output grid (as `Geography::makeCoordinates`
makes them, to within 1e-3 degrees), and the units in the table above. It must be stored
unpacked: a variable with `scale_factor` or `add_offset` is refused. NaN, or the variable's
`_FillValue` or `missing_value`, marks pixels without a value. A text attribute `provenance` on
the variable, if there is one, is copied into the output files; otherwise the file's path is.
`--albedo-file` cannot be combined with `--albedo MODIS`.

Such files can be made with any tool, or with the `specmagic` package from sources on other
grids:

    import specmagic, xarray as xr
    cams = xr.open_zarr("cams.zarr")
    aod = specmagic.regrid_input(cams.aod550, "aerosol", "2026-07-15T12:20")
    specmagic.write_inputs("inputs/", aerosol=aod)        # inputs/aerosol.nc

### From Python

`specmagic.run` takes the image and the inputs as xarray DataArrays, on their own grids and
with a time dimension, picks the time slice for the image, regrids it, and returns GHI, DNI,
CAL and CSR as a Dataset:

    import xarray as xr
    import specmagic

    mtg = xr.open_zarr("vis.zarr")
    image = mtg.vis_06.sel(start_time="2025-01-01T12:00")
    cams = xr.open_zarr("cams.zarr")
    dini = xr.open_zarr("dini.zarr", decode_coords="all")     # keeps dini_projection with wvint
    lsasaf = xr.open_dataset("LSA-105_MTG_MTDAL_MTG-FD_202501010000.nc")

    out = specmagic.run(image, aerosol=cams.aod550, ozone=cams.gtco3, water_vapour=dini.wvint,
                        **specmagic.lsasaf_albedo(lsasaf), grid=(31, 72, -19, 34, 0.025))

To fill the DINI domain's surroundings with CAMS, regrid both and combine them before the run:

    when = image.start_time.values
    h2o = specmagic.regrid_input(dini.wvint, "water_vapour", when).fillna(
        specmagic.regrid_input(cams.tcwv, "water_vapour", when))
    out = specmagic.run(image, water_vapour=h2o.assign_attrs(provenance="DINI wvint, then CAMS tcwv"))

`uv sync` builds the package, including the C++ (`python/bindings.cpp`, with pybind11), and
builds it again when the C++ sources change (`[tool.uv] cache-keys` in `pyproject.toml`).

**The image** is the channel as satpy gives it: named after the channel, with the attributes
`resolution` and `calibration`, and a start time (attribute or coordinate `start_time`). It
may be counts, or for MTG reflectance, which is turned back into counts with the FCI
calibration in `python/specmagic/_image.py` (taken from the L1c files; satpy's
`radiance_unit_conversion_coefficient` attribute is checked against it). It may be the full
disc or a crop of it with its `x`/`y` coordinates in metres; a crop is put in its place on an
otherwise empty full disc. The output grid must then lie inside the crop. Counts packed with
the radiance calibration, as in the files `preprocessMTG.py` writes, must be opened with
`xr.open_dataset(path, mask_and_scale=False)`; otherwise xarray unpacks them to radiance, and
`specmagic.run` refuses the image.

**The output grid** is that of `magic-config.asc`, or `grid=(lat_min, lat_max, lon_min,
lon_max, resolution)`. Python gets its coordinates from the C++ code
(`specmagic.output_grid`), so both always agree.

**The grid of an input** is recognised from its dimensions: `latitude`/`longitude` (or
`lat`/`lon`) for a regular grid, possibly global and in 0..360 like CAMS; or `x`/`y` for a
projected grid, whose grid mapping variable must be one of its coordinates (open the data
with `decode_coords="all"`, or assign it), like DINI's Lambert conformal grid or LSA SAF's
geostationary one, with `x`/`y` in metres or, for the geostationary grid, in radians. A field
already on the output grid, with dimensions `(lat, lon)`, is used as it is.

**Units** are checked: ozone in kg m^-2 is converted to DU (dividing by 2.1415e-5), and
anything else than the units in the table above is refused.

**Time.** An input with a time dimension (any `datetime64` dimension) gives the slice nearest
to the image start. A later time is only used if it is already there, otherwise the latest
earlier one is, so 00:50 uses 01:00 if it exists and 00:00 if not. Up to 3 h before and
30 min after the image are accepted. Forecasts with the dimensions `forecast_reference_time`
and `step`, like DINI, give the most recent run for that valid time. The albedo and land
mask are daily: the image's own day, else the most recent day up to 10 days before it. If a
slice has no value on the output grid, the next one is tried. An input without a time
dimension is used as it is.

**Regridding.** The aerosol, ozone and water vapour are interpolated bilinearly, so a regional
field such as DINI's only covers its domain. The albedo and land mask are taken from the
nearest pixel. Only the window of the source that the output grid needs is read.

**LSA SAF.** `specmagic.lsasaf_albedo(ds)` turns an LSA SAF MTDAL dataset (an archive or the
daily files, e.g. `LSA-105_MTG_MTDAL_MTG-FD_202501010000.nc`) into the three albedo inputs:
`AL-VI-DH` and `AL-NI-DH` where `Z_Age` is 0 to 10 days, and the land mask from bit 0 of
`QFLAGS`.

### What magic does with them

**Albedo.** The black-sky (directional-hemispherical) albedos are used with no BRDF model:
visible (0.4-0.7 um) for the Kato bands up to band 16 and for satellite channels below 700 nm,
near-infrared (0.7-4.0 um) otherwise. A pixel uses the albedo input only where both are
available, so that it never mixes it with the land-use albedo across bands. Over water
(`albedo_land` 0) a missing albedo is expected and not reported; LSA SAF has none there.

The albedo is taken to be for the sun at local noon, as LSA SAF gives it. According to the
LSA SAF albedo ATBD (SAF/LAND/MF/ATBD_MTDAL 1.1, section 4.3.5), the directional-hemispherical
albedo is given for the solar zenith angle at local noon of each pixel and day, capped at 85
degrees, not for the time of the image. specMAGIC moves it to the current sun with the ratio
of `zenithCorrection` at the current and at the noon zenith angle, for the pixel's land-use
class, and caps the result at 1. The noon angle is |latitude - declination| for the image's
day. The ATBD also suggests blending with the bi-hemispherical albedo by the diffuse
fraction; that is not done, since LSA SAF provides the bi-hemispherical albedo only as
broadband.

**Range.** These fields do leave the LUT axes: CAMS ozone over Europe spans about 197-549 DU
against the 210-525 DU of `o3corr-l.lut`, and AOD reaches about 4.6 in dust events against the
2.0 of `magic-clear_spectral.lut`. Input values are therefore clipped to the axes. Climatology
values are not clipped. If a climatology gives water vapour or ozone outside its axis, the run
stops with an error; AOD beyond 2.0 is extrapolated linearly.

**Unusable inputs.** A given input that cannot be used at all (a file that cannot be read, a
missing variable, a field not on the output grid, or in Python a source without a usable time
or units) is reported, and the climatology (or land-use albedo) is used for it. The report is
repeated at the end of the run. With `--strict-inputs` (or `strict=True` in Python, which
then raises `RuntimeError`) the run stops instead. The CLI and Python record such an input
in the same way: its `*_source` attribute names the climatology, and `input_problems` says
why. Pixels without a value always fall back without stopping the run; magic warns
when more than 1% of the day pixels (of the land pixels, for albedo) do so.

### Provenance in the output files

When input fields are given, magic writes as global attributes of every output file (GHI,
DNI, CAL and CSR, all with the same attributes), and of the Dataset that `specmagic.run`
returns, where each input came from. They are metadata only: specMAGIC never reads them back.
They let you tell from the output alone, without the run log, which data it was made from,
and whether part of it fell back to the climatologies. Without inputs, there are none.

| Attribute | Written with | Content |
|---|---|---|
| `aod_source`, `aod_usage` | aerosol | where the field came from; shares of the day pixels, see below |
| `o3_source`, `o3_usage` | ozone | the same for ozone |
| `h2o_source`, `h2o_usage` | water vapour | the same for water vapour |
| `albedo_vis_source`, `albedo_nir_source`, `albedo_usage` | albedo | the same for the albedo |
| `input_problems` | any, only if there were problems | every problem reported during the run, separated by `; ` |

A source attribute is the field's provenance. From `specmagic`, it names the variable, its
time and the source, e.g. `wvint run 2026-08-23T13:00 +0 h, valid 2026-08-23T13:00 from
/path/to/dini.zarr`, the 13:00 run at step 0; CAMS gives its run and step too if it has the
`forecast_reference_time` and `forecast_period` coordinates along `valid_time`. The source is
the file a DataArray was opened from, which xarray records for NetCDF files but not for zarr
stores; otherwise it is the DataArray's `source` attribute if set (e.g.
`cams.aod550.assign_attrs(source="cams.zarr")`), else `memory`. A field passed already on the
output grid keeps its own `provenance` attribute. If the input
was given but not usable, the attribute names what was used instead, ending in
`(input missing)`, such as `aerosol climatology (input missing)`.

The usage attributes count the day pixels, since night and invalid pixels use no inputs.

* `aod_usage`, `o3_usage` and `h2o_usage`, e.g. `20.0% input, 80.0% climatology, 0.0% clipped
  to LUT`. "input" and "climatology" add up to 100%. "clipped to LUT" is the share whose
  input value lay outside the LUT axis and was clipped to it; these pixels are also counted
  under "input".
* `albedo_usage`, e.g. `99.9% of land day pixels from the input, 0.1% land-use albedo; 53.0% of
  all day pixels are water and use the land-use albedo`. The first two shares are of the land
  day pixels and add up to 100%. The last share is of all day pixels.

The same lines are printed at the end of the run, prefixed with `[WARNING]:` when the
climatology or land-use share exceeds 1%.

To read the attributes, use `ncdump -h`:

    ncdump -h out/GHI/GHI202608231300.nc

or, in Python:

    import xarray as xr
    attrs = xr.open_dataset("out/GHI/GHI202608231300.nc").attrs
    print(attrs["h2o_source"], attrs.get("input_problems", "no problems"))

## Unused files and dead config fields

Present in the repository, referenced by nothing, retained from the legacy DWD code:

* `climatologies/h2oclim.dat` — NCEP long-term monthly water vapour, 144 x 73, in **cm**.
  Superseded by `h2o_ecmwf.dat`; see the units warning above.
* `climatologies/landuse.dat` — superseded by `IGBPa_2006.map`.
* `climatologies/CVS/` — leftover CVS metadata.
* `climatologies/readme_aerosol_data.txt` — documentation of alternative aerosol
  climatologies and their grid dimensions.

Parsed from the config but never used: the recorded climatology dimensions `xadim`, `yadim`,
`xhdim`, `yhdim`, `xo3dim` and `yo3dim` (only `xadim` and `yadim` are sanity-checked for
positivity), plus `iconflag` and `iconres`.
