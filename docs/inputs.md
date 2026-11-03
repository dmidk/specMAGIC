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

## Inputs from zarr archives

Two options of the driver replace static tables with fields for the time of the satellite
image, read from zarr archives:

* `--cams-dini-atmosphere` replaces the aerosol, ozone and water vapour climatologies.
* `--albedo LSASAF` replaces the land-use albedo (and `reflec.lut`) by the LSA SAF albedo.

They are independent, so that their effects can be judged separately, and either can be
combined with the other albedo options. The archive paths are taken from the environment
variables `CAMS_ZARR`, `DINI_ZARR` and `LSASAF_ZARR`; there are no defaults. An archive that
is requested but not set is skipped with a warning, or stops the run with `--strict-inputs`.
`py_utils/extract_zarr_inputs.py` reads the archives,
converts and regrids the fields onto the output grid, and writes one NetCDF file per archive
and source time into `zarr_handled/`. Their paths are listed in `zarr_inputs.txt`, which
`magic` reads through the `ZARR_INPUTS` environment variable (`ZarrInputs::Inputs` in
`headers/zarr_inputs.hpp`). A file that an earlier run wrote for the same source time and
output grid is reused, so the images of one hour share their CAMS and DINI extraction.

| Input | Archive and variable | Replaces | Conversion |
|---|---|---|---|
| Aerosol | CAMS `aod550` | `macc-ecwmf-clim03-06rich.dat` | none, dimensionless |
| Ozone | CAMS `gtco3` | `ERAint-o3-clim.asc` | kg m^-2 to DU, dividing by 2.1415e-5 |
| Water vapour | DINI `wvint`, then CAMS `tcwv` | `h2o_ecmwf.dat` | none, kg m^-2 |
| Albedo | LSA SAF `AL-VI-DH`, `AL-NI-DH` | land use and `reflec.lut` | none, the 1e-4 scale is applied on reading |

**Time.** For each archive, the nearest time to the image start is used. A later time is only
used if it is already in the archive, otherwise the latest earlier one is, so 00:50 uses
01:00 if it exists and 00:00 if not. CAMS and DINI accept up to 3 h before and 30 min after
the image; for DINI, the most recent run with that valid time is used. LSA SAF is daily: the
image's own day, else the most recent day up to 10 days before it, with a warning.

**Regridding.** CAMS is interpolated bilinearly from its 0.4 degree global grid, and DINI
bilinearly in its Lambert projection, so it is only available inside the DINI domain. LSA SAF
is taken from the nearest pixel of its geostationary grid.

**Albedo.** The black-sky (directional-hemispherical) albedos are used with no BRDF model:
visible (0.4-0.7 um) for the Kato bands up to band 16 and for satellite channels below 700 nm,
near-infrared (0.7-4.0 um) otherwise. A pixel uses LSA SAF only where both are available, so
that it never mixes the two sources across bands. Pixels whose albedo is more than 10 days old
(`Z_Age`) are dropped, and the pixel reverts to the land-use albedo. LSA SAF has no albedo 
over water (bit 0 of `QFLAGS` unset), so water uses the land-use albedo, which is expected and 
not reported.

According to the LSA SAF albedo ATBD (SAF/LAND/MF/ATBD_MTDAL 1.1, section 4.3.5), the
directional-hemispherical albedo is given for the solar zenith angle at local noon of each
pixel and day, capped at 85 degrees, not for the time of the image. specMAGIC moves it to the
current sun with the ratio of `zenithCorrection` at the current and at the noon zenith angle,
for the pixel's land-use class, and caps the result at 1. The noon angle is |latitude -
declination| for the image's day, also when the albedo comes from an earlier day. The ATBD also
suggests blending with the bi-hemispherical albedo by the diffuse fraction; that is not done,
since LSA SAF provides the bi-hemispherical albedo only as broadband.

**Range.** These fields do leave the LUT axes: CAMS ozone over Europe spans about 197-549 DU 
against the 210-525 DU of `o3corr-l.lut`, and AOD reaches about 4.6 in dust events against the 
2.0 of `magic-clear_spectral.lut`. Values from the zarr archives are therefore clipped to the 
axes.

**Missing inputs.** A requested input that is not available at all, because an archive cannot
be read, does not match the format below, no time lies within the window, the archive does
not cover the output grid, or a file is not on the output grid, is reported, and the
climatology (or land-use albedo) is used for it. The report is repeated at the end of the run.
With `--strict-inputs` the run stops instead. Pixels that an available field does not cover,
such as those outside the DINI domain, always fall back without stopping the run; `magic`
warns when more than 1% of the day pixels (of the land pixels, for albedo) do so.

**Provenance.** `magic` prints which share of the day pixels came from each source and how many
were clipped, and writes, as global attributes of the output files, the source of each
requested input (archive, variable, source time, DINI run), those shares, and any problems.

### Format of the zarr archives

The archives are not supplied with specMAGIC; you must provide them yourself. They must open
with `xarray.open_zarr` and have the layout below.
Names are case-sensitive, and dimensions must be in the order given. Time coordinates must
decode to `datetime64`.

**CAMS** (`CAMS_ZARR`), hourly fields on a global regular latitude/longitude grid:

| Name | Dimensions | Requirements |
|---|---|---|
| `valid_time` | `valid_time` | the time each field is valid at |
| `latitude` | `latitude` | degrees, strictly monotonic, either direction |
| `longitude` | `longitude` | degrees in 0..360, strictly ascending, evenly spaced, covering the globe |
| `aod550` | `valid_time, latitude, longitude` | units `~`, `1` or empty |
| `gtco3` | `valid_time, latitude, longitude` | units `kg m**-2` or `kg m-2` |
| `tcwv` | `valid_time, latitude, longitude` | units `kg m**-2` or `kg m-2` |

A CAMS variable that is missing or has the wrong units is reported on its own; the others are
still used.

**DINI** (`DINI_ZARR`), forecasts on a projected grid:

| Name | Dimensions | Requirements |
|---|---|---|
| `forecast_reference_time` | `forecast_reference_time` | the start of each run |
| `step` | `step` | forecast step, as `timedelta64` or in whole hours |
| `x`, `y` | `x`, `y` | projection coordinates, strictly ascending |
| `latitude`, `longitude` | `y, x` | degrees, used to find the domain |
| `dini_projection` | none | attribute `crs_wkt` with the projection as WKT |
| `wvint` | `forecast_reference_time, step, y, x` | units `kg m**-2` or `kg m-2` |

**LSA SAF** (`LSASAF_ZARR`), the daily MTDAL product on the geostationary grid:

| Name | Dimensions | Requirements |
|---|---|---|
| `time` | `time` | 00:00 of each product day |
| `x`, `y` | `x`, `y` | scan angles in radians, evenly spaced, either direction |
| `geostationary` | none | CF grid mapping attributes, including `perspective_point_height` |
| `AL-VI-DH`, `AL-NI-DH` | `time, y, x` | must decode to albedo in 0-1, e.g. int16 with `scale_factor` 1e-4 |
| `QFLAGS` | `time, y, x` | integer, bit 0 set on land |
| `Z_Age` | `time, y, x` | integer, age of the albedo in days, negative where there is none |

`py_utils/verify_zarr_inputs.py` checks an archive against this layout without reading any
data:

    uv run python py_utils/verify_zarr_inputs.py --cams $CAMS_ZARR --dini $DINI_ZARR --lsasaf $LSASAF_ZARR

`extract_zarr_inputs.py` runs the same checks on every archive it opens, so the pipeline needs
no separate step for it. An archive that fails is treated like a missing one.

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
