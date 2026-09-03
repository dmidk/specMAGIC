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
