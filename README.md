# Clean Magic

## Attributions 

The SpecMAGIC code has been developed by Richard Müller, Uwe Pfeifroth and Elloise Fangel-Lloyd. Input data and climatologies are provided by Richard Müller.

Annette Hammer and Axel Kemper contributed during their time at the University of Oldenburg to the development of the concepts of the methods and code. Tanja Behrendt contributed to the evaluation of SPECMAGIC and supported the project with her knowledge about the physics of spectral resolved irradiance.

The 2012 version of SPECMAGIC is described in
Mueller, R.; Behrendt, T.; Hammer, A.; Kemper, A. [A New Algorithm for the Satellite-Based Retrieval of Solar Surface Irradiance in Spectral Bands](https://doi.org/10.3390/rs4030622).

## Preparing to run

SpecMAGIC is using `uv` for package management. It can be installed with 

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

[Please see the uv docs for more information.](https://docs.astral.sh/uv/)

## Input data

This code is assuming the presence of input satellite data. Currently, only Meteosat Third Generation (MTG) is supported, but support for other satellites is forthcoming. The user must provide a file path to available satellite data. 

## To run the code in demo mode

```bash
uv sync
```
will set up the python environment and necessary plugins. 

The driver script, which pre-processes the MTG data, calls the SpecMAGIC code, and does some small postprocessing, is called `mtg.sh`. By default, it runs in demo mode, using the data that comes pre-packaged in the `test_data` directory. These test data originate from EUMETSAT Meteosat products published in 2026. The single `.zip` file corresponds to a single level 1c satellite image taken by MTG/FCI on the 23rd August 2026 at 13:00 UTC. These exist under a CC-BY-4.0 free license; for details, please [the EUMETSAT website.](https://user.eumetsat.int/catalogue/EO:EUM:DAT:0989)

```bash
./scripts/mtg.sh
```

The default MTG channel is `vis_06`, which is 640 nm. Other channels can be specified using `--channel`, e.g. 

```bash
./scripts/mtg.sh --channel "nir_13"
```

To see a full list of input options, 

```bash
./scripts/mtg.sh -h
```
Figures will appear in the `/figs` directory. By default, plots of GHI, DNI, CAL (cloud index) and CSR (clear-sky radiance) are produced. 

## To run the code with other MTG data 

As above, 
```bash
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


## Configuring SpecMAGIC

The SpecMAGIC config is found in `root/magic-config.asc`. Here also is a short overview of the available lookup tables and input data. Alterations here are at the user's own risk.



## MODIS BRDF albedo maps

This section explains how to upgrade from coarse land-use maps to MODIS hi-res land albedo maps for the clearsky calculation.

The MODIS BRDF NetCDF files are large and are not stored in Git.

Download them with:
```
bash uv run python py_utils/download_MODIS_maps.py climatologies/modis-brdf
``` 

If you have access to a pre-downloaded copy of the MODIS BRDF files, you can instead create a local symlink from the project root pointing to that location:

```
ln -s /path/to/shared/modis climatologies/modis-brdf
```