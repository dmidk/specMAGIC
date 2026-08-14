# Clean Magic

## Attributions 

The SpecMAGIC code has been developed by Richard Müller, Uwe Pfeifroth and Elloise Fangel-Lloyd. Input data and climatologies are provided by Richard Müller.

Annette Hammer and Axel Kemper contributed during their time at the University of Oldenburg
to the development of the concepts of the methods and code.

## Preparing to run

SpecMAGIC is using `uv` for package management. It can be installed with 

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

[Please see the uv docs for more information.](https://docs.astral.sh/uv/)

## Input data

This code is assuming the presence of input satellite data. Currently, only Meteosat Third Generation (MTG) is supported, but support for other satellites is forthcoming. The user must provide a file path to available satellite data. 

## To run the code with MTG data

```bash
uv sync
```
will set up the python environment and necessary plugins. 

The driver script, which pre-processes the MTG data, calls the SpecMAGIC code, and does some small postprocessing, is called `mtg.sh`.

First, export the path to the raw MTG satellite data 

```bash
export MTG_RAW_DIR='path/to/data/archive'
```

To run the code for one timestep

```bash
./scripts/mtg.sh --ref-time "Yesterday 08:00"
```

The default MTG channel is `vis_06`, which is 640 nm. Other channels can be specified using `--channel`, e.g. 

```bash
./scripts/mtg.sh --ref-time "Today 09:00" --channel "nir_13"
```

To see a full list of input options, 

```bash
./scripts/mtg.sh -h
```

If no time is specified, "Yesterday 08:00" will be used as a default. 

Figures will appear in the `/figs` directory. By default, plots of GHI, DNI, CAL (cloud index) and CSR (clear-sky GHI) are produced. 

## Configuring SpecMAGIC

The SpecMAGIC config is found in `root/magic-config.asc`. Here also is a short overview of the available lookup tables and input data. Alterations here are at the user's own risk.



