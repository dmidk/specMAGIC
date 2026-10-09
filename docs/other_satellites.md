## To run the code with other satellite data 

As above, it is necessary to first set up the python environment and necessary plugins. 

As above, 
```bash
git clone [url]
uv sync
```
will set up the python environment and necessary plugins. 

There are several other steps first needed in order to run the code on data from other satellites.


### 1. The driver script

A new driver script must be created. The current one, `mtg.sh`, contains lots of MTG-exclusive things. It would be beneficial to copy this to a new script e.g. `msg.sh` and give the directories new names accordingly. 

The driver script also accepts the viewing channel as an input arg, and will reject those channels that do not correspond to MTG channels. The list of accepted channels must be updated. There may be other things also in `mtg.sh` that do not translate -- the intention was never to be able to reuse this script for other satellites, but to create a new script for each.

The driver script sets an env variable which specifies the satellite being used. 

```bash
export SATELLITE=mtg
```

This tells the specMAGIC to search for a config file `mtg.toml`. A new `.toml` will be needed, see below, but remember also to change this export statement.

The driver script defines the `$EXTENT`, which is the plotting extent used by the pre- and post-processing routines to make maps of the data. When modifying the driver script, make sure that the desired plotting `$EXTENT` is set to the region of interest.

### 2. Input data

SpecMAGIC is expecting a single `*.nc` file containing the entire image for a given timestep. This file is expected to only contain a single channel. The preprocessing step in the demo version (which calls `extract.py`) will extract the `.zip` test data and produce this single `.nc` file. For differing data storage or file formats, additional preprocessing may be required.

The input to specMAGIC needs to be a single file containing a single channel. 

### 3. Satellite information 

specMAGIC will extract the image size and horizontal resolution directly from the file itself. No need to supply the number of pixels, etc. 

The user needs to gather the following information:
- Satellite distance in kilometres (`satellite_radius_km`). This is the **geostationary radius**, i.e. the height of the satellite above ground plus the equatorial Earth radius.
- Subsatellite longitude in degrees (`subsatellite_longitude_deg`) in the `[navigation]` section of the satellite TOML.
- The list of channel names and the wavelengths they correspond to (`[channels]`)
- Number of minutes it takes to scan the entire image (`full_disk_scan_min`)
- The grid sampling angle in radians. For MTG this was found in the [EUMETSAT FCI level 1C user guide documentation.](https://user.eumetsat.int/resources/user-guides/mtg-fci-level-1c-data-guide) 

All of the above are hard requirements.

The code **also** assumes that the origin (phrased in the EUMETSAT documentation as line and column offset) is in the exact centre of the image. This may or may not be the case for other satellites. This assumption is made in `parseImage()`, where the exact centre is chosen.

```cpp
info.column_offset = static_cast<int>(nx / 2);
info.line_offset = static_cast<int>(ny / 2);
```

If the origin is for some reason not the centre of the image, the user should exercise caution.

### 4. Setting up to run

Under `satellites/`, there exists a file `mtg.toml`. To use a new satellite, a new `.toml` must be created containing all of this above information. See the existing `mtg.toml` for inspiration. 

Any channel which specMAGIC is intended to run on and its wavelength should be included in the list. If the satellite has channels at varying resolutions, all of these resolutions and their accompanying grid sampling angular distance must be listed.

For example,

```bash
[grids."1000"]
angular_sampling_rad = 2.7943576e-5
```

refers to 1km resolution at nadir and its corresponding grid sampling angle. Note that the spatial resolution in use -- in this case 1 km -- will be extracted directly from the `netcdf` file. This bit here is to translate the distance-resolution 1 km to the angular grid sampling distance shown above, which varies between satellites.

Follow the `readme.md` for instructions to specify the satellite data file path and run in non-demo mode.

### 5. Hardware

The driver script specifies that specMAGIC should use 8 threads/cores in parallel, which are available on generally any laptop computer. If running on a larger setup, there may be more available. Change the number of tweaks by altering 

```bash
export OMP_NUM_THREADS=8
```

Additional speedup should not be expected beyond 64 threads.

### 6. Geography 

See the [readme](https://github.com/dmidk/specMAGIC/blob/main/README.md) for notes on changing the simulated geographic area. 
