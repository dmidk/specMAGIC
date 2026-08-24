import os
## Do not remove this
## Otherwise it might try multithreaded writes and will segfault
os.environ["DASK_NUM_WORKERS"] = "1"
chatty = os.getenv("CHATTY", "0").lower() in ("1", "true", "yes", "on")
import hdf5plugin
import sys
from satpy import Scene
from glob import glob
import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np
from helpers import plot_handled_MTG, check_required_chunks, suppress_known_warnings

# Expecting: program.py <date> <cycle> <DATA_DIR> <OUT_DIR> <FIGS_DIR>
if len(sys.argv) < 7:
    print(
        " ERROR: Not enough command-line arguments.\n"
        " Usage: python program <date> <cycle> <DATA_DIR> <OUT_DIR> <GLOBAL_FIGS_DIR>\n"
        f" Got argv={sys.argv}"
    )
    sys.exit(2)

date = sys.argv[1]
cycle = sys.argv[2]
channel = sys.argv[3]
data_dir = str(sys.argv[4].rstrip("/"))  # avoid trailing slash issues
out_dir = str(sys.argv[5].rstrip("/"))  # avoid trailing slash issues
figs_dir = str(sys.argv[6].rstrip("/"))  # avoid trailing slash issues


if (chatty):
    print(" I am looking for files with date", date, "and cycle number", cycle, "...")

# files always have this funky name structure
pattern = (data_dir + "/W_XX-EUMETSAT-Darmstadt,IMG+SAT,MTI1+FCI-1C-RRAD-*-FD--CHK-BODY-*_EUMT_*_OPE_"
    + str(date) + "*N_*_O_0*" + str(cycle) + "_*.nc")

lat_min = 40 
lat_max = 65
lon_min = -15
lon_max = 30
extent = [lat_min, lon_min, lat_max, lon_max, 1]

filenames = glob(pattern)

if (chatty or len(filenames) < 20):
    print(" I am expecting 20 files. I found", len(filenames), "MTG files for the specified date/time")

# Didn't find any files :(
if not filenames:
    print(
        " ERROR: No input files found (yet).\n"
        f" Looked for files matching:\n  {pattern}\n"
        f" Arguments were:\n  date={date}\n  cycle={cycle}\n  DATA_DIR={data_dir}\n"
        " This can happen if the data have not arrived yet, the slot/time is mismatched, "
        "or the directory/pattern is incorrect."
    )
    sys.exit(1)

# Load the data into a satpy Scene
sc = Scene(filenames=filenames, reader="fci_l1c_nc", reader_kwargs={"engine":"h5netcdf"})

# Safety check for data arrival
check_required_chunks(filenames)

# We only care about this one band
sc.load([channel], calibration="counts")


# Make a filename
outfile = os.path.join(out_dir,
    "{area.area_id}-{name}-{start_time:%Y%m%d%H%M%S}-{end_time:%Y%m%d%H%M%S}-"+cycle+".nc"
)

# Save the data
# This is now all files belonging to this cycle but only for the single channel
with suppress_known_warnings():
    sc.save_dataset(
        channel,
        outfile,
        engine="h5netcdf",
        include_lonlats = False,
        encoding={channel: {"dtype": "int16"}},
    )

# Plot data, only if verbose
if not chatty: 
    exit()

figs_folder = figs_dir+"/MTG_preparation"
if not os.path.exists(figs_folder):
    os.makedirs(figs_folder+"/")

# Looks for data
pattern = os.path.join(out_dir, f"*-{channel}-*.nc")
matches = glob(pattern)
if not matches:
    raise FileNotFoundError(f" No output file found matching: {pattern}")
out_path = max(matches, key=os.path.getmtime)

# Opens the same file that was just open
# this is for idiot checking
ds = xr.open_dataset(out_path, engine="h5netcdf")
da = ds[channel]

plot_handled_MTG(ds, channel, figs_folder, extent)
print("Wrote plot of MTG counts to directory ", figs_folder)