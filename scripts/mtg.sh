#!/usr/bin/env bash
set -euo pipefail

## Project root
ROOT="$PWD"

## Defaults
REF_TIME="" # UTC time!
NOCLEAN="${NOCLEAN:-1}"      # set to 1 to keep existing outputs
REBUILD="${REBUILD:-1}"      # set to 0 to skip cmake rebuild
VERBOSE="${VERBOSE:-1}"
TIMER="${TIMER:-0}"          # whether or not to make output of timing the c code
HI_PRECISION=ON
CHANNEL="vis_06"
DEMO="${DEMO:-0}"              # 1 = use local demo data, 0 = user supplies MTG_RAW_DIR
MTG_RAW_DIR="${MTG_RAW_DIR:-}"
ALBEDO="${ALBEDO:-LANDMAP}"
STRICT_INPUTS="${STRICT_INPUTS:-0}"   # 1 = stop if a given input file is not usable, instead of warning
# Optional input fields for the time of the image, one NetCDF file per table on the output grid
# (see docs/inputs.md). Pixels without a value in them use the climatology.
AEROSOL_FILE="${AEROSOL_FILE:-}"
OZONE_FILE="${OZONE_FILE:-}"
WATER_VAPOUR_FILE="${WATER_VAPOUR_FILE:-}"
ALBEDO_FILE="${ALBEDO_FILE:-}"
EXTENT="40.0 -15.0 65.0 30.0 1.0"

valid_channels=(vis_04 vis_05 vis_06 vis_08 vis_09 nir_13 nir_16 nir_22 ir_38 wv_63 wv_73 ir_87 ir_97 ir_105 ir_123 ir_133)

usage() {
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  --ref-time TIME       Reference time (not available in demo mode)"
    echo "  --channel CHANNEL     Satellite channel (default: $CHANNEL). Options: vis_04 vis_05 vis_06 vis_08 vis_09 nir_13 nir_16 nir_22 ir_38 wv_63 wv_73 ir_87 ir_97 ir_105 ir_123 ir_133"
    echo "  --albedo ALBEDO       Albedo type (default: $ALBEDO). Options: LANDMAP MODIS"
    echo "  --aerosol-file FILE   Aerosol optical depth for the image's time, instead of the climatology"
    echo "  --ozone-file FILE     Ozone for the image's time, instead of the climatology"
    echo "  --water-vapour-file FILE"
    echo "                        Water vapour for the image's time, instead of the climatology"
    echo "  --albedo-file FILE    Black-sky albedo at local noon, e.g. from LSA SAF, instead of the"
    echo "                        land-use albedo. The input files must be on the output grid; see"
    echo "                        docs/inputs.md"
    echo "  --strict-inputs       Stop if a given input file is not usable, instead of warning and"
    echo "                        falling back to the climatologies"
    echo "  --demo                Use demo data (default)"
    echo "  --no-demo             Use user-specified MTG_RAW_DIR instead of demo data"
    echo "  --clean               Delete existing outputs before running"
    echo "  --noclean             Keep existing outputs"
    echo "  --no-rebuild          Skip magic rebuild"
    echo "  --whole-domain        Plot the whole domain, rather than the DK cutout"
    echo "  --timer               Enable timing output"
    echo "  --verbose             Enable verbose output"
    echo "  -h, --help            Show this help"
}

is_valid_channel() {
    local ch="$1"

    for valid in "${valid_channels[@]}"; do
        if [[ "$ch" == "$valid" ]]; then
            return 0
        fi
    done

    return 1
}

export SATELLITE=mtg 

while [[ $# -gt 0 ]]; do
    case "$1" in

        --channel)
            CHANNEL="$2"

            if ! is_valid_channel "$CHANNEL"; then
                echo "Invalid channel: $CHANNEL"
                echo "Available channels: ${valid_channels[*]}"
                exit 1
            fi

            shift 2
            ;;

        --albedo)
            ALBEDO="$2"
            shift 2
            ;;

        --aerosol-file)
            AEROSOL_FILE="$2"
            shift 2
            ;;

        --ozone-file)
            OZONE_FILE="$2"
            shift 2
            ;;

        --water-vapour-file)
            WATER_VAPOUR_FILE="$2"
            shift 2
            ;;

        --albedo-file)
            ALBEDO_FILE="$2"
            shift 2
            ;;

        --strict-inputs)
            STRICT_INPUTS=1
            shift
            ;;

        --ref-time)
            REF_TIME="$2"
            shift 2
            ;;

        --demo)
            DEMO=1
            shift
            ;;

        --no-demo)
            DEMO=0
            shift
            ;;

        --clean)
            NOCLEAN=0
            shift
            ;;

        --noclean)
            NOCLEAN=1
            shift
            ;;

        --no-rebuild)
            REBUILD=0
            shift
            ;;

        --whole-domain)
            EXTENT="0.0 -60.0 65.0 60.0 1.0"
            shift
            ;;

        --timer)
            TIMER=1
            shift
            ;;

        --verbose)
            VERBOSE=1
            shift
            ;;

        -h|--help)
            usage
            exit 0
            ;;

        *)
            echo "Unknown argument: $1"
            usage
            exit 1
            ;;

    esac
done

case "$ALBEDO" in
    LANDMAP|MODIS)
        ;;
    *)
        echo "Error: invalid albedo type '$ALBEDO'."
        echo "Valid options are: LANDMAP, MODIS."
        exit 1
        ;;
esac

if [[ -n "$ALBEDO_FILE" && "$ALBEDO" == "MODIS" ]]; then
    echo "Error: --albedo MODIS and --albedo-file cannot be used together."
    exit 1
fi

# Check for ref time
# In demo mode, the time is fixed, so the user should not pass a time
if [[ "$DEMO" -eq 1 ]]; then

    if [[ -n "$REF_TIME" ]]; then
        echo "Error: --ref-time cannot be used in demo mode."
        echo "The demo data always uses 23/08/2026 13:00 UTC."
        echo "Do not pass a --ref-time in demo mode."
        exit 1
    fi

    REF_TIME="2026-08-23 13:00"

else

    if [[ -z "$REF_TIME" ]]; then
        echo "Error: --ref-time must be specified when not in demo mode."
        exit 1
    fi

fi

# MTG raw directory must be supplied if not in demo mode
if [[ "$DEMO" -eq 1 ]]; then
    MTG_RAW_DIR="${ROOT}/test_data"
else
    if [[ -z "$MTG_RAW_DIR" ]]; then
        echo "Error: MTG_RAW_DIR is not set."
        echo "Set it with, for example:"
        echo "  export MTG_RAW_DIR=/path/to/MTG_data"
        echo " or else run specMAGIC in demo mode."
        exit 1
    fi
fi

## Parse time of interest
ymd=$(date +"%Y%m%d" --date="$REF_TIME")
hh=$(date +"%H" --date="$REF_TIME")
mm=$(date +"%M" --date="$REF_TIME")



if [[ "$DEMO" -eq 1 ]]; then

    echo "Extracting demo data..."
    uv run python "${ROOT}/py_utils/extract.py" "$MTG_RAW_DIR" || { echo "Failed to extract demo data!!" >&2; exit 1;}

fi

## Output locations
MTG_READY_DIR="${ROOT}/MTG_handled"
FIGS_DIR="${ROOT}/figs"
OUTPATH="${ROOT}/out"

mkdir -p "$MTG_READY_DIR" "$FIGS_DIR" "$OUTPATH"

if [[ "$NOCLEAN" != "1" ]]; then
    rm -rf "${MTG_READY_DIR:?}"/*
    rm -rf "${FIGS_DIR:?}"/*
    rm -rf "${OUTPATH:?}"/*
fi

mkdir -p "$OUTPATH/CAL" "$OUTPATH/DNI" "$OUTPATH/GHI" "$OUTPATH/CSR"

## Cycle number (10-minute resolution)
# Add +1 to get correct time
cycle=$(( 10#$hh * 6 + 10#$mm / 10 ))
((cycle++))


## HDF5 env
export HDF5_USE_FILE_LOCKING=FALSE
export HDF5_PLUGIN_PATH="$(uv run python -c 'import hdf5plugin; print(hdf5plugin.PLUGINS_PATH)')"
export CHATTY="$VERBOSE"

echo " ------------------------------------------------------------------ "
echo "Preparing MTG data for UTC time $(date +"%d/%m/%Y %H:%M" --date="$REF_TIME") ($REF_TIME)..."

## Preprocess
uv run python "$ROOT/py_utils/preprocessMTG.py" "$ymd" "$cycle" "$CHANNEL" "$MTG_RAW_DIR" "$MTG_READY_DIR" "$FIGS_DIR" $EXTENT\
    || { echo "Failed to preprocess MTG data!!" >&2; exit 1; }

echo "Finished pre-processing MTG data."

## Build list of handled data files
cd "$MTG_READY_DIR"
filename="/fill-${CHANNEL}-${ymd}*${cycle}.nc"
rm -f "$ROOT/image_list.txt"
ls -t $MTG_READY_DIR$filename > "$ROOT/image_list.txt"
cd "$ROOT"

## Build C++ only if requested
if [[ "$REBUILD" == "1" ]]; then
    mkdir -p build
    rm -rf build/*
    cd build
    cmake -DMAGIC_HI_PRECISION=$HI_PRECISION .. > cmake_configure.log
    cmake --build . > cmake_configure.log
    cd "$ROOT"
fi

if [[ "$ALBEDO" == "MODIS" ]]; then
    echo "Downloading MODIS data..."
    uv run python py_utils/download_MODIS_maps.py climatologies/modis-brdf
fi



export OMP_NUM_THREADS=8
export AEROSOL_FILE OZONE_FILE WATER_VAPOUR_FILE ALBEDO_FILE STRICT_INPUTS
export ALBEDO=$ALBEDO
echo "Calling magic..."

## Run C++ driver
"$ROOT/build/magic" "$ROOT/" $CHANNEL $TIMER


## Post-process plots
uv run python "$ROOT/py_utils/post.py" "$OUTPATH" "$FIGS_DIR" $EXTENT || { echo "Failed to post-process MTG data!!" >&2; exit 1; }

echo "Done! Thanks for now."
echo " ------------------------------------------------------------------ "