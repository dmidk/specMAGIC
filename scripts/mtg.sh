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
CAMS_DINI_ATMOSPHERE="${CAMS_DINI_ATMOSPHERE:-0}"   # 1 = take aerosol, ozone and water vapour from the zarr archives
STRICT_INPUTS="${STRICT_INPUTS:-0}"                 # 1 = stop if a requested zarr input is missing, instead of warning
# Paths of the zarr archives. An unset archive is skipped, and what it would provide
# falls back to the climatologies (or stops the run with --strict-inputs)
CAMS_ZARR="${CAMS_ZARR:-}"
DINI_ZARR="${DINI_ZARR:-}"
LSASAF_ZARR="${LSASAF_ZARR:-}"
EXTENT="40.0 -15.0 65.0 30.0 1.0"

valid_channels=(vis_04 vis_05 vis_06 vis_08 vis_09 nir_13 nir_16 nir_22 ir_38 wv_63 wv_73 ir_87 ir_97 ir_105 ir_123 ir_133)

usage() {
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  --ref-time TIME       Reference time (not available in demo mode)"
    echo "  --channel CHANNEL     Satellite channel (default: $CHANNEL). Options: vis_04 vis_05 vis_06 vis_08 vis_09 nir_13 nir_16 nir_22 ir_38 wv_63 wv_73 ir_87 ir_97 ir_105 ir_123 ir_133"
    echo "  --albedo ALBEDO       Albedo type (default: $ALBEDO). Options: LANDMAP MODIS LSASAF"
    echo "                        LSASAF reads the LSA SAF zarr archive (LSASAF_ZARR)"
    echo "  --cams-dini-atmosphere"
    echo "                        Take aerosol and ozone (CAMS) and water vapour (DINI, then CAMS) from the"
    echo "                        zarr archives (CAMS_ZARR, DINI_ZARR) instead of the climatologies"
    echo "  --strict-inputs       Stop if a requested zarr input is missing, instead of warning and"
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

        --cams-dini-atmosphere)
            CAMS_DINI_ATMOSPHERE=1
            shift
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
    LANDMAP|MODIS|LSASAF)
        ;;
    *)
        echo "Error: invalid albedo type '$ALBEDO'."
        echo "Valid options are: LANDMAP, MODIS, LSASAF."
        exit 1
        ;;
esac

# The zarr archives that were asked for must have a path
required_zarr=()
if [[ "$CAMS_DINI_ATMOSPHERE" == "1" ]]; then
    required_zarr+=(CAMS_ZARR DINI_ZARR)
fi
if [[ "$ALBEDO" == "LSASAF" ]]; then
    required_zarr+=(LSASAF_ZARR)
fi
for var in "${required_zarr[@]}"; do
    if [[ -z "${!var}" ]]; then
        if [[ "$STRICT_INPUTS" == "1" ]]; then
            echo "Error: $var is not set, and --strict-inputs was given."
            echo "Set it with, for example:"
            echo "  export $var=/path/to/archive.zarr"
            exit 1
        fi
        echo "[WARNING] $var is not set, so that archive is skipped. magic reports what it falls back to." >&2
    fi
done

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
ZARR_READY_DIR="${ROOT}/zarr_handled"
FIGS_DIR="${ROOT}/figs"
OUTPATH="${ROOT}/out"

mkdir -p "$MTG_READY_DIR" "$FIGS_DIR" "$OUTPATH"

if [[ "$NOCLEAN" != "1" ]]; then
    rm -rf "${MTG_READY_DIR:?}"/*
    rm -rf "${ZARR_READY_DIR:?}"/*
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



## Extract the zarr slices for this image. Unless --strict-inputs is given, failing here
## is not fatal: magic reports what is missing and uses the climatologies for it
ZARR_INPUTS="$ROOT/zarr_inputs.txt"
rm -f "$ZARR_INPUTS"
if [[ "$CAMS_DINI_ATMOSPHERE" == "1" || "$ALBEDO" == "LSASAF" ]]; then
    extract_args=()
    if [[ "$CAMS_DINI_ATMOSPHERE" == "1" ]]; then
        extract_args+=(--cams "$CAMS_ZARR" --dini "$DINI_ZARR")
    fi
    if [[ "$ALBEDO" == "LSASAF" ]]; then
        extract_args+=(--lsasaf "$LSASAF_ZARR")
    fi
    if [[ "$STRICT_INPUTS" == "1" ]]; then
        extract_args+=(--strict)
    fi

    echo "Extracting zarr inputs..."
    if ! uv run python "$ROOT/py_utils/extract_zarr_inputs.py" "$ROOT/image_list.txt" "$ROOT/magic-config.asc" \
        "$ROOT/satellites/${SATELLITE}.toml" "$ZARR_READY_DIR" "$ZARR_INPUTS" "${extract_args[@]}"; then
        if [[ "$STRICT_INPUTS" == "1" ]]; then
            echo "[ERROR] Could not extract the requested zarr inputs, and --strict-inputs was given." >&2
            exit 1
        fi
        echo "[WARNING] Could not extract zarr inputs. Using the climatologies." >&2
    fi
fi

export OMP_NUM_THREADS=8
export ZARR_INPUTS CAMS_DINI_ATMOSPHERE STRICT_INPUTS
export ALBEDO=$ALBEDO
echo "Calling magic..."

## Run C++ driver
"$ROOT/build/magic" "$ROOT/" $CHANNEL $TIMER


## Post-process plots
uv run python "$ROOT/py_utils/post.py" "$OUTPATH" "$FIGS_DIR" $EXTENT || { echo "Failed to post-process MTG data!!" >&2; exit 1; }

echo "Done! Thanks for now."
echo " ------------------------------------------------------------------ "