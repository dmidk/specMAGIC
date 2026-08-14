#!/usr/bin/env bash
set -euo pipefail

ROOT="$PWD"

END_REF="${1:-Today 08:00}"   # when the gif should END (UTC)
N_FRAMES="${N_FRAMES:-10}"
STEP_MIN="${STEP_MIN:-10}"
PRODUCT="${PRODUCT:-GHI}"     # choose GHI, DNI, or CAL
FPS="${FPS:-2}"               # playback speed

FRAMES_DIR="${ROOT}/gif_frames/${PRODUCT}"
GIF_OUT="${ROOT}/${PRODUCT,,}_last_${N_FRAMES}_frames.gif"

mkdir -p "$FRAMES_DIR"
rm -f "$FRAMES_DIR"/*.png

echo "Compiling C++ once before looping..."
mkdir -p "$ROOT/build"
cd "$ROOT/build"
cmake .. > cmake_configure.log
cmake --build .
cd "$ROOT"

## Parse end reference time, in UTC
end_epoch=$(TZ=UTC date --date="$END_REF" +%s) || {
    echo "[ERROR] Could not parse END_REF='$END_REF'"
    exit 1
}

echo "Generating ${N_FRAMES} frames for ${PRODUCT} ending at '${END_REF}'..."

for ((i=N_FRAMES-1; i>=0; i--)); do
    offset_sec=$(( i * STEP_MIN * 60 ))
    this_epoch=$(( end_epoch - offset_sec ))

    this_ref=$(TZ=UTC date --date="@$this_epoch" +"%Y-%m-%d %H:%M")
    stamp=$(TZ=UTC date --date="@$this_epoch" +"%Y%m%d%H%M")

    echo
    echo "=== Frame $((N_FRAMES - i)) / ${N_FRAMES} : ${this_ref} UTC ==="

    ## Run the timestep
    ## No need to rebuild inside it
    REBUILD=0 NOCLEAN=0 bash "$ROOT/scripts/mtg.sh" --ref-time "$this_ref" || {
        echo "[WARN] Failed for ${this_ref}; skipping frame."
        continue
    }

    ## Find newest output plot 
    latest_png=$(ls -1t "$ROOT/figs/${PRODUCT,,}_"*.png 2>/dev/null | head -n 1 || true)

    if [[ -z "${latest_png}" ]]; then
        echo "[WARN] No ${PRODUCT} PNG found for ${this_ref}; skipping frame."
        continue
    fi

    cp "$latest_png" "$FRAMES_DIR/${PRODUCT,,}_${stamp}.png"
    echo "Saved frame: $FRAMES_DIR/${PRODUCT,,}_${stamp}.png"
done

## Build GIF
uv run python "$ROOT/py_utils/make_gif.py" "$FRAMES_DIR" "$GIF_OUT" --fps "$FPS"

echo
echo "GIF written to: $GIF_OUT"
