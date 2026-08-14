#!/usr/bin/env bash
set -euo pipefail

## Inputs
DATE_INPUT="${1:?Provide date (e.g. 2026-06-14)}"
START_HOUR="${2:-00}"
END_HOUR="${3:-23}"



echo "-------------------------------------------------------------"
echo "Processing date: $DATE_INPUT"
echo "Hours: ${START_HOUR}:00 → ${END_HOUR}:59"
echo "-------------------------------------------------------------"

## Force read in base 10
START_HOUR=$((10#$START_HOUR))
END_HOUR=$((10#$END_HOUR))

## Loop 
for hour in $(seq -w "$START_HOUR" "$END_HOUR"); do
    for minute in 00 10 20 30 40 50; do

        TIME_STR="$DATE_INPUT ${hour}:${minute}"

        echo ""
        echo ">>> Running for $TIME_STR"
        echo "-------------------------------------------------------------"

        ./scripts/mtg.sh --ref-time "$TIME_STR"

    done
done

echo ""
echo "Finished processing $DATE_INPUT (${START_HOUR}-${END_HOUR})"