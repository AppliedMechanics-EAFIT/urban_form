#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/private/tmp/osmx_matplotlib}"
mkdir -p "$MPLCONFIGDIR"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="$SCRIPT_DIR/corcoran_node_filtration_pipeline.py"
RESULTS_DIR="$SCRIPT_DIR/results/experiment1"

NORTH_AMERICA_INPUT="/Users/gregorioperez/Documents/Maestria/Investigacion/mde_dijsktra/experiment_bounding_city/graphs/osmnx/norte_/America"
SOUTH_AMERICA_INPUT="/Users/gregorioperez/Documents/Maestria/Investigacion/mde_dijsktra/experiment_bounding_city/graphs/osmnx/sur_/America"
SOUTH_ASIA_INPUT="/Users/gregorioperez/Documents/Maestria/Investigacion/mde_dijsktra/experiment_bounding_city/graphs/osmnx/sur_/Asia"

"$PYTHON_BIN" "$PIPELINE" \
  "$NORTH_AMERICA_INPUT" \
  --cities "$SCRIPT_DIR/north_america_cities.txt" \
  --out "$RESULTS_DIR/north_america"

"$PYTHON_BIN" "$PIPELINE" \
  "$SOUTH_AMERICA_INPUT" \
  --cities "$SCRIPT_DIR/south_america_cities.txt" \
  --out "$RESULTS_DIR/south_america"

"$PYTHON_BIN" "$PIPELINE" \
  "$SOUTH_ASIA_INPUT" \
  --cities "$SCRIPT_DIR/south_asia_cities.txt" \
  --out "$RESULTS_DIR/south_asia"
