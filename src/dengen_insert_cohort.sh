#!/bin/bash
set -euo pipefail

# =====================================
# DenGen: Insert JSON Metadata into Beacon MongoDB
# =====================================

BASE_DIR=$(dirname "$0")/..
RESULTS_DIR="$BASE_DIR/results"

# -----------------------------
# Configurable parameters
# -----------------------------
MONGO_CONTAINER="${MONGO_CONTAINER:-mongoprod}"
BEACON_CONTAINER="${BEACON_CONTAINER:-beaconprod}"
MONGO_URI="${MONGO_URI:-mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin}"

# JSON files and their target collections
declare -A JSON_COLLECTIONS=(
    ["datasets.json"]="datasets"
    ["analyses.json"]="analyses"
    ["biosamples.json"]="biosamples"
    ["cohorts.json"]="cohorts"
    ["individuals.json"]="individuals"
    ["runs.json"]="runs"
    ["genomicVariations.json"]="genomicVariations"
)

# -----------------------------
# Copy JSON files to Mongo container
# -----------------------------
echo "[INFO] Copying JSON files to Mongo container ($MONGO_CONTAINER)..."
for file in "${!JSON_COLLECTIONS[@]}"; do
    src="$RESULTS_DIR/$file"
    if [ -f "$src" ]; then
        docker cp "$src" "$MONGO_CONTAINER:/tmp/$file"
        echo "[INFO] Copied $file"
    else
        echo "[WARNING] File not found: $src, skipping."
    fi
done

# -----------------------------
# Import JSON files into MongoDB
# -----------------------------
echo "[INFO] Importing JSON files into MongoDB..."
for file in "${!JSON_COLLECTIONS[@]}"; do
    src="/tmp/$file"
    collection="${JSON_COLLECTIONS[$file]}"
    docker exec "$MONGO_CONTAINER" mongoimport \
        --jsonArray --uri "$MONGO_URI" \
        --file "$src" --collection "$collection"
    echo "[INFO] Imported $file -> $collection"
done

# -----------------------------
# Reindex MongoDB in Beacon container
# -----------------------------
echo "[INFO] Reindexing Beacon MongoDB..."
docker exec "$BEACON_CONTAINER" python /beacon/connections/mongo/reindex.py

# -----------------------------
# Extract filtering terms
# -----------------------------
echo "[INFO] Extracting filtering terms..."
docker exec "$BEACON_CONTAINER" python beacon/connections/mongo/extract_filtering_terms.py

# -----------------------------
# Done
# -----------------------------
echo "[INFO] JSON data insertion pipeline complete."

