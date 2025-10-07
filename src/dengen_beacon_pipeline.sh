#!/bin/bash
set -euo pipefail

# =====================================
# DenGen Beacon Pipeline
# =====================================

BASE_DIR=$(dirname "$0")/..
RESULTS_DIR="$BASE_DIR/results"
SRC_DIR="$BASE_DIR/src"

# Input VCF: already merged, SNP filtered, decomposed, and renamed
INPUT_VCF="$RESULTS_DIR/dengen_SNP_filtered.decomposed.renamed.vcf.gz"

# Check that input exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "[ERROR] Input VCF not found: $INPUT_VCF"
    exit 1
fi

# -----------------------------
# Optional: Copy to ri-tools mountpoint
# -----------------------------
RI_TOOLS_DIR="$BASE_DIR/beacon2-pi-api/ri-tools/files/vcf/files_to_read"
mkdir -p "$RI_TOOLS_DIR"

echo "[INFO] Copying VCF to ri-tools mountpoint..."
cp "$INPUT_VCF" "$RI_TOOLS_DIR/"

# -----------------------------
# Run ri-tools Docker container
# -----------------------------
BEACON_API_DIR="$BASE_DIR/beacon2-pi-api"

echo "[INFO] Restarting beacon-ri-tools container..."
cd "$BEACON_API_DIR"
docker compose restart beacon-ri-tools

echo "[INFO] Executing genomicVariations_vcf.py in ri-tools..."
docker exec -it ri-tools python genomicVariations_vcf.py

# -----------------------------
# Cleanup
# -----------------------------
echo "[INFO] Cleaning up ri-tools mountpoint..."
rm "$RI_TOOLS_DIR/$(basename "$INPUT_VCF")"

echo "[INFO] Cleaning up temporary VCF files..."
DATA_DIR="$BASE_DIR/beacon2-pi-api/ri-tools/files/vcf/data"
rm -f "$DATA_DIR"/*.SNP_filtered.* "$DATA_DIR"/*.head "$DATA_DIR"/*.variants

echo "[INFO] Beacon pipeline complete."

