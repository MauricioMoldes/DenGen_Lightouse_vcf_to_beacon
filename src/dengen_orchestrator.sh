#!/bin/bash
set -euo pipefail

# =========================================
# DenGen in a Box: Full Beacon Workflow
# =========================================

BASE_DIR=$(dirname "$0")/..
RESULTS_DIR="$BASE_DIR/results"
SRC_DIR="$BASE_DIR/src"
DATA_DIR="$BASE_DIR/data"
VCF_DIR="$RESULTS_DIR"

mkdir -p "$RESULTS_DIR"

# -----------------------------
# Step 1: Prepare Beacon-ready VCF
# -----------------------------
VCF_INPUT="$VCF_DIR/dengen_2211_merged.v2.tags.vcf.gz"
VCF_FINAL="$VCF_DIR/dengen_SNP_filtered.decomposed.renamed.vcf.gz"

if [ ! -f "$VCF_FINAL" ]; then
    echo "[STEP 1] Running DenGen VCF pipeline..."
    bash "$SRC_DIR/dengen_merged_pipeline.sh"
else
    echo "[STEP 1] Final VCF already exists: $VCF_FINAL, skipping pipeline."
fi

# -----------------------------
# Step 2: Generate JSON cohort metadata
# -----------------------------
echo "[STEP 2] Generating cohort JSONs..."
bash "$SRC_DIR/dengen_create_cohort.sh"

# -----------------------------
# Step 3: Insert JSON into Beacon MongoDB
# -----------------------------
echo "[STEP 3] Inserting cohort JSONs into MongoDB..."
bash "$SRC_DIR/dengen_insert_cohort.sh"

