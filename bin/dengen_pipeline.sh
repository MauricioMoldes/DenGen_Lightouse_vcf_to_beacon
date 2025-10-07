#!/bin/bash
set -euo pipefail   # stop on errors, undefined vars

# ========================
# Paths (portable)
# ========================
BASE_DIR=$(dirname "$0")/..
DATA_DIR="$BASE_DIR/data"
SRC_DIR="$BASE_DIR/src"
RESULTS_DIR="$BASE_DIR/results"
NORMAL_DIR="$RESULTS_DIR/normalization"

mkdir -p "$NORMAL_DIR"

FASTA_REF="$DATA_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta"
VCF_LIST="$DATA_DIR/dengen_annonymous_list_samples.txt"

# ========================
# Step 1: Normalize each VCF
# ========================
echo "=== NORMALIZATION ==="
while IFS= read -r line; do
    echo "[NORMALIZE] $line"
    filename=$(basename -- "$line")
    filename="${filename%.*.*}"
    
    bcftools norm \
        --fasta-ref "$FASTA_REF" \
        --multiallelics - \
        --threads 8 \
        -Oz -o "$NORMAL_DIR/${filename}.normalized.vcf.gz" \
        "$line"
done < "$VCF_LIST"

# ========================
# Step 2: Merge normalized VCFs
# ========================
echo "=== MERGE ==="
NORMALIZED_LIST="$RESULTS_DIR/normalized.paths"
ls "$NORMAL_DIR"/*.vcf.gz > "$NORMALIZED_LIST"

MERGED_VCF="$RESULTS_DIR/dengen_merged.vcf.gz"
bcftools merge -l "$NORMALIZED_LIST" --missing-to-ref -m none -Oz -o "$MERGED_VCF"

# ========================
# Step 3: Beacon pipeline
# ========================
echo "=== BEACON PIPELINE ==="
FILTERED="$RESULTS_DIR/dengen_SNP_filtered.vcf.gz"
bcftools filter -i 'TYPE="snp"' "$MERGED_VCF" -Oz -o "$FILTERED"

DECOMPOSED="$RESULTS_DIR/dengen_SNP_filtered.decomposed.vcf.gz"
vt decompose "$FILTERED" -o "$DECOMPOSED"

RENAMED="$RESULTS_DIR/dengen_SNP_filtered.decomposed.renamed.vcf.gz"
bcftools annotate --rename-chrs "$SRC_DIR/chr_name_conv.txt" "$DECOMPOSED" -Oz -o "$RENAMED"
tabix -p vcf "$RENAMED"

echo "=== PIPELINE COMPLETE ==="
echo "Final Beacon-ready VCF: $RENAMED"

