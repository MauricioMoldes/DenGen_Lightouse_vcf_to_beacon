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
# Step 3: Annotate with population fields
# ========================
echo "=== POPULATION ANNOTATION ==="

ANNOTATED="$RESULTS_DIR/dengen_SNP_filtered.decomposed.renamed.annotated.vcf.gz"

# Fill standard tags: AC, AF, AN, plus AC_HOM, AC_HET
bcftools +fill-tags "$RENAMED" \
    -- -t AC,AF,AN,AC_HOM,AC_HET -Oz -o "$ANNOTATED"

tabix -p vcf "$ANNOTATED"

# Optional: if you have per-population groups defined in a file, you can annotate per-population AF
# For example, a file pop_file.txt with: sample_name<TAB>population
# bcftools +fill-tags "$RENAMED" -- -t AC,AF,AN,AC_HOM,AC_HET -Oz -o "$ANNOTATED"

echo "=== ANNOTATION COMPLETE ==="
echo "Final annotated VCF: $ANNOTATED"


# ========================
# Step 4: Beacon pipeline
# ========================
echo "=== BEACON PIPELINE ==="

# Read from the annotated VCF
BEACON_INPUT="$ANNOTATED"

# Optional: filter to SNPs only (if needed)
BEACON_FILTERED="$RESULTS_DIR/dengen_SNP_beacon.vcf.gz"
bcftools filter -i 'TYPE="snp"' "$BEACON_INPUT" -Oz -o "$BEACON_FILTERED"

# Decompose multi-allelics (if needed)
BEACON_DECOMPOSED="$RESULTS_DIR/dengen_SNP_beacon.decomposed.vcf.gz"
vt decompose "$BEACON_FILTERED" -o "$BEACON_DECOMPOSED"

# Rename chromosomes
BEACON_RENAMED="$RESULTS_DIR/dengen_SNP_beacon.decomposed.renamed.vcf.gz"
bcftools annotate --rename-chrs "$SRC_DIR/chr_name_conv.txt" "$BEACON_DECOMPOSED" -Oz -o "$BEACON_RENAMED"
tabix -p vcf "$BEACON_RENAMED"

echo "=== BEACON PIPELINE COMPLETE ==="
echo "Final Beacon-ready VCF: $BEACON_RENAMED"

