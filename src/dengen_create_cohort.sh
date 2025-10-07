#!/bin/bash
set -euo pipefail

# =====================================
# DenGen: Create Beacon JSON Metadata
# =====================================

# -----------------------------
# Configurable inputs
# -----------------------------
BASE_DIR=$(dirname "$0")/..
DATA_DIR="$BASE_DIR/data"
RESULTS_DIR="$BASE_DIR/results"
mkdir -p "$RESULTS_DIR"

INPUT_FILE="${DATA_DIR}/dengen_list_of_samples.txt"
COHORT_NAME="DenGen"
DATASET_NAME="DenGen Dataset"
DATASET_ID="DenGen"
CREATED_DATE=$(date +"%Y-%m-%d")
CURRENT_YEAR=$(date +"%Y")
CHR_CONV="$BASE_DIR/src/chr_name_conv.txt"

# JSON output files
INDIVIDUALS_JSON="$RESULTS_DIR/individuals.json"
DATASETS_JSON="$RESULTS_DIR/datasets.json"
COHORTS_JSON="$RESULTS_DIR/cohorts.json"
RUNS_JSON="$RESULTS_DIR/runs.json"
ANALYSES_JSON="$RESULTS_DIR/analyses.json"
BIOSAMPLES_JSON="$RESULTS_DIR/biosamples.json"

# Remove previous outputs safely
rm -f "$INDIVIDUALS_JSON" "$DATASETS_JSON" "$COHORTS_JSON" "$RUNS_JSON" "$ANALYSES_JSON" "$BIOSAMPLES_JSON"

# -----------------------------
# Function to determine gender
# -----------------------------
get_gender() {
    local gender_char="${1: -1}"
    case "$gender_char" in
        f) echo "female" ;;
        m) echo "male" ;;
        *) echo "unknown" ;;
    esac
}

# -----------------------------
# Individuals JSON
# -----------------------------
echo "[" > "$INDIVIDUALS_JSON"
FIRST_ENTRY=true
while IFS= read -r line; do
    [[ $FIRST_ENTRY == false ]] && echo "," >> "$INDIVIDUALS_JSON"
    FIRST_ENTRY=false

    # Parse ID and birth year
    IFS='-' read -ra FIELDS <<< "$line"
    INDIVIDUAL_ID="${FIELDS[0]}"
    BIRTH_YEAR_PREFIX="${INDIVIDUAL_ID:0:2}"
    if (( 10#$BIRTH_YEAR_PREFIX >= 30 )); then
        BIRTH_YEAR="19$BIRTH_YEAR_PREFIX"
    else
        BIRTH_YEAR="20$BIRTH_YEAR_PREFIX"
    fi
    AGE=$((CURRENT_YEAR - BIRTH_YEAR))
    SEX=$(get_gender "$INDIVIDUAL_ID")

    cat <<EOF >> "$INDIVIDUALS_JSON"
    {
        "id": "$line",
        "sex": {"id": "NCIT:C20197","label": "$SEX"},
        "ethnicity": {"id": "NCIT:C42331","label": "European"},
        "geographicOrigin": {"id":"NCIT:C16496","label":"Denmark"},
        "measures": [
            {"assayCode":{"id":"LOINC:35925-4","label":"BMI"},"date":"2021-09-24","measurementValue":{"unit":{"id":"NCIT:C49671","label":"Kilogram per Square Meter"},"value":26.638}},
            {"assayCode":{"id":"LOINC:3141-9","label":"Weight"},"date":"2021-09-24","measurementValue":{"unit":{"id":"NCIT:C28252","label":"Kilogram"},"value":85.6358}},
            {"assayCode":{"id":"LOINC:8308-9","label":"Height-standing"},"date":"2021-09-24","measurementValue":{"unit":{"id":"NCIT:C49668","label":"Centimeter"},"value":179.2973}}
        ]
    }
EOF
done < "$INPUT_FILE"
echo "]" >> "$INDIVIDUALS_JSON"
echo "[INFO] Individuals JSON saved to $INDIVIDUALS_JSON"

# -----------------------------
# Dataset JSON
# -----------------------------
cat <<EOF > "$DATASETS_JSON"
[{
    "id": "$DATASET_ID",
    "name": "$DATASET_NAME",
    "type": "Whole Genome Sequencing",
    "description": "Genomic data from the Danish population",
    "version": "v1.0",
    "created": "$CREATED_DATE"
}]
EOF
echo "[INFO] Dataset JSON saved to $DATASETS_JSON"

# -----------------------------
# Cohorts JSON
# -----------------------------
cat <<EOF > "$COHORTS_JSON"
[
    {
        "cohortType": "study-defined",
        "id": "$COHORT_NAME",
        "name": "$COHORT_NAME",
        "collectionEvents": [
            {
                "eventEthnicities": {"availability": true,"availabilityCount": 2211,"distribution":{"ethnicities":{"European":2211}}},
                "eventGenders": {"availability": true,"availabilityCount": 2211,"distribution":{"genders":{"female":1165,"male":1046}}}
            }
        ],
        "inclusionCriteria": {
            "ageRange": {"start":{"iso8601duration":"P0Y"},"end":{"iso8601duration":"P94Y"}},
            "genders":[{"id":"NCIT:C16576","label":"female"},{"id":"NCIT:C20197","label":"male"}],
            "locations":[{"id":"GAZ:00002635","label":"Denmark"}]
        }
    }
]
EOF
echo "[INFO] Cohorts JSON saved to $COHORTS_JSON"

# -----------------------------
# Runs JSON
# -----------------------------
echo "[" > "$RUNS_JSON"
FIRST_ENTRY=true
while read -r SAMPLE; do
    [[ $FIRST_ENTRY == false ]] && echo "," >> "$RUNS_JSON"
    FIRST_ENTRY=false
    cat <<EOF >> "$RUNS_JSON"
    {
        "datasetId": "$DATASET_ID",
        "biosampleId": "$SAMPLE",
        "id": "$SAMPLE",
        "individualId": "$SAMPLE",
        "libraryLayout": "PAIRED",
        "librarySelection": "RANDOM",
        "librarySource": {"id":"GENEPIO:0001966","label":"genomic Source"},
        "libraryStrategy": "Short read Illumina WGS sequencing",
        "platform": "Illumina",
        "platformModel": {"id":"OBI:0002048","label":"Illumina NovaSeq6000"},
        "runDate": "$CREATED_DATE"
    }
EOF
done < "$INPUT_FILE"
echo "]" >> "$RUNS_JSON"
echo "[INFO] Runs JSON saved to $RUNS_JSON"

# -----------------------------
# Analyses JSON
# -----------------------------
echo "[" > "$ANALYSES_JSON"
FIRST_ENTRY=true
while read -r SAMPLE; do
    [[ $FIRST_ENTRY == false ]] && echo "," >> "$ANALYSES_JSON"
    FIRST_ENTRY=false
    cat <<EOF >> "$ANALYSES_JSON"
    {
        "datasetId": "$DATASET_ID",
        "aligner": "bwa-0.5.9",
        "analysisDate": "$CREATED_DATE",
        "biosampleId": "$SAMPLE",
        "id": "$SAMPLE",
        "individualId": "$SAMPLE",
        "pipelineName": "1000G-low-coverage-WGS",
        "pipelineRef": "https://www.nature.com/articles/nature15393",
        "runId": "$SAMPLE",
        "variantCaller": "GATK4.0"
    }
EOF
done < "$INPUT_FILE"
echo "]" >> "$ANALYSES_JSON"
echo "[INFO] Analyses JSON saved to $ANALYSES_JSON"

# -----------------------------
# Biosamples JSON
# -----------------------------
echo "[" > "$BIOSAMPLES_JSON"
FIRST_ENTRY=true
while read -r SAMPLE; do
    [[ $FIRST_ENTRY == false ]] && echo "," >> "$BIOSAMPLES_JSON"
    FIRST_ENTRY=false
    cat <<EOF >> "$BIOSAMPLES_JSON"
    {
        "datasetId": "$DATASET_ID",
        "biosampleStatus":{"id":"EFO:0009654","label":"reference sample"},
        "collectionDate": "2020-09-11",
        "collectionMoment": "P40Y1M1D",
        "id": "$SAMPLE",
        "individualId": "$SAMPLE",
        "info": {"sampleName":"$SAMPLE","taxId":9606,"characteristics":[{"organism":[{"ontologyTerms":["http://purl.obolibrary.org/obo/NCBITaxon_9606"],"text":"Homosapiens"}]}]},
        "sampleOriginType":{"id":"UBERON:0000178","label":"blood"}
    }
EOF
done < "$INPUT_FILE"
echo "]" >> "$BIOSAMPLES_JSON"
echo "[INFO] Biosamples JSON saved to $BIOSAMPLES_JSON"

echo "[INFO] Cohort creation pipeline complete."

