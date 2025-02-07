######################
## Creates the beacon friendly format files with Cohort Information 
######################

###############
## rm previous files
##############


rm analysis.json
rm runs.json
rm biosamples.json
rm datasets.json
rm cohort.json
rm individuals.json


#############
## Individual 
############

#!/bin/bash

# Input file containing the 20K entries
INPUT_FILE="has_vcfs.txt"
# Output JSON file
OUTPUT_FILE="individuals.json"

# Function to determine gender
get_gender() {
    local gender_char="${1: -1}"  # Extract last character

    case "$gender_char" in
        f) echo "female" ;;
        m) echo "male" ;;
        *) echo "unknown" ;;
    esac
}

# Get current year dynamically
CURRENT_YEAR=$(date +"%Y")

# Start JSON array
echo "[" > "$OUTPUT_FILE"

# Process each line
FIRST_ENTRY=true
while IFS= read -r line; do
    IFS='-' read -ra FIELDS <<< "$line"  # Split line by '-'

    INDIVIDUAL_ID="${FIELDS[0]}"           # Extract the individual ID (e.g., 99whagref)
    BIRTH_YEAR_PREFIX="${INDIVIDUAL_ID:0:2}" # Extract first two digits (e.g., 99 or 09)

    # Ensure proper numeric conversion
    if (( 10#$BIRTH_YEAR_PREFIX >= 30 )); then
        BIRTH_YEAR="19$BIRTH_YEAR_PREFIX"  # Assume 1930-1999
    else
        BIRTH_YEAR="20$BIRTH_YEAR_PREFIX"  # Assume 2000-2029
    fi

    AGE=$((CURRENT_YEAR - BIRTH_YEAR))    # Calculate age
    SEX=$(get_gender "$INDIVIDUAL_ID")    # Determine gender

    # Append JSON entry
    [[ $FIRST_ENTRY == false ]] && echo "," >> "$OUTPUT_FILE"
    FIRST_ENTRY=false

    cat <<EOF >> "$OUTPUT_FILE"
    {
        "id": "$line",
        "age": $AGE,
        "sex": "$SEX",
        "ethnicity": "European",
        "geographicOrigin": "Denmark",
        "datasets": ["DenGen"]
    }
EOF

done < "$INPUT_FILE"

# Close JSON array
echo "]" >> "$OUTPUT_FILE"

echo "JSON output saved to $OUTPUT_FILE"


################
## Dataset
################


# Dataset information
DATASET_ID="DenGen"
DATASET_NAME="DenGen Dataset"
DATASET_TYPE="Whole Genome Sequencing"
DATASET_DESCRIPTION="Genomic data from the Danish population"
DATASET_VERSION="v1.0"
CREATED_DATE="2025-01-01"

# Output file
OUTPUT_FILE="datasets.json"

# Create the dataset JSON
cat <<EOF > "$OUTPUT_FILE"
{
    "id": "$DATASET_ID",
    "name": "$DATASET_NAME",
    "type": "$DATASET_TYPE",
    "description": "$DATASET_DESCRIPTION",
    "version": "$DATASET_VERSION",
    "created": "$CREATED_DATE"
}
EOF

echo "Dataset model saved to $OUTPUT_FILE"


#####################
## Cohort 
######################

# Cohort information
COHORT_ID="DenGen"
COHORT_NAME="DenGen Cohort"
COHORT_DESCRIPTION="A cohort of individuals from the DenGen dataset"
CREATED_DATE="2025-02-07"

# Input file containing sample IDs (e.g., the list of individuals for this cohort)
INPUT_FILE="has_vcfs.txt"
# Output file
OUTPUT_FILE="cohort.json"

# Start the JSON object
echo "{" > "$OUTPUT_FILE"
echo "  \"id\": \"$COHORT_ID\"," >> "$OUTPUT_FILE"
echo "  \"name\": \"$COHORT_NAME\"," >> "$OUTPUT_FILE"
echo "  \"description\": \"$COHORT_DESCRIPTION\"," >> "$OUTPUT_FILE"
echo "  \"individuals\": [" >> "$OUTPUT_FILE"

# Flag for formatting JSON correctly
FIRST_ENTRY=true

# Process each line in the input file and add individual sample IDs to the JSON
while IFS= read -r line; do
    [[ $FIRST_ENTRY == false ]] && echo "," >> "$OUTPUT_FILE"
    FIRST_ENTRY=false
    echo "    \"$line\"" >> "$OUTPUT_FILE"
done < "$INPUT_FILE"

# Finish the JSON structure
echo "  ]," >> "$OUTPUT_FILE"
echo "  \"datasets\": [\"DenGen\"]," >> "$OUTPUT_FILE"
echo "  \"created\": \"$CREATED_DATE\"" >> "$OUTPUT_FILE"
echo "}" >> "$OUTPUT_FILE"

echo "Cohort model saved to $OUTPUT_FILE"

#############
## Runs
##############



#!/bin/bash

# Define the input file that contains the sample names
SAMPLES_FILE="has_vcfs.txt"

# Output file for the run models
OUTPUT_FILE="runs.json"

# Initialize the JSON array
echo "[" > "$OUTPUT_FILE"

# Read each sample name from the samples file
while read -r SAMPLE_NAME; do
    # Define the static values for the run model
    DATASET_ID="DenGen"
    BIOSAMPLE_ID="$SAMPLE_NAME"
    RUN_ID="$SAMPLE_NAME"  # Use the same as sample name
    INDIVIDUAL_ID="$SAMPLE_NAME"  # Use the same as sample name
    LIBRARY_LAYOUT="PAIRED"
    LIBRARY_SELECTION="RANDOM"
    LIBRARY_SOURCE_ID="GENEPIO:0001966"
    LIBRARY_SOURCE_LABEL="genomic Source"
    LIBRARY_STRATEGY="Short read Illumina WGS sequencing"
    PLATFORM="Illumina"
    PLATFORM_MODEL_ID="OBI:0002048"
    PLATFORM_MODEL_LABEL="Illumina NovaSeq6000"
    RUN_DATE="2025-01-01"  # Adjust date if needed

    # Construct the JSON object for the current run
    echo "  {" >> "$OUTPUT_FILE"
    echo "    \"datasetId\": \"$DATASET_ID\"," >> "$OUTPUT_FILE"
    echo "    \"biosampleId\": \"$BIOSAMPLE_ID\"," >> "$OUTPUT_FILE"
    echo "    \"id\": \"$RUN_ID\"," >> "$OUTPUT_FILE"
    echo "    \"individualId\": \"$INDIVIDUAL_ID\"," >> "$OUTPUT_FILE"
    echo "    \"libraryLayout\": \"$LIBRARY_LAYOUT\"," >> "$OUTPUT_FILE"
    echo "    \"librarySelection\": \"$LIBRARY_SELECTION\"," >> "$OUTPUT_FILE"
    echo "    \"librarySource\": {" >> "$OUTPUT_FILE"
    echo "      \"id\": \"$LIBRARY_SOURCE_ID\"," >> "$OUTPUT_FILE"
    echo "      \"label\": \"$LIBRARY_SOURCE_LABEL\"" >> "$OUTPUT_FILE"
    echo "    }," >> "$OUTPUT_FILE"
    echo "    \"libraryStrategy\": \"$LIBRARY_STRATEGY\"," >> "$OUTPUT_FILE"
    echo "    \"platform\": \"$PLATFORM\"," >> "$OUTPUT_FILE"
    echo "    \"platformModel\": {" >> "$OUTPUT_FILE"
    echo "      \"id\": \"$PLATFORM_MODEL_ID\"," >> "$OUTPUT_FILE"
    echo "      \"label\": \"$PLATFORM_MODEL_LABEL\"" >> "$OUTPUT_FILE"
    echo "    }," >> "$OUTPUT_FILE"
    echo "    \"runDate\": \"$RUN_DATE\"" >> "$OUTPUT_FILE"
    echo "  }," >> "$OUTPUT_FILE"
done < "$SAMPLES_FILE"

# Remove the last comma to make the JSON valid
sed -i '$ s/,$//' "$OUTPUT_FILE"

# Close the JSON array
echo "]" >> "$OUTPUT_FILE"

echo "Run models saved to $OUTPUT_FILE"

###############
## Analysis
###############


#!/bin/bash

# Define the input file that contains the sample names
SAMPLES_FILE="has_vcfs.txt"

# Output file for the analysis models
OUTPUT_FILE="analysis.json"

# Initialize the JSON array
echo "[" > "$OUTPUT_FILE"

# Read each sample name from the samples file
while read -r SAMPLE_NAME; do
    # Define the static values for the analysis model
    DATASET_ID="DenGen"
    ALIGNER="bwa-0.5.9"
    ANALYSIS_DATE="2025-01-01"  # Adjust date if needed
    PIPELINE_NAME="1000G-low-coverage-WGS"
    PIPELINE_REF="https://www.nature.com/articles/nature15393"
    VARIANT_CALLER="GATK4.0"

    # Construct the JSON object for the current analysis
    echo "  {" >> "$OUTPUT_FILE"
    echo "    \"datasetId\": \"$DATASET_ID\"," >> "$OUTPUT_FILE"
    echo "    \"aligner\": \"$ALIGNER\"," >> "$OUTPUT_FILE"
    echo "    \"analysisDate\": \"$ANALYSIS_DATE\"," >> "$OUTPUT_FILE"
    echo "    \"biosampleId\": \"$SAMPLE_NAME\"," >> "$OUTPUT_FILE"
    echo "    \"id\": \"$SAMPLE_NAME\"," >> "$OUTPUT_FILE"
    echo "    \"individualId\": \"$SAMPLE_NAME\"," >> "$OUTPUT_FILE"
    echo "    \"pipelineName\": \"$PIPELINE_NAME\"," >> "$OUTPUT_FILE"
    echo "    \"pipelineRef\": \"$PIPELINE_REF\"," >> "$OUTPUT_FILE"
    echo "    \"runId\": \"$SAMPLE_NAME\"," >> "$OUTPUT_FILE"
    echo "    \"variantCaller\": \"$VARIANT_CALLER\"" >> "$OUTPUT_FILE"
    echo "  }," >> "$OUTPUT_FILE"
done < "$SAMPLES_FILE"

# Remove the last comma to make the JSON valid
sed -i '$ s/,$//' "$OUTPUT_FILE"

# Close the JSON array
echo "]" >> "$OUTPUT_FILE"

echo "Analysis models saved to $OUTPUT_FILE"

############
## Biosamples
############

#!/bin/bash

# Define the input file that contains the sample names
SAMPLES_FILE="has_vcfs.txt"

# Output file for the biosample models
OUTPUT_FILE="biosamples.json"

# Initialize the JSON array
echo "[" > "$OUTPUT_FILE"

# Define the static values for the biosample model
DATASET_ID="DenGen"
BIOSAMPLE_STATUS_ID="EFO:0009654"
BIOSAMPLE_STATUS_LABEL="reference sample"
COLLECTION_DATE="2020-09-11"  # Adjust the date if needed
COLLECTION_MOMENT="P40Y1M1D"
TAX_ID=9606
SAMPLE_ORIGIN_ID="UBERON:0000178"
SAMPLE_ORIGIN_LABEL="blood"

# Read each sample name from the samples file
while read -r SAMPLE_NAME; do
    # Construct the JSON object for the current biosample
    echo "  {" >> "$OUTPUT_FILE"
    echo "    \"datasetId\": \"$DATASET_ID\"," >> "$OUTPUT_FILE"
    echo "    \"biosampleStatus\": {" >> "$OUTPUT_FILE"
    echo "      \"id\": \"$BIOSAMPLE_STATUS_ID\"," >> "$OUTPUT_FILE"
    echo "      \"label\": \"$BIOSAMPLE_STATUS_LABEL\"" >> "$OUTPUT_FILE"
    echo "    }," >> "$OUTPUT_FILE"
    echo "    \"collectionDate\": \"$COLLECTION_DATE\"," >> "$OUTPUT_FILE"
    echo "    \"collectionMoment\": \"$COLLECTION_MOMENT\"," >> "$OUTPUT_FILE"
    echo "    \"id\": \"$SAMPLE_NAME\"," >> "$OUTPUT_FILE"
    echo "    \"individualId\": \"$SAMPLE_NAME\"," >> "$OUTPUT_FILE"
    echo "    \"info\": {" >> "$OUTPUT_FILE"
    echo "      \"characteristics\": [" >> "$OUTPUT_FILE"
    echo "        {" >> "$OUTPUT_FILE"
    echo "          \"organism\": [" >> "$OUTPUT_FILE"
    echo "            {" >> "$OUTPUT_FILE"
    echo "              \"ontologyTerms\": [" >> "$OUTPUT_FILE"
    echo "                \"http://purl.obolibrary.org/obo/NCBITaxon_9606\"" >> "$OUTPUT_FILE"
    echo "              ]," >> "$OUTPUT_FILE"
    echo "              \"text\": \"Homosapiens\"" >> "$OUTPUT_FILE"
    echo "            }" >> "$OUTPUT_FILE"
    echo "          ]" >> "$OUTPUT_FILE"
    echo "        }" >> "$OUTPUT_FILE"
    echo "      ]," >> "$OUTPUT_FILE"
    echo "      \"sampleName\": \"$SAMPLE_NAME\"," >> "$OUTPUT_FILE"
    echo "      \"taxId\": \"$TAX_ID\"" >> "$OUTPUT_FILE"
    echo "    }," >> "$OUTPUT_FILE"
    echo "    \"sampleOriginType\": {" >> "$OUTPUT_FILE"
    echo "      \"id\": \"$SAMPLE_ORIGIN_ID\"," >> "$OUTPUT_FILE"
    echo "      \"label\": \"$SAMPLE_ORIGIN_LABEL\"" >> "$OUTPUT_FILE"
    echo "    }" >> "$OUTPUT_FILE"
    echo "  }," >> "$OUTPUT_FILE"
done < "$SAMPLES_FILE"

# Remove the last comma to make the JSON valid
sed -i '$ s/,$//' "$OUTPUT_FILE"

# Close the JSON array
echo "]" >> "$OUTPUT_FILE"

echo "Biosample models saved to $OUTPUT_FILE"


