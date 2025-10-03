# DenGen Beacon Data Generator

This repository contains scripts and resources to generate Beacon-ready data from the DenGen dataset.
It takes raw VCF files and metadata from DenGen, processes them, and produces metadata in formats that can be ingested into a Beacon MongoDB database.

Note: This pipeline can be adapted for other datasets, making it reusable for any Beacon implementation.

---

## Table of Contents

1. Features
2. Repository Structure
3. Prerequisites
4. Usage
5. Scripts Overview
6. Data Flow
7. License
8. Contact

---

## Features

- Converts DenGen VCF and metadata into Beacon-compatible formats
- Generates cohort and variant metadata for MongoDB ingestion
- Bash-based pipeline for easy automation
- Reusable for other Beacon datasets with minor configuration changes

---

## Repository Structure

.
├── bin                     # Optional: helper executables or scripts
├── data                    # Raw input data
│   └── dengen_annonymous_list_samples.txt
├── docs                    # Documentation
├── README.md
├── results                 # Generated Beacon-ready output
└── src                     # Core scripts
    ├── annonymous_dengen.sh
    ├── chr_name_conv.txt
    ├── create_cohort.sh
    ├── dengen_merged_pipeline.sh
    └── insert_cohort.sh

---

## Prerequisites

- Unix/Linux environment
- Bash shell
- Access to DenGen VCF files and metadata
- MongoDB (for ingestion into Beacon)

Optional:

- Python if any scripts rely on it (check src/ scripts)

---

## Usage

### 1. Prepare your data

- Place your VCF files and metadata in the data/ folder.
- Ensure dengen_annonymous_list_samples.txt contains the list of anonymized sample IDs.

### 2. Run the pipeline

From the root directory:

bash src/dengen_merged_pipeline.sh

This will:

1. Anonymize metadata (annonymous_dengen.sh)
2. Convert chromosome names (chr_name_conv.txt)
3. Create the cohort (create_cohort.sh)
4. Insert the cohort into the Beacon database (insert_cohort.sh)

### 3. Inspect results

- Output files will be placed in the results/ directory.
- These files are ready for ingestion into a Beacon MongoDB instance.

---

## Scripts Overview

Script | Description
-------|-------------
annonymous_dengen.sh | Anonymizes DenGen sample metadata
chr_name_conv.txt | Chromosome name conversion mapping
create_cohort.sh | Generates cohort JSON/TSV files for Beacon
insert_cohort.sh | Inserts cohort into MongoDB for Beacon
dengen_merged_pipeline.sh | Master script to run the Beacon RI tools for the genomic_variants pipeline

---

## Data Flow

Raw VCF + Metadata --> src/ scripts --> results/ (Beacon-ready files) --> MongoDB (Beacon)

- The pipeline ensures that metadata is anonymized, chromosome names are standardized, and cohort/variant files conform to the Beacon schema.

---

## License

MIT

---

## Contact

For questions, suggestions, or issues:

- Mauricio Moldes
- GitHub Issues: https://github.com/MauricioMoldes/dengen_beacon_data_generator/issues

