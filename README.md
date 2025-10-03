# DenGen Beacon Data Generator

This repository contains scripts and resources to generate Beacon-ready data from the DenGen dataset.
It takes raw VCF files and metadata from DenGen, processes them, and produces metadata in formats that can be ingested into a Beacon MongoDB database.

Note: This pipeline can be adapted for other datasets, making it reusable for any Beacon implementation.

--------------------------------------------------------------------------------

## Table of Contents

1. Features
2. Repository Structure
3. Prerequisites
4. Usage
5. Scripts Overview
6. Data Flow
7. License
8. Contact

--------------------------------------------------------------------------------

## Features

- Converts DenGen VCF and metadata into Beacon-compatible formats
- Generates cohort and variant metadata for MongoDB ingestion
- Bash-based pipeline for easy automation
- Includes cluster-side utilities for large-scale VCF manipulation
- Reusable for other Beacon datasets with minor configuration changes

--------------------------------------------------------------------------------

## Repository Structure

.
├── bin                     # Cluster-side scripts for VCF manipulation
│   ├── beacon_pipeline_anon.sh
│   ├── merge_dengen_vcf_anon.sh
│   └── normalize_vcf_snps_anon.sh
├── data                    # Raw input data
│   └── dengen_annonymous_list_samples.txt
├── docs                    # Documentation
├── README.md
├── results                 # Generated Beacon-ready output
└── src                     # Core scripts for metadata and cohort generation
    ├── chr_name_conv.txt
    ├── create_cohort.sh
    ├── dengen_merged_pipeline.sh
    └── insert_cohort.sh

--------------------------------------------------------------------------------

## Prerequisites

- Unix/Linux environment
- Bash shell
- Access to DenGen VCF files and metadata
- MongoDB (for ingestion into Beacon)

Optional:

- HPC/cluster environment for large-scale VCF manipulation (bin/ scripts)

--------------------------------------------------------------------------------

## Usage

### 1. Prepare your data

- Place your VCF files and metadata in the data/ folder.
- Ensure dengen_annonymous_list_samples.txt contains the list of anonymized sample IDs.

### 2. Run the pipeline

From the root directory:

bash src/dengen_merged_pipeline.sh

This will:

1. Normalize and merge VCFs (bin/ scripts)  
2. Convert chromosome names (chr_name_conv.txt)  
3. Create the cohort (create_cohort.sh)  
4. Insert the cohort into the Beacon database (insert_cohort.sh)  

### 3. Inspect results

- Output files will be placed in the results/ directory.  
- These files are ready for ingestion into a Beacon MongoDB instance.  

--------------------------------------------------------------------------------

## Scripts Overview

### Cluster-side (bin/)

Script | Description
-------|-------------
beacon_pipeline_anon.sh | Wrapper script to run the full VCF anonymization + normalization pipeline on cluster
merge_dengen_vcf_anon.sh | Merges multiple anonymized VCFs into a single dataset
normalize_vcf_snps_anon.sh | Normalizes SNPs in VCFs for Beacon compatibility

### Metadata & Cohort (src/)

Script | Description
-------|-------------
chr_name_conv.txt | Chromosome name conversion mapping
create_cohort.sh | Generates cohort JSON/TSV files for Beacon
insert_cohort.sh | Inserts cohort into MongoDB for Beacon
dengen_merged_pipeline.sh | Master script to run the full metadata/ingestion pipeline

--------------------------------------------------------------------------------

## Data Flow

Raw VCF + Metadata  
↓ (bin/ scripts on cluster)  
Normalized + Merged VCF  
↓ (src/ scripts)  
Beacon-ready metadata in results/  
↓  
MongoDB (Beacon ingestion)  

--------------------------------------------------------------------------------

## License

MIT

--------------------------------------------------------------------------------

## Contact

For questions, suggestions, or issues:

- Mauricio Moldes  
- GitHub Issues: https://github.com/MauricioMoldes/dengen_beacon_data_generator/issues

