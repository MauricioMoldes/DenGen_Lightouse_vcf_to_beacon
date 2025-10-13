# DenGen Beacon Data Generator

This project provides a fully automated pipeline to generate Beacon-ready data from the DenGen dataset. It takes raw VCF files and metadata from DenGen, processes them, and produces JSON and VCF files that can be ingested into a Beacon MongoDB database.

The pipeline is reusable and can be adapted for other genomic datasets with minimal configuration.

---

## Table of Contents

1. Features
2. Prerequisites
3. Usage
4. Scripts Overview
5. Data Flow
6. License
7. Contact

---

## Features

- Converts DenGen VCFs and metadata into Beacon-compatible formats
- Generates cohort and variant metadata for MongoDB ingestion
- Fully automated “DenGen in a Box” workflow
- Supports large-scale VCF normalization, merging, decomposition, and chromosome renaming
- Portable and reusable for other datasets
- Provides logging and safe execution with idempotency (skips steps if output already exists)

---

## Prerequisites

- Unix/Linux environment with Bash
- Access to DenGen VCF files and metadata
- MongoDB for ingestion into Beacon
- Optional: HPC/cluster environment for large-scale VCF processing

Required software/tools (modules or system installs):

- bcftools
- vt (variant tools)
- tabix
- docker and docker-compose

---

## Usage

### 1. Prepare your data

- Place your VCF files and metadata in the data/ folder.
- Ensure dengen_annonymous_list_samples.txt contains the list of anonymized sample IDs.

### 2. Run the full pipeline

From the root directory:

bash src/dengen_beacon_pipeline.sh


This orchestrator will:

1. Normalize and merge individual VCFs into a merged DenGen VCF
2. Filter SNPs, decompose multiallelics, and rename chromosomes
3. Generate cohort and metadata JSONs (individuals.json, biosamples.json, datasets.json, etc.)
4. Insert the metadata into the Beacon MongoDB
5. Push the final Beacon-ready VCF to the ri-tools container for ingestion

### 3. Inspect results

- All output JSON and VCF files are placed in the results/ directory.
- These files are ready for ingestion into a Beacon MongoDB instance.

---

## Scripts Overview

| Script | Description |
|--------|-------------|
| dengen_beacon_pipeline.sh | Orchestrator that runs the full workflow from raw VCFs to Beacon ingestion |
| dengen_create_cohort.sh | Generates cohort JSONs and metadata for Beacon |
| dengen_insert_cohort.sh | Inserts cohort and metadata JSONs into Beacon MongoDB |
| dengen_merged_pipeline.sh | Processes the VCF: SNP filtering, decomposition, renaming, and preparation for Beacon ingestion |
| chr_name_conv.txt | Chromosome name conversion mapping used in VCF processing |

---

## Data Flow

Raw VCF + Metadata
↓ (VCF normalization & merging)
Merged & Beacon-ready VCF
↓ (JSON generation)
Cohort and variant metadata
↓ (MongoDB ingestion)
Beacon MongoDB
↓
ri-tools container ingestion

---

## License

MIT License

---

## Contact

For questions, suggestions, or issues:

- Mauricio Moldes
- GitHub Issues: https://github.com/MauricioMoldes/DenGen_Lightouse/issues


