# Processed FastQ QC Workflow

This repository contains a Snakemake workflow for processing and quality control (QC) of FastQ files from sequencing experiments. The workflow is modular, highly configurable, and supports both single-end (SE) and paired-end (PE) data. It integrates several bioinformatics tools and custom scripts for preprocessing, QC, species detection, and reporting.

## Features
- **Support for UMI:** Handles Unique Molecular Identifiers (UMIs) if present.
- **Species Detection:** Optional module for species detection using BioBloom and BLAST.
- **MultiQC Reporting:** Generates a comprehensive HTML report summarizing QC metrics.
- **File Sending:** Optionally packages results for sending via FileSender.

## Workflow Overview
1. **FastQ Preparation:**
   - Prepares raw FastQ files (SE or PE) and processes UMIs if present.
2. **Trimming:**
   - Trims adapters and low-quality bases using parameters from the config file.
3. **Quality Control:**
   - Runs FastQC on processed FastQ files and merges QC results.
4. **Species Detection (Optional):**
   - Detects species using BioBloom and BLAST, summarizes results.
5. **Reporting:**
   - Generates a MultiQC HTML report.
6. **File Sending (Optional):**
   - Packages results for transfer.

## Directory Structure
- `Snakefile` — Main workflow file.
- `workflow.config.json` — Configuration file.
- `rules/` — Snakemake rule files for each workflow step.
- `wrappers/` — Scripts and conda environments for each step.
- `qc_reports/` — Output directory for QC results and reports.
- `processed_fastq/` — Output directory for processed FastQ files.
- `logs/` — Log files for each step.

## Usage
1. **Configure the workflow:**
   - Edit `workflow.config.json` to set paths, parameters, and options.
2. **Run the workflow:**
   ```bash
   snakemake --use-conda --cores <N>
   ```
   Replace `<N>` with the number of CPU cores to use.
3. **View results:**
   - Final QC report: `qc_reports/processed_fastq_multiqc.html`
   - (Optional) Packaged results: `sequencing_results.tar.gz`

## Requirements
- [Snakemake](https://snakemake.readthedocs.io/)
- [Conda](https://docs.conda.io/)
- Python 3.6+

## Customization
- Modify rule files in `rules/` to add or change workflow steps.
- Update wrapper scripts in `wrappers/` for custom processing.
- Adjust conda environment YAMLs for tool versions.

## Contact
For questions or contributions, please contact the BioIT-CEITEC team.
