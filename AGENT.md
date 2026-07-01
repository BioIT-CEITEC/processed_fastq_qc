# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a Snakemake workflow for processing and quality control of FastQ files from sequencing experiments. It supports both single-end (SE) and paired-end (PE) data, handles various UMI (Unique Molecular Identifier) formats, and integrates bioinformatics tools including FastQC, Cutadapt, BioBloom, and BLAST.

## Key Commands

```bash
# Run the full workflow
snakemake --use-conda --cores <N>

# Run with specific sample
snakemake --use-conda --cores <N> -s Snakefile --configfile workflow.config.json

# View available rules
snakemake --list-rules

# Dry run to check workflow without executing
snakemake --dryrun --use-conda --cores <N>

# Generate DAG visualization
snakemake --dag | dot -Tpdf > dag.pdf
```

## Architecture

### Main Components

1. **Snakefile** - Entry point that:
   - Imports the `bioroots_utilities` module from GitHub (BioIT-CEITEC/bioroots_utilities)
   - Loads sample configuration via `BR.load_sample()`
   - Sets up read pair tags for SE/PE handling
   - Loads UMI configuration from `workflow.config.json`
   - Includes rule modules from `rules/` directory

2. **rules/** - Modular Snakemake rule files:
   - `fastq_prepare.smk` - UMI extraction and fastq preparation (SE/PE rules)
   - `trimming.smk` - Adapter trimming via Cutadapt, FileSender integration
   - `fastqc.smk` - Quality control with FastQC and MultiQC merging
   - `species_detection.smk` - BioBloom and BLAST-based species detection

3. **wrappers/** - Custom scripts for each workflow step:
   - Python scripts (`.py`) for most operations
   - R scripts (`.R`) for species detection and merging
   - Each wrapper has its own conda environment (`env.yaml`)

4. **workflow.config.json** - Central configuration containing:
   - Workflow metadata and I/O patterns
   - UMI parameter presets for different library prep kits (CS_UMI, LYNX, BRONCO, Qiaseq, TruSight_Oncology, etc.)
   - Trimming parameters (adapters, quality thresholds, length filters)
   - Optional module toggles (species_detector, biobloom, filesender)

### Data Flow

```
raw_fastq/ → fastq_prepare (UMI extraction) → umi_fastq/
                                         ↓
                              preprocess (Cutadapt trimming) → processed_fastq/
                                         ↓
                          processed_fastq_qc (FastQC) → qc_reports/
                                         ↓
                    [optional] biobloom / species_detector → qc_reports/
                                         ↓
                              merge_fastq_qc → multiqc HTML report
```

### Key Configuration Parameters

- `UMI` - UMI type: `no_umi`, `custom_umi`, `CS_UMI`, `LYNX`, `BRONCO`, `Qiaseq`, `TruSight_Oncology`, etc.
- `is_paired` - Boolean for PE vs SE mode
- `trim_adapters` - Enable adapter trimming
- `quality_trim` - Quality trimming cutoffs (5',3')
- `biobloom` / `species_detector` - Enable optional detection modules
- `filesender` - Enable automated result packaging/sending

### Global Resources

The workflow references external resources defined in config:
- `globalResources` - Path to reference data (BioBloom databases, BLAST nt_db, taxdb)
- `globalTmpdPath` - Temporary directory for intermediate files

### Wrapper Script Patterns

- All wrapper scripts use `snakemake.shell` for command execution
- Logs are written to specified log files with command history
- Conda environments are declared per-rule via `conda:` directive
- Scripts receive configuration via `snakemake.params` and wildcards via `snakemake.wildcards`

## Development Notes

- The workflow depends on the `bioroots_utilities` Snakemake module - changes to sample loading or UMI configuration should account for this dependency
- UMI handling logic is complex and kit-specific; test thoroughly when modifying `fastq_prepare.smk` or its wrapper scripts
- Species detection requires pre-built BLAST databases at paths configured in `GLOBAL_REF_PATH`
- Output directories (`qc_reports/`, `processed_fastq/`, `logs/`) are created automatically by the workflow
