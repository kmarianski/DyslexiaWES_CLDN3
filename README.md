# Whole-Exome Sequencing (WES) Pipeline

This repository contains a modular and reproducible pipeline for Whole-Exome Sequencing (WES) analysis built using Snakemake.

## Overview

The pipeline includes the following steps:
1. **Quality Control**: Using Trimmomatic to trim raw sequencing reads.
2. **Alignment**: Mapping reads to a reference genome using BWA.
3. **Post-alignment Processing**:
   - Merging BAM files with Samtools.
   - Replacing read groups using Picard.
   - Marking duplicates.
   - Base recalibration with GATK.
   - BAM indexing.
4. **Variant Calling**: Using DeepVariant.
5. **Annotation**: Annotating variants with Annovar.
6. **Post-processing ANNOVAR Outputs**: Parsing, filtering, and merging per-sample ANNOVAR outputs.

## Installation

1. Clone this repository:
    ```bash
    git clone <repository_url>
    cd <repository_folder>
    ```

2. Install Snakemake and required tools:
    ```bash
    conda install -c bioconda -c conda-forge snakemake
    conda env create -f envs/biotools.yml
    ```

3. Ensure all required tools (e.g., BWA, Samtools, GATK, Annovar) are installed.

## Configuration

1. Update `config.yaml` with the following:
   - Path to the reference genome (`reference`).
   - Known sites for base recalibration (`known_sites`).
   - Adapter sequences for Trimmomatic (`adapter`).
   - Annovar database path (`annovar_db`).

2. Place raw sequencing data in the `raw/` directory. Files should be named as `<sample>_1.fq.gz` and `<sample>_2.fq.gz`.

## Processing ANNOVAR Outputs

The script `parse_ANNOVAR_refactored.py` parses, filters, and combines per-sample ANNOVAR outputs into a single DataFrame for downstream analysis.

### Usage

1. Ensure that all ANNOVAR output files (e.g., `multianno.csv`) are stored in the directory specified in the script (default: `results/annovar/multianno_csv/`).
2. Run the script:
    ```bash
    python process_ANNOVAR_files_refactored.py
    ```
3. Output will be saved as a single file (`data/by_variant.csv`) containing:
    - Sample ID
    - Genomic position (`chr:pos`)
    - Filtered variants based on pathogenicity and minor allele frequency (MAF).

### Key Features
- **Pathogenicity Filtering**: Filters variants based on annotation scores (e.g., SIFT, PolyPhen2, MutationTaster).
- **MAF Filtering**: Retains variants with MAF < 0.01.
- **Genomic Region Filtering**: Excludes intronic, intergenic, and synonymous variants.
- **Quality Filtering**: Retains variants with QS ≥ 30 and RD ≥ 10.

## Running the Pipeline

1. Dry-run to check configuration:
    ```bash
    snakemake -np
    ```

2. Execute the pipeline:
    ```bash
    snakemake --cores <number_of_threads>
    ```

3. Run on a cluster (SLURM example):
    ```bash
    snakemake --profile slurm
    ```

## Output

All results are stored in the `results/` directory. Key outputs include:
- Quality-controlled reads: `results/trimmomatic/`
- Aligned BAM files: `results/bwa/`
- Final recalibrated BAM files: `results/recal/`
- VCF files from DeepVariant: `results/deepvariant/`
- Annotated variants: `results/annovar/multianno_csv/`
- Processed variants: `data/by_variant.csv`

## Logs

Log files for each step are stored in the `logs/` directory for debugging and reproducibility.

## Contributing

Feel free to contribute by submitting issues or pull requests.

## License

This pipeline is released under the MIT License.
