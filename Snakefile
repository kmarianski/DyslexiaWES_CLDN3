
# Full Snakefile: WES Pipeline
# Description: Comprehensive, modular pipeline for Whole-Exome Sequencing (WES) analysis.

import os
from pathlib import Path

# Configuration
configfile: "config.yaml"

# Read samples from the raw directory
RAW_DIR = Path("raw")
samples = [sample.split("_")[0] for sample in os.listdir(RAW_DIR)]
SAMPLES = list(set(samples))

# Rule all specifies the final targets of the pipeline
rule all:
    input:
        expand("results/trimmomatic/{sample}_{lane}P.fq.gz", sample=samples, lane=[1, 2]),
        expand("results/bwa/{sample}.bam", sample=samples),
        expand("results/merged_bams/{s}.bam", s=SAMPLES),
        expand("results/recal/{s}.bam.bai", s=SAMPLES),
        expand("results/deepvariant/{s}.vcf.gz", s=SAMPLES),
        expand("results/annovar/multianno_csv/{s}.hg38_multianno.csv", s=SAMPLES)

# Directory creation rule
rule setup_directories:
    output:
        touch("logs/setup_directories.done")
    run:
        for subdir in config["directories"]:
            Path("results/" + subdir).mkdir(parents=True, exist_ok=True)
        for subdir in config["log_dirs"]:
            Path("logs/" + subdir).mkdir(parents=True, exist_ok=True)
        shell("touch {output}")

# Trimmomatic
rule trimmomatic:
    input:
        r1="raw/{sample}_1.fq.gz",
        r2="raw/{sample}_2.fq.gz"
    output:
        r1_paired="results/trimmomatic/{sample}_1P.fq.gz",
        r2_paired="results/trimmomatic/{sample}_2P.fq.gz",
        r1_unpaired="results/trimmomatic/{sample}_1U.fq.gz",
        r2_unpaired="results/trimmomatic/{sample}_2U.fq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    threads:
        config["threads"]
    params:
        adapter=config["adapter"]
    conda:
        "envs/biotools.yml"
    shell:
        "trimmomatic PE -threads {threads} {input.r1} {input.r2} "
        "{output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} "
        "ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 "
        "&> {log}"

# BWA Mapping
rule bwa_map:
    input:
        r1=rules.trimmomatic.output.r1_paired,
        r2=rules.trimmomatic.output.r2_paired
    output:
        bam="results/bwa/{sample}.bam"
    params:
        index=config["reference"]
    threads:
        config["threads"]
    log:
        "logs/bwa/{sample}.log"
    conda:
        "envs/biotools.yml"
    shell:
        "bwa mem -t {threads} {params.index} {input.r1} {input.r2} | "
        "samtools sort -@{threads} -o {output.bam}"

# Samtools Merge
rule samtools_merge:
    input:
        bam_files=expand("results/bwa/{sample}.bam", sample=samples)
    output:
        "results/merged_bams/{s}.bam"
    conda:
        "envs/biotools.yml"
    log:
        "logs/merged_bams/{s}.log"
    shell:
        "samtools merge {output} {input}"

# Replace Read Groups
rule replace_rg:
    input:
        "results/merged_bams/{s}.bam"
    output:
        "results/replace_rg/{s}.bam"
    params:
        extra="--RGLB lib1 --RGPL illumina --RGPU {s} --RGSM {s}"
    log:
        "logs/picard/replace_rg/{s}.log"
    wrapper:
        "v1.14.1/bio/picard/addorreplacereadgroups"

# Mark Duplicates
rule markduplicates_bam:
    input:
        bams="results/replace_rg/{s}.bam"
    output:
        bam="results/dedup/{s}.bam",
        metrics="results/dedup/{s}.metrics.txt"
    log:
        "logs/picard/dedup/{s}.log"
    params:
        extra="--REMOVE_DUPLICATES true"
    wrapper:
        "v3.3.3/bio/picard/markduplicates"

# Base Recalibrator
rule base_recalibrator:
    input:
        bam="results/dedup/{s}.bam",
        ref=config["reference"],
        known_sites=config["known_sites"]
    output:
        recal_table="results/recal/{s}.table"
    log:
        "logs/picard/recal_table/{s}.log"
    conda:
        "envs/biotools.yml"
    shell:
        "gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known_sites} "
        "-O {output.recal_table}"

# Apply BQSR
rule apply_bqsr:
    input:
        bam="results/dedup/{s}.bam",
        recal_table="results/recal/{s}.table"
    output:
        bam="results/recal/{s}.bam"
    log:
        "logs/picard/recal_apply/{s}.log"
    conda:
        "envs/biotools.yml"
    shell:
        "gatk ApplyBQSR -I {input.bam} --bqsr-recal-file {input.recal_table} -O {output.bam}"

# Index BAM
rule samtools_index:
    input:
        bam="results/recal/{s}.bam"
    output:
        bai="results/recal/{s}.bam.bai"
    conda:
        "envs/biotools.yml"
    log:
        "logs/samtools/index/{s}.log"
    shell:
        "samtools index {input.bam}"

# DeepVariant
rule deepvariant:
    input:
        bam="results/recal/{s}.bam",
        ref=config["reference"]
    output:
        vcf="results/deepvariant/{s}.vcf.gz"
    params:
        model="WES"
    threads:
        64
    log:
        "logs/deepvariant/{s}.log"
    singularity:
        "singularity/deepvariant_latest.sif"
    shell:
        "/opt/deepvariant/bin/run_deepvariant --model_type {params.model} "
        "--ref {input.ref} --reads {input.bam} --output_vcf {output.vcf} --num_shards {threads}"

# Annovar Recode
rule annovar_recode:
    input:
        "results/deepvariant/{s}.vcf.gz"
    output:
        recode="results/annovar/{s}.Info.recode.vcf"
    log:
        "logs/annovar/{s}_recode.log"
    conda:
        "envs/biotools.yml"
    shell:
        "vcftools --gzvcf {input} --recode --out {output.recode}"

# Convert to Annovar
rule convert2annovar:
    input:
        recode="results/annovar/{s}.Info.recode.vcf"
    output:
        avinput="results/annovar/{s}.avinput"
    shell:
        "convert2annovar.pl -format vcf4 {input.recode} > {output.avinput}"

# Annotate Variation
rule annotate_variation:
    input:
        avinput="results/annovar/{s}.avinput"
    output:
        csv="results/annovar/multianno_csv/{s}.hg38_multianno.csv"
    params:
        annovar_db=config["annovar_db"]
    threads:
        config["threads"]
    shell:
        "table_annovar.pl {input.avinput} {params.annovar_db} -buildver hg38 "
        "-out {output.csv} -remove -protocol refGene,dbnsfp31a_interpro,gnomad30_genome "
        "-operation g,f,f -nastring . -csvout"
