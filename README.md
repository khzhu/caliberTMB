# Tumor Mutational Burden (TMB) Project

## Overview
Tumor Mutational Burden (TMB) - the total number of mutations (changes) found in the DNA of cancer cells. TMB varies from different sequencing platforms to different type of cancer. Hence, it is important to harmonize variant selection methods and choose reasonable thresholds for the alternative allele count and the read depth. Accurate estimate of the total number of mutations in the targeted regions may help doctors plan the best treatment for each individual patient.

We present a nextflow workflow that calculates tumour mutational burden based on the criteria suggested by the [TMB Harmonization Consortium](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7174078/).

## Contents
This workflow is designed to run genomic analysis on Illumina targeted exome sequencing data, in support of the tumor cancer diagnostic panel for Yale's Molecular Pathology Department.

This pipeline starts from paired-end fastq data (.fastq.gz), and is meant to accompany the output from the [YCGA Illumina demultiplexing pipeline](https://medicine.yale.edu/genetics/research/ycga/faq/). It includes adapter trimming, QC assessment, alignment to the genome, variant calling, and TMB estimation, along with many other steps.

## Containers
The containers directory contains instructions and recipes for building Singularity containers used in the pipeline. Singularity containers are used on the Yale McCleary HPC cluster. The current pipeline configuration for McCleary uses .simg files stored in a shared location on the file system.
## Directory structure
```
TMB-estimation/
├── bin
│   └── calculate_tmb.R
├── modules
│   ├── bcftools
│   │   └── norm
│   │       ├── main.nf
│   │       └── meta.yml
│   ├── bwa
│   │   ├── index
│   │   │   ├── main.nf
│   │   │   └── meta.yml
│   │   └── mem
│   │       ├── main.nf
│   │       └── meta.yml
│   ...
│   ├── vcf2maf
│   │   ├── main.nf
│   │   └── meta.yml
│   └── vep
│       ├── main.nf
│       └── meta.yml
└── subworkflows
    ├── fastq/trim/fastqc
    │   ├── main.nf
    │   └── meta.yml
    ├── index genome
    │   ├── main.nf
    │   ├── meta.yml
    │   └── nextflow.config
    ├── alignment
    │   ├── main.nf
    │   ├── meta.yml
    │   └── nextflow.config
    ├── mutect2/strelka2
    │   ├── main.nf
    │   ├── meta.yml
    │   └── nextflow.config
    ├── tmb calibration
    │   ├── main.nf
    │   ├── meta.yml
    │   └── nextflow.config
├── nextflow.config
├── samples.json
├── README.md
├── main.nf
```
## Set up and run a workflow
1. This repository should first be cloned from GitHub:
```
git clone https://github.com/
```
2. You will need to load JAVA runtime and install a copy of nextflow in your work environment on McCleary prior to running the workflow.
```
module load ANTLR/2.7.7-GCCcore-12.2.0-Java-11
curl -s https://get.nextflow.io | bash
chmod +x nextflow
```
3. The input sample JSON file should have at least five columns: specimen_id, patient_id, tissue (tumor source site), purity (estimated tumor cell percentages in tissue samples), and read1/2 (paired end raw reads).
```
[
    {
        "specimen_id": "CPCT1",
        "patient_id":"CPCT1",
        "tissue": "other",
        "purity": 100,
        "read1":"~/workspace/dsl2/data/CPCT12345678R/CPCT1/CPCT12345678R_AHHKYHDSXX_S13_L001_R1_001.fastq.gz",
        "read2":"~/workspace/dsl2/data/CPCT12345678R/CPCT1/CPCT12345678R_AHHKYHDSXX_S13_L001_R2_001.fastq.gz"
    }
]
```

4. Finally, you can use the following commands to start an interactive session in the McCleary cluster and submit a job. It is recommended that you export the TMPDIR variable to somewhere other than the default directory on McCleary which is /tmp. Alternatively,
you can use the -w option to set a temp directory outside your home directory for the large number of intermediate files produced by each process.
```
ssh user123@mccleary.ycrc.yale.edu
salloc -p ycga -c 1 -t 00-3:00 --mem=8000
```
```
module load ANTLR/2.7.7-GCCcore-12.2.0-Java-11
```
```
nextflow run main.nf -profile hg19 \
--output_dir /SolidTumor/v2_0_validation \
--input_json samples.json \
-w /SolidTumor/v2_0_validation/tmp \
-bg
```
5. For multiple sample processing on HPC, you will need to use slurm job arrays to submit the workflow to avoid hitting [the limit of job submission per hour](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/common-job-failures/#rate-limits)
```
#!/bin/bash
#SBATCH -J CPCT1
#SBATCH --array=1-2%3
#SBATCH --output=CPCT1_%A_%a.out
#SBATCH --error=CPCT1_%A_%a.err
#SBATCH -p ycga
#SBATCH -t 8:0:0
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=32G

module load ANTLR/2.7.7-GCCcore-12.2.0-Java-11
source ~/.bashrc
nextflow run somatic_variant.nf --output_dir /SolidTumor/v2_0_validation --input_json samples.json \
-w /SolidTumor/v2_0_validation/tmp -profile hg19 -resume
```
