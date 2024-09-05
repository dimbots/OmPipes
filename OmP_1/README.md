
# ChIPipe - A ChIP-seq Data Processing Pipeline

## Description
ChIPipe is a bioinformatics pipeline designed for the analysis of ChIP-seq data. It automates the quality control, read trimming, alignment, and peak calling processes, leveraging powerful tools such as FastQC, Trimmomatic, HISAT2, and MACS2.

## Features
- Quality control with FastQC.
- Read trimming using Trimmomatic.
- Flexible alignment options with HISAT2.
- Peak calling with MACS2.

## Prerequisites
- [Nextflow](https://www.nextflow.io/) installed.
- [Conda](https://docs.conda.io/en/latest/) for managing dependencies.
- Reference genome files (FASTA format and indexed).

## Installation
Install Nextflow:
```bash
curl -s https://get.nextflow.io | bash
```

Install Conda dependencies:
```bash
bash <conda-installer-name>-latest-Linux-x86_64.sh
```

## Usage
Run the pipeline with:
```bash
nextflow run ChIPipe.1.0.9.SE.nf -profile conda -c merge.config
```

## Parameters
- `--merge.config`: Path to the merge.config file, used to specify how to merge replicates.

