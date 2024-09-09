
# ChIP, ATAC, and RNA-seq Data Processing Pipelines

## Description
This suite of bioinformatics pipelines is designed for ChIP-seq, ATAC-seq, and RNA-seq data analysis. It streamlines quality control, read trimming, alignment, peak calling, visualization, and quantification using state-of-the-art tools such as FastQC, Trimmomatic, HISAT2, MACS2, deepTools, and featureCounts.

## Prerequisites
- [Nextflow](https://www.nextflow.io/): A workflow tool to run tasks across multiple compute infrastructures in a very portable manner.
- [Conda](https://docs.conda.io/en/latest/): Package, dependency, and environment management for any language—Python, R, Ruby, Lua, Scala, Java, JavaScript, C/ C++, FORTRAN.

## Installation
1. **Install Nextflow**:
   ```bash
   curl -s https://get.nextflow.io | bash
   ```
2. **Install Dependencies with Conda**:
   ```bash
   bash <conda-installer-name>-latest-Linux-x86_64.sh
   ```

## Usage

- **For ChIP-seq**:
  ```bash
  nextflow run chip-seq-pipeline-v1.0.0.SE.nf --profile conda -c merge.config
  ```
- **For ATAC-seq**:
  ```bash
  nextflow run atac-seq-pipeline-v1.0.0.SE.nf --profile conda -c merge.config
  ```
- **For RNA-seq**:
  ```bash
  nextflow run rna-seq-pipeline-v1.2.0.SE.nf --profile conda
  ```

Replace `.SE.nf` with `.PE.nf` for paired-end data processing.

## Output
The output files are organized within the `results` directory, with subdirectories corresponding to each analysis step.
