
# ChIP and ATAC Data Processing Pipelines

## Description
This suite of bioinformatics pipelines is designed for ChIP-seq and ATAC-seq data analysis. It streamlines quality control, read trimming, alignment, peak calling, and visualization, utilizing tools like FastQC, Trimmomatic, HISAT2, MACS2, and deepTools.

## Features
- **Quality Control**: Uses FastQC for assessing the quality of raw and trimmed reads.
- **Read Trimming**: Employs Trimmomatic for adapter and quality trimming.
- **Alignment**: Utilizes HISAT2, set to avoid spliced alignments, suitable for ChIP-seq and ATAC-seq data.
- **Peak Calling**: Applies MACS2 for detecting enriched regions effectively.
- **Visualization**: Converts BAM files to BigWig format with deepTools for easy visualization of read coverage.

## ChIP-seq Specific Steps
- **Duplicate Removal**: Removes duplicates using Samtools to avoid overrepresentation of any region.
- **Quality Filtering**: Filters reads based on quality metrics to ensure data integrity.
- **Merging of BAM Files**: Combines multiple BAM files for comprehensive analysis.
- **Downsampling**: Ensures uniform read depth across samples for consistent peak calling.

## ATAC-seq Specific Steps
- **Mitochondrial Read Filtering**: Excludes mitochondrial DNA to prevent biases in the analysis.
- **Read Shifting**: Adjusts reads to represent transposase cut sites accurately, crucial for ATAC-seq.
- **Downsampling**: Standardizes read depth to minimize batch effects during peak calling.

## Prerequisites
- [Nextflow](https://www.nextflow.io/): A workflow tool to run tasks across multiple compute infrastructures in a very portable manner.
- [Conda](https://docs.conda.io/en/latest/): Package, dependency and environment management for any language—Python, R, Ruby, Lua, Scala, Java, JavaScript, C/ C++, FORTRAN.

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
- Replace `.SE.nf` with `.PE.nf` for paired-end data processing.

