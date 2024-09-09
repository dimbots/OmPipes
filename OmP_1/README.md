
# ChIP, ATAC, and RNA-seq Data Processing Pipelines

## Description
This suite of bioinformatics pipelines is designed for ChIP-seq, ATAC-seq, and RNA-seq data analysis. It streamlines quality control, read trimming, alignment, peak calling, visualization, and quantification using state-of-the-art tools such as FastQC, Trimmomatic, HISAT2, MACS2, deepTools, and featureCounts.

## Features

### ChIP-seq and ATAC-seq Pipelines
- **Quality Control**: Uses FastQC for assessing the quality of raw and trimmed reads.
- **Read Trimming**: Employs Trimmomatic for adapter and quality trimming.
- **Alignment**: Utilizes HISAT2, set to avoid spliced alignments, suitable for ChIP-seq and ATAC-seq data.
- **Peak Calling**: Applies MACS2 for detecting enriched regions effectively.
- **Visualization**: Converts BAM files to BigWig format with deepTools for easy visualization of read coverage.
- **Read Counting and Plotting**: Compares read counts between original FASTQ files and filtered BAM files, generating informative plots.

#### ChIP-seq Specific Steps
- **Duplicate Removal**: Removes duplicates using Samtools to avoid overrepresentation of any region.
- **Quality Filtering**: Filters reads based on quality metrics to ensure data integrity.
- **Merging of BAM Files**: Combines multiple BAM files for comprehensive analysis.
- **Downsampling**: Ensures uniform read depth across samples for consistent peak calling.

#### ATAC-seq Specific Steps
- **Mitochondrial Read Filtering**: Excludes mitochondrial DNA to prevent biases in the analysis.
- **Read Shifting**: Adjusts reads to represent transposase cut sites accurately, crucial for ATAC-seq.
- **Downsampling**: Standardizes read depth to minimize batch effects during peak calling.

### RNA-seq Pipeline
- **Quality Control**: Conducts initial and post-trimming quality assessments using FastQC.
- **Read Trimming**: Utilizes Trimmomatic to remove adapters and low-quality bases, enhancing read quality.
- **Alignment**: Aligns reads to the reference genome using HISAT2, followed by conversion to sorted BAM files.
- **Filtering**: Applies Samtools to filter reads by quality, excluding low-quality and unmapped reads.
- **Merging BAM Files**: Merges BAM files for comprehensive analysis across replicates or conditions.
- **Read Counting**: Counts reads per sample and normalizes read depth by downsampling to the minimum count across samples.
- **Indexing**: Indexes the downsampled BAM files to prepare for visualization.
- **Conversion to BigWig**: Generates BigWig files using deepTools' bamCoverage for visualization of read coverage.
- **Gene Quantification**: Employs featureCounts to generate gene-level read count matrices for downstream differential expression analysis.

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
The pipelines generate various output files stored in respective directories:
- **Quality Control Reports**: FastQC reports (`*.html`, `*.zip`) for both raw and trimmed reads.
- **Trimmed Reads**: Adapter-trimmed FASTQ files (`*_trimmed.fastq.gz`).
- **Aligned BAM Files**: HISAT2-aligned reads (`*.bam`).
- **Filtered BAM Files**: Quality-filtered BAM files (`*.filtered.bam`).
- **Merged BAM Files**: Combined BAM files across samples (`*.bam`).
- **Downsampled BAM Files**: BAM files with normalized read depth (`*_downsampled.bam`).
- **Indexed BAM Files**: Indexed BAM files for visualization (`*.bam`, `*.bam.bai`).
- **BigWig Files**: BigWig files for read coverage visualization (`*.bw`).
- **Read Counts**: CSV files and bar plots comparing read counts between raw and filtered files.
- **Gene Quantification**: Gene-level count matrices (`*.counts.txt`).

The output files are organized within the `results` directory, with subdirectories corresponding to each analysis step.
