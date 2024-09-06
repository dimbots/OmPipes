#!/usr/bin/env Rscript

# List of required packages
required_packages <- c("ggplot2", "dplyr", "tidyr")

# Function to check and install packages
install_if_missing <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, repos = "http://cran.us.r-project.org")
      library(package, character.only = TRUE)
    }
  }
}

# Install necessary packages
install_if_missing(required_packages)

# Rest of your script
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to count reads in BAM files
count_bam_reads <- function(bam_file) {
  cmd <- paste("samtools view", bam_file, "| wc -l")
  as.numeric(system(cmd, intern = TRUE))
}

# Function to count reads in FASTQ files
count_fastq_reads <- function(fastq_file) {
  cmd <- paste("zcat", fastq_file, "| wc -l")
  as.numeric(system(cmd, intern = TRUE)) / 4  # Divide by 4 because each read in FASTQ is 4 lines
}

# Continue with the rest of your script...
