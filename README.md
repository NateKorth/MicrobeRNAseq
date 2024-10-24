# Maize Microbe RNAseq Analysis Pipeline
This pipeline is a series of unix commands built to be applied in a IBM bsub submission system as employed by the North Carolina State University High Performance Computing center. 
The Conda environment(s) can be found in the environments folder and can be built independently with the most current versions of the software by installing the Required Packages to your own Conda environment. Some scripts for visualization/stats are written in R and can be easily applied 

## Repo information
This pipeline is designed to analyze RNAseq data. It seperates reads from plants, bacteria, and fungi and attempts to assign functionality / pathway information to each read
To use it, first clone the directory by entering the following command in your terminal: 
```
git clone git@github.com:NateKorth/MicrobeRNAseq.git
```

## Required Packages
The following can be installed to a conda environment with the command conda install
* ncbi-genome-download
* hisat2
* ribodetector
* DEseq2
* Kraken2
* eggNog-mapper

## Step0 Download Required Genomes
The human genome can be found here: https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/

The Maize reference genome can be found here (Get both the genome and gff annotation file): https://www.maizegdb.org/assembly

## Step 1 Trim/QC with fastqc and trimmomatic

## Step 2 Align reads to human genome and remove from downstream analysis

## Step 3 Remove rRNA reads

## Step 4 Align reads to maize genome

## Step 4a Process maize reads
Calculate tpm

Conduct DEseq

## Step 5 Assign Bacterial / Fungal taxonomy to remaining reads

## Contact
For clarification on code missing annotation contact:
* Nate Korth: njkorth@ncsu.edu / nate.korth@gmail.com
* Joe Gage: jlgage@ncsu.edu
