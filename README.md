# NextFlow_RD_Genomic

## Description

A simple base Rare disease and germline genomics pipeline to test the effects of down-sampling on variant calling

## Basic Overview
Using the NextFlow workflow software to run the following pipeline

### Pipeline
Index genome > Fastqc analysis > Align reads > Downsample bam files > Sort bam > Mark duplicates > Index bam > 
Call variants > Hard filter

## Setup
To run the pipeline, we need to obtain 

- A genome build (GRCh38) - provided by the Broad institute
```bash
$ cd data/genome
$ wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
```
- FastQ sample (for workflow development)
```bash
$ cd ../samples
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518253/SRR1518253_1.fastq.gz
$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518253/SRR1518253_2.fastq.gz
```
- When scaling up, FastQ samples
```bash
$ wget https://genomics.viapath.co.uk/benchmark/files/FASTQ/NA12878_WES.zip
```

## Running the pipeline
See docker.md

## Validating the pipeline
See [https://genomics.viapath.co.uk/benchmark](https://genomics.viapath.co.uk/benchmark)