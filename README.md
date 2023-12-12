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
$ wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518253/SRR1518253_1.fastq.gz
$ wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518253/SRR1518253_2.fastq.gz
$ gunzip *.gz
```
- When scaling up, FastQ samples
```bash
$ wget https://genomics.viapath.co.uk/benchmark/files/FASTQ/NA12878_WES.zip
```

## Running the pipeline
```bash
$ nextflow run main.nf
```

## Validating the pipeline
See [https://genomics.viapath.co.uk/benchmark](https://genomics.viapath.co.uk/benchmark)

## DNANexus applet setup (A local applet for basic testing)
- DNANexus Python Bindings [Documentation](https://github.com/dnanexus/dx-toolkit) 
- [Install the app](https://documentation.dnanexus.com/downloads) 
```bash
pip install -r requirements.txt
```
- Routine maintenance
Periodically update dxpy
```bash
$ pip install --upgrade dxpy
```

### DNANexus Tutorial
- [Overview videos](https://documentation.dnanexus.com/getting-started)
- [Developer tutorial](https://documentation.dnanexus.com/getting-started/developer-quickstart)
```bash
$ dx select <your-project-name>
$ dx build --nextflow
```