# Outputs

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes the data using the steps presented in the main README file.  
Briefly, its goal is to generate ready-to-use annotations for downstream analysis pipelines.

The directories listed below will be created in the output directory after the pipeline has finished. 

## `genome`

The `genome` folder contains everything about the assembly itself, with :
- The genome in `fasta` format (.fa, .fas, .fasta)
- The genome `index` file (.fa.fai)
- The genome `dict` file (.dict) mainly for GATK usage
- The chromosome size file (.sizes)
- The effective genome size defined by the sum of non 'N' base on the genome (.effgsize)

## `indexes`

The current version of the pipeline allows to automatically generates indexes for the following mappers ;
- STAR
- bowtie2
- bwa
- hisat2

## `gtf`

The gtf folder contains the annotations files generated from `gtf` or `gff` files, with :
- The complete annotation in `gtf/gff` format (.gtf)
- The annotation in `bed12` format, usefull for some tools like `rseqc`
- The gene-based annotation in `bed` format
- If available, all these files are also created for `protein coding` annotation only, using the keyword `gene_type` or `gene_biotype`
