# Annotation Maker

**Institut Curie - annotationMaker**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![MultiQC](https://img.shields.io/badge/MultiQC-1.8-blue.svg)](https://multiqc.info/)
[![Install with](https://anaconda.org/anaconda/conda-build/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
<!--[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)-->

## !! UNDER CONSTRUCTION !! 

### Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It comes with conda / singularity containers making installation easier and results highly reproducible.

The goal of this pipeline is to generate annotation and indexes files in a standardized way for production analysis pipelines.

### Pipline summary

1. Process `.fasta` file and generate `.dict` and `.fai` files
2. Generate chromosome size file
3. Generate indexes for :
- [`BWA`](http://bio-bwa.sourceforge.net/)
- [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 
- [`STAR`](https://github.com/alexdobin/STAR)
- [`hisat2`](http://ccb.jhu.edu/software/hisat2/index.shtml)

### Quick help

```bash
N E X T F L O W  ~  version 20.01.0
Launching `main.nf` [furious_volta] - revision: ca21e76c70
======================================================================
annotationMaker v1.0.1dev
======================================================================

Usage:
nextflow run main.nf --fasta '*.fasta' -profile conda

Mandatory arguments:
  --fasta [file]                Path to input data (must be surrounded with quotes)
  --gtf [file]                  Path to GTF file with gene annotation
  -profile [str]                Configuration profile to use. test / conda / toolsPath / singularity / cluster (see below)
	  
Optional arguments:
  --build                       Build name to use for genome index
  --indexes                     List of indexes to build. Available: all,bwa,star,bowtie2,hisat2. Default: all
		  
Other options:
  --outdir [file]               The output directory where the results will be saved
  -w/--work-dir [file]          The temporary directory where intermediate data will be saved
  --email [str]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
  -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
				  
			  
=======================================================
Available Profiles
				  
  -profile test                Set up the test dataset
  -profile conda               Build a new conda environment before running the pipeline
  -profile toolsPath           Use the paths defined in configuration for each tool
  -profile singularity         Use the Singularity images for each process
  -profile cluster             Run the workflow on the cluster, instead of locally							
```

### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follow

#### Run the pipeline on a test dataset
See the conf/test.conf to set your test dataset.

```
nextflow run main.nf -profile test,conda

```

#### Run the pipeline on a computational cluster

```
echo "nextflow run main.nf --fasta '*.fasta' --outdir MY_OUTPUT_DIR -profile singularity,cluster" | qsub -N makeAnnot"

```

### Defining the '-profile'

By default (whithout any profile), Nextflow will excute the pipeline locally, expecting that all tools are available from your `PATH` variable.

In addition, we set up a few profiles that should allow you i/ to use containers instead of local installation, ii/ to run the pipeline on a cluster instead of on a local architecture.
The description of each profile is available on the help message (see above).

Here are a few examples of how to set the profile option.

```
## Run the pipeline locally, using the paths defined in the configuration for each tool (see conf.tool-path.config)
-profile toolsPath

## Run the pipeline on the cluster, using the Singularity containers
-profile cluster,singularity

## Run the pipeline on the cluster, building a new conda environment
-profile cluster,conda

```

### Sample Plan

A sample plan is a csv file (comma separated) that list all samples with their biological IDs.


Sample ID | Sample Name | Path R1 .fastq file | [Path R2 .fastq file]

### Full Documentation

1. [Installation](docs/installation.md)
2. [Reference genomes](docs/reference_genomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

#### Credits

This pipeline has been written by the bioinformatics platform of the Institut Curie (P. La Rosa, N. Servant)

#### Contacts

For any question, bug or suggestion, please use the issues system or contact the bioinformatics core facility.

