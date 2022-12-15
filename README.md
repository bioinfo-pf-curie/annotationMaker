# Annotation Maker

**Institut Curie - annotationMaker**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![Install with](https://anaconda.org/anaconda/conda-build/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)

### Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It comes with conda / singularity containers making installation easier and results highly reproducible.

The goal of this pipeline is to generate annotation and indexes files in a standardized way for production analysis pipelines.

### Pipline summary

1. Processes `.fasta` file and generates `.dict` and `.fai` files
2. Generates chromosome size file
3. Calculates the effective genome size as the sum of non 'N' base on the genome.
4. Generates indexes for :
- [`BWA`](http://bio-bwa.sourceforge.net/)
- [`BWA-mem2`](https://github.com/bwa-mem2/bwa-mem2)
- [`DRAGMAP`](https://github.com/Illumina/DRAGMAP)
- [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 
- [`STAR`](https://github.com/alexdobin/STAR)
- [`hisat2`](http://ccb.jhu.edu/software/hisat2/index.shtml)
- [`cellRanger`](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
- [`Kallisto`](https://pachterlab.github.io/kallisto/about)
- [`Salmon`](https://combine-lab.github.io/salmon/)
5. Processes GTF/GFF annotation file for downstream analysis tools

### Quick help

```bash
N E X T F L O W  ~  version 21.10.6
Launching `main.nf` [gloomy_khorana] - revision: 20613d4cd4
------------------------------------------------------------------------
                           _        _   _             __  __       _  
                          | |      | | (_)           |  \/  |     | |
    __ _ _ __  _ __   ___ | |_ __ _| |_ _  ___  _ __ | .  . | __ _| | _____ _ __ 
   / _` | '_ \| '_ \ / _ \| __/ _` | __| |/ _ \| '_ \| |\/| |/ _` | |/ / _ \ '__|
  | (_| | | | | | | | (_) | || (_| | |_| | (_) | | | | |  | | (_| |   <  __/ |   
   \__,_|_| |_|_| |_|\___/ \__\__,_|\__|_|\___/|_| |_\_|  |_/\__,_|_|\_\___|_| 

  v2.0.0
------------------------------------------------------------------------
Usage:
The typical command for running the pipeline is as follows:

nextflow run main.nf --profile STRING --genome STRING -profile PROFILES

MANDATORY ARGUMENTS:

--genome  STRING   Name of the reference genome.
--profile STRING [conda, cluster, docker, multiconda, conda, path, multipath, singularity]  Configuration profile to use. Can use multiple (comma separated).

REFERENCES:
--fasta                PATH   Path to genome fasta file
--genomeAnnotationPath PATH   Path to genome annotations folder
--gff                  PATH   Path to GFF annotation file
--gtf                  PATH   Path to GTF annotation file

INDEXES:

--build   STRING   Name of the genome build
--indexes STRING [all, bwa, bwamem2, dragmap, star, bowtie2, hisat2, cellranger, kallisto, salmon, none]  Genome indexes to generate

OTHER OPTIONS:

--cellRangerPath    PATH   CellRanger path
--outDir            PATH   The output directory where the results will be saved
--skipGtfProcessing        Skip GTF processing steps
--starVersion       STRING [2.7.6a, 2.7.8a]   Version of the STAR aligned to use
			  
=======================================================
Available Profiles
  -profile test                    Run the test dataset
  -profile conda                   Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
  -profile multiconda              Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
  -profile path                    Use the installation path defined for all tools. Use `--globalPath` to define the insallation path
  -profile multipath               Use the installation paths defined for each tool. Use `--globalPath` to define the insallation path
  -profile docker                  Use the Docker images for each process
  -profile singularity             Use the Singularity images for each process. Use `--singularityPath` to define the insallation path
  -profile cluster                 Run the workflow on the cluster, instead of locally
```

### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follow

#### Run the pipeline on a test dataset
See the conf/test.conf to set your test dataset.

```
nextflow run main.nf -profile test,multiconda
```

#### Generate all annotations for a given species from an existing conda cache

```
echo "nextflow run main.nf --genome 'hg38' \
      -profile multiconda,cluster --condaCacheDir MY_CONDA_CACHE \
      --outDir MY_OUTPUT_DIR -w MY_OUTPUT_DIR/work" | qsub -N makeAnnot"
```

#### Generate a few indexes from already downloaded genome fasta file

```
echo "nextflow run main.nf --genome 'mm39' \
      --starVersion 2.7.8a --indexes star,kallisto,salmon \
      --fasta GENOME_FASTA --skipGtfProcessing \
      -profile multiconda,cluster --condaCacheDir MY_CONDA_CACHE \
      --outDir MY_OUTPUT_DIR -w MY_OUTPUT_DIR/work -resume" | qsub -N hg19
```

### Defining the '-profile'

By default (whithout any profile), Nextflow will excute the pipeline locally, expecting that all tools are available from your `PATH` variable.

In addition, we set up a few profiles that should allow you i/ to use containers instead of local installation, ii/ to run the pipeline on a cluster instead of on a local architecture.
The description of each profile is available on the help message (see above).

Here are a few examples of how to set the profile option.

```
## Run the pipeline locally, using the paths defined in the configuration for each tool (see conf.tool-path.config)
-profile path --globalPath 'PATH_TO_BINARY'

## Run the pipeline on the cluster, using the Singularity containers
-profile cluster,singularity

## Run the pipeline on the cluster, building a new conda environment
-profile cluster,conda

```

Sample ID | Sample Name | Path R1 .fastq file | [Path R2 .fastq file]

### Full Documentation

1. [Installation](docs/installation.md)
2. [Reference genomes](docs/referenceGenomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

#### Credits

This pipeline has been written by the bioinformatics platform of the Institut Curie

#### Contacts

For any question, bug or suggestion, please use the issues system or contact the bioinformatics core facility.

