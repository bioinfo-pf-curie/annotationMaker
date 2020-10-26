# Usage

<!-- TODO - Update with the usage of your pipeline -->

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Main arguments](#main-arguments)
    * [`-profile`](#-profile)
    * [`--genome`](#-genome)
* [Nextflow profiles](#nextflow-profiles)
* [Job resources](#job-resources)
* [Other command line parameters](#other-command-line-parameters)
    * [`--skip*`](#-skip)
    * [`--outDir`](#-outDir)
    * [`-name`](#-name)
    * [`-resume`](#-resume)
    * [`-c`](#-c)
    * [`--maxMemory`](#-maxmemory)
    * [`--maxTime`](#-maxtime)
    * [`--maxCpus`](#-maxcpus)
    * [`--multiqcConfig`](#-multiqcconfig)


## General Nextflow info

Nextflow handles job submissions on SLURM or other environments, and supervises the job execution. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the pipeline is as follows:
```bash
nextflow run main.nf --genome 'hg38' -profile 'conda'
```

This will launch the pipeline with the `conda` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

You can change the output director using the `--outDir/-w` options.

## Main arguments

### `-profile`

Use this option to set the [Nextflow profiles](profiles.md). For example:

```bash
-profile path,cluster
```

### `--genome`

Use this to specify the genome annotation to download and process. For example :

```bash
-genome 'hg38'
```

The list of genomes (and their sources) is available in [`conf/genomes.conf`](../conf/genomes.config)

The pipeline expects a link to a `fasta` file, and an annotation file (`gtf` or `gff`).  
If your genome is not available in the config file, the options `--fasta`, `--gtf` or `--gff` can be used on command line.


## Other arguments

### `--build`

Name of the genome build for output prefix. By default the build is extracted from the name of the fasta file.

### `--indexes`

Specify which indexes you would like to generate (bwa,star,bowtie2,hisat2), as well as `none` or `all` (default).
Several tools can be specify (comma separated).

## Nextflow profiles

Different Nextflow profiles can be used. See [Profiles](profiles.md) for details.

## Job resources

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time (see the [`conf/process.conf`](../conf/process.config) file). 
For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

## Other command line parameters

### `--skip*`

The pipeline is made with a few *skip* options that allow to skip optional steps in the workflow.
The following options can be used:
* `--skipGtfProcessing`
				
### `--outDir`
The output directory where the results will be saved.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML.

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--maxMemory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--maxMemory '8.GB'`

### `--maxTime`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--maxTime '2.h'`

### `--maxCpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--maxCpus 1`


