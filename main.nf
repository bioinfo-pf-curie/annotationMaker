#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2020
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/


/*
========================================================================================
                         Annotation Maker
========================================================================================
 Annotation Maker Pipeline.
 #### Homepage / Documentation
 https://gitlab.curie.fr/data-analysis/annotationMaker
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    if ("${workflow.manifest.version}" =~ /dev/ ){
       dev_mess = file("$baseDir/assets/dev_message.txt")
       log.info dev_mess.text
    }

    log.info """
    as-mapping v${workflow.manifest.version}
    ======================================================================

    Usage:
    nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome mm9 -profile conda
    nextflow run main.nf --samplePlan sample_plan --genome mm9 -profile conda


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

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}


def defineAligners() {
  return [
        'bwa',
        'star',
        'bowtie2',
        'hisat2']
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

alignersList = defineAligners()
aligners = params.indexes ? params.indexes == 'all' ? alignersList : params.indexes.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(aligners, alignersList)) exit 1, 'Unknown Aligner(s), see --help for more information'



/*
 * CHANNELS
 */

// Reference Genome
Channel.fromPath("${params.fasta}")
       .ifEmpty { exit 1, "Reference Genome not found: ${params.fasta}" }
       .into { chFasta; chFastaBwa; chFastaStar; chFastaBowtie2; chFastaHisat2 }

// GTF
if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtfHisat2Splicesites; gtfHisat2Index; chGtf }
}

if ( params.build ){
  build = params.build
}else{
  build = params.fasta.baseName - ~/(.fa)?(.fasta)?/
}


// Header log info
if ("${workflow.manifest.version}" =~ /dev/ ){
   dev_mess = file("$baseDir/assets/dev_message.txt")
   log.info dev_mess.text
}

log.info """=======================================================

asmapping : Annotation Maker workflow v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Command Line'] = workflow.commandLine
summary['Fasta']          = params.fasta
summary['Build']          = build

summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Container Engine'] = workflow.containerEngine
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Config Profile'] = workflow.profile

if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 * PREPROCESSING - Prepare Genome
 */

process indexFasta {
  label 'process_small'
  publishDir "${params.outdir}/genome", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  input:
  file(fasta) from chFasta

  output:
  set file(fasta), file("*.fai") into chFastaDict, chFastaSize 

  script:
  """
  samtools faidx ${fasta}
  """
}

process makeDict {
  label 'process_small'

  publishDir "${params.outdir}/genome", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  input:
  set file(fasta), file(faidx) from chFastaDict

  output:
  file(".dict") into chDict

  script:
  pfix = fasta.toString() - /(.fa)?(.fasta)?/
  """
  picard CreateSequenceDictionary REFERENCE=${fasta} OUTPUT=${pfix}.dict
  """
}

process makeChromSizes {
  label 'process_small'

  publishDir "${params.outdir}/genome", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  input:
  set file(fasta), file(faidx) from chFastaSize

  output:
  file("*sizes") into chChromSize

  script:
  pfix = fasta.toString() - /(.fa)?(.fasta)?/
  """
  cut -f1,2 ${faidx} > chrom_${pfix}.sizes
  """
}

/*
 * PREPROCESSING - Build Index
 */ 

process makeBwaIndex {
  label 'process_medium'
  publishDir "${params.outdir}/indexes/bwa", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  when: 
  !('bwa' in aligners)
   
  input:
  file(fasta) from chFastaBwa

  output:
  file "bwa" into chBwaIdx

  script:
  """
  mkdir bwa
  bwa index -p bwa/${build} ${fasta} > bwa.log 2>&1
  """
}

process makeStarIndex {
  label 'process_high'
  publishDir "${params.outdir}/indexes/star", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  when:
  !('star' in aligners)

  input:
  file(fasta) from chFastaStar

  output:
  file "STAR" into chStarIdx

  script:
  """
  STAR --runMode genomeGenerate --limitGenomeGenerateRAM 33524399488 --runThreadN ${task.cpus} --genomeDir STAR --genomeFastaFiles $fasta
  """
}

process makeBowtie2Index {
  label 'process_low'
  publishDir "${params.outdir}/indexes/bowtie2", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  when:
  !('bowtie2' in aligners)

  input:
  file(fasta) from chFastaBowtie2

  output:
  file("bowtie2") into chBowtie2Idx

  script:
  """
  mkdir -p bowtie2
  bowtie2-build ${fasta} bowtie2/${build}
  """
}

process makeHisat2Splicesites {
  label 'process_low'
  publishDir "${params.outdir}/reference_genome/indexes", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  when:
  !('hisat2' in aligners)

  input:
  file gtf from gtfHisat2Splicesites.collect()

  output:
  file "${gtf.baseName}.hisat2_splice_sites.txt" into indexingSplicesites

  script:
  """
  hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
  """
}

process makeHisat2Index {
  label 'process_extra'
  publishDir "${params.outdir}/reference_genome/indexes", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  when:
  !('hisat2' in aligners)

  input:
  file fasta from chFastaHisat2
  file indexing_splicesites from indexingSplicesites.collect()
  file gtf from gtfHisat2Index.collect()

  output:
  file("hisat2") into chHisat2Idx

  script:
  """
  mkdir -p hisat2
  hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt
  hisat2-build -p ${task.cpus} --ss $indexing_splicesites --exon ${gtf.baseName}.hisat2_exons.txt $fasta hisat2/${build}
  """
}


/*
 * Sub-routine
 */

workflow.onComplete {

    /*pipeline_report.html*/

    def report_fields = [:]
    report_fields['version'] = workflow.manifest.version
    report_fields['runName'] = custom_runName ?: workflow.runName
    report_fields['success'] = workflow.success
    report_fields['dateComplete'] = workflow.complete
    report_fields['duration'] = workflow.duration
    report_fields['exitStatus'] = workflow.exitStatus
    report_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    report_fields['errorReport'] = (workflow.errorReport ?: 'None')
    report_fields['commandLine'] = workflow.commandLine
    report_fields['projectDir'] = workflow.projectDir
    report_fields['summary'] = summary
    report_fields['summary']['Date Started'] = workflow.start
    report_fields['summary']['Date Completed'] = workflow.complete
    report_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    report_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) report_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) report_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) report_fields['summary']['Pipeline Git branch/tag'] = workflow.revision

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/oncomplete_template.txt")
    def txt_template = engine.createTemplate(tf).make(report_fields)
    def report_txt = txt_template.toString()
    
    // Render the HTML template
    def hf = new File("$baseDir/assets/oncomplete_template.html")
    def html_template = engine.createTemplate(hf).make(report_fields)
    def report_html = html_template.toString()

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << report_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << report_txt }

    /*oncomplete file*/

    File woc = new File("${params.outdir}/workflow.oncomplete.txt")
    Map endSummary = [:]
    endSummary['Completed on'] = workflow.complete
    endSummary['Duration']     = workflow.duration
    endSummary['Success']      = workflow.success
    endSummary['exit status']  = workflow.exitStatus
    endSummary['Error report'] = workflow.errorReport ?: '-'
    String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
    println endWfSummary
    String execInfo = "${fullSum}\nExecution summary\n${logSep}\n${endWfSummary}\n${logSep}\n"
    woc.write(execInfo)

    /*final logs*/
    if(workflow.success){
        log.info "[rnaseq] Pipeline Complete"
    }else{
        log.info "[rnaseq] FAILED: $workflow.runName"
        if( workflow.profile == 'test'){
            log.error "====================================================\n" +
                    "  WARNING! You are running with the profile 'test' only\n" +
                    "  pipeline config profile, which runs on the head node\n" +
                    "  and assumes all software is on the PATH.\n" +
                    "  This is probably why everything broke.\n" +
                    "  Please use `-profile test,conda` or `-profile test,singularity` to run on local.\n" +
                    "  Please use `-profile test,conda,cluster` or `-profile test,singularity,cluster` to run on your cluster.\n" +
                    "============================================================"
        }
    }
 
}
