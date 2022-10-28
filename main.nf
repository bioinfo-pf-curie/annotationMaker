#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019-2022
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/

/*
========================================================================================
                         ANNOTATION MAKER
========================================================================================
Analysis Pipeline DSL2 template.
https://patorjk.com/software/taag/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Initialize lintedParams and paramsWithUsage
NFTools.welcome(workflow, params)

// Use lintedParams as default params object
paramsWithUsage = NFTools.readParamsFromJsonSettings("${projectDir}/parameters.settings.json")
params.putAll(NFTools.lint(params, paramsWithUsage))

// Run name
customRunName = NFTools.checkRunName(workflow.runName, params.name)

// Custom functions/variables
mqcReport = []
include {checkAlignmentPercent} from './lib/functions'

/*
===================================
  SET UP CONFIGURATION VARIABLES
===================================
*/

// Genome-based variables
if (!params.genome){
  exit 1, "No genome provided. The --genome option is mandatory"
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Initialize variable from the genome.conf file
//params.bowtie2Index = NFTools.getGenomeAttribute(params, 'bowtie2')

// Stage config files
outputDocsCh = Channel.fromPath("$projectDir/docs/output.md")
//outputDocsImagesCh = file("$projectDir/docs/images/", checkIfExists: true)

/*
==========================
 VALIDATE INPUTS
==========================
*/

if (params.starVersion != "2.7.8a" && params.starVersion != "2.7.6a" && params.starVersion != "2.6.1b"){
  exit 1, "The provided STAR version is not available. Use either 2.6.1b, 2.7.6a or 2.7.8a"
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
   exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

def defineAligners() {
  return [
        'bwa',
        'star',
        'bowtie2',
        'hisat2',
        'cellranger',
        'kallisto',
	'salmon',
        'none']
}

// Compare each parameter with a list of parameters
static def checkParameterExistence(it, list) {
  if (!list.contains(it)) {
    println("Unknown parameter: ${it}")
    return false
  }
  return true
}

def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

alignersList = defineAligners()
def aligners = params.indexes || params.indexes == 'none' ? params.indexes == 'all' ? alignersList : params.indexes.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(aligners, alignersList)) exit 1, 'Unknown Aligner(s), see --help for more information'



/*
==========================
 BUILD CHANNELS
==========================
*/

// Reference Genome
chFasta = Channel.empty()
chFastaURL = Channel.empty()
wasFastaURL = false
if (params.fasta){
  chFasta = Channel.fromPath("${params.fasta}").ifEmpty { exit 1, "Reference Genome not found: ${params.fasta}" }
}else if (params.genome){
  fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
  if (fasta){
    if (fasta.startsWith("http") || fasta.startsWith("ftp")){
      chFastaURL = Channel.from(fasta)
      wasFastaURL = true
    }else{
      chFasta = Channel.fromPath(fasta).ifEmpty { exit 1, "Reference Genome not found: ${fasta}" }
    }
  }else{
    exit 1, "Fasta file not found for ${params.genome} : ${fasta}"
  }
}

// Transcriptome
chTrsFasta = Channel.empty()
chTrsURL = Channel.empty()
wasTrsURL = false
if (params.transcripts){
  chTrsFasta = Channel.fromPath("${params.transcripts}").ifEmpty { exit 1, "Reference Genome not found: ${params.transcripts}" }
}else if (params.genome){
  transcriptome = params.genome ? params.genomes[ params.genome ].transcripts ?: false : false
  if (transcriptome){
    if (transcriptome.startsWith("http") || transcriptome.startsWith("ftp")){
      wasTrsURL = true
      chTrsURL = Channel.from(transcriptome)
    }else{
      chTrsFasta = Channel.fromPath(transcriptome).ifEmpty { exit 1, "Reference annotation file not found: ${transcriptome}" }
    }
  }else{
    log.warn("No transcripts file detected for ${params.genome}")
  }
}

// GTF
gff=false
gtf=false
chGtf = Channel.empty()
chGff = Channel.empty()
chGtfURL = Channel.empty()
chGffURL = Channel.empty()
wasGtfURL = false
wasGffURL = false
if (params.gtf){
  chGtf = Channel.fromPath(params.gtf).ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
}else if (params.gff){
  chGff = Channel.fromPath(params.gff).ifEmpty { exit 1, "GFF annotation file not found: ${params.gtf}" }
}else if (params.genome){
  gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
  gff = params.genome ? params.genomes[ params.genome ].gff ?: false : false

  if ( gtf && gff ){
    log.info "Both GTF and GFF have been provided: Using GTF."
  }else if ( gtf ){
    if (gtf.startsWith("http") || gtf.startsWith("ftp")){
      wasGtfURL = true
      chGtfURL = Channel.from(gtf)
    }else{
      chGtf = Channel.fromPath(gtf).ifEmpty { exit 1, "GTF annotation file not found: ${gtf}" }
    }
  }else if ( gff ){
    if ( gff.startsWith("http") || gff.startsWith("ftp")){
      wasGffURL = true
      chGffURL = Channel.from(gff)
    }else{
      chGff = Channel.fromPath(gff).ifEmpty { exit 1, "GFF annotation file not found: ${gff}" }
    }
  }else{
    log.warn("No GTF/GFF information detected for ${params.genome}")
  }
}

if ( params.build ){
  build = params.build
}else if (params.genome){
  build = params.genome
}else{
  fafile = file(params.fasta)
  build = fafile.baseName - ~/(\.fa)?(\.fasta)?(\.gz)?$/
}
chBuild = Channel.of(build)

/*
===========================
   SUMMARY
===========================
*/

summary = [
  'Pipeline Release': workflow.revision ?: null,
  'Build': build,
  'Fasta' : params.fasta ?: params.genome ? fasta : "",
  'Transcripts' : params.transcripts ?: params.genome ? transcriptome : null,
  'Gtf' : params.gtf ?: params.genome ? gtf : null,
  'Gff' : params.gff ?: params.genome ? gff : null,
  'Indexes' : aligners,
  'Max Resources': "${params.maxMemory} memory, ${params.maxCpus} cpus, ${params.maxTime} time per job",
  'Container': workflow.containerEngine && workflow.container ? "${workflow.containerEngine} - ${workflow.container}" : null,
  'Profile' : workflow.profile,
  'OutDir' : params.outDir,
  'WorkDir': workflow.workDir,
  'CommandLine': workflow.commandLine
].findAll{ it.value != null }

workflowSummaryCh = NFTools.summarize(summary, workflow, params)

/*
==================================
           INCLUDE
==================================
*/ 

// Workflows
include { fastaProcessingFlow } from './nf-modules/local/subworkflow/fastaProcessing'
include { gtfProcessingFlow } from './nf-modules/local/subworkflow/gtfProcessing'
include { hisat2IndexFlow } from './nf-modules/local/subworkflow/hisat2Index' 
include { cellrangerIndexFlow } from './nf-modules/local/subworkflow/cellrangerIndex'
include { kallistoIndexFlow } from './nf-modules/local/subworkflow/kallistoIndex'

// Processes
include { getSoftwareVersions } from './nf-modules/common/process/utils/getSoftwareVersions'
include { outputDocumentation } from './nf-modules/common/process/utils/outputDocumentation'
include { getDataFromURL as getTranscriptome } from './nf-modules/common/process/custom/getDataFromURL'
include { getDataFromURL as getAnnotation } from './nf-modules/common/process/custom/getDataFromURL'
include { getFastaFromURL } from './nf-modules/common/process/custom/getFastaFromURL'
include { gffreads } from './nf-modules/common/process/gffreads/gffreads'
include { bowtie2Index } from './nf-modules/common/process/bowtie2/bowtie2Index'
include { starIndex } from './nf-modules/common/process/star/starIndex'
include { bwaIndex } from './nf-modules/common/process/bwa/bwaIndex'
include { salmonIndex } from './nf-modules/common/process/salmon/salmonIndex'

/*
=====================================
            WORKFLOW 
=====================================
*/

workflow {
  versionsCh = Channel.empty()

  main:

    //************************************
    // Download

    getFastaFromURL(
      chFastaURL
    )
    chFasta = wasFastaURL ? getFastaFromURL.out.fasta : chFasta
    
    getTranscriptome(
      chTrsURL,
      Channel.of('fasta')
    )
    chTrsFasta = wasTrsURL ? getTranscriptome.out.output : chTrsFasta

    getAnnotation(
      chGffURL.concat(chGtfURL),
      Channel.of('gtf,gff')
    )
    if (wasGtfURL){
      chGtf = getAnnotation.out.output
    }else if (wasGffUrl){
      chGff = getAnnotation.out.output
    }

    if ((params.gff || gff)  && (!params.gtf && !gtf)){
      gffreads(
        getAnnotation.out.output
      )
      chGtf = gffreads.out.gtf
    }

    //**********************************
    // Fasta processing

    fastaProcessingFlow(
      chFasta,
      chBuild
    )

    //**********************************
    // Generate genome indexes

    bwaIndex(
      chFasta
    )

    starIndex(
      chFasta,
      fastaProcessingFlow.out.chromSizes,
      Channel.value([])
    )

    bowtie2Index(
      chFasta
    )

    hisat2IndexFlow(
      chFasta,
      chGtf
    )

    cellrangerIndexFlow(
      chFasta,
      chGtf,
      chBuild
    )

    kallistoIndexFlow(
      chTrsFasta
    )

    salmonIndex(
      chFasta,
      chTrsFasta
    )

    //**********************************
    // Process GTF file

    gtfProcessingFlow(
      chGtf
    )

    // subroutines
    //outputDocumentation(
    //  outputDocsCh,
    //  outputDocsImagesCh
    //)
}

workflow.onComplete {
  NFTools.makeReports(workflow, params, summary, customRunName, mqcReport)
}
