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
                         AnnotationMaker
========================================================================================
 AnnotationMaker Pipeline.
 #### Homepage / Documentation
 https://gitlab.curie.fr/data-analysis/annotationMaker
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    if ("${workflow.manifest.version}" =~ /dev/ ){
       devMess = file("$baseDir/assets/devMessage.txt")
       log.info devMess.text
    }

    log.info """
    annotationMaker v${workflow.manifest.version}
    ======================================================================

    Usage:
    nextflow run main.nf --genome 'hg19' -profile conda
    nextflow run main.nf --fasta '*.fasta' --gtf '*.gtf' -profile conda

    Mandatory arguments:
      --genome [str]                Reference genome name and annotations to use
      --genomeAnnotationPath [dir]  Path to genome annotation folder
      -profile [str]                Configuration profile to use. test / conda / toolsPath / singularity / cluster (see below)

    Optional arguments: If --genome is not specified
      --fasta [file]                Path to input data (must be surrounded with quotes)
      --gtf [file]                  Path to GTF file with gene annotation
      --gff [file]                  Path to GFF file with gene annotation    

    Optional arguments:
      --build [str]                 Build name to use for genome index
      --indexes [file]              List of indexes to build. Available: all,bwa,star,bowtie2,hisat2,none. Default: all

    Other options:
      --skipGtfProcessing [bool]    Skip the GTF file processing
      --outDir [file]               The output directory where the results will be saved
      -w/--work-dir [file]          The temporary directory where intermediate data will be saved
      -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    =======================================================
    Available Profiles
      -profile test                 Run the test dataset
      -profile conda                Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
      -profile multiconda           Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
      -profile path                 Use the installation path defined for all tools. Use `--globalPath` to define the insallation path
      -profile multipath            Use the installation paths defined for each tool. Use `--globalPath` to define the insallation path
      -profile docker               Use the Docker images for each process
      -profile singularity          Use the Singularity images for each process. Use `--singularityPath` to define the insallation path
      -profile cluster              Run the workflow on the cluster, instead of locally
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
        'hisat2',
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
aligners = params.indexes || params.indexes == 'none' ? params.indexes == 'all' ? alignersList : params.indexes.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(aligners, alignersList)) exit 1, 'Unknown Aligner(s), see --help for more information'


/*
 * CHANNELS
 */

// Reference Genome
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
   exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

if (params.genome){
  fastaURL = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
  if (fastaURL){
    Channel.from(fastaURL)
      .ifEmpty { exit 1, "Reference Genome not found: ${fastaURL}" }
      .set { chFastaLink }
  }else{
    exit 1, "Fasta file not found for ${params.genome} : ${fastaURL}"
  }
}else if (params.fasta){
  Channel.fromPath("${params.fasta}")
    .ifEmpty { exit 1, "Reference Genome not found: ${params.fasta}" }
    .into { chFasta; chFastaBwa; chFastaStar; chFastaBowtie2; chFastaHisat2 }
  chFastaLink = Channel.empty()
}

// GTF
if (params.genome){
  params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
  params.gff = params.genome ? params.genomes[ params.genome ].gff ?: false : false

  if (params.gtf && params.gff){
    log.info "Both GTF and GFF have been provided: Using GTF as priority."
  }
  if (params.gtf){
    Channel.from(params.gtf)
      .ifEmpty { exit 1, "Reference annotation not found: ${params.gtf}" }
      .set { chGtfLink }
    chGffLink = Channel.empty()
  }else if (params.gff){
    Channel.from(params.gff)
      .ifEmpty { exit 1, "Reference annotation not found: ${params.gff}" }
      .set { chGffLink }
    chGtfLink = Channel.empty()
  }else{
    log.warn("No GTF/GFF information detected for ${params.genome}")
    Channel.empty().into{ chGtf; chGtfBed12; chGtfGene; chGtfHisat2Splicesites; chGtfHisat2Index; chGffLink; chGtfLink }
  }
}else if (params.gtf){
  Channel
    .fromPath(params.gtf)
    .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
    .into { chGtfHisat2Splicesites; chGtfHisat2Index; chGtfBed12; chGtfGene }
  Channel.empty().into{ chGtfLink; chGffLink }
}else if (params.gff){
  Channel
    .fromPath(params.gff)
    .ifEmpty { exit 1, "GFF annotation file not found: ${params.gtf}" }
    .into { chGff }
  Channel.empty().into{ chGtfLink; chGffLink }
}else{
  Channel.empty().into{ chGtf; chGtfBed12; chGtfGene; chGtfHisat2Splicesites; chGtfHisat2Index; chGffLink; chGtfLink }
}

if ( params.build ){
  build = params.build
}else{
  if (params.genome){
    build = params.genome
  }else{
    fafile = file(params.fasta)
    build = fafile.baseName - ~/(.fa)?(.fasta)?(.gz)?/
  }
}

// Header log info
if ("${workflow.manifest.version}" =~ /dev/ ){
   devMess = file("$baseDir/assets/devMessage.txt")
   log.info devMess.text
}

log.info """=======================================================

Annotation Maker workflow v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Command Line'] = workflow.commandLine
if (params.genome){
summary['Fasta']          = fastaURL
}else{
summary['Fasta']          = params.fasta 
}
summary['Build']          = build
if (params.gtf){
summary['Gtf']            = params.gtf
}
if (params.gff){
summary['Gff']            = params.gff
}
summary['Indexes']        = aligners
summary['Max Memory']     = params.maxMemory
summary['Max CPUs']       = params.maxCpus
summary['Max Time']       = params.maxTime
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outDir
summary['Config Profile'] = workflow.profile

log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 * DOWNLOAD
 */

process getFasta {
  label 'unix'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outDir}/genome", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  input:
  val(url) from chFastaLink

  output:
  file("*.fa") into chFastaURL

  script:
  if (url.endsWith(".tar.gz")){
  """
  wget ${url} -O chromFa.tar.gz
  
  mkdir ./tmp
  tar zxvf chromFa.tar.gz -C ./tmp
  for i in \$(ls tmp/*.fa | grep -v "_" | sort -V); do cat \$i >> ${build}.fa; done
  for i in \$(ls tmp/*.fa | grep "_" | sort -V); do cat \$i >> ${build}.fa; done
  rm -rf ./tmp
  """
  }else if (url.endsWith(".gz")){
  """
  wget ${url} -O ${build}.fa.gz
  gunzip ${build}.fa.gz
  """
  }else{
  """
  wget ${url} -O ${build}.fa
  """
  }
}

process getAnnotation {
  label 'unix'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outDir}/gtf", mode: 'copy'

  input:
  val(url) from chGffLink.concat(chGtfLink).dump(tag:'annot')

  output:
  file("*.{gtf,gff}") into chAnnotURL

  script:
  if (url.endsWith(".gz")){
  """
  wget ${url} 
  gunzip *.gz
  """
  }else{
  """
  wget ${url} 
  """
  }
}

if (params.genome){
  chFastaURL.into{chFasta; chFastaBwa; chFastaStar; chFastaBowtie2; chFastaHisat2}
  if (params.gtf){
    chAnnotURL.into{chGtfHisat2Splicesites; chGtfHisat2Index; chGtf; chGtfBed12; chGtfGene}
  }else if (params.gff){
    chAnnotURL.set{chGff}
  }
}


/*
 * GTF/GFF processing
 */

if (params.gff && !params.gtf) {
  process convertGFFtoGTF {
    label 'gffread'
    label 'lowCpu'
    label 'lowMem'

    publishDir "${params.outDir}/gtf", mode: 'copy'

    input:
    file gff from chGff

    output:
    file "${gff.baseName}.gtf" into chGtfHisat2Splicesites, chGtfHisat2Index, chGtf, chGtfBed12, chGtfGene

    script:
    """
    gffread $gff --keep-exon-attrs -F -T -o ${gff.baseName}.gtf
    """
  }
}

/*
 * FASTA PROCESSING
 */

process indexFasta {
  label 'samtools'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outDir}/genome", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  input:
  file(fasta) from chFasta

  output:
  set file(fasta), file("*.fai") into chFastaDict, chFastaSize, chFastaEffgsize 

  script:
  """
  samtools faidx ${fasta}
  """
}

process makeDict {
  label 'picard'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outDir}/genome", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  input:
  set file(fasta), file(faidx) from chFastaDict

  output:
  file("*.dict") into chDict

  script:
  pfix = fasta.toString() - /(.fa)?(.fasta)?/
  """
  picard CreateSequenceDictionary REFERENCE=${fasta} OUTPUT=${pfix}.dict
  """
}

process makeChromSizes {
  label 'unix'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outDir}/genome", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  input:
  set file(fasta), file(faidx) from chFastaSize

  output:
  file("*sizes") into (chChromSize, chChromSizeStar)

  script:
  pfix = fasta.toString() - /(.fa)?(.fasta)?/
  """
  cut -f1,2 ${faidx} > chrom_${pfix}.sizes
  """
}

process effectiveGenomeSize {
  label 'ucsctools'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outDir}/genome", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  input:
  set file(fasta), file(faidx) from chFastaEffgsize

  output:
  file("*effgsize") into chEffSize

  script:
  pfix = fasta.toString() - /(.fa)?(.fasta)?/
  """
  faSize ${fasta} > ${pfix}_fasize.log
  awk -v genome=${pfix} 'NR==1{print genome"\t"\$5}' ${pfix}_fasize.log > ${pfix}.effgsize
  """
}

/*
 * PREPROCESSING - Build Index
 */ 

process makeBwaIndex {
  label 'bwamem'
  label 'medCpu'
  label 'medMem'

  publishDir "${params.outDir}/indexes/", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  when: 
  'bwa' in aligners
   
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
  label 'star'
  label 'medCpu'
  label 'highMem'

  publishDir "${params.outDir}/indexes/", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  when:
  'star' in aligners

  input:
  file(fasta) from chFastaStar
  file(chrSize) from chChromSizeStar

  output:
  file "STAR" into chStarIdx

  script:
  """
  odir=\$(STAR --version | cut -d_ -f1,2)
  mkdir -p \${odir}
  chrBinNbits=\$(awk -F"\t" '{s+=\$2;l+=1}END{p=log(s/l)/log(2); printf("%.0f", (p<18 ? p:18))}' ${chrSize})
  STAR --runMode genomeGenerate --limitGenomeGenerateRAM 33524399488 --genomeChrBinNbits \${chrBinNbits} --runThreadN ${task.cpus} --genomeDir ${odir} --genomeFastaFiles $fasta
  """
}

process makeBowtie2Index {
  label 'bowtie2'
  label 'medCpu'
  label 'highMem' 

  publishDir "${params.outDir}/indexes/", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  when:
  'bowtie2' in aligners

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
  label 'hisat2'
  label 'lowCpu'
  label 'medMem'

  publishDir "${params.outDir}/indexes/hisat2", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  when:
  'hisat2' in aligners

  input:
  file gtf from chGtfHisat2Splicesites.collect()

  output:
  file "${gtf.baseName}.hisat2_splice_sites.txt" into indexingSplicesites

  script:
  """
  mkdir -p hisat2
  hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2_splice_sites.txt
  """
}

process makeHisat2Index {
  label 'hisat2'
  label 'highCpu'
  label 'hugeMem'

  publishDir "${params.outDir}/indexes/", mode: 'copy',
    saveAs: {filename -> if (filename.indexOf(".log") > 0) "logs/$filename" else filename}

  when:
  'hisat2' in aligners

  input:
  file fasta from chFastaHisat2
  file indexing_splicesites from indexingSplicesites.collect()
  file gtf from chGtfHisat2Index.collect()

  output:
  file("hisat2") into chHisat2Idx

  script:
  """
  mkdir -p hisat2
  hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt
  hisat2-build -p ${task.cpus} --ss $indexing_splicesites --exon ${gtf.baseName}.hisat2_exons.txt $fasta hisat2/${build}
  """
}


/***********************
 * GTF Annotation
 */

// Extract protein-coding genes only
process reduceGtf {
  label 'unix'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outDir}/gtf", mode: 'copy',

  when:
  !params.skipGtfProcessing

  input:
  file(gtf) from chGtf

  output:
  file("*proteinCoding.gtf") optional true into (chGtfReducedBed12, chGtfReducedGene)

  script:
  """
  nbPC=\$(head -n 1000 | awk '\$0~"gene_biotype \\"protein_coding\\"" || \$0~"gene_type \\"protein_coding\\"" {print}' ${gtf} | wc -l)
  if [[ \$nbPC -gt 0 ]]; then
    awk '\$0~"gene_biotype \\"protein_coding\\"" || \$0~"gene_type \\"protein_coding\\"" {print}' ${gtf} > ${gtf.baseName}_proteinCoding.gtf
  fi
  """
}

// bed12 file are transcripts-based annotation files
process gtf2bed12 {
  label 'ucsctools'
  label 'lowCpu'
  label 'lowMem'

  publishDir "${params.outDir}/gtf", mode: 'copy',

  when:
  !params.skipGtfProcessing

  input:
  file(gtf) from chGtfBed12.concat(chGtfReducedBed12)

  output:
  file("*.bed12") into chBed12

  script:
  """
  gtfToGenePred -genePredExt -geneNameAsName2 -allErrors -ignoreGroupsWithoutExons ${gtf} ${gtf.baseName}.genepred 2> genepred.log
  genePredToBed ${gtf.baseName}.genepred ${gtf.baseName}.bed12
  """
}

// Extract gene coordinates
process gtf2genes {
  label 'r'
  label 'lowCpu'
  label 'medMem'

  publishDir "${params.outDir}/gtf", mode: 'copy', 

  when:
  !params.skipGtfProcessing

  input:
  file(gtf) from chGtfGene.concat(chGtfReducedGene)

  output:
  file("*_gene.bed") into chGeneBed

  script:
  """
  extractGeneFromGTF.r ${gtf} ${gtf.baseName}_gene.bed
  sort -k1,1V -k2,2n ${gtf.baseName}_gene.bed > ${gtf.baseName}_gene_sorted.bed
  mv ${gtf.baseName}_gene_sorted.bed ${gtf.baseName}_gene.bed
  """  
}


/*
 * Sub-routine
 */

workflow.onComplete {

    /*pipeline_report.html*/

    def report_fields = [:]
    report_fields['pipeline'] = workflow.manifest.name
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
    def tf = new File("$baseDir/assets/onCompleteTemplate.txt")
    def txt_template = engine.createTemplate(tf).make(report_fields)
    def report_txt = txt_template.toString()
    
    // Render the HTML template
    def hf = new File("$baseDir/assets/onCompleteTemplate.html")
    def html_template = engine.createTemplate(hf).make(report_fields)
    def report_html = html_template.toString()

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outDir}/pipelineInfo/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipelineReport.html" )
    output_hf.withWriter { w -> w << report_html }
    def output_tf = new File( output_d, "pipelineReport.txt" )
    output_tf.withWriter { w -> w << report_txt }

    /*oncomplete file*/
    File woc = new File("${params.outDir}/workflowOnComplete.txt")
    Map endSummary = [:]
    endSummary['Completed on'] = workflow.complete
    endSummary['Duration']     = workflow.duration
    endSummary['Success']      = workflow.success
    endSummary['exit status']  = workflow.exitStatus
    endSummary['Error report'] = workflow.errorReport ?: '-'
    String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
    println endWfSummary
    String execInfo = "Execution summary\n${endWfSummary}\n"
    woc.write(execInfo)

    /*final logs*/
    if(workflow.success){
      log.info "[annotationMaker] Pipeline Complete"
    }else{
      log.info "[annotationMaker] FAILED: $workflow.runName"
    } 
}
