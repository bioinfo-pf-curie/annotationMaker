/* 
 * Fasta processing workflow
 */

include { samtoolsFaidx } from '../../common/process/samtools/samtoolsFaidx'
include { createSequenceDictionary } from '../../common/process/picard/createSequenceDictionary'
include { getChromSizes } from '../../common/process/custom/getChromSizes'
include { faSize } from '../../common/process/ucsctools/faSize'

workflow fastaProcessingFlow {

  take:
  fasta
  build

  main:
  
  chVersions = Channel.empty()

  samtoolsFaidx(
    fasta
  )
  chFastaFai = fasta.combine(samtoolsFaidx.out.fai)
  chVersions = chVersions.mix(samtoolsFaidx.out.versions)

  createSequenceDictionary(
    chFastaFai
  )
  chVersions = chVersions.mix(createSequenceDictionary.out.versions)
  
  getChromSizes(
    chFastaFai,
    build
  )

  faSize(
    chFastaFai,
    build
  )

  emit:
  fasta = chFastaFai
  dict = createSequenceDictionary.out.dict
  chromSizes = getChromSizes.out.sizes
  effgsize = faSize.out.effgsize
  versions = chVersions
}
