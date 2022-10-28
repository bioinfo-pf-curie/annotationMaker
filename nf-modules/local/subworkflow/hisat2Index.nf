/* 
 * Fasta processing workflow
 */

include { makeHisat2Splicesites } from '../../common/process/hisat2/makeHisat2Splicesites'
include { hisat2Index } from '../../common/process/hisat2/hisat2Index'

workflow hisat2IndexFlow {

  take:
  fasta
  gtf

  main:
  chVersions = Channel.empty()

  makeHisat2Splicesites(
    gtf
  )
  chVersions = chVersions.mix(makeHisat2Splicesites.out.versions)

  hisat2Index(
    fasta,
    makeHisat2Splicesites.out.alignmentSplicesites,
    gtf
  )
  chVersions = chVersions.mix(hisat2Index.out.versions)

  emit:
  index = hisat2Index.out.index
  versions = chVersions
}
