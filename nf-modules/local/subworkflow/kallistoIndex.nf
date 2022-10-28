/* 
 * Build indexes for kallisto
 */

include { kallistoCorrectFasta } from '../../local/process/kallistoCorrectFasta'
include { kallistoIndex } from '../../common/process/kallisto/kallistoIndex'

workflow kallistoIndexFlow {

  take:
  fasta

  main:
  chVersions = Channel.empty()

  kallistoCorrectFasta(
    fasta
  )

  kallistoIndex(
    kallistoCorrectFasta.out.fasta,
  )
  chVersions = chVersions.mix(kallistoIndex.out.versions)

  emit:
  index = kallistoIndex.out.index
  versions = chVersions
}
