/* 
 * Build indexes for cellRanger
 */

include { cellrangerMkGtf } from '../../common/process/cellranger/cellrangerMkGtf'
include { cellrangerMkRef } from '../../common/process/cellranger/cellrangerMkRef'

workflow cellrangerIndexFlow {

  take:
  fasta
  gtf
  build

  main:
  chVersions = Channel.empty()

  cellrangerMkGtf(
    gtf
  )
  chVersions = chVersions.mix(cellrangerMkGtf.out.versions)

  cellrangerMkRef(
    fasta,
    cellrangerMkGtf.out.gtf,
    build
  )
  chVersions = chVersions.mix(cellrangerMkRef.out.versions)

  emit:
  gtf = cellrangerMkGtf.out.gtf
  index = cellrangerMkRef.out.index
  versions = chVersions
}
