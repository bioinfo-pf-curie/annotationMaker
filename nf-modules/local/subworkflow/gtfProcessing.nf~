/* 
 * Process GTF file
 */

include { reduceGtf } from '../../local/process/reduceGtf'
include { gtf2bed12 } from '../../local/process/gtf2bed12'
include { gtf2genes } from '../../local/process/gtf2genes'

workflow gtfProcessingFlow {

  take:
  gtf

  main:
  chVersions = Channel.empty()

  reduceGtf(
    gtf
  )

  gtf2bed12(
    gtf.concat(reduceGtf.out.gtf)
  )

  gtf2genes(
    gtf.concat(reduceGtf.out.gtf)
  )

  emit:
  reducedGtf = reduceGtf.out.gtf
  bed12 = gtf2bed12.out.bed12
  genes = gtf2genes.out.bed
  versions = chVersions
}
