/*
 * get chromosomes sizes
 */

process getChromSizes {
  label 'unix'
  label 'lowCpu'
  label 'lowMem'

  input:
  tuple path(fasta), path(faidx)
  val(build)

  output:
  path("*sizes"), emit: sizes

  script:
  """
  cut -f1,2 ${faidx} > chrom_${build}.sizes
  """
}
