/*
 * Calculate effective size
 */

process faSize {
  label 'ucsctools'
  label 'lowCpu'
  label 'lowMem'

  input:
  tuple path(fasta), path(faidx)
  val(build)

  output:
  path("*effgsize"), emit: effgsize

  script:
  """
  faSize ${fasta} > ${build}_fasize.log
  awk -v genome=${build} 'NR==1{print genome"\t"\$5}' ${build}_fasize.log > ${build}.effgsize
  """
}
