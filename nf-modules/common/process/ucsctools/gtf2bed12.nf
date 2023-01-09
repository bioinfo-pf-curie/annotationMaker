/*
 * GTF to BED12 convertion
 */

process gtf2bed12 {
  label 'ucsctools'
  label 'lowCpu'
  label 'lowMem'

  input:
  path(gtf)

  output:
  path("*.bed12"), emit: bed12

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: ${gtf.baseName}
  def args = task.ext.args ?: ''
  """
  gtfToGenePred ${args} ${gtf} ${prefix}.genepred 2> genepred.log
  genePredToBed ${prefix}.genepred ${prefix}.bed12
  """
}

