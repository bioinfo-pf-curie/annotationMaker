/*
 * Convert gff to gtf with gffreads
 */

process gffreads {
  label 'gffread'
  label 'lowCpu'
  label 'lowMem'

  when:
  !params.skipGtfProcessing

  input:
  path(gff)

  output:
  path("*.gtf"), emit: gtf

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${gff.baseName}"  
  """
  gffread $gff --keep-exon-attrs -F -T -o ${prefix}.gtf
  """
}
