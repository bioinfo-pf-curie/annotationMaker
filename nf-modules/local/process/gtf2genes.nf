process gtf2genes {
  label 'r'
  label 'lowCpu'
  label 'medMem'

  input:
  path(gtf)

  output:
  path("*_gene.bed"), emit: bed

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${gtf.baseName}"
  """
  extractGeneFromGTF.r ${gtf} ${prefix}_gene.bed
  sort -k1,1V -k2,2n ${prefix}_gene.bed > ${prefix}_gene_sorted.bed
  mv ${prefix}_gene_sorted.bed ${prefix}_gene.bed
  """
}
