process reduceGtf {
  label 'unix'
  label 'lowCpu'
  label 'lowMem'

  input:
  path(gtf)

  output:
  path("*proteinCoding.gtf") optional true, emit: gtf

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: ${gtf.baseName}
  """
  nbPC=\$(head -n 1000 | awk '\$0~"gene_biotype \\"protein_coding\\"" || \$0~"gene_type \\"protein_coding\\"" {print}' ${gtf} | wc -l)
  if [[ \$nbPC -gt 0 ]]; then
    awk '\$0~"gene_biotype \\"protein_coding\\"" || \$0~"gene_type \\"protein_coding\\"" {print}' ${gtf} > ${prefix}_proteinCoding.gtf
  fi
  """
}
