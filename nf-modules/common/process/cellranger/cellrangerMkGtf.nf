/*
 * cellranger mkgtf
 */

process cellrangerMkGtf {
  label 'unix'
  label 'medCpu'
  label 'medMem'

  input:
  path(gtf)

  output:
  path("cellranger_filtered.gtf"), emit: gtf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "cellranger_filtered"
  """
  cellranger mkgtf $gtf ${prefix}.gtf \
    --attribute=gene_biotype:protein_coding --attribute=gene_biotype:lincRNA \
    --attribute=gene_biotype:antisense --attribute=gene_biotype:IG_LV_gene \
    --attribute=gene_biotype:IG_V_gene --attribute=gene_biotype:IG_V_pseudogene \
    --attribute=gene_biotype:IG_D_gene --attribute=gene_biotype:IG_J_gene \
    --attribute=gene_biotype:IG_J_pseudogene --attribute=gene_biotype:IG_C_gene \
    --attribute=gene_biotype:IG_C_pseudogene --attribute=gene_biotype:TR_V_gene \
    --attribute=gene_biotype:TR_V_pseudogene --attribute=gene_biotype:TR_D_gene \
    --attribute=gene_biotype:TR_J_gene --attribute=gene_biotype:TR_J_pseudogene \
    --attribute=gene_biotype:TR_C_gene
  echo \$(cellranger --version | grep version | sed -e 's/--version (//' -e 's/)//' ) > versions.txt
  """
}
