/*
 * Build kallisto index for pseudoAlignment
 */

process kallistoIndex {
  label 'kallisto'
  label 'medCpu'
  label 'highMem'

  input:
  path(corrTranscrpitsFasta)

  output:
  path("kallisto_${suffix}"), emit:index
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  suffix=corrTranscrpitsFasta.toString() - ~/([_.])?(corrTranscripts.fa)?(.gz)?$/
  """
  mkdir -p kallisto_${suffix}
  kallisto index -i kallisto_${suffix}/transcriptome.idx $corrTranscrpitsFasta
  echo \$(kallisto version) > versions.txt
  """
}
