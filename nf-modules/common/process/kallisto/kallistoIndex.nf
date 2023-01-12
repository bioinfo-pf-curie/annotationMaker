/*
 * Build kallisto index for pseudoAlignment
 */

process kallistoIndex {
  label 'kallisto'
  label 'medCpu'
  label 'highMem'

  input:
  path(transcrpitsFasta)

  output:
  path("kallisto_${suffix}"), emit:index
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  suffix=transcrpitsFasta.toString() - ~/([_.])?(_corrected.fa)?(.gz)?$/
  """
  mkdir -p kallisto_${suffix}
  kallisto index -i kallisto_${suffix}/transcriptome.idx $transcrpitsFasta
  echo \$(kallisto version | sed -e 's/, version//') > versions.txt
  """
}
