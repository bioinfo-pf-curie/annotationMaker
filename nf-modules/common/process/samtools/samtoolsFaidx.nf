/*
 * Samtools - Index
 */

process samtoolsFaidx {
  label 'samtools'
  label 'minCpu'
  label 'lowMem'
 
  input:
  path(fasta)

  output:
  path("*fai"), emit: fai
  path("versions.txt") , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  echo \$(samtools --version | head -1) > versions.txt
  samtools faidx ${fasta}
  """
}
    