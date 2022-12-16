
/*
 * Build genome index for BWA mapping
 */

process bwamem2Index{
  label 'bwamem2'
  label 'highCpu'
  label 'highMem'

  input:
  path(fasta)

  output:
  path("bwamem2"), emit: index
  path("*.log"), emit: logs
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: fasta.toString() - ~/(\.fa)?(\.fasta)?$/
  """
  mkdir -p bwamem2
  bwa-mem2 index -p bwamem2/${prefix} ${fasta} > bwamem2.log 2>&1
  echo "Bwa-mem2 "\$(bwa-mem2 version 2>&1 | sed -n 3p) &> versions.txt
  """
}