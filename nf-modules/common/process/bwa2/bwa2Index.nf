
/*
 * Build genome index for BWA mapping
 */

process bwa2Index{
  label 'bwa2'
  label 'highCpu'
  label 'highMem'

  input:
  path(fasta)

  output:
  path("bwa2"), emit: index
  path("*.log"), emit: logs
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: fasta.toString() - ~/(\.fa)?(\.fasta)?$/
  """
  mkdir -p bwa2
  bwa-mem2 index -p bwa2/${prefix} ${fasta} > bwa2.log 2>&1
  echo "Bwa-mem2 "\$(bwa2 2>&1 | grep Version | cut -d" " -f2) &> versions.txt
  """
}