/*
 * Build genome index for dragmap mapping
 */

process dragmapIndex{
  label 'dragmap'
  label 'highCpu'
  label 'highMem'

  input:
  path(fasta)

  output:
  path("dragmap"), emit: index
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
  echo "bwa-mem2 "\$(bwa 2>&1 | grep Version | cut -d" " -f2) &> versions.txt
  """
// modifier commande
}
