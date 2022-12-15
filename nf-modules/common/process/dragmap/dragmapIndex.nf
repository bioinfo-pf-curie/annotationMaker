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
  mkdir -p dragmap
  dragen-os --build-hash-table true --ht-reference ${fasta} --output-directory dragmap > dragmap.log 2>&1
  echo "Dragmap "\$(dragen-os -V) &> versions.txt
  """
} 