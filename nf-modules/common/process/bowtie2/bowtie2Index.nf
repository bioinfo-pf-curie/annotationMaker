/*
 * Alignment on reference genome with Bowtie2
 */

process bowtie2Index{
  label 'bowtie2'
  label 'highCpu'
  label 'highMem'

  input:
  path(fasta)

  output:
  path("bowtie2"), emit: index
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: fasta.toString() - ~/(\.fa)?(\.fasta)?$/
  """
  mkdir -p bowtie2
  bowtie2-build ${fasta} bowtie2/${prefix}
  echo \$(bowtie2 --version | awk 'NR==1{print "bowtie2 "\$3}') > versions.txt
  """
}


