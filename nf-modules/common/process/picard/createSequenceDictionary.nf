/*
 * Picard CreateSequenceDictionary
 */

process createSequenceDictionary {
  label 'picard'
  label 'minCpu'
  label 'medMem'

  input:
  tuple path(fasta), path(fai)

  output:
  path('*.dict'), emit: dict
  path('versions.txt'), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: fasta.toString() - ~/(\.fa)?(\.fasta)?$/
  def javaArgs = task.ext.args ?: ''
  markdupMemOption = "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
  """
  echo \$(picard CreateSequenceDictionary --version 2>&1 | sed -e 's/Version:/picard /') > versions.txt
  picard ${markdupMemOption} ${javaArgs} CreateSequenceDictionary \\
    REFERENCE=${fasta} \\
    OUTPUT=${prefix}.dict
  """
}
