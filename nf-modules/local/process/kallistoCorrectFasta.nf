/*
 * Correct fasta file for Kallisto
 */

process kallistoCorrectFasta {
  label 'unix'
  label 'medCpu'
  label 'lowMem'

  input:
  path(fasta)

  output:
  path("*_corrected.fa"), emit: fasta
  

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = fasta.toString() - ~/(\.fa)?(\.fasta)?$/
  """
  cat ${fasta} | \
  awk -F  "|" '{if(\$0 ~ /^>.*/){print \$1;} else {print \$0;}}' | \
  sed 's/\\_mRNA//' | \
  grep -vE  ".*_PAR_Y.*" > \
  ${prefix}_corrected.fa
  """
}
