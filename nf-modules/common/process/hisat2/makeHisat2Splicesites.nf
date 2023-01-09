/*
 * Hisat2 / makeHisat2Splicesites
 */

process makeHisat2Splicesites {
  label 'hisat2'
  label 'minCpu'
  label 'lowMem'

  input:
  path gtf

  output:
  path "${gtf.baseName}.hisat2SpliceSites.txt", emit: alignmentSplicesites
  path ("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  echo \$(hisat2 --version | awk 'NR==1{print "hisat2 "\$3}') > versions.txt
  hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.hisat2SpliceSites.txt
  """
}