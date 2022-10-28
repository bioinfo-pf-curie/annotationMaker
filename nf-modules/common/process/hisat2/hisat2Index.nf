/*
 * Build Hisat2 mapping indexes
 */

process hisat2Index {
  label 'hisat2'
  label 'highCpu'
  label 'extraMem'

  input:
  path(fasta)
  path(indexing_splicesites)
  path(gtf)

  output:
  path("hisat2"), emit: index
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: fasta.toString() - ~/(\.fa)?(\.fasta)?$/
  """
  mkdir -p hisat2
  hisat2_extract_exons.py $gtf > ${gtf.baseName}.hisat2_exons.txt
  hisat2-build -p ${task.cpus} --ss $indexing_splicesites --exon ${gtf.baseName}.hisat2_exons.txt $fasta hisat2/${prefix}
  echo \$(hisat2 --version | awk 'NR==1{print "hisat2 "\$3}') > versions.txt
  """
}
  
  
