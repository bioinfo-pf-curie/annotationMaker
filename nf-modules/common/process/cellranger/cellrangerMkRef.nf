/*
 * Make cellranger index for reads aligments
 */

process cellrangerMkRef {
  label 'unix'
  label 'highCpu'
  label 'extraMem'

  input:
  path(fasta)
  path(gtf)
  val(build)

  output:
  path("cellranger"), emit: index
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  (full, mem) = (task.memory =~ /(\d+)\s*[a-z]*/)[0]
  """
  cellranger mkref \
    --genome=cellranger \
    --fasta=$fasta \
    --genes=$gtf \
    --nthreads ${task.cpus} --memgb ${mem} --ref-version ${build}
  echo \$(cellranger --version | grep version | sed -e 's/--version (//' -e 's/)//' ) > versions.txt
  """
}
