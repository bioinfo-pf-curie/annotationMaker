/*
 * Build index for Salmon pseudoAlignment
 */

process salmonIndex {
  label "salmon"
  label "medCpu"
  label "extraMem"

  input:
  path(genomeFasta)
  path(transcrpitsFasta)

  output:
  path("salmon_${suffix}/"), emit: index
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  suffix=transcrpitsFasta.toString() - ~/([_.])?(transcripts.fa)?(.gz)?$/
  """
  grep "^>" ${genomeFasta} | cut -d " " -f 1 > decoys.txt
  sed -i.bak -e 's/>//g' decoys.txt
  cat ${transcrpitsFasta} ${genomeFasta} > gentrome.fa

  salmon index \
    -t gentrome.fa \
    --decoy decoys.txt \
    -i salmon_${suffix} \
    -p ${task.cpus} \
    --gencode \
    ${args}

  mv gentrome.fa salmon_${suffix}
  gzip salmon_${suffix}/*.fa
  mv decoys.txt salmon_${suffix}
  echo \$(salmon --version 2>&1) > versions.txt
  """
}

