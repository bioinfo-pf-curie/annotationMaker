/*
 * Download fasta annotation file
 */

process getFastaFromURL {
  label 'unix'
  label 'lowCpu'
  label 'lowMem'

  input:
  val(url)

  output:
  path("*.{fa,fasta}"), emit: fasta

  script:
  if (url.endsWith(".tar.gz")){
  """
  wget --no-check-certificate ${url} -O chromFa.tar.gz

  mkdir ./tmp
  tar zxvf chromFa.tar.gz -C ./tmp
  for i in \$(find tmp/ -name "*.fa" | grep -v "_" | sort -V); do cat \$i >> ${build}.fa; done
  for i in \$(find tmp/ -name "*.fa" | grep "_" | sort -V); do cat \$i >> ${build}.fa; done
  rm -rf ./tmp
  """
  }else if (url.endsWith(".gz")){
  """
  wget --no-check-certificate ${url}
  gunzip *.gz
  """
  }else{
  """
  wget --no-check-certificate ${url}
  """
  }
}
