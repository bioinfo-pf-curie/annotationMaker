/*
 * Download data via wget
 */

process getDataFromURL {
  label 'unix'
  label 'lowCpu'
  label 'lowMem'

  input:
  val(url)
  val(extension)

  output:
  path("*.{${extension}}"), emit: output

  script:
  if (url.endsWith(".gz")){
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
