#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
      stop("Usage: extracGeneFromGTF <input> <format <output>", call.=FALSE)
  }
input <- args[1]
format <- args[2]
output <- args[3]

library(GenomicFeatures)
library(rtracklayer)

## Make TxDB
txdb <- makeTxDbFromGFF(input, format=format)
g <- genes(txdb)

## Sort genes
g <- sort(g)

## Export in BED format
export(g, con=output, format="BED")

