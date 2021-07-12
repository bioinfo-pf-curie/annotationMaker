#!/usr/bin/env Rscript
.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
      stop("Usage: extracGeneFromGTF <input> <output>", call.=FALSE)
  }
input <- args[1]
output <- args[2]

library(GenomicFeatures)
library(rtracklayer)

## Make TxDB
x <- import(input)
## Remove genes without gene_id which make MakeTxDb crashing
x <- x[which(!is.na(x$gene_id))]
## Bild db
txdb <- makeTxDbFromGRanges(x)
g <- genes(txdb)

## Sort genes
g <- sort(g)

## Export in BED format
export(g, con=output, format="BED")

## Update annotation
d.gtf <- rtracklayer::import(input)
if ( is.element("gene_name", colnames(elementMetadata(d.gtf))) ){
    annot <- unique(elementMetadata(d.gtf)[,c("gene_id", "gene_name")])
    names(g) <- annot$gene_name[match(names(g), annot$gene_id)]
}

## Export in BED format
export(g, con=output, format="BED")
