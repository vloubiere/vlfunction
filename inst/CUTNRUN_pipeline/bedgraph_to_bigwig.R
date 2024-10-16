#!/usr/bin/env Rscript

# Check arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript script.R input_file.bdg output_file.bw genome")
}

bdg_file <- args[[1]]
bw_file <- args[[2]]
genome <- args[[3]]

# Load required packages
suppressMessages(library(rtracklayer, warn.conflicts = FALSE))
suppressMessages(library(GenomeInfoDb, warn.conflicts = FALSE))
suppressMessages(library(BSgenome, warn.conflicts = FALSE))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10, warn.conflicts = FALSE))

# Check file extensions
if (!grepl("\\.bdg$", bdg_file) && grepl("\\.bw$", bw_file)) 
  stop("Please provide a .bdg file as input and a .bw file as output.")

# Import bedgraph
gr <- rtracklayer::import(bdg_file, format = "bedGraph")
seqlengths(gr) <- if(genome=="mm10")
  GenomeInfoDb::seqlengths(BSgenome.Mmusculus.UCSC.mm10)[seqlevels(gr)]else if(genome=="hg38")
    GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(gr)] else
      stop("Genome not supported")

# Export to .bw file
rtracklayer::export(gr, bw_file, format = "bigWig")