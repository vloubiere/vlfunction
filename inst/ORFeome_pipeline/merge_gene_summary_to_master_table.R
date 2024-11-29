#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Please specify:\n
       [required] 1/ Gene summary output file from mageck\n
       [required] 2/ output table .csv\n")
}

require(data.table)

# parse ----
input <- args[1]
output <- args[2]

# Source function
vl_merge_GS_MS <- function(master= "/groups/stark/pachano/projects/eORFeome/Rdata/Master_eORFeome_Mar24.csv",
                           gene_summary,
                           master.columns= c("id", "Effector class", "RefSeq protein", "UniProt", "Gene", "Description", "Species", "DNA sequence", "size", "Amino acid sequence","Size (aa)", "Uniprot", "Protein names", "Function"),
                           FC.columns= c("id", "num", "pos|rank", "pos|goodsgrna", "pos|score", "pos|lfc", "pos|fdr"),
                           output_path_csv= gsub("gene_summary.txt$", "gene_summary_master.csv", gene_summary))
{
  if(!grepl(".csv$", master))
    stop("master should be a path to a .csv file")
  if(!grepl(".txt$", gene_summary))
    stop("gene_summary should be a path to a .txt file")
  if(!grepl(".csv$", output_path_csv))
    stop("output_path_csv should be a path to a .csv path")
  
  # Import master ----
  master <- fread(master)
  master <- as.data.table(master)
  master[, id:= as.character(id)]
  master <- master[, (master.columns), with= F]
  
  # Import FC table ----
  gene_summary <- fread(gene_summary)
  gene_summary <- gene_summary[, (FC.columns), with= FALSE]
  if(!"id" %in% names(gene_summary))
    stop("gene_summary should contain an id column")
  gene_summary[, first_id:= tstrsplit(id, "_", keep= 1)]
  setnames(gene_summary, function(x) gsub("pos|", "", fixed = TRUE, x))
  
  # Merge
  res <- merge(gene_summary,
               master,
               by.x= "first_id",
               by.y= "id",
               all.x= TRUE)
  # Return
  fwrite(res,
         output_path_csv,
         col.names = T,
         row.names = F,
         na= NA,
         sep = ",",
         quote=F)
}

vl_merge_GS_MS(gene_summary = input,
               output_path_csv= output)