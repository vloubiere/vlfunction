% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAseq_pipeline.R
\name{vl_RNAseq_DESeq2}
\alias{vl_RNAseq_DESeq2}
\title{DESeq2 analysis RNAseq pipeline}
\usage{
vl_RNAseq_DESeq2(
  processed_metadata,
  FC_metadata_output,
  Rpath = "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
  dds_output_folder = "db/dds/RNAseq/",
  FC_output_folder = "db/FC_tables/RNAseq/",
  PDF_output_folder = "pdf/RNAseq/",
  cores = 2,
  mem = 8,
  submit = FALSE,
  wdir = getwd(),
  logs = "db/logs/RNAseq/DESeq2",
  time = "01:00:00"
)
}
\arguments{
\item{processed_metadata}{Path to the metadata file generated by vl_PROseq_processing (in .rds or .txt format), or the corresponding data.table.}

\item{FC_metadata_output}{An .rds path where to save the metadata file, which contains the directories containing .dds (DESeq2 objects) and FC table files and will be used to manage them.
By default, when the processed_metadata is a path to a processed_metadata files, "_FC_tables.rds" will be appended to the processed_metadata file path.}

\item{Rpath}{Path to an Rscript executable. Default= "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript".}

\item{dds_output_folder}{Output folder for .dds files (DESeq2 objects). Default= "db/dds/RNAseq/".}

\item{FC_output_folder}{Output folder for FC tables. Default= "db/FC_tables/RNAseq/".}

\item{PDF_output_folder}{Output folder for .pdf files of MA plots and statistics. Default= "pdf/RNAseq/".}

\item{cores}{Number of cores per job. Default= 4.}

\item{mem}{Memory per job (in Go). Default= 16.}

\item{submit}{Should the command be submitted? Default= FALSE.}

\item{wdir}{The working directory to use. Default= getwd(), meaning current working directory will be used.}

\item{logs}{Output folder for log files. Default= "db/logs/CUTNRUN/peak_calling".}

\item{time}{The time required for the SLURM scheduler. Default= '1-00:00:00'.}
}
\value{
Differential analyses using DESeq2.
}
\description{
This function uses the provided metadata file to run all DESeq2 comparisons. By modifying this input file, one can test several experiment designs is needed.
}
\examples{
# Processed metadata ----
processed <- readRDS("Rdata/metadata_RNAseq_processed.rds")

# Before starting, one can modify the processed metadata's design
processed$experiment
processed$DESeq2_control

# Differential analysis ----
vl_RNAseq_DESeq2(processed_metadata = processed,
                 FC_metadata_output = "Rdata/metadata_RNAseq_FC_tables.rds",
                 submit = TRUE)

# FC tables ----
FC <- readRDS("Rdata/metadata_RNAseq_FC_tables.rds")

# Import 6h depletion ----
IAA_6h <- fread(FC[DESeq2_condition=="IAA_6h", fcTable])

# There should be 23 down-regulated and 1 up-regulated genes ----
table(IAA_6h$diff)

}
