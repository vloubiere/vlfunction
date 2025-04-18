% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChIP_CUT&RUN_pipeline.R
\name{vl_CUTNRUN_peakCalling}
\alias{vl_CUTNRUN_peakCalling}
\title{Function to call peaks and compute bw tracks from CUTNRUN/CHIP data}
\usage{
vl_CUTNRUN_peakCalling(
  processed_metadata,
  peaks_metadata_output,
  peaks_output_folder = "db/peaks/CUTNRUN/",
  bw_output_folder = "db/bw/CUTNRUN/",
  Rpath = "/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript",
  extsize = 300,
  cores = 8,
  mem = 64,
  overwrite = FALSE,
  submit = FALSE,
  wdir = getwd(),
  logs = "db/logs/CUTNRUN/peak_calling",
  time = "1-00:00:00"
)
}
\arguments{
\item{processed_metadata}{Path to the metadata file generated by vl_CUTNRUN_processing, or the corresponding data.table.}

\item{peaks_metadata_output}{An .rds path where to save the metadata file, which contains the directories containing peaks filesand will be used to manage them. By default, when the processed_metadata is a path to a processed_metadata files, "_peaks_files.rds" will be appended to the processed_metadata file path.}

\item{peaks_output_folder}{Output folder for peak files. Default= "db/peaks/CUTNRUN/".}

\item{bw_output_folder}{Output folder for bw files. Default= "db/bw/CUTNRUN/".}

\item{Rpath}{The path to the Rscript executable to use, on which the latest version of the vlfunction package should be installed. Default= /software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript}

\item{cores}{Number of cores per job. Default= 8.}

\item{mem}{Memory per job (in Go). Default= 32.}

\item{overwrite}{Should existing files be overwritten?}

\item{submit}{Should the command be submitted? default= FALSE.}

\item{wdir}{The working directory to use. Default= getwd(), meaning current working directory will be used.}

\item{logs}{Output folder for log files. Default= "db/logs/CUTNRUN/peak_calling".}

\item{time}{The time required for the SLURM scheduler. Default= '1-00:00:00'.}

\item{extsise}{The extsize to be used, meaning that the building of the peak model will be skipped (--nomodel). To build the model, set extsize= NA. Default= 300.}
}
\value{
Command lines for peak calling.
}
\description{
Function to call peaks and compute bw tracks from CUTNRUN/CHIP data
}
\examples{
# Peak calling example dataset ----
library(vlfunctions)

# Processed metadata ----
processed <- readRDS("Rdata/metadata_CutNRun_processed.rds")

# Peak calling ----
vl_CUTNRUN_peakCalling(processed_metadata = processed,
                       extsize = 300,
                       cores = 8,
                       mem = 64,
                       overwrite = FALSE,
                       submit = TRUE)
# Peaks metadata ----
peaks <- readRDS("Rdata/metadata_CutNRun_peaks_files.rds")

}
