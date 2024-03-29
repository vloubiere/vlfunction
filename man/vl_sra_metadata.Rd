% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SRRTools.R
\name{vl_sra_metadata}
\alias{vl_sra_metadata}
\title{get SRA metadata}
\format{
metadata data.table object
\describe{
  \item{study_name}{GSE id}
  \item{run}{SRR ID}
  \item{library_layout}{"PAIRED - " or "SINGLE - "}
  \item{experiment_title}{experiment_title}
  \item{ftp}{List of ftp download link}
}
}
\usage{
vl_sra_metadata(
  GSE,
  SRAdb = "/groups/stark/vloubiere/exp_data/SRAmetadb.sqlite"
)
}
\arguments{
\item{GSE}{GSE number of interest.}
}
\value{
metadata data.table object
}
\description{
Donwload SRA metadata info from GSE
}
\examples{
vl_sra_metadata(GSE= "GSE119708", SRAdb= "/groups/stark/vloubiere/exp_data/SRAmetadb.sqlite")
}
