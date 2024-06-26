% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bamTools.R
\name{vl_importBamRaw}
\alias{vl_importBamRaw}
\title{Title}
\usage{
vl_importBamRaw(file, extra_arg, headN)
}
\arguments{
\item{file}{bam file path}

\item{extra_arg}{Extra arg to be passed to samtools view}

\item{headN}{Number of starting lines to import}
}
\value{
Imported bam file
}
\description{
Uses samtools to import a bam file
}
\examples{
vl_importBamRaw("path/to/bam/file.bam")
}
