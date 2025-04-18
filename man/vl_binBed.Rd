% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bedTools.R
\name{vl_binBed}
\alias{vl_binBed}
\title{Bin Genomic Regions}
\usage{
vl_binBed(
  bed,
  nbins = NULL,
  bins.width = NULL,
  steps.width = bins.width,
  ignore.strand = FALSE
)
}
\arguments{
\item{bed}{A data.table containing genomic ranges with columns: seqnames, start, end.}

\item{nbins}{An integer specifying the number of bins to create for each genomic range. If specified, this takes precedence over bins.width.}

\item{bins.width}{An integer specifying the width of each bin. This parameter is used if nbins is not specified.}

\item{steps.width}{An integer specifying the step size between the start positions of consecutive bins. Default= bins.width.}

\item{ignore.strand}{Should the strand be ignored? If TRUE, binning will always start from the leftmost coordinates. Default= FALSE}
}
\value{
A data.table containing the binned genomic regions, with an extra column (binIDX) containing bin indexes.
}
\description{
This function bins genomic regions from a BED file based on either a specified number 
of bins or a specified bin width and step size. It supports both fixed-width bins 
and sliding window approaches.
}
\details{
The function uses either the nbins parameter to divide each genomic region into a fixed number of bins, or the bins.width and steps.width parameters to create bins of a specific width with a defined step size between them. An extra binIDX column contains unique indexes for each bin created, depending on its stran if ignore.strand is set to FALSE (default). If set to TRUE, then the binIDX will increase with start coordinates.
}
\examples{
# Example BED data
library(data.table)
bed <- data.table(seqnames= "chr2L",
                  start= 101,
                  end= c(200, 210))

# Bin using a specified number of bins
vl_binBed(bed, nbins = 5)

# Bin using a specified bin width and step size
vl_binBed(bed, bins.width = 50, steps.width = 25)

}
