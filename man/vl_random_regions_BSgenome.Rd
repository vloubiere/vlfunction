% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BSgenomeTools.R
\name{vl_random_regions_BSgenome}
\alias{vl_random_regions_BSgenome}
\title{random region}
\usage{
vl_random_regions_BSgenome(genome, n, width = 1, restrict.seqnames = NULL)
}
\arguments{
\item{n}{Number of regions to sample}

\item{width}{Widths of regions to sample (length should either be 1 or equal n)}

\item{restrict.seqnames}{If specified, only the provided seqnames will be used}
}
\value{
data.table containing randomly sampled regions
}
\description{
Sample random regions from BSgenome
}
\examples{
vl_control_regions_BSgenome("dm3", 100, 1)

}
