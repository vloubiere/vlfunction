% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bedTools.R
\name{vl_subtractBed}
\alias{vl_subtractBed}
\title{Subtract bed coverage}
\usage{
vl_subtractBed(a, b, ignore.strand = TRUE)
}
\arguments{
\item{a}{Ranges for which overlaps with b have to be removed.}

\item{b}{Regions to subtract from a.}

\item{ingore.strand}{Should the strand be ignored? If set to FALSE, only the regions in b with the same strand then a will be sutracted. Default= TRUE}
}
\value{
For each range in 'a', reports the number of overlapping features in 'b'
}
\description{
Substracts regions in b to regions in a
}
\examples{
a <- data.table(seqnames= "chr3R", start= 1, end= 1000)
b <- data.table(seqnames= c("chr3R", "chr3R", "chrX"),
                start= c(100, 500, 100),
                end= c(200, 600, 200))

vl_subtractBed(a, b)

}
