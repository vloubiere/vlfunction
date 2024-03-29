% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MolBioTools.R
\name{vl_digest}
\alias{vl_digest}
\title{Digest sequence}
\usage{
vl_digest(seq, enzyme, keepsite = F)
}
\arguments{
\item{seq}{character sequence to digest}

\item{enzyme}{Character vector containing the enzymes to use}

\item{keepsite}{If TRUE, then the enzymmatic consensu will be paste at cutting sites}
}
\value{
A vector containing the digested pieces
}
\description{
Simulate DNA digest from thermofisher "simple" enzymes
}
\examples{
vl_digest(seq= "AAAAAAAAGGTACCTTTTTTTTTTTTGCGGCCGCAAAAAAAAAA",
enzyme= c("KpnI", "NotI"))

}
