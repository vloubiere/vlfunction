% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deepLearning_tools.R
\name{vl_plot_contrib_logo}
\alias{vl_plot_contrib_logo}
\title{Plot contribution scores matrix}
\usage{
vl_plot_contrib_logo(
  bed,
  h5,
  h5.bed = list.files(dirname(h5), ".bed$", full.names = TRUE),
  h5.fa = list.files(dirname(h5), ".fa$", full.names = TRUE),
  genome,
  agg.FUN = function(x) mean(x),
  mot,
  mot.name.column = "motif_ID",
  xlab = "nt",
  ylab = "Contribution",
  xlim,
  ylim
)
}
\arguments{
\item{bed}{A bed file containing a unique region for which contrib scores will be plotted.}

\item{h5}{Path(s) to h5 files containing the contribution scores.}

\item{h5.bed}{Bed files containing the coordinates of the regions corresponding to provided h5 files.}

\item{genome}{The genome to be used.}

\item{agg.FUN}{In the case were several contribution scores would be found for a single nt, how should they be aggregated? Default= function(x) mean(x)}

\item{mot}{An optional bed file containing motifs to be added.}

\item{mot.name.column}{Name of the column containing the motif name.}

\item{xlab}{Default= "nt"}

\item{ylab}{Default= "Contribution"}

\item{xlim}{Default= sequence length}

\item{ylim}{Default= range(contrib)}
}
\value{
contrib plot
}
\description{
Plot contribution scores matrix
}
