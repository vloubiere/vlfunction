% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MotifAnalysisTools.R
\name{vl_pwm_perc_to_log2}
\alias{vl_pwm_perc_to_log2}
\title{PWM perc to log2}
\usage{
vl_pwm_perc_to_log2(perc_pwm)
}
\arguments{
\item{perc_pwm}{A percentage PWM, where all columns sum to 1}
}
\value{
Returns a log2 odd ratio pwm
}
\description{
PWM perc to log2
}
\examples{
vl_Dmel_motifs_DB_full$pwms_perc[[5]]@profileMatrix
vl_pwm_perc_to_log2(vl_Dmel_motifs_DB_full$pwms_perc[[5]]@profileMatrix)
}
