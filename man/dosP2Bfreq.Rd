% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_fun.R
\name{dosP2Bfreq}
\alias{dosP2Bfreq}
\title{Transform dosage probabilities into B allele probabilities}
\usage{
dosP2Bfreq(dosP)
}
\arguments{
\item{dosP}{A numeric matrix with dosage probabilities in columns. Columns
are expected to be ordered from P0 (prob. of dosage 0) to P(ploidy).}
}
\value{
A vector where each element corresponds to the weighted mean
of each row of dosP.
}
\description{
Transform dosage probabilities into B allele probabilities
}
\examples{
## Create an example matrix of dosage probabilities for a tetraploid individual
dosP <- matrix(c(1,0,0,0,0,
                 0,1,0,0,0,
                 0,0,0,0,1,
                 0,0,0.5,0.5,0,
                 0,0,0.2,0.8,0,
                 0,0.95,0.05,0,0), ncol = 5, byrow = TRUE)
colnames(dosP) <- paste0("D",0:4)

## Convert dosage probabilities into B allele probabilities
dosP2Bfreq(dosP)
}
