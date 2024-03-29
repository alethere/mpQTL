% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapQTL_fun.R
\name{impute.knn}
\alias{impute.knn}
\title{KNN genotype imputator}
\usage{
impute.knn(geno, ploidy, map = NULL, kneighbors = 30, K = NULL)
}
\arguments{
\item{geno}{A matrix of SNP genotypes (continuous or discrete) or
haplotypes (multiallelic markers), where rows are markers and
columns may be either (1) individuals (for biallelic markers) or (2)
individual homologues (for multiallelic markers, such as haplotypes).}

\item{ploidy}{Numeric indicating the ploidy level.}

\item{map}{Optional. A data frame of three columns (marker, chromosome,
position) with map information. To be used in \code{\link{sample.marker}}
to sample a subset of evenly distributed markers.}

\item{kneighbors}{Number of nearest neighbours to use in the imputation. Default is 20.}

\item{K}{Optionally, provide the similarity matrix (K) directly.}
}
\value{
A complete matrix of haplotypes, SNP dosages or continuous SNP genotypes.
}
\description{
K-Nearest Neighbour imputation algorithms are useful for
imputing missing data. In this function, we use the kinship matrix
calculated with \code{\link{calc.K}} to determine the most similar
individuals, and substitute the missing values of one individual with
the most common dosage/haplotype of its nearest neighbours.
For continuous SNP genotypes, the nearest neighbours mean is used instead of
the most common dosage/haplotype.
}
\examples{
\dontrun{
## Get simulated genotypes for tetraploid individuals
data("mpsnpdose") # SNP dosages
data("mphapdose") # haplotypes
mpsnpdose[sample(1:length(mpsnpdose),500)] <- NA  # add NAs randomly
mphapdose[sample(1:length(mphapdose),500)] <- NA  # add NAs randomly

## Imputation
mpsnpdose <- impute.knn(mpsnpdose, ploidy = 4, kneighbors = 50)
mphapdose <- impute.knn(mphapdose, ploidy = 4, kneighbors = 50)
}

}
