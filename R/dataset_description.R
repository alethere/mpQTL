#' Example dataset of multiparental population
#'
#' List with simulated tetraploid data of a multiparental population
#' demonstrating the capabilities of \code{mpQTL}.
#' It includes:
#' * phenotype: matrix with two phenotypes, the first without a cofactor effect and
#' the second with a cofactor effect.
#' * map: a genetic map data.frame
#' * dosage: a dosage genotype matrix
#' * founder: a founder-allele matrix, that can be regarded as a "true haplotype" matrix
#' * cofactor: a numeric vector identifying the cofactor applied to "phenotype2"
#' in the phenotype matrix.
#' * result: a result using mixed model and dosages for this dataset.
#'
#' @docType data
#'
#' @usage data(data)
#'
#' @format A list with multiple entries, describing the phenotypes and genotypes of a
#' multiparental population.
#'
#' @keywords datasets
#'
#'
