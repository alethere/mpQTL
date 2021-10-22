
#' Phenotypes
#'
#' A vector of simulated phenotypes for a tetraploid multiparental population.
#'
#' @docType data
#' @format A vector of 460 elements.
#' @keywords datasets
"mppheno"

#' Pedigree
#'
#' A table of pedigrees (offspring, parent 1 and parent 2) for a simulated
#' tetraploid multiparental population.
#'
#' @docType data
#' @format A character matrix with three columns (for offspring, parent 1
#' and parent 2) and 460 rows.
#' @keywords datasets
"mpped"


#' Genetic map
#'
#' A simulated genetic map of 1013 loci, indicating marker name, chromosome
#' and position. Markers in chromosome 0 have unknown position.
#'
#' @docType data
#' @format A data frame with three columns for locus name, chromosome
#' and position.
#' @keywords datasets
"mpmap"


#' SNP dosages
#'
#' A matrix of simulated SNP dosages for a tetraploid multiparental population
#' of 460 individuals.
#'
#' @docType data
#' @format An integer matrix of SNP dosages, with markers in rows and
#' individuals in columns.
#' @keywords datasets
"mpsnpdose"


#' Haplotype dosages
#'
#' A matrix of simulated haplotype dosages for a tetraploid multiparental population
#' of 460 individuals. Since we simulated tetraploid individuals, the genotype of
#' each individual is represented by four columns, reporting the haplotype present at
#' each homologue.
#'
#' @docType data
#' @format A matrix (integer or character) indicating the haplotype present at each
#' individual homologue. Each row represents a locus and each column an individual
#' homologue.
#' @keywords datasets
"mphapdose"
