

# Data format conversion -------------------------------

#' Convert haplotype dosages into haplotype names (and viceversa)
#'
#' @description A genotype of an individual can be represented either by a
#' numerical vector of dosages of each haplotype or by a numerical vector of
#' haplotype names (numbers) for each homologue. In the former case, vector
#' length will correspond to the number of haplotypes, while in the latter
#' case it will be the number of homologues, that is the individual ploidy.
#' The function \code{HapdoseToHapname} converts a list of haplotype dosages
#' into a matrix of haplotypes names, while \code{HapnameToHapdose} does the
#' opposite.
#' @describeIn HapdoseToHapname This function converts a list of haplotype
#' dosages into a matrix of haplotypes names.
#'
#' @param hapdose A list where each element is a numeric matrix of haplotype
#' dosages for a certain locus and each element name is a locus name. Locus
#' names must be unique. In each matrix, haplotypes are in rows and
#' individuals in columns. Row names must contain haplotype names and
#' column names contain individual names.
#' @param ploidy numeric indicating the ploidy level (identical for all the
#' individuals).
#'
#' @return \code{HapdoseToHapname}: a matrix of multiallelic genotypes,
#' with loci (blocks) in rows and individual homologues in columns. Genotypes
#' are indicated by using haplotype names.
#' @export
#' @examples
#' ## Create random haplotype dosage data for ten tetraploid
#' ## individuals and three loci.
#' hapdose <- lapply(1:3, function(l) {
#'   sapply(1:10, function(x) {
#'     tabulate(sample(1:3, 4, replace = TRUE), nbins = 3)
#'   })
#' })
#' names(hapdose) <- paste0("locus",1:3)
#'
#' for(i in 1:3){ #add names to markers and individuals
#'   d <- dim(hapdose[[i]])
#'   rownames(hapdose[[i]]) <- paste0("hap",1:d[1])
#'   colnames(hapdose[[i]]) <- paste0("ind",1:d[2])
#' }
#'
#' colSums(hapdose$locus1) #each column sums to ploidy
#'
#' ## Convert to Hapname format
#' hapname <- HapdoseToHapname(hapdose, ploidy = 4)
#'
#' ## and viceversa
#' HapnameToHapdose(hapname)
#'
HapdoseToHapname <- function(hapdose, ploidy) {

  # to use as input directly the results from PolyHaplotyper
  if (inherits(hapdose[[1]], "list")) {
    PolyH <- TRUE
  } else {
    PolyH <- FALSE
  }

  if (PolyH) {
    hapdose <- lapply(hapdose, function(x) {
      x[[1]]
    })
  }


  # check each matrix has the same number of columns
  nind <- sapply(hapdose, ncol)
  if (!all(nind[1]==nind)) stop("hapdose matrices have different number of columns")

  # check each matrix has the same individual order
  indnames <- sapply(hapdose, function(x) {
    colnames(x)
  })
  if (!all(indnames[,1]==indnames)) stop("different individual names in hapdose matrices")

  # hapdose to hapname
  hapname <- t(sapply(hapdose, function(m) {
    if (all(dim(m) > 0)) {
      # print(rownames(m)[1])
      if (PolyH) {
        hapname <- strsplit(rownames(m), "_")
        # hapname <- as.numeric(unlist(lapply(hapname, function(n) {n[length(n)]})))
        hapname <- unlist(lapply(hapname, function(n) {n[length(n)]}))
      } else {
        hapname <- as.character(rownames(m))
      }

      c(sapply(1:ncol(m), function(i) {
        if (any(is.na(m[,i]))) {
          rep(NA, ploidy)
        } else {
          rep(hapname, m[,i])
        }
      }))
    } else {
      rep(NA, nind[1]*ploidy)
    }
  }))
  colnames(hapname) <- c(sapply(indnames[,1], function(x) paste(x, 1:ploidy, sep="_")))
  return(hapname)
}


#' Convert haplotype names into haplotype dosages
#'
#' @describeIn HapdoseToHapname This function converts a matrix of haplotypes
#' names into a list of haplotype dosages.
#'
#' @param hapname A matrix of haplotype names, with markers in rows and individual
#' homologues in columns. Column names for homologues of the same individual
#' are distinguished by adding the suffices "_1","_2", etc., depending on
#' the ploidy.
#'
#' @return \code{HapnameToHapdose}: a list of matrices. Each matrix contains
#' haplotype dosages for one marker (or block), with haplotypes in rows and
#' individuals in columns.
#' @export
#'
HapnameToHapdose <- function(hapname) {

  # delete the homologue identifier to get individual names
  indname2 <- substr(colnames(hapname),1,nchar(colnames(hapname))-2)
  indname <- unique(indname2)
  # list with column indices of individuals
  indid <- lapply(1:length(indname), function(x){
    which(indname2 %in% indname[x])
  })
  names(indid) <- indname
  # list with haplotype names per block
  hapname_list <- lapply(1:nrow(hapname), function(x){
    unique(hapname[x,])
  })
  names(hapname_list) <- rownames(hapname)
  hapname_list[1:5]

  # generate the output list
  output <- lapply(1:nrow(hapname), function(x){ #
    oneblock <- sapply(indid, function(n){
      sapply(hapname_list[[x]], function(h){
        sum(hapname[x,n] == h)
      })
    })
    rownames(oneblock) <- hapname_list[[x]]
    return(oneblock)
  })
  names(output) <- rownames(hapname)
  return(output)
}




#' Convert SNP dosages to haplotype names
#'
#' Missing values are not allowed yet
#'
#' @param snpdose A matrix of SNP dosages, with markers in row and individuals
#' in columns.
#' @param ploidy Numeric indicating the ploidy level.
#'
#' @return A matrix of haplotype names, with markers in rows and individual
#' homologues in columns. Column names for homologues of the same individual
#' are distinguished by adding the suffices "_1","_2", etc., depending on
#' the ploidy.
#' @noRd
SNPdoseToHapname <- function(snpdose, ploidy) {
  hapname <- t(sapply(1:nrow(snpdose), function(i) {
    do.call("c",lapply(snpdose[i,], function(x) {
      if (!is.na(x)) {
        c(rep("B",x),rep("A",ploidy-x))
      } else {
        rep(NA, ploidy)
      }
    }))
  }))
  colnames(hapname) <- sapply(colnames(snpdose), function(x) {
    paste(x, 1:ploidy, sep = "_")
  })
  rownames(hapname) <- rownames(snpdose)
  return(hapname)
}





# Dosage probability -------------

#' Transform dosage probabilities into B allele probabilities
#'
#' @param dosP A numeric matrix with dosage probabilities in columns. Columns
#' are expected to be ordered from P0 (prob. of dosage 0) to P(ploidy).
#'
#' @return A vector where each element corresponds to the weighted mean
#' of each row of dosP.
#' @export
#' @examples
#' ## Create an example matrix of dosage probabilities for a tetraploid individual
#' dosP <- matrix(c(1,0,0,0,0,
#'                  0,1,0,0,0,
#'                  0,0,0,0,1,
#'                  0,0,0.5,0.5,0,
#'                  0,0,0.2,0.8,0,
#'                  0,0.95,0.05,0,0), ncol = 5, byrow = TRUE)
#' colnames(dosP) <- paste0("D",0:4)
#'
#' ## Convert dosage probabilities into B allele probabilities
#' dosP2Bfreq(dosP)
dosP2Bfreq <- function(dosP) {
  ploidy <- ncol(dosP)-1
  BfreqP <- as.matrix(dosP) %*% matrix(0:ploidy/ploidy, ncol=1)
  return(BfreqP)
}



