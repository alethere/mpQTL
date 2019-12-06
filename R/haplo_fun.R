# Haplotyping (not really) functions -----------------

#' This script contains functions to generate haplotypes based on the genotype
#' data of PedigreeSim. IT IS NOT A HAPLOTYPE ESTIMATOR, but rather a set of functions
#' to transform SNP phase data into haplotype data.

#First we create the indeces of the SNPs we must concatenate
.hap_indexer <- function(nmark,l,method = "adjacent"){
  if(method == "adjacent"){
    index <- t(sapply(1:floor(nmark/l),function(n) 1:l+(n-1)*l ))
  }else if(method == "sliding"){
    index <- t(sapply(1:(nmark-l+1),function(i) i:(i+l-1)))
  }
  return(index)
}


#Then we calculate the position as the average of the position of each marker
.hap_pos <- function(hap_index,map){
  position <- apply(hap_index,1,function(n) mean(map$position[n]))
  chromosome <- apply(hap_index,1,function(n) mean(map$chromosome[n]))
  marker <- paste0("hap",sprintf("%04.f",1:length(position)))
  return(data.frame(marker,position,chromosome))
}
#Now we have an updated haplotype map

#Now we must create a set of haplotypes
.happer <- function(gen,hap_index){
  t(apply(hap_index,1,function(h_i){
    Reduce(paste0,data.frame(t(gen[h_i,])))
  }))

}

#' Haplotype corrector
#'
#' Sometimes haplotypes might be built accross chromosomes. In order
#' to avoid this (which results in some haplotypes belonging to
#' chromosome 7.4), we can remove those haplotypes that bridge multiple
#' chromosomes. This will, at most, be as many haplotypes as chromosomes.
#'
#' @param hap_index
#' @param hap_map
#'
#' @return
#' @export
#'
#' @examples
.hap_correct <- function(hap_index,hap_map){
  wrong <- hap_map$chromosome%%1 != 0
  return(list(hap_index = hap_index[!wrong,],
              map = hap_map[!wrong,]))
}

# Wrapper function that creates haplotypes based on a map, length and genotype
haplotyper <- function(
  gen,
  map,
  l = NULL,
  index = NULL,
  method = "adjacent"
){

  if(!is.null(l)){
    hap_i <- .hap_indexer(nrow(gen),l,method = method)
  }else if(!is.null){
    hap_i <- index
  }else{ stop("Either l or index must be specified")}

  new_map <- .hap_pos(hap_i,map)
  corrected <- .hap_correct(hap_i,new_map)
  haps <- .happer(gen,corrected$hap_index)
  return(list(haplotypes = haps,
              map = corrected$map))
}


# Haploblocking -----------------------------------

#' Group SNPs into haploblocks based on a genetic map
#'
#' @param map a data frame with three columns named "marker", "chromosome"
#' and "position".
#' @param winsize numeric specifing the maximum block size (in cM).
#' @param sldpace numeric indicating the minimum distance between consecutive.
#' blocks
#'
#' @return a list where each element contains the names of markers grouped
#' in a block
#' @export
map2blocks <- function(map,
                       winsize,    # window size
                       sldpace=winsize)  { # sliding pace

  chrs <- unique(map$chromosome)
  hb_list <- list()
  for(n in chrs){
    tempmap <- map[map$chromosome==n,]
    mpos <- tempmap$position[1]
    i <- length(hb_list) + 1

    while(mpos <= max(tempmap$position)) {
      midx <- which(tempmap$position >= mpos &
                      tempmap$position < (mpos + winsize))
      mnames <- tempmap$marker[midx]
      hbname <- paste0("chr",n,"_",i)
      hb_list[[i]] <- mnames
      names(hb_list)[i] <- hbname

      i <- i + 1
      mpos <- mpos + sldpace
      if(mpos < max(tempmap$position)) {
        mpos <- tempmap$position[min(which(tempmap$position >= mpos))]
      }
    }
  }
  return(hb_list)
}


#' Re-define haploblocks exceeding a maximum number of markers
#'
#' @param hb_list haploblock list created by \code{map2blocks}
#' @param nmrk numeric vector indicating the allowed numbers of markers per
#' block
#' @param mrkDosage not used yet
#' @param method character indicating the method. Choose between c("split",
#' "random"):
#'
#' * `split` will subset the block into blocks with maximum x SNP markers
#' (where x is max(nmrk)), according to the marker order provided.
#' * `random` will randomly select x SNP markers (where x is max(nmrk)).
#'
#' @param nrand numeric. If \code{(method = "random")}, \code{nrand} indicates
#' the number of subsets (blocks) randomly extracted.
#'
#' @return a list where each element contains the names of markers grouped
#' in a block.
#' @export
refineBlocks <- function(hb_list,
                         nmrk=2:8,
                         mrkDosage=NULL,
                         method="split", # c("split","random")
                         nrand=3) {  # number of random samplings for each window

  maxnmrk <- max(nmrk) # maximum number of markers per block

  hb_names <- names(hb_list) # store haploblock names


  # strategy for blocks with more than maxnmrk

  ## 1) split (it will produce more blocks)
  if(method=="split") {
    hb_list <- lapply(hb_list, function(x) {
      if(length(x) > maxnmrk) {
        idx <- 1:length(x)
        binvect <- idx %/% maxnmrk
        binvect <- c(0, binvect[-length(x)]) + 1
        # binvect <- (idx %/% (maxnmrk + 1e-5)) + 1
        split(x, as.factor(binvect))
      } else {
        list(x)
      }
    })
    hb_list <- unlist(hb_list, recursive = F)
  }


  ## 2) select markers up to maxnmrk
  ### less missingness

  ### less dependancy

  ### random
  if(method=="random") {
    hb_list <- lapply(hb_list, function(x) {
      if(length(x) > maxnmrk) {
        # x[sort(sample(1:length(x), maxnmrk))]
        a <- lapply(1:nrand, function(i) {
          x[sort(sample(1:length(x), maxnmrk))]
        })
        names(a) <- as.character(1:nrand)
        return(a)
      } else {
        # x
        list(x)
      }
    })
    hb_list <- unlist(hb_list, recursive = F)
    #names(hb_list) <- hb_names
  }


  # select haploblocks based on number of markers
  cond <- sapply(hb_list, length) %in% nmrk
  hb_list <- hb_list[cond]
  print(table(sapply(hb_list, length)))

  return(hb_list)
}


# Data format conversion -------------------------------

#' Convert haplotype names into haplotype dosages
#'
#' @description A genotype of an individual can be represented either by a
#' numerical vector of dosages of each haplotype or by a numerical vector of
#' haplotype names (numbers) for each homologue. In the former case, vector
#' length will correspond to the number of haplotypes, while in the latter
#' case it will be the number of homologues, that is the individual ploidy.
#' This function converts a matrix of haplotypes names into a list of
#' haplotype dosages.
#'
#' @param x a matrix of haplotype names, with markers in rows and individual
#' homologues in columns. Column names for homologues of the same individual
#' are distinguished by adding the suffices "_1","_2", etc., depending on
#' the ploidy.
#'
#' @return a list of matrices. Each matrix contains haplotype dosages for one
#' marker (or block), with haplotypes in rows and individuals in columns.
#' @export
HapnameToHapdose <- function(m) {

  # delete the homologue identifier to get individual names
  indname2 <- substr(colnames(m),1,nchar(colnames(m))-2)
  indname <- unique(indname2)
  # list with column indices of individuals
  indid <- lapply(1:length(indname), function(x){
    which(indname2 %in% indname[x])
  })
  names(indid) <- indname
  # list with haplotype names per block
  hapname <- lapply(1:nrow(m), function(x){
    unique(m[x,])
  })
  names(hapname) <- rownames(m)
  hapname[1:5]

  # generate the output list
  output <- lapply(1:nrow(m), function(x){ #
    oneblock <- sapply(indid, function(n){
      sapply(hapname[[x]], function(h){
        sum(m[x,n] == h)
      })
    })
    rownames(oneblock) <- hapname[[x]]
    return(oneblock)
  })
  names(output) <- rownames(m)
  return(output)
}


#' Convert haplotype dosages into haplotype names
#'
#' @description A genotype of an individual can be represented either by a
#' numerical vector of dosages of each haplotype or by a numerical vector of
#' haplotype names (numbers) for each homologue. In the former case, vector
#' length will correspond to the number of haplotypes, while in the latter
#' case it will be the number of homologues, that is the individual ploidy.
#' This function converts a list of haplotype dosages into a matrix of
#' haplotypes names.
#'
#' @param hapdose Either the output of PolyHaplotyper or a list where each element
#' is a numeric matrix of haplotype dosages for a certain locus and each element name
#' is a locus name. Locus names must be unique.
#' In each matrix, haplotypes are in rows and individuals in columns.
#' Row names must contain haplotype
#' names.
#' Column names contain individual names.
#' @param ploidy numeric indicating the ploidy level (identical for all the
#' individuals).
#'
#' @return a matrix of multiallelic genotypes, with loci (blocks) in rows and
#' individual homologues in columns. Genotypes are indicated by using haplotype
#' names.
#' @export
HapdoseToHapname <- function(hapdose, ploidy) {

  # to use as input directly the results from PolyHaplotyper
  if (inherits(hapdose[[1]], "list")) {
    PolyH <- T
  } else {
    PolyH <- F
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


# missing values are not allowed yet
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



# create a matrix of true haplotypes based on phased genotypes and blocks
TrueGeno2TrueHap <- function(geno, hb_list) {

  truegeno_list <- lapply(1:length(hb_list), function(i) {
    geno[rownames(geno) %in% hb_list[[i]],]
  })
  names(truegeno_list) <- names(hb_list)

  truehap <- lapply(truegeno_list, function(x) {
    x <- as.matrix(x)
    if (ncol(x) == 1) {
      return(t(x))
    } else {
      haps <- unique(x, MARGIN=2)
      # if (!is.null(dim(x)))
      colnames(haps) <- sprintf("H_%02d", 1:ncol(haps))
      return(haps)
    }
  })
  names(truehap) <- names(truegeno_list)

  truehap_list <- lapply(1:length(truegeno_list), function(i) {
    colnames(truehap[[i]])[match(as.data.frame(truegeno_list[[i]]),
                                 as.data.frame(truehap[[i]], stringsAsFactors = F))]
  })
  names(truehap_list) <- names(truehap)

  truehap_tab <- do.call("rbind", truehap_list)
  colnames(truehap_tab) <- colnames(geno)
  return(truehap_tab)
}

# Inferred haplotypes curation ----------------------
#' Prepare haplotypes for QTL mapping
#'
#' @param haplo It can be one of the following options:
#' \itemize{
#'   \item A matrix of multi-allelic genotypes, with markers in rows and
#'   individual homologues in columns. Marker names are provided as rownames
#'   and individual homologues names as column names. Genotypes can be both
#'   numbers and characters.
#'   \item A list where each element
#'   is a numeric matrix of haplotype dosages for a certain locus and each element name
#'   is a locus name. Locus names must be unique.
#'   In each matrix, haplotypes are in rows and individuals in columns.
#'   Row names must contain haplotype
#'   names.
#'   Column names contain individual names.
#'   \item The output of PolyHaplotyper.
#' }
#'
#' @param hb_list
#' @param map
#' @param ploidy
#' @param na.rate
#' @param use.SNPs
#' @param snpdose
#'
#' @return
#' @export
HapCurate <- function(haplo, #results of PolyHaplotyper
                      hb_list, #containing also 1 SNP blocks
                      map = NULL,  #snp map
                      ploidy,
                      na.rate = 1, #from 0 to 1. 1 means no filtration
                      use.SNPs = F,
                      snpdose) {

  if(inherits(haplo, "list")) {
    # from haplotype dose to haplotype name
    haplo <- HapdoseToHapname(haplo, ploidy)
  }


  # NA screen for markers
  namrk <- apply(haplo, 1, function(x) sum(is.na(x))/length(x))
  na_hb <- namrk > na.rate

  # use SNPs (optional)
  haplosnp2 <- NULL
  haplosnp1 <- NULL
  if (use.SNPs) {
    # add back SNPs for blocks with many NAs
    snp_nahb <- unlist(hb_list[na_hb])
    if (sum(na_hb) > 0) {
      # snp_nahb <- do.call("c", lapply(results[na_hb], function(x) {
      #   x$markers
      # }))
      haplosnp2 <- SNPdoseToHapname(snpdose[rownames(snpdose) %in% snp_nahb,,drop=F],
                                    ploidy)
    }

    # add back SNPs of 1 SNP blocks
    snp_1snphb <- unlist(hb_list[sapply(hb_list, length) == 1])
    if (length(snp_1snphb) > 0) {
      haplosnp1 <- SNPdoseToHapname(snpdose[rownames(snpdose) %in% snp_1snphb,,drop=F],
                                    ploidy)
    }
  }


  # apply filters (and merge snps)
  haplo <- rbind(haplo[!na_hb,], haplosnp2, haplosnp1)

  # map
  newmap <- NULL
  if (!is.null(map)) {

    ## haplotype map
    names(hb_list) %in% rownames(haplo)

    haplomap.temp <-
      t(sapply(hb_list[names(hb_list) %in% rownames(haplo)], function(x) {
        tmean <- mean(map$position[map$marker %in% as.character(x)])
        tchr <- unique(map$chromosome[map$marker %in% as.character(x)])
        return(c(tchr, tmean))
      }))

    haplomap <-
      data.frame(marker = rownames(haplomap.temp),
                 chromosome = haplomap.temp[,1],
                 position = haplomap.temp[,2])

    ## add snp map?
    snpmap <- NULL
    if (!is.null(haplosnp1) | !is.null(haplosnp2)) {
      snpmap <- map[map$marker %in% c(rownames(haplosnp2),rownames(haplosnp1)),]
    }
    newmap <- rbind(haplomap, snpmap) #combine maps

    ## order map and haplo (same order)
    newmap <- newmap[order(newmap$chr,newmap$pos),]
    haplo <- haplo[match(newmap[,1], rownames(haplo)),]
    rownames(newmap) <- NULL
    if (any(newmap[,1] != rownames(haplo))) stop("map and genotypes have markers in different order")
  }

  return(list(genotypes = haplo,
              map = newmap))
}
