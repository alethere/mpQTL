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


