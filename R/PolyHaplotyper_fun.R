


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







# Inferred haplotypes curation ----------------------

#' Prepare haplotypes for QTL mapping
#'
#' @param haplo It can be one of the following options:
#' \itemize{
#'   \item A matrix of multi-allelic genotypes ('haplotype names'), with
#'   markers in rows and individual homologues in columns. Marker names
#'   are provided as rownames and individual homologues names as column names.
#'   Genotypes can be both numbers and characters.
#'   \item A list where each element
#'   is a numeric matrix of haplotype dosages for a certain locus and each element name
#'   is a locus name. Locus names must be unique.
#'   In each matrix, haplotypes are in rows and individuals in columns.
#'   Row names contain haplotype names.
#'   Column names contain individual names.
#'   \item The output of PolyHaplotyper.
#' }
#'
#' @param hb_list A list where each element is a vector of SNP names contained
#' in a block.
#' @param map A data frame with three columns named "marker", "chromosome"
#' and "position".
#' @param ploidy Numeric indicating the ploidy level.
#' @param na.rate Numeric, from 0 to 1, indicating the rate of missingness
#' allowed per marker (or block). Using 1 no filtration is applied.
#' @param use.SNPs Logical value. If TRUE, SNP markers of discarded blocks will be
#' included in the final haplotype table. Default is FALSE.
#' @param snpdose A numeric matrix of SNP dosages, with markers in rows and
#' individuals in columns. Row names are marker names and colmn names are
#' individual names.
#'
#' @return A list with two elements
#' \itemize{
#'   \item *$genotypes* Curated matrix of multi-allelic genotypes
#'   ('haplotype names'), with markers in rows and individual homologues in
#'   columns. Marker names are provided as rownames and individual homologues
#'   names as column names. Genotypes can be both numbers and characters.
#'   \item *$map* A map for the haplotypes, where the position of each marker
#'   is the average position of the SNP markers contained in a block.
#' }
#' @export
HapCurate <- function(haplo,
                      hb_list, #may contain also 1-SNP blocks
                      map = NULL,  #snp map
                      ploidy,
                      na.rate = 1, #from 0 to 1, where 1 means no filtration
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
        # tmean <- mean(map$position[map$marker %in% as.character(x)]) #mean
        tmean <- mean(range(map$position[map$marker %in% as.character(x)])) #mid point of range
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






inspect_haps <- function(results, hb_list, ind=NULL, fld, outname) {


  fld_hap <- fld
  if(is.null(ind)){
    ind <- colnames(results[[1]][[1]])
  }


  # + NA rate per individual -------------
  na_ind <- sapply(results, function(x){
    is.na(x[[1]][1,ind])
  })
  narate_ind <- rowSums(na_ind)/ncol(na_ind)


  # + NA rate per marker --------------------
  narate <- sapply(results, function(x){
    sum(is.na(x[[1]][1,ind]))/length(x[[1]][1,ind])
  })
  # narate[1:5]

  png(paste0(fld_hap,"/NArate_sorted_",outname,".png"), width = 480, height = 480)
  plot(sort(narate))
  dev.off()

  # plot(narate)


  # summarize block results per window
  # prepare plot data
  window_fact <- sapply(strsplit(names(results),"\\."), function(x) x[[1]])
  wfc <- factor(window_fact, levels=names(hb_list)) # factor including all the bins in hb_list
  # class(window_fact)
  # factor(window_fact, levels=unique(window_fact))
  # boxplot(narate ~ factor(window_fact, levels=unique(window_fact)),
  #         las=2, xlab = "")
  #
  # barplot(sapply(hb_list_chr13_25to50, length), las=2)


  # # plot 1
  # # make labels and margins smaller
  # par(cex=0.7, mai=c(0.1,0.1,0.2,0.1))
  # # define area for the histogram
  # par(fig=c(0.1,1,0.7,1))
  # barplot(sapply(hb_list_chr13_25to50, length), las=2)
  # # define area for the boxplot
  # par(fig=c(0.1,1,0.1,0.6), new=TRUE)
  # boxplot(narate ~ factor(window_fact, levels=unique(window_fact)),
  #         las=2, xlab = "")

  # plot 2
  # names(hb_list)[1:5]
  # make labels and margins smaller
  png(paste0(fld_hap,"/NArate_",outname,".png"), width = 700, height = 480)
  par(cex=0.7, mai=c(0.2,0.7,0.2,0.1))
  # define area for the histogram
  par(fig=c(0,1,0.7,1))
  barplot(sapply(hb_list, length), #xaxt='n',
          names.arg = table(wfc), las = 2)
  # define area for the boxplot
  par(fig=c(0,1,0.1,0.75), new=TRUE)
  boxplot(narate ~ wfc, #unique(window_fact)
          ylim=c(0,1), las=2, xlab = "", ylab = "NA rate")
  dev.off()



  # # plot 3 (incorrect, scale of upper plot is different since some windows have no marker)
  # # make labels and margins smaller
  # par(cex=0.7, mai=c(0.1,0.1,0.2,0.1))
  # # define area for the histogram
  # par(fig=c(0.1,1,0.7,1))
  # hist(map_biglip_locus$position, breaks = seq(25,50,0.5))
  # # define area for the boxplot
  # par(fig=c(0.1,1,0.1,0.7), new=TRUE)
  # boxplot(narate ~ factor(window_fact, levels=unique(window_fact)),
  #         las=2, xlab = "")



  # + number of haplotypes --------------------

  # str(results[1])
  # results$chr13_1.1$hapdos[,1:5]

  ## Haplotype frequency
  hap_freq <- lapply(results, function(x) {
    rowSums(x[[1]][,ind], na.rm = T)/(ncol(x[[1]][,ind])*4)
  })
  # hap_freq[1:3]
  # hap_freq_check <- sapply(hap_freq, sum)

  ## number of haplotypes
  hap_num <- sapply(hap_freq, length)
  # plot(sort(hap_num))

  ## number of haplotypes above a certain frequency
  hap_num_freq <- sapply(hap_freq, function(x){
    sum(x >= 0.02)
  })
  # plot(sort(hap_num_freq))

  # plot
  # make labels and margins smaller
  png(paste0(fld_hap,"/HapNum_MAF0.02_",outname,".png"), width = 700, height = 480)
  # par(cex=0.7, mai=c(0.2,0.7,0.2,0.1))
  # # define area for the histogram
  # par(fig=c(0,1,0.7,1))
  # barplot(sapply(hb_list_chr13_25to50, length), xaxt='n')
  # # define area for the boxplot
  # par(fig=c(0,1,0.1,0.74), new=TRUE)
  boxplot(hap_num_freq ~ wfc,
          las=2, xlab = "", ylab = "hap number")
  dev.off()


  ##
  return(list(NArate = narate,
              NArate_ind = narate_ind,
              hapfreq = hap_freq))
}



# function to filter NA per window
filter_haps <- function(results, insp, n=10, na.thr=0.25){

  window_fact <- sapply(strsplit(names(results),"\\."), function(x) x[[1]])
  window_fact <- factor(window_fact, levels=unique(window_fact))
  # table(window_fact)
  # length(insp$NArate) == length(window_fact)
  out <- NULL
  for(i in unique(window_fact)) {
    x <- sort(insp$NArate[window_fact == i])
    if(length(x) > n) x <- x[1:n]
    x <- x[x < na.thr]

    out <- c(out,x)
  }
  return(out)
}



# For simulated data --------

#' create a matrix of true haplotypes based on phased genotypes and blocks
#'
#' @param geno
#' @param hb_list
#'
#' @return
#' @export
#'
#' @noRd
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

