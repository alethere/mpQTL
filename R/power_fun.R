### Functions for Power Calculations ###

# In this file there's functions stored that are used for generating
# power calculations based on the following parameters:
# 1. A distribution of p-values along a genetic map.
# 2. A set of true QTL positions, defined by position and chromosome in the genetic map.
# 3. A threshold of significance

lag_diff <- function(x)  x[-1] - x[-length(x)]

#' QTL reader
#'
#' Given a p-value vector, a genetic map and a significance threshold,
#' this function obtains sets of p-values that all belong to the same
#' QTL interval. These are defined by the parameter `space`, a linking
#' distance parameter. When two significant markers are at a distance of
#' `space` or less they are considered to be part of the same QTL interval.
#'
#' @param pv numeric vector of p-values
#' @param map data.frame containing at least chromosome and position columns.
#' It is assumed that the pvalue order in pv and the row order in map are the same.
#' @param threshold numeric threshold for p-value significance (values below the threshold)
#' are considered significant.
#' @param space numeric indicating the linking distance to be used to define a QTL
#' interval. That is, if two significant p-values are at a distance of space or less
#' they are considered to be part of the same QTL interval.
#' @param log10 logical indicating whether p-values and threshold are expressed in
#' -log10(pval) or not.
#'
#' @return A list of QTLs  (class QTL_list) where each element is a data.frame containing
#' all the significant markers in a QTL interval, their position, chromosome and
#' pvalue.
#' @export
#'
#' @examples
read_QTL <- function(pv,map,threshold = 10^-3, space = 3,log10 = F){
  #First we add the pvalue on the map data.frame, that makes it easier to work with the related
  #data of position, chromosome, and pvalue
  if(log10){pv <- 10^-pv; threshold <- 10^-threshold}
  mapv <- cbind(map,pval = pv)

  #We split the map per chromosome
  sig <- pv < threshold
  mapv <- split(mapv[sig,],mapv$chromosome[sig])

  #We split the mapv into QTL groups, per chromosome
  mapv <- lapply(mapv,function(m){
    distance <- lag_diff(m$position)

    #We calculate a vector, peak, that tells us
    #to which QTL peak of that chromosome does each marker belong
    peak <- 1; q <- 1
    for(d in distance){
      if(d >= space) q <- q+1
      peak <- c(peak,q)
    }
    peak <- paste0("QTL",peak)


    return(split(m,peak))
  })

  mapv <- unlist(mapv,recursive = F)
  n_mark <- sapply(mapv,nrow)
  mapv <- mapv[n_mark > 1 ] #All QTLs should have more than one marker
  if(length(mapv) == 0){
    mapv <- list(data.frame(marker = NA,chromosome = NA,
                       position = NA, pval = NA))
  }

  #We must rename the list
  chroms <- sapply(mapv,function(m) unique(m$chromosome))
  qtl_num <- unlist(sapply(table(chroms),function(t) 1:t))
  qtl_names <- paste0(chroms,".QTL",qtl_num)

  names(mapv) <- qtl_names
  class(mapv) <- "QTL_list"
  return(mapv)
}


#' Summary method for QTL list
#'
#' Summarizes a QTL list to a data.frame containing information on
#' number of markers supporting a QTL, positional interval, p-value peak
#' (most significant p-value in the interval), peak position and average p-value.
#'
#' @param QTL A QTL list as obtained from `read_QTL`
#'
#' @return A summary data.frame with columns `n`, number of markers
#' supporting a QTL interval; `chr`, chromosome in which the QTL interval
#' is found; `peak`, position of the most significant p-value; `min`,
#' left boundary of the QTL interval; `max`, right boundary of the QTL
#' interval; `mean_pv`, average p-value in the QTL interval; `peak_pv`,
#' p-value of the most significant marker.
#' @export
#'
#' @examples
summary.QTL_list <- function(QTL){

  #We summarize the results of each QTL peak
  res <- sapply(QTL,function(q){
    n <- nrow(q)
    peak <- q[which.min(q$pval),"position"]
    if(length(peak) == 0){
      #no peaks are present
      res <- cbind(n = 0,
                   chr = NA,
                   peak = NA,
                   min = NA,
                   max = NA,
                   mean_pv = NA,
                   peak_pv = NA)
    }else{
      #At least one peak found
      res <- cbind(n = n,
                   chr = unique(q$chromosome),
                   peak = peak,
                   min = min(q$position,na.rm = T),
                   max = max(q$position,na.rm=T),
                   mean_pv = mean(q$pval,na.rm=T),
                   peak_pv = min(q$pval,na.rm=T))
    }
    return(res)


  })

  res <- t(res)
  res <- as.data.frame(res)
  colnames(res) <- c("n","chr","peak","min","max","mean_pv","peak_pv")

  return(res)
}


#' QTL power calculator
#'
#' Given a `QTL_list` (list of QTL dataframes) and the positions
#' of true QTLs, it calculates QTL precision, QTL sensitivity, marker precision
#' and accuracy (mean distance). It does so by checking whether the
#' detected QTL intervals contain or not the true positions.
#'
#' @param QTL a `QTL_list` object as returned from `readQTL`
#' @param trueQTL a data.frame containing columns chromosome and position,
#' indicating the chromosome and position of true QTLs.
#'
#' @return
#' @export
#'
#' @examples
power_QTL <- function(QTL,trueQTL){
  #QTL is a QTL_list object (return from read_QTL)
  #trueQTL is a data.frame with at least position and chromosome columns

  q <- summary(QTL)

  res <- t(sapply(1:nrow(trueQTL),function(i){
    p <- trueQTL[i,]
    #this evaluates, are found QTLs on the same chromosome and within the interval?
    in_chr <- q$chr == p$chromosome
    in_int <- (p$position >= q$min & p$position <= q$max)
    true_q <- in_chr & in_int

    if(all(is.na(true_q))){
      true_q <- F
      peak_dist <- NA
      n <- NA
      mean_pv <- NA
    }else if(any(true_q)){
      peak_dist <- abs(q[true_q,]$peak - p$position) #distance from peak to position
      n <- q[true_q,"n"] #Number of markers supporting that
      mean_pv <- q[true_q,"mean_pv"]
    }else{
      peak_dist <- NA
      n <- NA
      mean_pv <- NA
    }
    return(c(found = any(true_q),
             peak_dist = peak_dist,
             markers = n,
             mean_pv = mean_pv))
  }))
  rownames(res) <- rownames(trueQTL)

  #Here we calculate the different statistics
  detected_QTL <- nrow(q) - sum(is.na(q$chr)) #number of detected QTLs
  nt_QTL <- nrow(trueQTL) #number of true QTLs
  tp_QTL <- sum(res[,"found"]) #number of detected true QTLs

  detected_mark <- sum(q$n) #only within QTL intervals
  tp_mark <- sum(res[,"markers"],na.rm = T) #Number of detected true markers

  #Here we calculate all the different statistics
  QTL_precision <- tp_QTL/detected_QTL
  QTL_sensitivity <- tp_QTL/nt_QTL
  mark_precision <- tp_mark/detected_mark
  mean_dist <- mean(res[,"peak_dist"],na.rm=T)

  power <- c(QTL_precision,QTL_sensitivity,mark_precision,mean_dist)
  power[is.nan(power)] <- 0
  names(power) <- c("QTL_precision","QTL_sensitivity","mark_precision","mean_dist")
  power_summary <- list(data = res,
                        power = power)
  return(power_summary)
}


