# Functions ---------------------------------------------------------------
source("R/viz_fun.R")
# Data loading ------------------------------------------------------------
dos <- data.table::fread(
  "research/PedigreeSim/NAM_crosses/3_ancestral/cross200_alleledose.dat",
  header = T)[,-1]
map <- read.table("research/PedigreeSim/Potato.map",header=T)


# LD decay study ------------
LD_decay <- function(dos,map,win_size = 0.1,max_dist = NULL){

  #LD calculation
  LD <- cor(t(dos))
  not_dup <- lower.tri(LD,diag=F)
  LD <- LD[not_dup]

  #centimorgan distance between markers
  cm_d <- sapply(map$position,function(p) abs(p-map$position))
  cm_d <- cm_d[not_dup]

  #but distances between markers at different chromosomes are meaningless
  same_c <- sapply(map$chromosome, function(p) p == map$chromosome)
  same_c <- same_c[not_dup]

  cm_d <- cm_d[same_c]
  LD <- LD[same_c]

  if(is.null(max_dist)) max_dist <- max(cm_d)

  #This is inefficient, first, markers at max_dist should
  #be selected and then LD calculated in order to minimize the number
  #of operations. Optimization would be marginal, though
  windows <- seq(0,max_dist,win_size)
  ld_estimates <- t(sapply(seq_along(windows),function(w){
    sel <- cm_d > windows[w] & cm_d < windows[w] + win_size
    return(quantile(LD[sel]^2,c(0.5,0.8,0.9,0.95)))
  }))

  res <- data.frame(ld_estimates, distance = windows)
  attr(res,"class") <- c("LD","data.frame")
  attr(res,"max_dist") <- max_dist

  return(res)
}

LD_decay <- function(dos,map,win_size = 0.1,max_dist = NULL,per_chr = F){

  #First we split dosages per chromosome and calculate correlations
  #only within chromosomes
  dos_per_chr <- split(dos,map$chromosome)
  LD <- lapply(dos_per_chr,function(g){
    ld <- cor(t(g))
    not_dup <- lower.tri(ld,diag=F)
    return(ld[not_dup]^2)
  })

  #then we calculate the distance between markers
  pos_per_chr <- split(map$position,map$chromosome)
  dis <- lapply(pos_per_chr,function(d){
    res <- sapply(d,function(p) abs(p-d))
    not_dup <- lower.tri(res,diag=F)
    return(res[not_dup])
  })

  if(!per_chr){

    #Per window we calculate the percentile estimates
    LD <- do.call(c,LD)
    dis <- do.call(c,dis)

    if(is.null(max_dist)) max_dist <- max(dis)

    windows <- seq(0,max_dist,win_size)
    ld_estimates <- t(sapply(seq_along(windows),function(w){
      sel <- dis >= windows[w] & dis < windows[w] + win_size
      return(quantile(LD[sel],c(0.5,0.8,0.9,0.95)))
    }))

    #We calculate the background correlation between chromosomes
    dos_sample <- lapply(dos_per_chr,function(d) d[sample(1:nrow(d),100),] )
    dos_sample <- do.call(rbind,dos_sample)
    back_ld <- cor(t(dos_sample))^2
    back_ld <- back_ld[lower.tri(back_ld,diag=F)]
    back_ld <- quantile(back_ld,c(0.5,0.8,0.9,0.95))

    #We add some extra features
    res <- list(LD = data.frame(ld_estimates,
                                distance = windows),
                background = back_ld)
    attr(res,"max_dist") <- max(dis,na.rm=T)
    attr(res,"class") <- c("LD","list")

  }else{
    res <- lapply(1:length(unique(map$chromosome)),function(i){
      ld <- LD[[i]]
      d <- dis[[i]]

      windows <- seq(0,max(d,na.rm=T),win_size)
      ld_estimates <- t(sapply(seq_along(windows),function(w){
        sel <- d > windows[w] & d < windows[w] + win_size
        return(quantile(ld[sel],c(0.5,0.8,0.9,0.95)))
      }))

      back_ld <- quantile(ld[d > max(d,na.rm=T)*0.9],c(0.5,0.8,0.9,0.95))

      #We add some extra features
      res <- list(LD = data.frame(ld_estimates,distance = windows),
                  background = back_ld)
      attr(res,"max_dist") <- max(d,na.rm=T)
      attr(res,"class") <- c("LD","list")
      return(res)
    })

  }

  return(res)
}

LD <- LD_decay(dos,map)

plot.LD <- function(LD,max_dist = NULL,main=NULL){
  if(is.null(max_dist)) max_dist <- attr(LD,"max_dist")
  if(max_dist > attr(LD,"max_dist"))
    stop("Cannot choose a larger max_dist than ",attr(LD,"max_dist"))
  if(!"LD" %in% class(LD)){
    stop("Object provided is not an LD list")
  }

  linkage <- LD$LD

  plot(0,type="n",
       xlab = "cM",ylab = "r2",main = main,
       xlim = c(0,max_dist),ylim=c(0,0.5))

  cols <- colorspace::qualitative_hcl(4)
  halflife <- c()
  for(i in 1:4){
    ld <- linkage[,i][linkage$distance <= max_dist]
    dis <- linkage$distance[linkage$distance <= max_dist]

    #We take out missing percentile values because spline cannot handle it
    if(any(is.na(ld))){
      ld <- na.omit(ld)
      dis <- dis[-attr(ld,"na.action")]
    }

    points(dis,ld,
           pch=19,cex=0.85,col=cols[i])

    sp <- smooth.spline(dis,ld,spar = 0.6)
    lines(sp,lwd=2,col= lighten(cols[i]))

    top <- max(ld,na.rm=T)
    halfpoint <- (top - LD$background[i])/2 + LD$background[i]
    abline(v = sp$x[which.min(abs(sp$y - halfpoint))],
           col=cols[i],
           lty=2,lwd=1.5)
    halflife[i] <- round(sp$x[which.min(abs(sp$y - halfpoint))],1)

  }
  legend("topright",pch=19,
         legend=paste(c("50th","80th","90th","95th")," ld1/2=",halflife)
         ,col=cols,bty="n")

}

plot(LD,main="NAM3 LD decay",max_dist = 50)

LD$background


