#Visualisation functions ----------------------------


#' This file contains functions to visualize different types of results
#' To do so, a series of helper functions have been created to use colorspace
#' in an elegant and sensible manner
#' - skyplot to generate manhattan plots based on a map
#' - comp.skyplot to compare multiple manhattan plots based on a map
#' - QQ.plot for plotting quantile-quantile plots
#' - comp.QQ.plot for comparing multiple p-value distributions
#' - pcoa.plot to generate pcoa distribution of distance matrices
#'

# Color -------------------------------------------

#' Colour lightening
#'
#' @param col One or more colours.
#' @param lighten Lightening factor as a value between 0 and 1 (default 0.55).
#' If negatives are used, it darkens the colour.
#'
#' @return Lightened (or darkened) colour(s)
#' @export
#'
#' @examples
lighten <- function(
  col,
  factor=0.55
){

  light<-col2rgb(col)

  if(factor>0){ result<-light+(255-light)*factor
  }else{result<-light+light*factor}

  result <- rgb(t(result), maxColorValue=255)
  return(result)

}

#' Color selection based on HCL model
#'
#' @description
#' Choice of color palettes can have a relevant impact not only
#' in the outlook of data visualizations, but also on their interpretation.
#' There exist many models of colour description, such as rgb, where colours
#' are expressed as combinations of red, green and blue; or cymk, where colours
#' are expressed as combinations of cyan, yellow, magenta and black. Alternatively
#' the hcl model proposes the use of Hue, a value between 0 and 360 in a colour
#' wheel; Chroma, the intensity of the pigment (similar to saturation); and
#' Lightness, the amount of white/black in a colour. Such organization of colour
#' aligns very well with human colour perception and thus is more intuitive to
#' use than rgb or cymk models, where the colour outcome is not so intuitive.
#'
#' For this reason, the package \code{colorspace} has been developed, which helps
#' in generating pleasant and adequate palettes for data visualization. For more
#' information on how to choose colour palettes be sure to visit
#' \link[http://colorspace.r-forge.r-project.org/articles/hcl_palettes.html]{their blog}. Or
#' if you whish to learn about colour spaces and the HCL colour space visit
#' \link[http://hclwizard.org/why-hcl/]{"why HCL"}.
#'
#' @param n Integer. Number of colours to return.
#' @param coltype Either "sequential", "qualitative" (default), "divergent" or "rainbow".
#' For a few categories, such as different treatments, choose "qualitative" or "rainbow".
#' For ordered categories, such as increasingly high levels of a compound, use
#' "sequential". For a gradient between two opposites, chose "divergent".
#' @param h One or two values between 0 and 360 (degrees in the colour wheel) to represent hue.
#' Default is c(120,240).#' If "qualitative" or "divergent" is used, the two hues will correspond
#' to each end of the colour
#' gradient. Otherwise only the first hue will be used. For reference, 0 is red, 140 is green, 240
#' is blue 300 is purple and 360 is back to red.
#' @param c Value between 0 and 100 (default 100). Represents chroma, or colour intensity and is similar to saturation.
#' 0 means grey, 100 means completely saturated.
#' @param l Value between 0 and 100 (default 60). Represents brightness, or the amount of white
#'
#' @return A colour palette of n colours.
#' @export
#'
#' @examples
select.col <- function(
  n,
  coltype = "qualitative",
  h = c(120, 240),
  c = 100,
  l = 60
){
  if(coltype=="divergent"){
    #creates a divergent palette from a central neutral color
    cols<-colorspace::diverge_hcl(n+1,h=h,c=c,power=0.7)
    half<-(n+1)%/%2+(n+1)%%2
    cols<-c(rev(cols[1:(half-1)]),cols[-1:-(half-1)])
    cols<-cols[-c(half)]#to take out the central neutral color

  }else if(coltype=="sequential"){

    cols<-colorspace::sequential_hcl(n,h=h,c=c)

  }else if(coltype=="rainbow"){
    cols<-colorspace::rainbow_hcl(n,c=c,l=l)

  }else if(coltype=="qualitative"){
    cols<-colorspace::qualitative_hcl(n,h=h,c=c)
  }
  return(cols)
}

# QQ plots -------------------

#' QQ plot calculation
#'
#' @describeIn plot.QQ Used to calculate the expected p-values and the
#' sorted observed p-values. Automatically applies -log10 on the p-values.
#'
#' @param pvals a vector of p-values.
#'
#' @return observed and expected p-values
#'
#' @examples
QQcalc<-function(pvals){
  o <- -log10(sort(pvals,decreasing=F))
  e <- -log10(1:length(o)/length(o))
  e <- e/max(e)
  return(data.frame(exp=e,obs=o))
}

#' P-value Quantile-Quantile plot
#'
#' @description Quantile-Quantile plots are useful for determining
#' whether the p-value distribution of a QTL analysis follows the expected
#' distribution. Moreover, they allow us to compare models based on their
#' p-value distribution. When the p-values follow the distribution,
#' they will be around the trend line (in red). Observed ignificant p-values
#' should have an expected \eqn{-log10(pval) > 1}, and will deviate from the
#' trend line. If the whole p-value distribution (including the region below
#' expected \eqn{-log10(pval)}) is over or under the red line,
#' we can say that there is inflation or deflation of p-values.
#'
#' For more information visit \link[https://en.wikipedia.org/wiki/Q%E2%80%93Q_plot]{Q-Q plot}
#'
#' @param pvals Vector, matrix or list of p-values.
#' Colnames/list names will be used for the legend.
#' @param main character string for title
#' @param coltype either "sequential", "rainbow", "qualitative" or "divergent".
#' Used for determining the color palette to use, for more information check
#' \code\link[select.col]{select.col}
#' @param legtitle
#' @param legend
#' @param lim
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
comp.QQ <- function(
  pvals,
  main=NULL,
  coltype="qualitative",
  legtitle="Models",
  legnames=NULL,
  ylim = NULL,
  h=240,
  ...
){
  if(is.numeric(pvals)) pvals <- list(pvals)
  if(is.matrix(pvals)) pvals <- mat2list(pvals)

  qqval <- lapply(pvals,function(p) QQcalc(p))


  if(is.null(ylim)) ylim <- c(0,max(unlist(sapply(qqval,'[', ,2))))

  plot(0,type="l",col="red",main=main,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,1),ylim=ylim)

  m <- length(pvals)
  cols <- select.col(m,coltype,h=h)

  for(i in 1:length(qqval)){
    qvals <- qqval[[i]]
    points(qvals$exp,qvals$obs,col=cols[i],pch=19,cex=0.7)
  }

  abline(a=0,b=1,col="red",lty=2)

  if(is.null(legnames)) legnames <- names(pvals)
  if(is.null(legnames)) legnames <- paste("pval",1:length(pvals))
  legend("topleft",legend = legnames,pch=19,col=cols,bty="n")
}

# Manhattan plots -------------------------------

#' Manhattan Skyline Plot
#'
#' @description Takes a vector of values and a dataframe of a genetic map containing the columns
#' "chromosome" and "position", and produces a plot where values are mapped onto the genome.
#'
#' @param pval Numerical vector. Usually -log10(p-values), but other values are accepted.
#' @param map Dataframe containing columns "chromosome" and "position"
#' @param col Base colour to use for plotting. Odd chromosomes will be plotted with this colour,
#' even chromosomes will be plotted with a lighter version of the same colour.
#' @param threshold A threshold value to draw the threshold line.
#' @param chrom A vector of chromosome names to be included in the plot
#' @param ... Other parameters to be passed to plot()
#'
#' @return
#' @export
#'
#' @examples
skyplot<-function(
  pval,
  map,
  col="navyblue",
  threshold=NULL,
  ylab=NULL,
  ylim=NULL,
  chrom = NULL,
  ...
){
  #In case the markers are not in order
  map <- map[with(map,order(map$chromosome,map$position)),]

  #We filter per chromosome
  if(is.null(chrom)) chrom <- as.character(unique(map$chromosome))
  chrom_filter <- sapply(map$chromosome,function(x) any(x == chrom))
  if(!any(chrom_filter)){
    stop("The specified chrom: ",paste(chrom,collapse=", "),
         "; is not found in the chromosomes of map: ",
         paste(as.character(unique(map$chromosome)), collapse=", "))}
  map <- map[chrom_filter, ]
  pval <- pval[chrom_filter]

  tot_length <- sum(sapply(split(map$position,map$chromosome),max))
  new_map <- map_axis(map,space = 0.05*tot_length)

  #Calculate a lighter colour of the "col"
  lightcol <- lighten(col)
  chrom_col <- sapply(map$chromosome,function(m) which(m == unique(map$chromosome)))
  col <- c(lightcol,col)[chrom_col%%2+1]


  if(is.null(ylim)) ylim <- c(min(pval,na.rm = T),
                              ceiling(max(pval,na.rm = T)/10)*10)

  if(is.null(ylab)) ylab <- expression(-log[10](italic(p)))

  plot(new_map$axis,pval,xlab="",
       ylim=ylim,axes=F,ylab=ylab,
       pch=19,col=col)

  #Y axis
  at <- axisTicks(round(ylim),log = F); axis(2,at)

  #big X axis
  #This line allows us to calculate the geometric position at which
  #each chromosome starts and ends
  ch_pos <- sapply(split(new_map$axis,map$chromosome),range)[,as.character(chrom)]
  #For each chromosome we draw a "big axis" specifying which chromosome it is
  if(length(chrom) > 2) small <- T
  else small <- F
  draw_chrom_axis(ch_pos[1,],ch_pos[2,],chrom,small)

  mtext("Chromosome",1,line = 3)

  if(!is.null(threshold)) abline(h=threshold,col="red",lwd=1.3,lty=2)
}


#' Axis calculator
#'
#' Using a map dataframe, it calculates the
#' position of each marker onto a plot, including
#' spaces between chromosomes, if desired.
#'
#' @param map map dataframe (with position, chromosome and marker)
#' @param space space to leave between chromosomes (in cM)
#' @param maxes the maximum length of each chromosome
#'
#' @return
#'
#' @examples
map_axis <- function(map,space = 0,maxes = NULL){
  sp_map <- split(map,map$chromosome)[as.character(unique(map$chromosome))]
  if(is.null(maxes)){
    maxes <- sapply(sp_map,function(x) max(x$position))
  }else{
    maxes <- maxes[as.character(unique(map$chromosome))]
  }
  maxes <- c(0,maxes)

  for(i in 1:length(sp_map)){
    sp_map[[i]]$axis <- sp_map[[i]]$position + sum(maxes[1:i]) + space*(i-1)
  }
  return(do.call(rbind,sp_map))
}

#' Comparative Skyline Manhattan plot
#'
#' @describeIn skyplot
#' @description Generates overlapped manhattan plots and adds a legend.
#' @inheritParams skyplot
#'
#' @param pvals list of pvalues to be plotted together.
#' @param map list of map dataframes. If there is only a data.frame, it will be assumed
#' that all p-value sets have the same underlying genetic map.
#' @param coltype Argument to be passed to \code{|link{select.col}}.Either "sequential",
#' "qualitative", "divergent" or "rainbow".For a few categories, such as different treatments,
#' choose "qualitative" or "rainbow". For ordered categories, such as increasingly high levels
#' of a compound, use "sequential". For a gradient between two opposites, chose "divergent".
#' @param ...
#'
#' @return A comparative manhattan plot
#' @export
#'
#' @examples
comp.skyplot<-function(
  pvals,
  map,
  coltype="qualitative",
  h=240,
  threshold=NULL,
  chrom = NULL,
  ...
){
  n <- length(pvals)
  if(is.data.frame(map)){ map <- lapply(1:n,function(x) map)
  }else{
    if(length(map) != n) stop("The number of maps provided and pvalues does not coincide")
  }

  #We filter per chromosome
  if(is.null(chrom)){
    chrom <- unique(unlist(lapply(new_maps,function(nmp){
      as.character(unique(nmp$chromosome))})))
    chrom <- sort(chrom)
  }
  for(i in 1:n){
    mp <- map[[i]]
    chrom_filter <- sapply(mp$chromosome,function(x) any(x == chrom))

    if(!any(chrom_filter)){
      stop("The specified chrom: ",paste(chrom,collapse=", "),
           "; is not found in the chromosomes of map: ",
           paste(as.character(unique(mp$chromosome)), collapse=", "))}

    map[[i]] <- mp[chrom_filter, ]
    pvals[[i]] <- pvals[[i]][chrom_filter]
  }


  #Then, we create the colour palette we are going to use
  cols <- select.col(n,coltype,h=h)

  #First we calculate the maximum positions of each map
  maxes <- sapply(1:n,function(i){
    sapply(split(map[[i]]$position,map[[i]]$chromosome),max)
  })
  maxes <- apply(maxes,1,max)
  #This allows us to calculate unified mapping axes
  new_maps <- lapply(map,map_axis,space = sum(maxes)*0.05,maxes = maxes)

  #then we calculate the maximum pvalue
  if(is.null(ylim)){
    ylim <- sapply(pvals,function(pv){
      c(min(pv,na.rm = T),
        ceiling(max(pv,na.rm = T)/10)*10)
    })
    ylim <- c(min(ylim[1,]),max(ylim[2,]))
  }

  if(is.null(ylab)) ylab <- expression(-log[10](italic(p)))

  xlim <- c(0,max(sapply(new_maps,function(x) max(x$axis,na.rm=T))))

  plot(0,type="n",ylim=ylim,ylab=ylab,xlim=xlim,axes = F,xlab="Chromosomes")

  for(i in 1:n){
    pv <- pvals[[i]]
    mp <- new_maps[[i]]
    col <- cols[i]

    lightcol <- lighten(col)
    chrom_col <- sapply(mp$chromosome,function(m) which(m == unique(mp$chromosome)))
    col <- c(lightcol,col)[chrom_col%%2+1]

    points(mp$axis,pv,pch=19,col = col)
  }

  #Plotting the legend
  leg_names <- names(pvals)
  if(is.null(leg_names)) leg_names <- paste("pvals",1:n)
  legend("topright",legend = leg_names,pch=19,col=cols,bty="n")

  #Y axis
  at <- axisTicks(round(ylim),log = F); axis(2,at)

  #We calculate the start and end of each chromosome for all maps
  ch_pos <- lapply(new_maps,function(nmp){
    sapply(split(nmp$axis,nmp$chromosome),range)[,as.character(chrom)]
  })

  #We take the max/min of each end/start for each chromosome
  ch_pos <- matrix(sapply(1:length(ch_pos[[1]]),function(i){
    val <- sapply(ch_pos,'[',i)
    if(i%%2 == 0) return(min(val))
    else return(max(val))
  }),nrow=2)

  if(length(chrom) < 2) small <- T
  else small <- F
  draw_chrom_axis(ch_pos[1,],ch_pos[2,],chrom,small)

  if(!is.null(threshold)) abline(h=threshold,col="red",lwd=1.3,lty=2)
}

#' Chromosome axis drawer
#'
#' @param ch_start vector of chromosome starts
#' @param ch_end vector of chromosome ends
#' @param chrom names of chromosomes
#' @param small logical, should small axis be also drawn?
#'
#' @return
#'
#' @examples
draw_chrom_axis <- function(ch_start,ch_end,chrom,small=F){
  if(!identical(length(ch_start),length(ch_end),length(chrom))){
    stop("Not all arguments provided of equal length")
  }

  for(i in 1:length(chrom)){
    ch_s <- ch_start[i]
    ch_e <- ch_end[i]
    axis(1,c(ch_s,ch_e),labels=c("",""),tck=-0.03)

    #Small axis only if there's one or two chromosomes
    if(small){
      at <- round(seq(ch_s,ch_e,length.out = 50))
      axis(1,at,cex.axis=0.7,padj = -1)
    }

    at <- (ch_e - ch_s)/2 + ch_s
    mtext(chrom[i],1,at=at,line=1.8,cex=1.2)
  }

}

# PCoA Plots ------------------

#' PCoA plotter
#'
#' @description It produces a principal component plot for the K
#' distance matrix based on the prcomp() function.
#'
#' @param K distance matrix
#' @param col vector for each individual in the matrix, indicating
#' a grouping. For instance, c("parent1","pop1",...). Each element will be
#' used to produce a colour.
#' @param coltype Argument to be passed to \code{\link{select.col}}.Either "sequential",
#' "qualitative", "divergent" or "rainbow".For a few categories, such as different treatments,
#' choose "qualitative" or "rainbow". For ordered categories, such as increasingly high levels
#' of a compound, use "sequential". For a gradient between two opposites, chose "divergent".
#' @param h Argument to be passed to \code{\link{select.col}}. It indicates the hue of
#' the colour palette to use.
#' @param comp Vector of two integers, indicating which principal comonents to plot.
#' 1 and 2 by default.
#' @param add Logical, whether to add the points to the current plot (using points()),
#' or make a new plot.
#' @param legend Logical, whether legend should be plotted.
#' @param legpos Text string indicating the position of the legend. "topright",
#' "topleft", "centre"... See \code{legend} for more information.
#' @param ... Additional parameters to be passed to the \code{plot()} or \code{points()}
#' method.
#'
#' @return
#' @export
#'
#' @examples
pcoa.plot <- function(
  K,
  comp=c(1,2),
  legend=T,
  col=NULL,
  coltype="sequential",
  h=140,
  add=F,
  legpos="topright",
  ...
){
  pc <- prcomp(K)

  if(is.null(col)){
    cols <- "black"
  }else{
    if(length(col) != dim(K)[1]){
      stop("col length and number of individuals do not match")}
    n <- length(unique(col))
    cols <- select.col(n,coltype=coltype,h=h,c=100)
    id <- sapply(unique(col),function(i) col==i)
    cols <- cols[id%*%1:length(cols)]
  }

  var.pc <- paste0("PC",1:length(pc$sdev),"  ",
                 round(pc$sdev^2/sum(pc$sdev^2)*100,2),"% variance")
  if(add){
    points(pc$rotation[,comp],col=cols,...)
  }else{
    plot(pc$rotation[,comp],col=cols,...,
         xlab=var.pc[comp[1]],ylab=var.pc[comp[2]])
  }

  if(legend){
    if(length(col)==dim(K)[1]) legend(legpos,legend=unique(col),
                                      col=unique(cols),bty="n",pch=19)
  }
}

# LD plots ---------------------------
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

# Miscellaneous ---------------
mat2list <- function(matrix,rowise=F){
  if(rowise){
    #if rowise we transpose the matrix so that rows are columns
    matrix<-t(matrix)
  }

  n <- ncol(matrix)
  res <- lapply(1:n,function(i) matrix[,i])

  return(res)
}
