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
#' @export
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
#' Colnames will be used for the legend.
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
  if(is.vector(pvals)) pvals <- list(pvals)
  if(is.matrix(pvals)) pvals <- mat2list(pvals)
  
  opar <- par(no.readonly = T)
  
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
  
  if(is.null(legnames)) legnames <- names(pvals)
  if(is.null(legnames)) legnames <- paste("pval",1:length(pvals))
  legend("topleft",legend = legnames,pch=19,col=cols,bty="n")
}

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
  add=F,
  ...
){
  #First we calculate the positions based on cumulative cM 
  chroms<-unique(map$chromosome)
  ch.length<-sapply(chroms,function(i) max(map$position[map$chromosome==i]))
  if(length(chroms)>1){
    space<-sum(ch.length)*0.1/(length(chroms)-1)
  }else{
    space <- sum(ch.length)*0.1/(length(chroms))
  }
  spaces<-sapply(1:length(ch.length),function(i) space*(i-1))
  #absolute chromosome length (cumulative)
  abschlen<-sapply(1:length(ch.length),function(i) sum(ch.length[1:i]))
  abschlen<-c(0,abschlen)
  ch.start<-abschlen[-length(abschlen)]+spaces
  ch.end<-abschlen[-1]+spaces
  abspos<-map$position+ch.start[factor(map$chromosome)]
  
  #Calculate a lighter colour of the "col"
  lightcol<-lighten(col)
  col<-c(lightcol,col)[map$chromosome%%2+1]
  
  
  if(is.null(ylim)) ylim <- c(min(pval,na.rm = T),
                              ceiling(max(pval,na.rm = T)/10)*10)
  
  if(is.null(ylab)) ylab <- expression(-log[10](italic(p)))
  
  if(add){
    points(abspos,pval,pch=19,cex=1.5,col=col,
           ...)
  }else{
    plot(abspos,pval,xlab="",ylim=ylim,axes=F,
         ylab=ylab,pch=19,cex=1.5,col=col,
         ...)
  }
  
  
  #Y axis
  at<-axisTicks(round(ylim),log = F)
  axis(2,at)
  
  #big X axis
  for(i in 1:length(ch.length)){
    axis(1,c(ch.start[i],ch.end[i]),labels=c("",""),tck=-0.03)
    at<-(ch.end[i]-ch.start[i])/2+ch.start[i]
    mtext(chroms[i],1,at=at,line=1.8,cex=1.2)
  }
  
  #small X axis
  nint<-100/length(ch.length)
  for(i in 1:length(ch.length)){
    pos<-map$position[map$chromosome==chroms[i]]
    apos<-pos+ch.start[i]
    at<-axisTicks(range(pos),F,nint=nint)+ch.start[i]
    
    labs<-round(at-ch.start[i])
    labs[1]<-""
    labs[length(labs)]<-""
    axis(1,at,labels = labs,padj=-1,cex.axis=0.75)
  }
  
  mtext("Chromosome and cM position",1,line = 3)
  
  if(!is.null(threshold)) abline(h=threshold,col="red",lwd=1.3,lty=2)
}

#' Comparative Skyline Manhattan plot
#' 
#' @describeIn skyplot
#' @description Generates overlapped manhattan plots and adds a legend.
#' @inheritParams skyplot
#' 
#' @param pvals Matrix of values. Each column will be plotted in a different colour.
#' @param map 
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
  ...
){
  n <- ncol(pvals)
  
  #Then, we create the colour palette we are going to use
  cols<-select.col(n,coltype,h=h)
  
  #First we calculate the positions based on cumulative cM 
  chroms<-unique(map$chromosome)
  ch.length<-sapply(chroms,function(i) max(map$position[map$chromosome==i]))
  space<-sum(ch.length)*0.1/(length(chroms)-1)
  spaces<-sapply(1:length(ch.length),function(i) space*(i-1))
  #absolute chromosome length (cumulative)
  abschlen<-sapply(1:length(ch.length),function(i) sum(ch.length[1:i]))
  abschlen<-c(0,abschlen)
  ch.start<-abschlen[-length(abschlen)]+spaces
  ch.end<-abschlen[-1]+spaces
  abspos<-map$position+ch.start[map$chromosome]
  
  good<-which.max(apply(pvals,2,max))
  skyplot(pvals[,good],map,col=cols[good],...)
  for(i in 1:ncol(pvals)){
    lightcol<-lighten(cols[i])
    col<-c(lightcol,cols[i])[map$chromosome%%2+1]
    points(abspos,pvals[,i],pch=19,col=col,cex=0.7)
  }
  
  legend("topright",legend = colnames(pvals),pch=19,col=cols,bty="n")
}


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
  
  var.pc<-paste0("PC",1:length(pc$sdev),"  ",
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
