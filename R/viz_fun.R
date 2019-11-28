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
#'
#' @examples
lighten <- function(
  col,
  factor=0.55
){

  light<-col2rgb(col,alpha = T)
  alpha <- light["alpha",]
  if(factor>0){ result<-round(light+(255-light)*factor)
  }else{result<-round(light+light*factor)}

  result <- rgb(t(result), alpha=alpha, maxColorValue=255)

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
  coltype = NULL,
  h = NULL,
  c = 100,
  l = NULL,
  alpha = 1
){

  extra <- F
  if(is.null(coltype)) coltype <- "qualitative"

  if(all(is.null(h))){
    if(coltype == "qualitative"){
      if(n < 7){
        h <- c(140,500)
        extra <- T}
      else h <- c(120,240)
    }

    if(coltype == "sequential") h <- 240

    if(coltype == "divergent") h <- c(240,20)
  }

  if(length(h) >1) if(h[1] == abs(h[2]-360)) extra <- T

  if(all(is.null(l))){
    if(coltype == "qualitative") l <- 60
    if(coltype == "sequential") l <- c(30,90)
  }

  if(coltype=="divergent"){
    #creates a divergent palette from a central neutral color
    cols<-colorspace::diverge_hcl(n+1,h=h,c=c,power=0.7,alpha = alpha)
    half<-(n+1)%/%2+(n+1)%%2
    cols<-c(rev(cols[1:(half-1)]),cols[-1:-(half-1)])
    cols<-cols[-c(half)]#to take out the central neutral color

  }else if(coltype=="sequential"){

    cols<-colorspace::sequential_hcl(n,h=h,c=c,l=l,alpha=alpha)

  }else if(coltype=="rainbow"){
    cols<-colorspace::rainbow_hcl(n,c=c,l=l,h=h,alpha=alpha)

  }else if(coltype=="qualitative"){
    if(extra) n <- n+1
    cols <- colorspace::qualitative_hcl(n,h=h,c=c,l=l,alpha = alpha)
    if(extra) cols <- cols[-n]
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
  e <- e
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
#' @param legend
#' @param lim
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
QQ.plot <- function(
  pvals,
  ylim = NULL,
  plot_legend = T,
  legnames=NULL,
  coltype= NULL,
  h = NULL,
  l = NULL,
  ...
){
  if(is.numeric(pvals)) pvals <- list(pvals)
  if(is.matrix(pvals)) pvals <- mat2list(pvals)

  qqval <- lapply(pvals,function(p) QQcalc(p))


  if(is.null(ylim)) ylim <- c(0,max(-log10(unlist(pvals))))
  xlim <- range(unlist(sapply(qqval,'[',1)))

  plot(0,type="l",col="red",
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       ylim=ylim,xlim=xlim,...)

  m <- length(pvals)
  cols <- select.col(m,coltype,h=h,l=l)

  for(i in 1:length(qqval)){
    qvals <- qqval[[i]]
    points(qvals$exp,qvals$obs,col=cols[i],pch=19,cex=0.7)
  }

  abline(a=0,b=1,col="red",lty=2)

  if(plot_legend){
    if(is.null(legnames)) legnames <- names(pvals)
    if(is.null(legnames)) legnames <- paste("pval",1:length(pvals))
    legend("topleft",legend = legnames,pch=19,col=cols,bty="n")
  }
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
  threshold = NULL,
  ylab = NULL,
  xlab = NULL,
  ylim=NULL,
  chrom = NULL,
  small = NULL,
  col = NULL,
  h = NULL,
  l = NULL,
  pch = NULL,
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
  if(is.null(col)) col <- select.col(1,h = h,l=l)
  lightcol <- lighten(col)
  chrom_col <- sapply(map$chromosome,function(m) which(m == unique(map$chromosome)))
  col <- c(lightcol,col)[chrom_col%%2+1]


  if(is.null(ylim)) ylim <- c(min(pval,na.rm = T),
                              ceiling(max(pval,na.rm = T)/10)*10)

  if(is.null(ylab)) ylab <- expression(-log[10](italic(p)))
  if(is.null(xlab)) xlab <- "Chromosome"

  if(is.null(pch)) pch <- 19
  plot(new_map$axis,pval,
       ylim=ylim,axes=F,
       ylab=ylab,xlab=xlab,
       pch=pch,col=col,...)

  #Y axis
  at <- axisTicks(round(ylim),log = F); axis(2,at)

  #big X axis
  #This line allows us to calculate the geometric position at which
  #each chromosome starts and ends
  ch_pos <- sapply(split(new_map$axis,map$chromosome),range)[,as.character(chrom),drop=F]
  if(is.vector(ch_pos)) ch_pos <- matrix(ch_pos,nrow=2)
  #For each chromosome we draw a "big axis" specifying which chromosome it is
  if(is.null(small)){
    if(length(chrom) <= 2) small <- T
    else small <- F
  }
  draw_chrom_axis(ch_pos[1,],ch_pos[2,],chrom,small)



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
#' @param pval list of pvalues to be plotted together.
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
  pval,
  map,
  threshold=NULL,
  chrom = NULL,
  ylim=NULL,
  ylab = NULL,
  xlab = NULL,
  legnames = NULL,
  coltype=NULL,
  h = NULL,
  l = NULL,
  pch = NULL,
  ...
){
  n <- length(pval)
  if(is.data.frame(map)){ map <- lapply(1:n,function(x) map)
  }else{
    if(length(map) != n) stop("The number of maps provided and pvalues does not coincide")
  }

  if(is.matrix(pval)) pval <- mat2list(pval)

  #We filter per chromosome
  if(is.null(chrom)){
    chrom <- unique(unlist(lapply(map,function(nmp){
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
    pval[[i]] <- pval[[i]][chrom_filter]
  }

  #First we calculate the maximum positions of each map
  maxes <- sapply(1:n,function(i){
    sapply(split(map[[i]]$position,map[[i]]$chromosome),max)
  })
  if(is.vector(maxes)) maxes <- matrix(maxes,nrow=2)
  maxes <- apply(maxes,1,max)
  #This allows us to calculate unified mapping axes
  new_maps <- lapply(map,map_axis,space = sum(maxes)*0.05,maxes = maxes)


  #Then, we create the colour palette we are going to use
  cols <- select.col(n,coltype,h=h,l=50,alpha=0.3)

  #then we calculate the maximum pvalue
  if(is.null(ylim)){
    ylim <- sapply(pval,function(pv){
      c(min(pv,na.rm = T),
        ceiling(max(pv,na.rm = T)/10)*10)
    })
    ylim <- c(min(ylim[1,]),max(ylim[2,]))
  }

  if(is.null(ylab)) ylab <- expression(-log[10](italic(p)))
  if(is.null(xlab)) xlab <- "Chromosomes"

  xlim <- c(0,max(sapply(new_maps,function(x) max(x$axis,na.rm=T))))

  plot(0,type="n",ylim=ylim,ylab=ylab,xlim=xlim,axes = F,xlab=xlab,...)

  if(is.null(pch)) pch <- 19
  for(i in 1:n){
    pv <- pval[[i]]
    mp <- new_maps[[i]]
    col <- cols[i]

    # lightcol <- lighten(col)
    # chrom_col <- sapply(mp$chromosome,function(m) which(m == unique(mp$chromosome)))
    # col <- c(lightcol,col)[chrom_col%%2+1]

    points(mp$axis,pv,pch=19,col = col)
  }

  #Plotting the legend
  if(is.null(legnames)) legnames <- names(pval)
  if(is.null(legnames)) legnames <- paste("pval",1:n)
  legend("topright",legend = legnames,pch=19,col=cols,bty="n")

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
      lab <- round(at - ch_s)
      axis(1,at,labels = lab,cex.axis=0.7,padj = -1)
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
#' @param legspace Numeric indicating the proportion of space to be left for plotting the legend
#' to the right of the plot. By default takes value 0.1 (10% of the x-value range). If legend
#' names are very long, increase this number to tweak the amount of space left to the right.
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
  plot_legend=T,
  col=NULL,
  coltype=NULL,
  h=NULL,
  l=NULL,
  legspace = 0.1,
  legname = NULL,
  pch=19,
  ...
){
  pc <- prcomp(K)

  if(is.null(col)){
    cols <- "black"
  }else{
    if(length(col) != ncol(K)){
      stop("col length and number of individuals do not match")}
    n <- length(unique(col))
    cols <- select.col(n,coltype=coltype,h=h,l=l)
    id <- sapply(sort(unique(col)),function(i) col == i)

    cols <- cols[id %*% 1:length(cols)]
  }

  var.pc <- paste0("PC",1:length(pc$sdev),"  ",
                 round(pc$sdev^2/sum(pc$sdev^2)*100,2),"% variance")

  if(length(pch) == 1){
    pch <- rep(pch,ncol(K))
  }else{
    pch <- rep(pch,ncol(id))
    pch <- pch[id %*% 1:ncol(id)]
  }

  if(!plot_legend) legspace <- 0

  xlim <- range(pc$rotation[,comp[1]])
  xlim[2] <- xlim[2] + legspace*(xlim[2] - xlim[1])

  plot(pc$rotation[,comp],col=cols,pch=pch,...,xlim = xlim,
         xlab=var.pc[comp[1]],ylab=var.pc[comp[2]])

  if(plot_legend){
    if(length(col)== ncol(K)){
      if(is.null(legname)) legname <- sort(unique(col))

      pch_1 <- pch[!duplicated(col)][order(unique(col))]

      legend("topright",
             legend = legname,
             col = unique(cols)[order(unique(col))],bty="n",
             pch = pch_1)

    }
  }
}

# Boxplots ----------------------------

pheno_box <- function(
  phe,
  gen,
  haplotype = F,
  ...
  ){
  if(haplotype){
    pheno_haplo(phe,gen,...)
  }else{
    pheno_dosage(phe,gen,...)
  }
}

#' Phenotype boxplot
#'
#' Creates a boxplot per dosage class
#'
#' @param phe phenotype vector
#' @param gen dosage vector
#' @param coltype
#' @param h
#' @param ... other arguments to be passed to boxplot/plot (main, xlab, ylab...).
#' Not ylim.
#'
#' @return
#' @export
#' @inheritParams select.col
#'
#' @examples
pheno_dosage <- function(
  phe,
  gen,
  coltype=NULL,
  h=NULL,
  l=NULL,
  draw.points = T,
  ...
  ){
  #Boxplot works the following way:
  #It transforms whatever data you give into a list in which each
  #element is a "box". It will draw as many elements as there are in the list,
  #taking as axis names the names of the list.
  #So, you can create an empty space in the boxplot by creating a list
  #that has all the levels you need, but with one of them empty.
  boxlist <- split(phe,gen)
  box_class <- min(gen):max(gen)
  new_boxlist <- boxlist[as.character(box_class)]
  names(new_boxlist) <- box_class

  col <- select.col(length(box_class),coltype = coltype,h = h,l=l)

  boxplot(new_boxlist,
          border = col,
          outline = F,
          ylim=range(phe),...)

  if(draw.points){
    set.seed(7)
    points(jitter(unlist(gen)+1,amount = 0.25)-min(gen),
           phe,col=col[gen-min(gen)+1],
           pch=19,cex=0.85)
  }

}

pheno_haplo <- function(
  phe,
  gen,
  ploidy,
  draw.points = T,
  hap.select = NULL,
  coltype = NULL,
  h = NULL,
  l= NULL,
  ...
){

  #Here we obtain the matrix of dosages per haplotype
  data <- dosage.X(gen,haplotype = T,ploidy=ploidy)
  data <- data[,as.character(colnames(data)),drop=F]

  #Filtering
  if(is.null(hap.select)) hap.select <- colnames(data)
  if(!all(as.character(hap.select) %in% colnames(data))){
    not_in <- !as.character(hap.select) %in% colnames(data)
    stop(paste(hap.select[not_in],collapse=" "),
         " not found in provided haplotypes:\n",
         paste(unique(gen),collapse = " "))
  }
  data <- data[,as.character(hap.select),drop=F]

  #Create some basic parameters
  n_hap <- ncol(data)
  #space should be relative to the number of elements
  nbox_hap <- apply(data,2,function(x) length(unique(x)))
  #Select some colour
  col <- rep(select.col(n_hap,coltype = coltype,h=h,l=l),nbox_hap)
  space <- sum(nbox_hap)*0.05

  #This function allows us to transform haplotype dosages
  #into x axis calculation, where i is the box group
  #we want to draw into.
  axis_calc <- function(dh,nbox_hap,space,i){
    dh - min(dh) + sum(nbox_hap[0:(i-1)]) + space*(i-1)
  }

  #dh is dosage haplotype
  #Here we calculate where should the boxes be plotted
  at <- lapply(1:ncol(data),function(i){
    dh <- data[,i]
    axis_calc(sort(unique(dh)),nbox_hap,space,i)
  })
  big_at <- sapply(at,mean) #This will be used later for the axis
  at <- unlist(at)

  #Here we create the list of boxes for the boxplot
  #We split phenotypes according to haplotype dosage
  boxlist <- apply(data,2,function(x) split(phe,x))
  boxlist <- unlist(boxlist,recursive = F)

  #We draw the boxplots
  boxplot(boxlist,at = at,outline = F,
          border = col,axes = F,
          ylim=range(phe),...)

  #We add the observation points
  if(draw.points){
    x <- sapply(1:n_hap,function(i){
      dh <- data[,i]
      x <- axis_calc(dh,nbox_hap,space,i)
      set.seed(10)
      x <- jitter(x,amount = 0.15)
      return(x)
    })
    x <- as.vector(x)
    y <- rep(phe,n_hap)
    col_points <- rep(unique(col),each=length(phe))
    points(x,y,pch=19,cex=0.5,col=col_points)

  }

  #Drawing the axis of the boxplot
  axis(2,pretty(range(phe)))
  lab <- unlist(apply(data,2,function(x) sort(unique(x))))
  axis(1,at,labels = lab,cex.axis = 0.5,padj = - 2)
  axis(1,at = big_at, lab = colnames(data),tick = F,padj=1)

}

# LD plots ---------------------------
plot.LD <- function(
  LD,
  max_dist = NULL,
  main=NULL
  ){

  if(is.null(max_dist)) max_dist <- attr(LD,"max_dist")
  if(max_dist > attr(LD,"max_dist"))
    stop("Cannot choose a larger max_dist than ",attr(LD,"max_dist"))
  if(!"LD" %in% class(LD)){
    stop("Object provided is not an LD list")
  }

  linkage <- LD$LD

  plot(0,type="n",
       xlab = "cM",ylab = "r2",main = main,
       xlim = c(0,max_dist),ylim=c(0,1))

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
         legend=paste(names(LD$background)," ld1/2=",halflife)
         ,col=cols,bty="n")

}

# Circular plotting ----------------

#' Angular coordinates
#'
#' Given a circle, it calculates where in the circumference will a line
#' at a certain angle cut through the circumference.
#'
#' @param angles angles in radians
#' @param cen numeric, x and y coordinates to locate the center of the circle
#' @param r numeric, radius of the circle
#'
#' @return
#' @export
#'
#' @examples
angle2coord <- function(angles, cen = c(0,0),r = 1){
  res <- data.frame(x = cos(angles)*r + cen[1],
                    y = sin(angles)*r + cen[2])
  return(res)
}

#' Colour wheel
#'
#' Draws a colour wheel with an indicated resolution
#'
#' @param res Number of sections of different colour to draw in the circle
#' @param cen numeric, x and y coordinates to locate the center of the circle
#' @param r numeric, radius of the circle
#' @param l numeric, values of luminance for the colours
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
colour_wheel <- function(res,cen = c(0,0), r = 1,l=NULL,...){
  col <- select.col(res+1,"qualitative",h=c(0,360),l=l)
  angles <- seq(0,2*pi,length.out = res + 1)
  loc <- angle2coord(angles,cen,r)

  lim <- c(-r*1.1+cen[1],r*1.1+cen[2])
  plot(cen,type="n",axes=F, xlim = lim,ylim = lim,...)

  for(i in 1:res){
    polygon(x = c(cen[1],loc$x[i],loc$x[i+1]),
            y = c(cen[2],loc$y[i],loc$y[i+1]),
            col = col[i],border=NA)
  }
}

#' Axis wheel
#'
#' Draws axis ticks (not circle)
#'
#' @param n resolution of the circle
#' @param cen numeric, x and y coordinates to locate the center of the circle
#' @param r numeric, radius of the circle
#' @param ... other parameters to pass to segments
#'
#' @return
#'
#' @examples
axis_wheel <- function(res, cen = c(0,0),r=1,...){
  angles <- seq(0,2*pi,length.out = res +1)[-(res+1)]
  loc <- angle2coord(angles,cen,r)

  text(loc*1.2,labels = paste0(angles*360/(2*pi)),xpd = T) # something wrong with the symbol for degree, when building the3 package
  ticks <- cbind(loc*1.05,loc*1.1)
  segments(ticks[,1],ticks[,2],ticks[,3],ticks[,4],...)
  circle(r = 1.05)
}

#' Title
#'
#' @param cen numeric, x and y coordinates to locate the center of the circle
#' @param r numeric, radius of the circle
#' @param res numeric, resolution of the cirlce, over 60 recommended
#' @param ... other parameters to pass to segment
#'
#' @return
#' @export
#'
#' @examples
circle <- function(cen = c(0,0),r=1,res = 1200,...){

  angles <- seq(0,2*pi,length.out = res + 1)
  loc <- angle2coord(angles,cen,r)

  segments(loc$x,loc$y,loc$x[c(2:res,1)],loc$y[c(2:res,1)],...)
}

#' Draws a tick wheel
#'
#'
#' @param n number of ticks
#' @param cen numeric, x and y coordinates to locate the center of the circle
#' @param r numeric, radius of the circle
#'
#' @return
#'
#' @examples
tick_wheel <- function(n,cen = c(0,0),r=1){
  angles <- seq(0,2*pi,length.out = n +1)[-1]
  loc <- angle2coord(angles,cen,r)

  segments(loc$x*1.08,loc$y*1.08,loc$x*1.05,loc$y*1.05)
}

#' Hue wheel
#'
#' Draws the hue wheel of colorspace in order to know what number corresponds to each
#' colour
#'
#' @param l luminance parameter, which can be modified to obtained lighter/darker
#' colours. Values between 20 and 100 are recommended (above and below
#' not all colours exist).
#'
#' @return
#' @export
#'
#' @examples
hue_wheel <- function(l = NULL){
  opar <- par(no.readonly = T)
  par(mar = c(3,3,3,3))
  if(is.null(l)) l<-60
  main <- paste("Hue wheel l =",l)
  colour_wheel(1200,xlab="",ylab="",main=main,l=l)
  axis_wheel(18)
  tick_wheel(360)
  par(opar)
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
