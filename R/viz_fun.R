# Visualisation functions ----------------------------
# This file contains functions to visualize different types of results
# To do so, a series of helper functions have been created to use colorspace
# in an elegant and sensible manner
# - skyplot to generate manhattan plots based on a map
# - comp.skyplot to compare multiple manhattan plots based on a map
# - QQ.plot for plotting quantile-quantile plots
# - comp.QQ.plot for comparing multiple p-value distributions
# - pcoa.plot to generate pcoa distribution of distance matrices

# Color -------------------------------------------

#' Colour lightening
#'
#' Returns a lighter (if lighten parameter is positiove) version of the colour(s) provided
#'
#' @param col One or more colours.
#' @param lighten Lightening factor as a value between 0 and 1 (default 0.55).
#' If negatives are used, it darkens the colour.
#'
#' @return Lightened (or darkened) colour(s)
#'
#' @examples
#'
#' lighten("black")
#' lighten("#000000")
#' lighten("red",-0.55) #this returns a darker colour
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
#' Choice of color palettes can have a relevant impact not only
#' in the outlook of data visualizations, but also on their interpretation.
#' There exist many models of colour description, such as rgb, or cymk, but in
#' this package we chose the hcl model (Hue, Chroma, Luminance). Hue is expressed as
#' an angle (usually between 0 and 360) in the colour wheel; Chroma, between0 and 100,
#' reflects the intensity of the pigment (similar to saturation); and
#' Lightness, the amount of white/black in a colour. Such organization of colour
#' aligns very well with human colour perception and thus is more intuitive to
#' use than rgb or cymk models, where the colour outcome is not so intuitive.
#' For this reason, the package \code{colorspace} has been developed, which helps
#' in generating pleasant and adequate palettes for data visualization. For more
#' information on how to choose colour palettes be sure to visit the
#' \href{http://colorspace.r-forge.r-project.org/articles/hcl_palettes.html}{colorspace blog}. Or
#' if you whish to learn about colour spaces and the HCL colour space visit
#' \href{http://hclwizard.org/why-hcl/}{"why HCL"}.
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
#' @param alpha Value between 0 (transparent) and 1 (opaque), to add some degree of
#' transparency to the colour palette.
#'
#' @return A colour palette of n colours.
#' @export
#'
#' @examples
#'
#' select.col(10,coltype = "qualitative")
#' select.col(10,coltype = "sequential", h = 280)
select.col <- function( n, coltype = NULL, h = NULL, c = 100, l = NULL, alpha = 1
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
#' @keywords internal
#'
QQcalc<-function(pvals){
  o <- -log10(sort(pvals,decreasing=F))
  e <- -log10(1:length(o)/length(o))
  e <- e
  return(data.frame(exp=e,obs=o))
}

#' P-value Quantile-Quantile plot
#'
#' Quantile-Quantile plots are useful for determining
#' whether the p-value distribution of a QTL analysis follows the expected
#' distribution. Moreover, they allow us to compare models based on their
#' p-value distribution. When the p-values follow the distribution,
#' they will be around the trend line (in red). Observed ignificant p-values
#' should have an expected \eqn{-log10(pval) > 1}, and will deviate from the
#' trend line. If the whole p-value distribution (including the region below
#' expected \eqn{-log10(pval)}) is over or under the red line,
#' we can say that there is inflation or deflation of p-values.
#'
#' @param pvals Vector, matrix or list of p-values.
#' Colnames/list names will be used for the legend.
#' @param ... extra parameters to be passed to plot() (not xlim or ylim)
#' @param ylim vector of two, the limits of y axis, defaults to range(pval)
#' @param plot_legend logical, should legend be plotted?
#' @param legnames vector of names for the legend, defaults to column/list names of
#' pvals. If empty "pval 1", "pval 2" etc.
#'
#' @inheritParams select.col
#'
#' @return Draws a QQ plot
#' @export
#'
#' @examples
#'
#' pvals <- pnorm(rnorm(100),lower.tail = T)
#' QQ.plot(pvals)
QQ.plot <- function( pvals, ylim = NULL, plot_legend = T, legnames=NULL, coltype= NULL,
                     h = NULL, l = NULL, legspace = 0 ,cex=0.7,...
){
  if(is.numeric(pvals)) pvals <- list(pvals)
  if(is.matrix(pvals)) pvals <- mat2list(pvals)

  qqval <- lapply(pvals,function(p) QQcalc(p))


  if(is.null(ylim)) ylim <- c(0,max(-log10(unlist(pvals)),na.rm=T))
  xlim <- range(unlist(sapply(qqval,'[',1)),na.rm = T)
  xlim[1] <- xlim[1] - legspace*(xlim[2] - xlim[1])

  plot(0,type="l",col="red",
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       ylim=ylim,xlim=xlim,...)

  m <- length(pvals)
  cols <- select.col(m,coltype,h=h,l=l)

  for(i in 1:length(qqval)){
    qvals <- qqval[[i]]
    points(qvals$exp,qvals$obs,col=cols[i],pch=19,cex=cex)
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
#' Takes a vector of values and a dataframe of a genetic map containing the columns
#' "chromosome" and "position", and produces a plot where values are mapped onto the genome.
#'
#' @param pval Numerical vector. Usually -log10(p-values), but other values are accepted.
#' @param map Dataframe containing columns "chromosome" and "position"
#' @param col Base colour to use for plotting. Odd chromosomes will be plotted with this colour,
#' even chromosomes will be plotted with a lighter version of the same colour. Alternatively,
#' a character vector specifying two colours can be provided.
#' @param threshold A threshold value to draw the threshold line.
#' @param chrom A vector of chromosome names to be included in the plot
#' @param ... Other parameters to be passed to plot()
#' @param ylab character, label to add on the y axis, defaults to "-log10(pval)"
#' @param xlab character, label to add on the x axis, defaults to "Chromosome"
#' @param ylim numeric vector of length two, minimum and maximum y axis values
#' @param small logical, should the small axis (per chromosome) be drawn? By default
#' T only if number of chromosomes <3.
#' @param pch numeric indicating the type of point to be passed to plot()
#' @param chromspace numeric determining the space between chromosomes on the
#' x axis. It is expressed as a proportion of the total map lenght. Values in the
#' range 0 - 0.05 are recommended. Default set to 0.05.
#'
#' @inheritParams select.col
#' @return creates a skyline plot
#' @export
#'
skyplot<-function( pval, map, threshold = NULL, ylab = NULL, xlab = NULL, ylim=NULL, chrom = NULL,
<<<<<<< HEAD
                   small = NULL, col = NULL, h = NULL, l = NULL, pch = NULL, cex.big = NULL,
                   cex.small = NULL, line.big = NULL, ...
){
  #In case the markers are not in order
  # map <- map[with(map,order(map$chromosome,map$position)),]
=======
                    small = NULL, col = NULL, h = NULL, l = NULL, pch = NULL, chromspace = 0.05, ...
){
  #In case the markers are not in order
  neworder <- order(map$chromosome,map$position)
  map <- map[neworder,]
  #any marker order change must be applied to pval accordingly,
  #since map and pval are expected to be provided in the same marker order
  pval <- pval[neworder]
>>>>>>> 250adb70a1ce95c7f585a3928dfc840b7af7186e

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
  axis_res <- map_axis(map,space = chromspace*tot_length)
  new_map <- axis_res[[1]][[1]]

  #Calculate a lighter colour of the "col"
  if(is.null(col)) col <- select.col(1,h = h,l=l)
  if(length(col)==1) lightcol <- lighten(col)
  if(length(col)==2) {
    lightcol <- col[2]
    col <- col[1]
  }
  palette(c(col, lightcol)) #gt
  # chrom_col <- sapply(map$chromosome,function(m) which(m == unique(map$chromosome)))
  # col <- c(lightcol,col)[chrom_col%%2+1]


  if(is.null(ylim)) ylim <- c(min(pval,na.rm = T),
                              ceiling(max(pval,na.rm = T)/10)*10)

  if(is.null(ylab)) ylab <- expression(-log[10](italic(p)))
  if(is.null(xlab)) xlab <- "Chromosome"

  if(is.null(pch)) pch <- 16 #gt: changed to 16 to reduce plotting time (when transparent)

  plot(new_map$axis,pval,
       ylim=ylim,axes=F,
       ylab=ylab,xlab=xlab,col=as.factor(map$chromosome),pch = pch,...)

  palette("default")

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
  draw_chrom_axis(axis_res$chr_edges,small = small,
                  bigcex = cex.big,
                  smallcex = cex.small)



  if(!is.null(threshold)) abline(h=threshold,col="red",lwd=1.3,lty=2)
}








#' Axis calculator
#'
#' Using a map dataframe, it calculates the
#' position of each marker onto a plot, including
#' spaces between chromosomes, if desired.
#'
#' @param maplist map dataframe (or list of map dataframes)
#' @param space space to leave between chromosomes (in map$position units)
#'
#' @return a maplist with the original maps and an extra "axis" column, and
#' a data.frame with the axis chromosomes start and end, for axis plotting.
#'
#' @keywords internal
map_axis <- function(maplist,space = 0){
  if(class(maplist) != "list") maplist <- list(maplist)

  #All chromosomes present in all maps
  chroms <- lapply(maplist,function(m){
    as.character(unique(m$chromosome))
  })
  chroms <- unique(unlist(chroms))

  #Per chromosome max per map
  ch_max <- lapply(chroms,function(ch){
    sapply(maplist,function(m){
      max_map(m,ch)
    })
  })

  #Total chromosome max
  ch_max <- unlist(lapply(ch_max,max,na.rm=T))
  ch_max <- c(0, ch_max)
  names(ch_max) <- c(chroms,"end")

  #Now we add the space between chromosomes
  ch_add <- sapply(1:length(ch_max),function(i) sum(ch_max[1:i]) + space*(i-1) )
  names(ch_add) <- names(ch_max)

  #Chromosome edges
  ch_start <- ch_add[-length(ch_add)]
  ch_end <- ch_add[-length(ch_add)] + ch_max[-1]

  #Now we create the axis columns on each map
  new_maps <- lapply(maplist,function(m){
    m$axis <- m$position + ch_add[as.character(m$chromosome)]
    return(m)
  })

  return(list(new_maps = new_maps,
              chr_edges = data.frame(ch_start,ch_end)))
}

#' Comparative Skyline Manhattan plot
#'
#' Similar to \code{skyplot}, it creates a skyline plot, but of multiple
#' p-value distributions. It takes a list or matrix of pvalues, and a single or
#' multiple maps identifying the position of each p-value set, and creates
#' a single plot that allows to compare multiple skyline plots at once.
#'
#' @inheritParams skyplot
#' @inheritParams select.col
#'
#' @param pval list or matrix (per column) of pvalues to be plotted together.
#' @param map list of map dataframes. If there is only a data.frame, it will be assumed
#' that all p-value sets have the same underlying genetic map.
#' #' @param legspace Numeric indicating the proportion of space to be left for plotting the legend
#' to the right of the plot. By default takes value 0.1 (10% of the x-value range). If legend
#' names are very long, increase this number to tweak the amount of space left to the right.
#' @param ... additional parameters to be passed to plot() (not xlim or ylim)
#' @param pch numeric vector. Each point type provided will be used for
#' each of the pvalue sets provided.
#'
#' @return A comparative manhattan plot
#' @export
#'
comp.skyplot <- function( pval, map, threshold=NULL, chrom = NULL, ylim=NULL, ylab = NULL,
                          xlab = NULL, legnames = NULL, coltype=NULL, h = NULL, l = NULL,
                          pch = NULL, chromspace = 0.05, legspace = 0.1,alpha = 0.3,small = F,...
){

  #Pvalues must be stored in a list, there needs to be as many maps
  #as p-values
  if(is.matrix(pval)) pval <- mat2list(pval)
  n <- length(pval)
  if(is.data.frame(map)){ map <- lapply(1:n,function(x) map)
  }else{
    if(length(map) != n) stop("The number of maps provided and pvalues does not coincide")
  }

  #We filter per chromosome
  if(!is.null(chrom)){
    pval <- lapply(1:length(pval),function(i){
      p <- pval[[i]]
      p <- p[map[[i]]$chromosome %in% chrom]
      if(length(p) == 0) return(NULL)
      return(p)
    })

    map <- lapply(map,function(m){
      m <- m[m$chromosome %in% chrom, ]
      if(nrow(m) == 0 ) return(NULL)
      return(m)
    })
  }
  pval <- pval[!sapply(pval,is.null)]
  map <- map[!sapply(map,is.null)]
  n <- length(pval)

  #space calculation
  all_chrom <- lapply(map,function(m) as.character(unique(m$chromosome)))
  all_chrom <- unique(unlist(all_chrom))
  tot <- sapply(all_chrom,function(ch){
    max(sapply(map,max_map,ch),na.rm=T)
  })
  #The sum of the max of each all_chromosome is the total range of the plot
  tot <- sum(tot)
  space <- tot*chromspace

  #Maps with axis column
  axis_maps <- map_axis(map,space = space)
  new_maps <- axis_maps[[1]]
  ch_edges <- axis_maps[[2]]

  #Then, we create the colour palette we are going to use
  cols <- select.col(n,coltype,h=h,l=50,alpha=alpha)

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

  #Create the xlim
  xlim <- c(0,max(sapply(new_maps,function(x) max(x$axis,na.rm=T))))
  xlim[2] <- xlim[2] + legspace*(xlim[2] - xlim[1])

  plot(0,type="n",ylim=ylim,ylab=ylab,xlim=xlim,axes = F,xlab=xlab,...)

  if(is.null(pch)) pch <- 16 #gt: changed to 16 to reduce plotting time (when transparent)
  if(length(pch) > n) warning("More pch values than pvals have been provided, only the first",n,"pch values will be used.")
  if(length(pch) < n) pch <- rep(pch,n)

  #We calculate and plot the points
  pv <- do.call(c,pval)
  ax <- do.call(c,lapply(new_maps,'[[',"axis"))
  p_i <- lapply(1:length(pval),function(i){
    rep(i,length(pval[[i]]))
  })
  p_i <- do.call(c,p_i)

  col <- cols[p_i]
  pc <- pch[p_i]

  points(ax,pv,pch = pc, col = col)

  #Plotting the legend
  if(is.null(legnames)) legnames <- names(pval)
  if(is.null(legnames)) legnames <- paste("pval",1:n)
  legend("topright",legend = legnames,pch=pch[1:n],col=cols,bty="n")

  #Y axis
  at <- axisTicks(round(ylim),log = F); axis(2,at)

  if(length(all_chrom) < 3) small <- T
  draw_chrom_axis(ch_edges,small)

  if(!is.null(threshold)) segments(0,threshold,xlim[2]*(1-legspace),
                                   threshold,lwd=1.3,lty=2,col="red")
}

#' Chromosome maximum position
#'
#' Given a genetic map (data.frame with columns chromosome
#' and position), returns the maximum position per chromosome.
#'
#' @param map data.frame containing a genetic map (at least a chromosome and
#' a position column)
#' @param chr optional, vector with chromosome(s) name(s) to return the
#' maximums of.
#'
#' @keywords internal
max_map <- function(map,chr = NULL){
  if(is.null(chr)) chr <- as.character(unique(map$chromosome))

  sapply(chr,function(ch){
    if(!any(map$chromosome == as.character(ch))) return(c(NA))
    else  return( max(map$position[map$chromosome == as.character(ch)]) )

  })
}

#' Chromosome axis drawer
#'
#' Helper function for \code{link\{skyplot}}
#'
#' @param ch_edges data.frame with first column having the axis start of each
#' chromosome, and the second column the axis end. Row names must have the
#' names of each chromosome. This data.frame can be obtained using \code{map_axis}
#' @param small logical, should small axis be also drawn?
#'
#' @return
#' @keywords internal
draw_chrom_axis <- function(ch_edges,small = F,bigcex = 1.2,smallcex=0.7,
                            bigspace = 1.8){
  for(i in 1:nrow(ch_edges)){
    #first we draw the main axis
    #the chromosome identity
    ch_edg <- unlist(ch_edges[i,])
    axis(1,ch_edg,labels = c("",""))
    #We calculate where we're gonna put them
    at <- (ch_edg[2] - ch_edg[1])/2 + ch_edg[1]
    mtext(rownames(ch_edges)[i],1,at=at,line= bigspace,cex=bigcex)

    #Small axis only if there's one or two chromosomes
    if(small){
      at <- round(seq(ch_edg[1],ch_edg[2],length.out = 50))
      lab <- round(seq(0,ch_edg[2]-ch_edg[1],length.out = length(at)),1)
      axis(1,at,labels = lab,cex.axis=smallcex,padj = -1)
    }
  }

}

# PCoA Plots ------------------

#' PCoA plotter
#'
#' It produces a principal component plot for the K
#' distance matrix based on the prcomp() function.
#' @inherit select.col
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
#' @param pch numeric vector containing 1 or as many values as unique values in "col"
#' @param colmode character, either "discrete" or "continuous", "discrete" by default. If
#' "discrete" col will be used to assign a colour to each unique value, using col as
#' a grouping variable. If "continous", col must be numeric and will be used in a gradient
#' of colour.
#' @param ... Additional parameters to be passed to the \code{plot} or \code{points}
#' method.
#'
#'
#' @inheritParams select.col
#' @return calculates and plots a pcoa plot
#' @export
#'
pcoa.plot <- function( K, comp=c(1,2), plot_legend=T, col=NULL, coltype=NULL, h=NULL,
                       l=NULL, alpha = 1, legspace = 0.1, legname = NULL, pch = 19, colmode = "discrete" ,
                       ...
){
  pc <- prcomp(K)

  if(is.null(col)){
    cols <- "black"
    plot_legend <- F
  }else{
    if(length(col) != ncol(K)){
      stop("col length and number of individuals do not match")}

    if(colmode == "discrete"){
      n <- length(unique(col))
      cols <- select.col(n,coltype=coltype,h=h,l=l,alpha=alpha)
      id <- sapply(sort(unique(col)),function(i) col == i)
      cols <- cols[id %*% 1:length(cols)]

    }else if(colmode == "continuous"){
      if(!is.numeric(col)) stop("If colmode == continuous, col must be numeric")
      #first we create a colour map on the colour space (if col starts at 0
      #and ends at 1, we put 200 values between 0 and 1)
      col_map <- seq(min(col),max(col),length.out = 200)
      #And then we assign to each value of col, which col_map is closer to it
      col_index <- sapply(col,function(co){
        which.min(abs(co - col_map))
      })
      #We obtain the colours from the colour space
      col_mat <- select.col(200,coltype=coltype,h=h,l=l,alpha=alpha)
      #And we apply the colour index to the colour space
      cols <- col_mat[col_index]

    }else stop("Wrong colmode specification")


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

  xlim <- range(pc$x[,comp[1]])
  xlim[2] <- xlim[2] + legspace*(xlim[2] - xlim[1])

  plot(pc$x[,comp],col=cols,pch=pch,xlim = xlim,
         xlab=var.pc[comp[1]],ylab=var.pc[comp[2]])

  if(plot_legend){
    if(colmode == "discrete"){

      if(is.null(legname)) legname <- sort(unique(col))
      pch_1 <- pch[!duplicated(col)][order(unique(col))]
      legend("topright",
             legend = legname,
             col = unique(cols)[order(unique(col))],bty="n",
             pch = pch_1)

    }else if(colmode == "continuous"){
      rasterImage(matrix(col_mat,ncol=1),
                  xleft = 0.6,xright = 1,ytop = 1,ybottom = -1)
    }

  }
}

# Boxplots ----------------------------

#' Phenotype boxplot
#'
#' Plots a boxplot per dosage of SNPs / per haplotype, with overlapped
#' points.
#'
#' @describeIn pheno_box Provided with a vector of phenotypes and a vector
#' of genotypes, it plots a boxplot grouping phenotypes per dosage of each genotype.
#'
#' @param phe numeric vector of phenotypes
#' @param gen if haplotype = F, numeric vector of same length as phe. If
#' haplotype = T, vector of length(phe)*ploidy.
#' @param haplotype logical, does gen contain haplotypes?
#' @param ploidy integer, ploidy of the individual
#' @param draw.points logical, should points be drawn? Defaults to T
#' @param hap.select vector, if haplotype = T, names of the haplotypes
#' to be drawn. All by default.
#' @inheritParams select.col
#'
#' @param ... Additonal parameters to "plot" (not xlim and ylim)
#'
#' @export
#'
pheno_box <- function( phe, gen, haplotype = F, ... ){
  if(haplotype){
    pheno_haplo(phe,gen,...)
  }else{
    pheno_dosage(phe,gen,...)
  }
}

#' @describeIn pheno_box method for when dosages are passed to pheno_box
pheno_dosage <- function( phe, gen, coltype=NULL, h=NULL, l=NULL, draw.points = T, ... ){
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

#' @describeIn pheno_box method for when haplotypes are passed to pheno_box
pheno_haplo <- function( phe, gen, ploidy, draw.points = T, hap.select = NULL, coltype = NULL,
                         h = NULL, l= NULL, ...
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
#' Method for LD plotting
#'
#' Given an LD list object (obtained via \code{LD_decay} function), it creates an LD
#' decay plot using all the percentiles calculated in \code{LD_decay}.
#'
#' @param LD LD list object as generated by \code{LD_decay}
#' @param max_dist numeric, maximum distance to plot
#' @param main character, title for the plot
#'
#' @return creates a plot
#' @export
#'
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
#' @return Angular coordinate matrix based on a set of angles
#' @keywords internal
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
#' @return plots a colour wheel based on triangles
#' @keywords internal
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
#' @return plots the axis of around a circle
#' @keywords internal
#'
axis_wheel <- function(res, cen = c(0,0),r=1,...){
  angles <- seq(0,2*pi,length.out = res +1)[-(res+1)]
  loc <- angle2coord(angles,cen,r)

  text(loc*1.2,labels = paste0(angles*360/(2*pi)),xpd = T) # something wrong with the symbol for degree, when building the3 package
  ticks <- cbind(loc*1.05,loc*1.1)
  segments(ticks[,1],ticks[,2],ticks[,3],ticks[,4],...)
  circle(r = 1.05)
}

#' Circle drawing
#'
#' @param cen numeric, x and y coordinates to locate the center of the circle
#' @param r numeric, radius of the circle
#' @param res numeric, resolution of the cirlce, over 60 recommended
#' @param ... other parameters to pass to segment
#'
#' @return plots a circle
#' @keywords internal
#'
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
#' @return plots ticks around a center
#' @keywords internal
#'
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
#' @param l luminance parameter, which can be modified to obtain lighter/darker
#' colours. Values between 20 and 100 are recommended (above and below
#' not all colours exist). Defaults to 60.
#'
#' @return plots a hue wheel of a specific luminance
#' @export
#'
#' @examples
#'
#' hue_wheel() #luminance = 60 by default
#' hue_wheel(90) #luminance = 90
hue_wheel <- function(l = NULL){
  if(is.null(l)) l<-60
  main <- paste("Hue wheel l =",l)
  colour_wheel(1200,xlab="",ylab="",main=main,l=l,asp = 1)
  axis_wheel(18)
  tick_wheel(360)
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
