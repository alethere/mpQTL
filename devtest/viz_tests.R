#Testing and development of new visualization functions

# Development of some new testing phenotypes
source("R/pheno_fun.R")
dos <- data.table::fread("research/PedigreeSim/HWE_NAM/3_ancestral/cross020_alleledose.dat")[,-1]
map <- data.table::fread("research/PedigreeSim/Potato.map")
anc <- link_NAM(crossfile = "research/PedigreeSim/HWE_NAM/3_ancestral/cross020_founderalleles.dat",
                totallele = "research/PedigreeSim/HWE_Parents/HWE_Totallele.txt",ploidy = 4)

#2866 is the last marker of chromosome 4
new_last <- max(which(map$chromosome == 4))
set.seed(10)
chr_0 <- sample((new_last+1):nrow(map),150)

#Dosages for four chromosomes and 150 unmapped markers
dos <- dos[c(1:new_last,chr_0),]
map <- rbind(map[1:new_last,],data.frame(marker = paste0("marker",1:150),
                                         position = 1:150,
                                         chromosome = 0))
#Ancestral alleles for the same markers
anc <- anc[c(1:new_last,chr_0),]

set.seed(17)
QTLpos <- sample(1:nrow(anc),4)

QTL <- anc[QTLpos,]
eff <- effect_gen(QTL,matrix(0:399,ncol=10),anc_sd = 2)
phe <- pheno(QTL,eff,4,herit=0.7,mu=50,return.effects = F)

source("R/mapQTL_fun.R")
res <- map.QTL(phe,dos,4,map,K=T)


res2 <- map.QTL(phe,dos,4,map,K=T,
                cofactor = unlist(dos[QTLpos[1],]),
                cofactor.type = "numerical")


data <- list(map = map,
             pheno = phe,
             dosage = dos,
             result = res,
             result2 = res2)

saveRDS(data, "devtest/test_data.RDS")
data <- readRDS("devtest/test_data.RDS")

source("R/viz_fun.R")
pval <- list(data$result$pheno1$pval, data$result2$pheno1$pval)
map <- data$map
map <- map[sample(1:nrow(map),nrow(map)),]
skyplot(-log10(pval),map,col="blue")

map$chromosome <- letters[1:5][map$chromosome+1]



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
  for(i in 1:length(unique(map$chromosome))){
    ch_start <- ch_pos[1,i]
    ch_end <- ch_pos[2,i]
    axis(1,c(ch_start,ch_end),labels=c("",""),tck=-0.03)

    #Small axis only if there's one or two chromosomes
    if(length(chrom) <= 2){
      at <- round(seq(ch_start,ch_end,length.out = 50))
      axis(1,at,cex.axis=0.7,padj = -1)
    }

    at <- (ch_end - ch_start)/2 + ch_start
    mtext(chrom[i],1,at=at,line=1.8,cex=1.2)
  }

  mtext("Chromosome",1,line = 3)

  if(!is.null(threshold)) abline(h=threshold,col="red",lwd=1.3,lty=2)
}


map_axis <- function(map,space = 0){
  sp_map <- split(map,map$chromosome)[as.character(unique(map$chromosome))]
  maxes <- sapply(sp_map,function(x) max(x$position))
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
#'

pvals <- list(-log10(data$result$pheno1$pval),
              -log10(data$result2$pheno1$pval))
map <- data$map
comp.skyplot(pvals,map,chrom=c(1,0))

# Boxploots
boxplot(phe~unlist(dos[9,]),col=select.col(5,coltype = "sequential"))
gen <- unlist(dos[9,])

pheno_box <- function(phe,gen){
  col <- select.col(length(unique(gen)))
  boxplot(phe~gen,border = col,outline = F,ylim=range(phe))
  points(jitter(gen+1),phe,col=col[gen+1],pch=19,cex=0.85)
}

pheno_box(phe,unlist(dos[1001,]))
