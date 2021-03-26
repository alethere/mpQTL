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
             hap = anc,
             result = res,
             result2 = res2)

saveRDS(data, "devtest/test_data.RDS")
data <- readRDS("devtest/test_data.RDS")

source("R/viz_fun.R")
pval <- list(-log10(data$result$pheno1$pval),-log10(data$result2$pheno1$pval))
map <- data$map
map <- map[sample(1:nrow(map),nrow(map)),]
skyplot(-log10(pval[[1]]),map,col="blue",chrom = c(0,4))
comp.skyplot(pval,map,chrom=1,threshold = 3.8)
map$chromosome <- letters[1:5][map$chromosome+1]

# Boxploots
phe <- data$pheno
dos <- data$dosage
pheno_box(phe,unlist(dos[4,]))


length(anc[1,])
gen <- anc[QTLpos[1],]

pheno_haplo(phe,anc[QTLpos[4],],4,
            main="NO",xlab="popo",ylab="more popo")








#' Calculates design matrix of X
#'
#' @param genotypes vector of dosages per individual, or matrix of genotypes per
#' chromosome at one single marker
#' @param ploidy integer indicating ploidy. Defaults to 4.
#'
#' @return if genotypes is vector, a matrix will be returned with first column having an
#' intercept (all 1's), second column having the dosages. If genotype is a matrix,
#' a design matrix is return with ncol=unique()
#' @export
#'
#' @examples
dosage.X <- function(genotypes,
                     haplotype=F,
                     ploidy=NULL,
                     normalize = F) {

  if(!haplotype){
    alcount <- matrix(genotypes,ncol=1)
    if(normalize) alcount <- (alcount-mean(alcount))/sd(alcount)
  }else{
    #we obtain the different alleles present
    unals <- unique(genotypes)
    #we obtain a design matrix indicating the allele of each chromosome
    #atn: changed this function to add a NA column
    match <- sapply(unals, function(x) {
      if (!is.na(x)) {
        m <- x == genotypes
        m[is.na(m)] <- F
      } else{
        m <- is.na(genotypes)
      }
      return(m)
    })


    #we count the number of each allele for each individual
    #For haplotypes we need to combine ploidy columns into one count
    n<-length(genotypes)/ploidy
    alcount <- t(sapply(1:n, function(x){
      colSums(match[1:ploidy + (x - 1) * ploidy, ,drop=F])
    }))

    #this line will not give the correct answer if we have a single individual
    if(nrow(alcount)==1) alcount <- t(alcount)

    if(normalize & ncol(alcount)!= 1) alcount <- apply(alcount,2,function(a) (a-mean(a))/sd(a))
    else if(normalize) alcount <- apply(alcount,2,function(a) (a-mean(a)))

    # #This part is also not nice
    # inds <- sapply(1:n,function(i){
    #   j <- 1:ploidy + (i - 1) * ploidy
    #   s <- names(genotypes)[j]
    # })

    inds <- unique(substr(names(genotypes), 1, nchar(names(genotypes)) - 2))
    rownames(alcount) <- inds
    colnames(alcount) <- unals
  }

  return(alcount)
}

#Circle of colour ---------------
n <- 1000
colour_wheel(120,xlab = "",ylab = "",r=1)




hue_wheel()




polygon(x = c(0,0.01,-0.01),
        y = c(-0.5,0.5,0.5))


for(i in c(seq(3,29),
           seq(30,by =5,length.out = 10),
           seq(60,by = 10, length.out = 10),
           seq(200,by = 50, length.out = 10)
           )){
  i <- round(i)
  jpeg(paste0("devtest/Plots/hue_wheel/",sprintf("%03d",i),".jpeg"),height = 1200,width = 1200)
  par(mar=c(3,3,3,3))
  colour_wheel(i,xlab="",ylab="")
  dev.off()
}

dir.create("devtest/Plots/lumi_wheel")

for(i in c(0:29,seq(30,170,3))){
  png(paste0("devtest/Plots/lumi_wheel/",sprintf("0%03d",i),".png"),height = 1000,width = 1000)
  hue_wheel(l = i)
  text(0,0,labels = i,cex=3,xpd=T)
  dev.off()
}

# Colour choice ----------
col <- select.col(5, alpha = 0.1)

col2rgb(col,alpha = T)

factor
light <- col2rgb(col,alpha = T)
alpha <- light["alpha",]
if(factor>0){ result<-round(light+(255-light)*factor)
}else{result<-round(light+light*factor)}

result <- rgb(t(result),alpha = alpha,maxColorValue=255)

return(result)

plot(1:5,col=result,pch=19,cex=5)

col2rgb(result,alpha=T)


### Continuous legend
plot(0,type = "n",xlim=c(0,1),ylim=c(0,1))

col_mat <- select.col(100)

xmax <- 1
xleft <- 0
ybottom <- 0
ytop <- 1

xright <- (xmax - xleft)*2 / 3


gradient_legend <- function(cols,xleft,ybottom,ytop,xright,lab_range = NULL,cex = 1){
  xmax <- xright
  xright <- xleft + (xmax - xleft)*2 / 3

  if(is.null(lab_range)) lab_range <- c(ybottom,ytop)

  rasterImage(matrix(cols,ncol=1),
              xleft = xleft,
              ybottom = ybottom,
              xright = xright,
              ytop = ytop)

  xrange <- (xmax - xleft)
  yrange <- (ytop - ybottom)
  if(xrange >= yrange) ran <- yrange
  else ran <- xrange

  inner_axis(xright + (xmax-xleft)*0.05,
             ybottom,ytop,
             tick = ran*0.02,
             offset = ran*0.07,
             lab_range = lab_range,
             cex = 0.6)
}


gradient_legend(cols= c("red","green","black"),
                xleft = 0.0,ybottom = 0, ytop = 01, xright = 0.8)

inner_axis <- function(x,ymin,ymax,lab_range = NULL,tick = 0.05,offset = 0.05,cex = 1){
  segments(x,ymin,x,ymax)
  at <- seq(ymin,ymax,length.out = 10)
  if(is.null(lab_range)) lab_range <- range(at)

  lab <- round(seq(lab_range[1],lab_range[2],length.out = 10),2)
  if(!all(lab%%1 == 0)){
    lab <- sprintf("%.2f",lab)
  }

  ticksize <- x + tick
  segments(x,at,ticksize,at)
  textpos <- ticksize + offset
  text(textpos,at,labels = lab,cex = cex)
}

