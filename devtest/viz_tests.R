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

pheno_haplo <- function(
  phe,
  gen,
  ploidy,
  draw.points = T,
  hap.select = NULL,
  coltype = "qualitative",
  h = c(120,240),
  ...
  ){

  #Here we obtain the matrix of dosages per haplotype
  data <- dosage.X(gen,haplotype = T,ploidy=4)
  data <- data[,as.character(sort(as.numeric(colnames(data)))),drop=F]

  #Filtering
  if(is.null(hap.select)) hap.select <- colnames(data)
  if(!all(as.character(hap.select) %in% colnames(data))){
    not_in <- !as.character(hap.select) %in% colnames(data)
    stop(paste(hap.select[not_in],collapse=" "),
         " not found in provided haplotypes")
  }
  data <- data[,as.character(hap.select),drop=F]

  #Create some basic parameters
  n_hap <- ncol(data)
  #space should be relative to the number of elements
  nbox_hap <- apply(data,2,function(x) length(unique(x)))
  #Select some colour
  col <- rep(select.col(n_hap,coltype = coltype,h=h),nbox_hap)
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
          ylim=range(pretty(range(phe))),...)

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
