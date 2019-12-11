#' Workshop December 2019 -------------------
#' Code and explanations on mpQTL v0.2.0
#' by Alejandro Thérèse Navarro & Giorgio Tumino


# FIRST PART: Data introduction ------------------------
# Set your working directory
# if you followed the instructions, it should be here:
setwd("C:/workshop/myanalyses")


# Installing and loading the package
install.packages("../mpQTL/mpQTL_0.2.0.tar.gz", repos=NULL)
library("mpQTL")

# the object "data" contains an example dataset
data <- readRDS("../mpQTL/new_workshop_data.RDS")
names(data)


# 1 Population design ----------------------------
## pedigree file
ped <- data$pedigree
head(ped,15)


# 2 Haplotyping ----------------------------------
## parents used in each cross
parents <- as.matrix(unique(ped[!is.na(ped[,2]) & !is.na(ped[,3]),2:3]))
parents

## merge reciprocal crosses
for(i in 1:nrow(parents)) {
  parents[i,] <- parents[i, order(parents[i,])]
}
parents <- as.matrix(unique(parents[,1:2]))
parents

F1 <- list()
for (p in seq_len(nrow(parents))) {
  F1[[p]] <- ped[,1][(!is.na(ped[,2]) & !is.na(ped[,3]) & ped[,2]==parents[p, 1] & ped[,3]==parents[p, 2]) |
                       (!is.na(ped[,2]) & !is.na(ped[,3]) & ped[,2]==parents[p, 2] & ped[,3]==parents[p, 1])]
}
str(F1)


# Make haploblocks
## genetic map
map <- data$map
head(map)
table(map$chromosome) # number of markers per chromosome
tapply(map$position, map$chromosome, max) # chromosome lenght

## create blocks by a sliding window
hb_list <- map2blocks(map=map[map$chromosome != 0,], #excluding chrom 0
                      winsize=1, # window size
                      sldpace=1) # window step (sliding pace)

length(hb_list)
table(sapply(hb_list, length))
names(hb_list)[1:5]

## refine blocks with more than 8 SNPs
set.seed(3)
hb_list <- refineBlocks(hb_list,   # haploblock list
                        nmrk=2:7,  # numeric vector indicating allowed block length
                        method="random",
                        nrand = 2) # number of random selections per block
names(hb_list)[1:5]

# Haplotype inference (DO NOT RUN, we get true haplotypes from our simulation)
# ## infer haplotypes by PolyHaplotyper (DO NOT RUN)
# results <- inferHaplotypes(mrkDosage=snp, ploidy=4,
#                            haploblock=hb_list,
#                            parents=parents, F1=F1,
#                            maxmrk=7)

## get true haplotypes from the simulation
PolyHap_res <- data$PolyHap_res
str(PolyHap_res[[1]])
## first block, 21 haplotypes, first 5 individuals
PolyHap_res[[1]]$hapdos[,1:5]
## the sum per column (individual) must be equal to ploidy
colSums(PolyHap_res[[1]]$hapdos[,1:5])


# 3 Haplotype data cuaration and format conversion -------------------------------------
hapres <- HapCurate(PolyHap_res,  # results of PolyHaplotyper or 'haplotype dosages' or 'haplotype names'
                    hb_list = hb_list, #may contain also 1-SNP blocks
                    map = map,
                    ploidy = 4,
                    na.rate = 1,  # allowed missingness per marker
                    use.SNPs = F, # whether you want to include SNPs of discarded blocks
                    snpdose = snp)

## haplotypes ('haplotype names' format)
hap <- hapres$genotypes
hap[1:5,1:8]
## haplotype map
hapmap <- hapres$map
head(hapmap)

## save hap and hapmap
saveRDS(hap, "hap.rds")
saveRDS(hapmap, "hapmap.rds")


# 4 Types of genetic data accepted by `mpQTL` --------------------------
## SNP dosages
snp <- as.matrix(data$snp)
snp[1:5,1:4]

## Haplotype names
hap[1:5,1:8]

## from (a list of) haplotype dosages to haplotype names
eg1 <- HapdoseToHapname(PolyHap_res[1:5], ploidy = 4)
eg1[,1:8]

## from SNP dosages to haplotype names
eg2 <- SNPdoseToHapname(snp[1:5,1:4], ploidy = 4)
eg2

## from haplotype names to haplotype dosages
eg3 <- HapnameToHapdose(eg1)
str(eg3)
eg3[[1]][,1:5]




# SECOND PART: mpQTL functionalities  -----------------------

# #Installing and loading the package
# install.packages("../mpQTL_0.2.0.tar.gz",repos=NULL)
# library("mpQTL")


# Set your working directory again in 'myanalyses'
setwd("C:/workshop/myanalyses")


#the object "data" contains an example dataset
data <- readRDS("../mpQTL/new_workshop_data.RDS")
names(data)

map <- data$map
phe <- data$pheno
snp <- as.matrix(data$snp)
anc <- as.matrix(data$founder)
sample_pv <- -log10(data$result[[1]]$pval)
sample_pv2 <- -log10(data$result[[2]]$pval)
cof <- data$cofactor

pval <- function(res,n = 1) -log10(res[[n]]$pval)
ploidy <- 4

#load hap and hapmap (which you created during the previuos part of the vignette)
hap <- readRDS("hap.rds")
hapmap <- readRDS("hapmap.rds")

## haplotypes ('haplotype names' format)
hap[1:5,1:8]
## haplotype map
head(hapmap)



# 1 Recall data structures -----------

#Genetic map has three columns, "marker", "position" and "chromosome"
map[1:20,]
#SNP dosages matrices are expressed with markers in rows and individuals in columns
snp[1:20,1:6]
#haplotypes are expressed with p columns per individual, where p is ploidy
hap[1:20,1:8] #8 columns of a tetraploid = 2 individuals
#phenotypes are expressed as a numeric matrix, with individuals in rows and each trait in a column
phe[1:20,]

#2 Genetic Structure (K) --------------------
pop <- substr(rownames(phe),1,2)
pop
table(pop)

Kd <- calc.K(t(snp))
Kh <- calc.K(t(hap),haplotypes = T, ploidy = 4)

#We can visualize the matrices using heatmaps
heatmap(Kd, Colv = NA, Rowv = NA, main = "Dosage matrix")
heatmap(Kh, Colv = NA, Rowv = NA, main = "Haplotype matrix")

#Or a bit more clearly using PCoA plots
pcoa.plot(Kd,col = pop, main = "Dosage matrix")
pcoa.plot(Kh,col = pop, main = "Haplotype matrix")


#3 Linear QTL models ---------------
#3.1 Linear -------------

#With dosages
result_dos <- map.QTL(phenotypes = phe,
                      genotypes = snp,
                      ploidy = 4,
                      map = map)
str(result_dos,max.level = 2, give.attr = F)
result_dos[[1]]$beta[1:3]

pv <- pval(result_dos)
skyplot(pv,map,main="QTL detection with linear model using dosages")

#With haplotypes
result_hap <- map.QTL(phenotypes = phe,
                      genotypes = hap,
                      ploidy = 4,
                      map = hapmap)

str(result_hap,max.level = 2, give.attr = F)
result_hap$phenotype1$beta[1:3]

pv <- pval(result_hap)
skyplot(pv,hapmap,main="QTL detection with linear model using haplotypes")

#3.2 Linear + Q ---------------
result_Q <- map.QTL(phenotypes = phe,
                    genotypes = snp,
                    ploidy = 4,
                    map = map,
                    Q = pop)

result_Q$phenotype1$beta[1]

pvQ <- pval(result_Q)
skyplot(pvQ,map,main="QTL detection with linear model and Q correction")

#3.3 Linear + Qpco ----------------
result_Qpco <- map.QTL(phenotypes = phe,
                       genotypes = hap,
                       ploidy = 4,
                       map = hapmap,
                       Q = T,
                       Qpco = 2)

result_Qpco$phenotype1$beta[1]

pvQpco <- pval(result_Qpco)
skyplot(pvQpco,hapmap,main="QTL detection with linear model and Qpco correction")

#We can check the similarity between pvalues with both methods
comp.skyplot(list(Q = pvQ,Qpco = pvQpco),map = list(map,hapmap),
             pch = c(4,1),main="Comparison between Q and Qpco corrections")

#4 Mixed QTL models ------------------

#4.1 Mixed -------------
result_mix <- map.QTL(phenotypes = phe,
                      genotypes = snp,
                      ploidy = 4,
                      map = map,
                      K = T,
                      cM = 1)

str(result_mix,max.level = 2, give.attr = F)

pv <- -log10(result_mix$phenotype1$pval)
skyplot(pv,map,main="QTL detection with mixed model")

result_mix$phenotype1$beta[1:3]


#5 Cofactors and covariates ----------------
table(cof)

result_cof <- map.QTL(phenotypes = phe,
                      genotypes = snp,
                      ploidy = 4,
                      map = map,
                      K = T,
                      cofactor = cof,
                      cofactor.type = "categorical")

without_cofactor <- -log10(result_mix$phenotype2$pval)
with_cofactor <- -log10(result_cof$phenotype2$pval)
comp.skyplot(list(without_cofactor,with_cofactor),map,legnames = c("No cofactor","Cofactor"),
             main = "QTL detection with and without cofactor")

# 6 Significance threshold ------------
result_perm <- map.QTL(phenotypes = phe,
                       genotypes = snp,
                       ploidy = 4,
                       map = map,
                       Q = T,
                       permutation = "pop",
                       nperm = 20,
                       alpha = 0.05)

str(result_perm,max.level = 2,give.attr = F)

## for SNPs only, you can use the Li & Ji method
result_liji <- thr.LiJi(snp,
                        map$chromosome,
                        alpha = 0.05,
                        ploidy = 4)
result_liji
-log10(result_liji$threshold)

# 7 Genotype imputation --------------
snp_NA <- snp
snp_NA[sample(1:length(snp),2000)] <- NA
result_NA <- map.QTL(phenotypes = phe,
                     genotypes = snp_NA,
                     ploidy = 4,
                     map = map,
                     Q = T)

str(result_NA,max.level = 2,give.attr = F)

# 8 Visualisation ------------------

#8.1 PCoA ---------------
K <- calc.K(t(snp))
#We must transform the result of dist into a matrix to use it with pcoa.plot()
euc <- as.matrix(dist(t(snp),method = "euclidean"))

pop <- substr(colnames(K),1,2)
table(pop)

#Interpretation examples
layout(matrix(1:4,byrow = T,ncol=2))
par(mar=c(4,4,2,1))
pcoa.plot(K,col = pop,
          pch=19,main="PCoA plot of K (similarity)", h = c(0,360))

pcoa.plot(1-K,col = pop,
          pch=19,main="PCoA plot of 1-K (distance)", h = c(0,360))

pcoa.plot(euc,col = pop,
          pch=19,main="PCoA plot of Euclidean distance", h = c(0,360))

#If we want to, we can choose different components of the PCoA, for instance the
#second and third
pcoa.plot(K,col = pop,comp = c(5,6),
          pch=19,main="PCoA plot of K (similarity)", h = c(0,360))

# Visualization examples
pcoa.plot(K,main="Default PCoA plot", h = c(0,360),plot_legend = F)

pcoa.plot(K,col = pop,
          pch=19,main="PCoA plot with colour mapping", h = c(0,360))

pcoa.plot(K,col = pop,pch = c(15,17,19),
          main="PCoA plot with different point types", h = c(0,360))

#In this case the legend does not work properly, it's better to turn it off
pcoa.plot(K,col = 1:nrow(K),
          pch=19,main="PCoA plot with  different color per individual")

dev.off()

#8.2 QQ-plot --------
#The p-values are from 500 random normal values
not_sig <- pnorm(rnorm(500),lower.tail = F)

#The p-values are from some random normal values and some non-random
some_sig <- pnorm(c(rnorm(450),rnorm(50,mean = 3)),lower.tail = F)

#The p-values are all too significant
high_sig <- pnorm(rnorm(1000,mean=3),lower.tail = F)

#The p-values are all too non-significant
low_sig <- pnorm(rnorm(300),sd = 3,lower.tail = F)

pvals <- list(not_sig,some_sig,high_sig,low_sig)
QQ.plot(pvals,main="Example QQ-plot",
        legnames = c("Not significant",
                     "Good significant",
                     "Overestimated",
                     "Underestimated"))

#We can also look at the results we have been producing
pvals <- list(lin_dosage = pval(result_dos),
              lin_hap = pval(result_hap),
              lin_Q = pval(result_Q),
              lin_Qpco = pval(result_Qpco),
              mix = pval(result_mix),
              lin_Q_NA = pval(result_NA))
#We need to undo the -log10 pval, because QQplot performs the transformation internally
pvals <- lapply(pvals,function(p) 10^-p)
QQ.plot(pvals, main = "Model comparison for p-values of pheno1",legspace = 0.22)
QQ.plot(pvals[-1:-2], main = "Model comparison for pvalues of pheno1")

#8.3 Skyline plots ------
pv1 <- pval(result_mix)

#Just a skyline plot
skyplot(pv1,map = map,main="Example Skyline plot")

#We want to focus on chromosome 2
skyplot(pv1,map,chrom = 2,
        main ="Skyline lot of chromosome 2",
        threshold = 3.6)

#The "chromosomes" can also be expressed as characters
map_letters <- map
map_letters$chromosome <- LETTERS[map_letters$chromosome + 1]

skyplot(pv1,map_letters,
        main="Skyline plot where chromosomes are characters",
        chrom = c("C","E"),small = T)

#We can also create "comparative" skyplots
pv2 <- pval(result_mix,2)
comp.skyplot(list(pv1,pv2),map,
             main = "Comparison of two p-value distributions",
             threshold = 4,pch = c(19,21))

# And compare the different models
#Remember skyplot will not tranform values into -log10
pvals <- lapply(pvals,function(i) -log10(i))
comp.skyplot(pvals,map = list(map,hapmap,map,hapmap,map,map),
             legspace = 0.22)
#Again the non-corrected models are too outlying
comp.skyplot(pvals[-1:-2],
             map= list(map,hapmap,map,map),
             pch = c(15,17,19),legspace = 0.22,
             main= "Comparison of models on example data")

# 8.4 Phenotype boxplots ----------------
#We choose the most significant marker in our QTL analysis
best_snp <- which.min(result_mix$phenotype1$pval)
best_hap <- which.min(result_Qpco$phenotype1$pval)
gen_snp <- unlist(snp[best_snp,])
gen_hap <- unlist(hap[best_hap,])

pheno_box(phe[,1],gen_snp,
          xlab="Dosage",ylab="phenotype",main="A boxplot of dosages")

#But now there are too many things plotted and I can't see anything
pheno_box(phe[,1],gen_hap,haplotype = T,ploidy = 4,
          xlab="Haplotypes",ylab="phenotype",main="A boxplot of haplotype dosages")

#This is better but still too many boxes
pheno_box(phe[,1],gen_hap,haplotype = T,ploidy = 4,draw.points = F,
          xlab="Haplotypes",ylab="phenotype",main="A boxplot of haplotype dosages (no points)")

#This is better
pheno_box(phe[,1],gen_hap,haplotype = T,ploidy = 4,
          hap.select = c("01","03", "06","08"),
          xlab="Haplotypes",ylab="phenotype",
          main="A boxplot of some haplotype dosages")

#If you'd like to know the underlying haplotype data (the haplotype dosages)
#you can use the dosage.X function
dosage.X(gen_hap,haplotype = T,ploidy = 4)

#8.5 Colour -------------
layout(matrix(1:4,byrow = T, ncol = 2))
for(l in c(40,60,80,100)) hue_wheel(l)
dev.off()



