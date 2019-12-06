#' Workshop December 2019 -------------------
#' Code and explanations on mpQTL v0.2.0
#' by Alejandro Thérèse Navarro & Giorgio Tumino

#Installing and loading the package
install.packages("../mpQTL_0.2.0.tar.gz",repos=NULL)
library("mpQTL")

#the object "data" contains an example dataset
data <- readRDS("new_workshop_data.RDS")
map <- data$map
phe <- data$pheno
snp <- as.matrix(data$snp)
anc <- as.matrix(data$founder)
sample_pv <- -log10(data$result[[1]]$pval)
sample_pv2 <- -log10(data$result[[2]]$pval)
cof <- data$cofactor
hapmap <- data$hapmap

ploidy <- 4
hap_names <- sapply(colnames(data$PolyHap_res$chr1_1.1$hapdos),function(n) paste0(n,"_",1:ploidy))
hap_list <- lapply(data$PolyHap_res,function(chrom){
  hap_mat <- chrom$hapdos
  haps <- rownames(hap_mat)
  hap_vec <- as.vector(apply(hap_mat,2,function(h) rep(haps,h)))
  return(hap_vec)
})

hap <- do.call(rbind,hap_list)
rownames(hap) <- names(hap_list)
colnames(hap) <- hap_names

pval <- function(res,n = 1) -log10(res[[n]]$pval)

#Data structures -----------

#Genetic map has three columns, "marker", "position" and "chromosome"
map[1:20,]
#SNP dosages matrices are expressed with markers in rows and individuals in columns
snp[1:20,1:6]
#haplotypes are expressed with p columns per individual, where p is ploidy
hap[1:20,1:8] #8 columns of a tetraploid = 2 individuals
#phenotypes are expressed as a numeric matrix, with individuals in rows and each trait in a column
phe[1:20,]

#1 Linear QTL models ---------------
#1.1 Linear -------------

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

pv <- pval(result$hap)
skyplot(pv,hapmap,main="QTL detection with linear model using haplotypes")

#1.2 Linear + Q ---------------
pop <- substr(rownames(phe),1,2)
pop
table(pop)

result_Q <- map.QTL(phenotypes = phe,
                    genotypes = snp,
                    ploidy = 4,
                    map = map,
                    Q = pop)

result_Q$phenotype1$beta[1]

pvQ <- pval(result_Q)
skyplot(pvQ,map,main="QTL detection with linear model and Q correction")

#1.3 Linear + Qpco ----------------
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

#2 Mixed QTL models ------------------

#2.1 Mixed -------------
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

#2.2 Kinship matrix ------------------
Kd <- calc.K(t(snp))
Kh <- calc.K(t(hap),haplotypes = T, ploidy = 4)

#We can visualize the matrices using heatmaps
heatmap(Kd, Colv = NA, Rowv = NA, main = "Dosage matrix")
heatmap(Kh, Colv = NA, Rowv = NA, main = "Haplotype matrix")

#Or a bit more clearly using PCoA plots
pcoa.plot(Kd,col = pop, main = "Dosage matrix")
pcoa.plot(Kh,col = pop, main = "Haplotype matrix")

#3 Cofactors and covariates ----------------
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

# 4 Permutation threshold ------------
result_perm <- map.QTL(phenotypes = phe,
                       genotypes = snp,
                       ploidy = 4,
                       map = map,
                       Q = T,
                       permutation = "pop",
                       nperm = 10,
                       alpha = 0.95)

str(result_perm,max.level = 2,give.attr = F)

# 5 Genotype imputation --------------
dos_NA <- snp
dos_NA[sample(1:length(snp),2000)] <- NA
result_NA <- map.QTL(phenotypes = phe,
                     genotypes = dos_NA,
                     ploidy = 4,
                     map = map,
                     Q = T)

str(result_NA,max.level = 2,give.attr = F)

# 6 Visualisation ------------------

#6.1 PCoA ---------------
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

gpdev.off()

#6.2 QQ-plot --------
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
QQ.plot(pvals, main = "Model comparison for p-values of pheno1",legspace = 0.22)
QQ.plot(pvals[-1:-2], main = "Model comparison for pvalues of pheno1")

#6.3 Skyline plots ------
pv1 <- -log10(result_mix$phenotype1$pval)

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
pv2 <- -log10(result_mix$phenotype2$pval)
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

# 6.4 Phenotype boxplots ----------------
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
          hap.select = c("H_01","H_03", "H_06","H_08"),
          xlab="Haplotypes",ylab="phenotype",
          main="A boxplot of some haplotype dosages")

#If you'd like to know the underlying haplotype data (the haplotype dosages)
#you can use the dosage.X function
dosage.X(gen_hap,haplotype = T,ploidy = 4)

#6.5 Colour -------------
layout(matrix(1:4,byrow = T, ncol = 2))
for(l in c(40,60,80,100)) hue_wheel(l)
dev.off()


