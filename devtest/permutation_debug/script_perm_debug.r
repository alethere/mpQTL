
source("R/mapQTL_fun.R")



# data ----------------------------------------
set.seed(3)
pheno <- matrix(c(sample(c(5:10,NA), 50, replace = T),
                  sample(c(5:10), 50, replace = T)), ncol = 2)
pheno <- cbind(pheno, pheno, pheno)


set.seed(5)
geno <- matrix(sample(0:4, 500*50, replace = T), ncol = 50)
geno[1:5,1:5]



map <- data.frame(marker=paste0("snp",1:500),
                  chromosome=c(rep(1,300),
                               rep(2,200)),
                  position=c(seq(0,100, length.out = 300),
                             seq(0,100, length.out = 200)))
head(map)





# bug cused by NA in phenotypes, when using the linear model ----------

## select a very small dataset
pheno2 <- pheno[10:19,1:2]
pheno3 <- cbind(pheno2[,c(1,2)], pheno2[c(6:10,1:5),1])
pheno4 <- cbind(pheno2[,c(2)], pheno2[c(6:10,1:5),2])
geno2 <- geno[1:3,10:19]


## run the analysis below comparing different phenotype datasets:
## 1) pheno2[,2] (a phenotype without NAs)
## 2) pheno4 (adding a phenotype without NAs)
## 3) pheno2 (adding a phenotype with NA)
## the p-values of pheno2[,2] when anylized together with a phenotype with NA (example 3)

res <- map.QTL(phenotypes = pheno2[,2],
               genotypes = geno2, #genotype matrix
               ploidy = 4,
               map = map[1:3,], #genetic map table
               K=NULL, #distance matrix
               Q=NULL, #population effect matrix
               Z=NULL,
               cofactor=NULL,
               cofactor.type=NULL,
               cM=1, #
               seed=NULL,
               Qpco=2, #number of axis used for pco decomposition
               no_cores=parallel::detectCores()-1,
               approximate = T,
               permutation = NULL, #permutation strategy: "pop" or "fam"
               nperm = 2, #number of permutations
               alpha = 0.95,
               impute=T,
               k=20,
               linear = NULL,
               K_identity = F)

res

res$pheno1
res$pheno2

plot(-log10(res$pheno1$pval))
all.equal(temp1, res$pheno1$pval)
all.equal(temp2, res$pheno2$pval)
temp1 <- res$pheno1$pval
temp2 <- res$pheno2$pval



# analysis with permutation --------------------
res <- map.QTL(phenotypes = pheno[,c(1,2,4)],
               genotypes = geno, #genotype matrix
               ploidy = 4,
               map = map, #genetic map table
               K=NULL, #distance matrix
               Q=NULL, #population effect matrix
               Z=NULL,
               cofactor=NULL,
               cofactor.type=NULL,
               cM=1, #
               seed=NULL,
               Qpco=2, #number of axis used for pco decomposition
               no_cores=parallel::detectCores()-1,
               approximate = T,
               permutation = NULL, #permutation strategy: "pop" or "fam"
               nperm = 3, #number of permutations
               alpha = 0.95,
               impute=T,
               k=20,
               linear = NULL,
               K_identity = F)



plot(-log10(res$pheno1$pval))
all.equal(temp1, res$pheno1$pval)
all.equal(temp2, res$pheno2$pval)
temp1 <- res$pheno1$pval
temp2 <- res$pheno2$pval




res0 <- map.QTL(phenotypes = pheno[,1],
               genotypes = geno, #genotype matrix
               ploidy = 4,
               map = map, #genetic map table
               K=T, #distance matrix
               Q=NULL, #population effect matrix
               Z=NULL,
               cofactor=NULL,
               cofactor.type=NULL,
               cM=1, #
               seed=NULL,
               Qpco=2, #number of axis used for pco decomposition
               no_cores=parallel::detectCores()-1,
               approximate = T,
               permutation = "pop", #permutation strategy: "pop" or "fam"
               nperm = 3, #number of permutations
               alpha = 0.95,
               impute=T,
               k=20,
               linear = NULL,
               K_identity = F)

res0$pheno1$pval[1:10]
plot(-log10(res0$pheno1$pval))














# covariate analysis --------------------
# set.seed(5)
# mycov <- matrix(sample(c(0:2,NA),50,replace = T), ncol = 1)
mycov <- geno[which.max(res0$pheno1$pval),]

res2 <- map.QTL(phenotypes = pheno,
               genotypes = geno, #genotype matrix
               ploidy = 4,
               map = map, #genetic map table
               K=NULL, #distance matrix
               Q=NULL, #population effect matrix
               Z=NULL,
               cofactor=mycov,
               cofactor.type="num",
               cM=1, #
               seed=NULL,
               Qpco=2, #number of axis used for pco decomposition
               no_cores=parallel::detectCores()-1,
               approximate = T,
               permutation = "pop", #permutation strategy: "pop" or "fam"
               nperm = 2, #number of permutations
               alpha = 0.95,
               impute=T,
               k=20,
               linear = NULL,
               K_identity = F)

res2$pheno1$perm.thr
res2$pheno1$pval[1:10]
plot(-log10(res2$pheno1$pval))


#
mycov <- geno[which.max(-log10(res0$pheno1$pval)),]

res3 <- map.QTL(phenotypes = pheno,
                genotypes = geno, #genotype matrix
                ploidy = 4,
                map = map, #genetic map table
                K=T, #distance matrix
                Q=NULL, #population effect matrix
                Z=NULL,
                cofactor=mycov,
                cofactor.type="num",
                cM=1, #
                seed=NULL,
                Qpco=2, #number of axis used for pco decomposition
                no_cores=parallel::detectCores()-1,
                approximate = T,
                permutation = "pop", #permutation strategy: "pop" or "fam"
                nperm = 2, #number of permutations
                alpha = 0.95,
                impute=T,
                k=20,
                linear = NULL,
                K_identity = F)


plot(-log10(res3$pheno1$pval))






