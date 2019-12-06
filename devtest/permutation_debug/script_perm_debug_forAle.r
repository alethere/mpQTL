
source("R/mapQTL_fun.R")



# data ----------------------------------------
## create a small dataset with 10 individuals, 5 markers and 3 phenotypes
set.seed(3)
pheno <- matrix(c(sample(c(5:10), 10, replace = T),
                  sample(c(5:10), 10, replace = T),
                  sample(c(5:10,NA), 10, replace = T)), ncol = 3)
colnames(pheno) <- c("A","B","C")

set.seed(5)
geno <- matrix(sample(0:4, 5*10, replace = T), ncol = 10)
geno



map <- data.frame(marker=paste0("snp",1:5),
                  chromosome=c(rep(1,3),
                               rep(2,2)),
                  position=c(seq(0,100, length.out = 3),
                             seq(0,100, length.out = 2)))
map





# bug cused by NA in phenotypes, when using the linear model ----------
## run the analysis below comparing different phenotype datasets:

## 1) only pheno "A" (a phenotype without NAs)
res1 <- map.QTL(phenotypes = pheno[,"A"],
               genotypes = geno, #genotype matrix
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
               nperm = NULL, #number of permutations
               alpha = 0.95,
               impute=T,
               k=20,
               linear = NULL,
               K_identity = F)


## 2) pheno "A" and "B" (adding a phenotype without NAs)
res2 <- map.QTL(phenotypes = pheno[,c("A","B")],
                genotypes = geno, #genotype matrix
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
                nperm = NULL, #number of permutations
                alpha = 0.95,
                impute=T,
                k=20,
                linear = NULL,
                K_identity = F)

## 3) now add pheno "C" (adding a phenotype with NA)
res3 <- map.QTL(phenotypes = pheno,
                genotypes = geno, #genotype matrix
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
                nperm = NULL, #number of permutations
                alpha = 0.95,
                impute=T,
                k=20,
                linear = NULL,
                K_identity = F)

data <- test.compatibility(pheno,dosage.X(geno[1,]),Z = diag(nrow(pheno)), K = diag(nrow(pheno)))
length(data)
## compare the p-values of pheno "A" across the three analyses (for the first marker only)
res1[[1]]$pval
res2[[1]]$pval
res3[[1]]$pval
## res1 and res returned the same pval (0.0706), but in res3 it is different (0.0851)
## the results for one phenotype should not be affected by other phenotypes
## this situation is exacerbated in permutation tests

sapply(res1,'[[',"pval")
sapply(res3,'[[',"pval")
