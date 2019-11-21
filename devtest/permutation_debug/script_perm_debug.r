
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


# analysis with permutation --------------------
res <- map.QTL(phenotypes = pheno,
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
               permutation = "pop", #permutation strategy: "pop" or "fam"
               nperm = 3, #number of permutations
               alpha = 0.95,
               impute=T,
               k=20,
               linear = NULL,
               K_identity = F)

res$pheno1$perm.thr
res$pheno1$pval[1:10]
plot(-log10(res$pheno1$pval))


res <- map.QTL(phenotypes = pheno,
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
               permutation = NULL, #permutation strategy: "pop" or "fam"
               nperm = NULL, #number of permutations
               alpha = 0.95,
               impute=T,
               k=20,
               linear = NULL,
               K_identity = F)

res$pheno1$pval[1:10]
plot(-log10(res$pheno1$pval))





# covariate analysis --------------------
set.seed(5)
# mycov <- matrix(sample(c(0:2,NA),50,replace = T), ncol = 1)
mycov <- geno[which.max(res$pheno1$pval),]

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
               nperm = 3, #number of permutations
               alpha = 0.95,
               impute=T,
               k=20,
               linear = NULL,
               K_identity = F)

res2$pheno1$perm.thr
res2$pheno1$pval[1:10]
plot(-log10(res2$pheno1$pval))

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
                nperm = 3, #number of permutations
                alpha = 0.95,
                impute=T,
                k=20,
                linear = NULL,
                K_identity = F)


plot(-log10(res3$pheno1$pval))






