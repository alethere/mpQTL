
source("../mpQTL_fun_clean.R")



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







# covariate analysis --------------------
set.seed(5)
mycov <- matrix(sample(c(0:2,NA),50,replace = T), ncol = 1)


res <- map.QTL(
  phenotypes = pheno, 
  genotypes = geno,  
  ploidy = 4,
  map = map, 
  K = T,
  cofactor = mycov,
  cofactor.type = "num",
  cM = 2,
  no_cores = 6, 
  approximate = T)


res$pheno1$pval[1:10]







