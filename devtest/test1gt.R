

founderalleles[1:5,1:5] #phased founder alleles
genotypes[1:5,1:5] #phased snp alleles

g[1:5,1:5] # dosage matrix
sum(is.na(g))



# pheno standardization ------------------
phenotypes1 <- matrix(rnorm(40,100,10), ncol = 4)
nas <- sample(1:length(phenotypes1), length(phenotypes1)/10)
phenotypes1[nas] <- NA
rownames(phenotypes1) <- paste0("ind", 1:nrow(phenotypes1))
colnames(phenotypes1) <- paste0("pheno", 1:ncol(phenotypes1))

phenotypes2 <- phenotypes1[,1]

phenotypes3 <- as.data.frame(phenotypes1)


#We need to standardize the phenotypes
#Put them in a matrix if they come in a vector
f1gt <- function(phenotypes) {
  if(is.vector(phenotypes)) {
    phenotypes<-(phenotypes-mean(phenotypes, na.rm = T))/sd(phenotypes, na.rm = T)
    inds<-names(phenotypes)
    phenotypes<-matrix(phenotypes,ncol=1)
    rownames(phenotypes)<-inds
    return(phenotypes)
  }else{
    inds<-rownames(phenotypes)
    colns<-colnames(phenotypes)
    phenotypes<-apply(phenotypes,2,function(p){
      (p-mean(p, na.rm = T))/sd(p, na.rm = T)
    })
    rownames(phenotypes)<-inds
    colnames(phenotypes)<-colns
    return(phenotypes)
  }
}

phenotypes1res <- f1gt(phenotypes1)
class(phenotypes1res)
phenotypes2res <- f1gt(phenotypes2)
class(phenotypes2res)
phenotypes3res <- f1gt(phenotypes3)
class(phenotypes3res)


f2gt <- function(phenotypes) {
  if(is.vector(phenotypes)) {
    phenotypes <- matrix(phenotypes, ncol=1, 
                         dimnames = list(names(phenotypes)))
  } else if (is.data.frame(phenotypes)) {
    phenotypes <- as.matrix(phenotypes)
  }
  phenotypes <- apply(phenotypes,2,function(p){
    (p-mean(p, na.rm = T))/sd(p, na.rm = T)
  })
}

phenotypes1resf2 <- f2gt(phenotypes1)
class(phenotypes1resf2)
phenotypes2resf2 <- f2gt(phenotypes2)
class(phenotypes2resf2)
phenotypes3resf2 <- f2gt(phenotypes3)
class(phenotypes3resf2)

identical(phenotypes1res,phenotypes1resf2)
identical(phenotypes2res,phenotypes2resf2)
identical(phenotypes3res,phenotypes3resf2)



# genotypes curation --------------------------

n <- 4
ploidy <- 4

## pheno
pheno <- matrix(rnorm(n,100,10), ncol = 1)
nas <- sample(1:length(pheno), length(pheno)/10)
pheno[nas] <- NA
rownames(pheno) <- paste0("ind", 1:nrow(pheno))
colnames(pheno) <- paste0("pheno", 1:ncol(pheno))
pheno

## haplotypes
gen <- matrix(sample(letters, 10*n*ploidy, replace = T), ncol=n*ploidy)
colnames(gen) <- sapply(rownames(pheno), function(x) paste(x, 1:ploidy, sep = "_"))
gen
gen2 <- gen[,-1]
gen3 <- gen
colnames(gen3)[1] <- "ind9_1"

## dosages
dos <- matrix(sample(0:ploidy, 10*n, replace = T), ncol=n)
colnames(dos) <- rownames(pheno)
dos

dos2 <- dos
colnames(dos2)[1] <- "ind9"


map.QTL(
  phenotypes = pheno,
  genotypes = dos, #genotype matrix
  ploidy = ploidy,
  map, #genetic map table
  K=NULL, #distance matrix
  Q=NULL, #population effect matrix
  Z=NULL,
  cofactor=NULL,
  cofactor.type=NULL,
  dosage = NULL, #dosage matrix
  cM=1, #
  Qpco=2, #number of axis used for pco decomposition
  no_cores=parallel::detectCores()-1,
  P3D=T,
  EMMAX=T
)




# kinship ----------------------------------
rownames(g) <- g[,1]
g <- g[,-1]
class(g)
K2 <- calc.K(t(g))
identical(K,K2)

g2 <- as.matrix(g)
dim(g2)
g2[1:10,1:5]
nas <- sample(1:length(g2), length(g2)/1000)
g2[nas] <- NA
K3 <- calc.K(t(g2))
K3[1:5,1:5]

g3 <- g2[1:10,1:5]
nas <- sample(1:length(g3), length(g3)/10)
g3[nas] <- NA
g3

imputeNA <- function(m,
                     miss = NULL) {
  
  if(!is.null(miss) & is.numeric(miss)) m[m==miss] <- NA
  
  meanDose <- rowMeans(m, na.rm=T)
  na <- is.na(m)
  means <- na * meanDose
  m[na] <- means[na]
  
  # if (any(is.na(m))) stop("There are still NAs!!")
  return(m)
}
g3noNA <- imputeNA(g3)

K3 <- calc.K(t(g2))
K3[1:5,1:5]



# Q matrix --------------------------------
Q <- sample(letters[1:3], 20, replace = T)
pop <- Q
Q<-Q.mat(Q)


















# permutation ----------------------------
phenotypes <- matrix(rnorm(20,100,10), ncol = 2)
nas <- sample(1:length(phenotypes), length(phenotypes)/10)
phenotypes[nas] <- NA
rownames(phenotypes) <- paste0("ind", 1:nrow(phenotypes))
colnames(phenotypes) <- paste0("pheno", 1:ncol(phenotypes))


nperm <- 3
fam <- c(rep("p",2),
         rep("f1",4),
         rep("f2",4))

##
sample_data <- readRDS("sample_data.rds")
names(sample_data)
sample_data$phenotype[1:5,]
sample_data$genotypes[1:5,1:5]
dim(sample_data$genotypes)
sample_data$genotypes.na005[1:5,1:5]
dim(sample_data$dosage)
sample_data$dosage[1:5,1:5]
rownames(sample_data$phenotype)[1:20]
colnames(sample_data$dosage)[1:20]



### individual names in phenotypes are wrong
rownames(sample_data$phenotype) <- colnames(sample_data$dosage)



### test with no missing values
res1 <- map.QTL(
  phenotypes = sample_data$phenotype,
  genotypes = sample_data$genotypes, #genotype matrix
  ploidy = 4,
  map = sample_data$map, #genetic map table
  K=T, #distance matrix
  Q=NULL, #population effect matrix
  Z=NULL,
  cofactor=NULL,
  cofactor.type=NULL,
  cM=1, #
  Qpco=2, #number of axis used for pco decomposition
  no_cores=6,
  P3D=T,
  EMMAX=T,
  permutation = NULL, #permutation strategy: "pop" or "fam"
  nperm = NULL, #number of permutations 
  impute=T,
  k=20
)

str(res1[[1]])
names(res1)
names(res1[[1]])
res1pval <- sapply(res1, function(x) x$pval)
res1pval[1:5,]
plot(-log10(res1pval[,1]))
plot(-log10(res1pval[,2]))
plot(-log10(res1pval[,3]))



# ### bug
# res2 <- myf(phenotypes = pheno1,
#             genotypes = geno1,
#             dosage = dosage1,
#             map = sample_data$map,
#             w = 1)
# res2[1:10]

sample_data$genotypes.na001[1:10,1:10]
sum(is.na(sample_data$genotypes.na001))

### test with missing values
res3 <- map.QTL(
  phenotypes = sample_data$phenotype,
  genotypes = sample_data$genotypes.na001, #genotype matrix
  ploidy = 4,
  map = sample_data$map, #genetic map table
  K=T, #distance matrix
  Q=NULL, #population effect matrix
  Z=NULL,
  cofactor=NULL,
  cofactor.type=NULL,
  cM=1, #
  Qpco=2, #number of axis used for pco decomposition
  no_cores=6,
  P3D=T,
  EMMAX=T,
  permutation = NULL, #permutation strategy: "pop" or "fam"
  nperm = NULL #number of permutations 
)

names(res3)
names(res3[[1]])
res3pval <- sapply(res3, function(x) x$pval)
res3pval[1:5,]
plot(-log10(res3pval[,1]))
plot(-log10(res3pval[,2]))
plot(-log10(res3pval[,3]))



### test permutation, mixed model
res4 <- map.QTL(
  phenotypes = sample_data$phenotype,
  genotypes = sample_data$genotypes, #genotype matrix
  ploidy = 4,
  map = sample_data$map, #genetic map table
  K=T, #distance matrix
  Q=NULL, #population effect matrix
  Z=NULL,
  cofactor=NULL,
  cofactor.type=NULL,
  cM=1, #
  Qpco=2, #number of axis used for pco decomposition
  no_cores=6,
  P3D=T,
  EMMAX=T,
  permutation = "pop", #permutation strategy: "pop" or "fam"
  nperm = 3, #number of permutations
  alpha = c(0.95,0.99)
)


names(res4)
res4$h08$perm.thr
names(res4[[1]])

res4pval <- sapply(res4, function(x) x$pval)
res4pval[1:5,]
plot(-log10(res4pval[,1]))
abline(h=res4[[1]]$perm.thr[1])
abline(h=res4[[1]]$perm.thr[2])
plot(-log10(res4pval[,2]))
abline(h=res4[[2]]$perm.thr[1])
abline(h=res4[[2]]$perm.thr[2])
plot(-log10(res4pval[,3]))
abline(h=res4[[3]]$perm.thr[1])
abline(h=res4[[3]]$perm.thr[2])


# ---------------------
tags <- colnames(sample_data$genotypes)
pmatch(tags[1], tags[1:8])
a <- strsplit(tags[1:12], "")
sapply(a, function(x) {
  intersect(a[[1]],a[[2]])
  a[[1]]==a[[2]]
})


