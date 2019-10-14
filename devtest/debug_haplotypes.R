source("mpQTL_fun_clean.R")
data<-readRDS("output_inferredHaps/data_debug_impute.knn.rds")
phenotypes<-readRDS("Phenotypes_for_Giorgio.RDS")





phenotypes = phenotypes$cross200
genotypes = data #genotype matrix
ploidy = 4
map = cross200_0.1$map #genetic map table
K=T #distance matrix
Q=NULL #population effect matrix
Z=NULL
cofactor=NULL
cofactor.type=NULL
cM=1 #
Qpco=2 #number of axis used for pco decomposition
no_cores=6
P3D=T
EMMAX=T
permutation = NULL #permutation strategy: "pop" or "fam"
nperm = NULL #number of permutations
impute=T
k=20


matrix<-apply(matrix,2,dosage.X,ploidy=ploidy)
matrix<-do.call(cbind,matrix)


genotypes<-data$genotypes[6,]


