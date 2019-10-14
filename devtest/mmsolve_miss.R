phenotypes<-data$phenotype
genotypes<-impgeno

impute.knn <- function(
  geno,
  ploidy,
  map = NULL,
  k = 20
){
  
  #We determine whether we are looking at dosages or haplotypes
    if (all(unique(geno) %in% c(0:ploidy, NA))){
      cat("Imputation performed on dosages")
      haplo<-F
    }else if (ncol(geno) %% ploidy == 0){
      cat("Imputation performed on haplotypes")
      haplo<-T
    }else stop(paste("Genotypes not recognized as dosage nor haplotypes.",
                    "Dosage values are not between 0 and ploidy or",
                    "ncol of genotypes is not a multiple of ploidy."))

  
  #Calculate K distance matrix. If possible using fewer markers
  #homogeneously distributed across the genome. 
  if (!is.null(map)) genoK <- sample.cM(geno, map)
  else genoK <- geno
  K <- calc.K(t(genoK), ploidy = ploidy, haplotypes = haplo)
  
  #Imputation on dosages
  #Each individual is a column
  if (!haplo) {
    result <- sapply(1:ncol(K), function(d) {
      #Find the closest individual to individual d
      best <- order(K[, d][-d], decreasing = T)[1:k]
      imputed <- geno[,d]
      
      #Get the genotypes of the best individuals
      #at the position where current individual has missing values
      knei <- geno[is.na(imputed), best]
      #Percentage of missing values per marker on the nearest neighbours
      naperc<-apply(knei,1,function(i) sum(is.na(i))/k)
      if(any(naperc>0.5)){
        markers<-names(naperc)[naperc>0.5]
        error<-paste("Nearest neighbours of individual",d,
                     "have over 50% missing values at:",
                     paste(markers,collapse=" "),
                     "\nNon confident imputation. Try increasing k? Remove markers?")
        warning(error)
      }
      
      #Are they expressed numerically? The following function coerces to char
      asnum <- is.numeric(knei)
      #For each marker, select the most common dosage
      nearest <- apply(knei, 1, function(i) {
        res<-names(sort(table(i), decreasing = T)[1])
        if(asnum) res<-as.numeric(res)
        return(res)
      })
      
      #the nearest neigbhours are put in. Careful, if too many missing markers
      #less neigbhours are being used to impute.
      imputed[is.na(imputed)] <- nearest
      return(imputed)
    })
  }
  
  #The same thing must be done for haplotypes.
  #The main difference is that 1) haplotypes can be characters
  #and 2) each individual is represented by multiple columns
  #in the genotype matrix
  if (haplo) {
    result <- lapply(1:ncol(K), function(d) {
      #Best individual selection is done the same
      best <- order(K[-d, d], decreasing = T)[1:k]
      #but we need to transform individual index into 
      #the chromosome column index
      best <- as.vector(sapply(best, function(b)
          (b * ploidy - (ploidy - 1)):(b * ploidy)))
      
      #Same goes for the individual of interest
      dG <- (d * ploidy - (ploidy - 1)):(d * ploidy)
      #Missing values are detected in ANY of the haplotypes
      miss <- is.na(geno[, dG])
      miss <- apply(miss, 1, any)
      
      #We get the genotypes of the neighbours
      knei <- geno[miss, best]
      naperc<-apply(knei,1,function(i) sum(is.na(i))/ncol(knei))
      if(any(naperc>0.5)){
        markers<-names(naperc)[naperc>0.5]
        error<-paste("Nearest neighbours of individual",d,
                     "have over 50% missing values at:",
                     paste(markers,collapse=" "),
                     "\nNon confident imputation. Try increasing k? Remove markers?")
        warning(error)
      }
      
      #Are they expressed numerically? The following function coerces to char
      asnum <- is.numeric(knei)
      
      nearest <- t(apply(knei, 1, function(alleles) {
        #We get alleles ordered by frequency
        common <- sort(table(alleles) / length(alleles), decreasing = T)
        
        #We need to select the four most common alleles
        #This probably needs to be improved
        res <- c()
        for (i in 1:ploidy) {
          res[i] <- names(common)[1]
          common[1] <- common[1] - 1 / ploidy
          common <- sort(common, decreasing = T)
        }
        #coerce back to numeric if it was numeric
        if (asnum) res <- as.numeric(res)
        return(res)
      }))
      
      imputed <- geno[, dG]
      imputed[miss, ] <- nearest
      return(imputed)
    })
    result <- do.call(cbind, result)
  }
  
  return(result)
  
  }

geno2<-impute.knn(data$dosage.na005,4,k=50)

check<-is.na(data$dosage.na005)
imputedcheck<-as.vector(geno2)[as.vector(check)]
realcheck<-as.vector(data$dosage)[as.vector(check)]
cor(imputedcheck,realcheck)
# plot(imputedcheck,realcheck)
sum(imputedcheck==realcheck)/sum(check)

apply(geno2[1:10,],1,table)[[1]]
apply(data$genotypes[1:10,],1,table)[[1]]
profvis::profvis(impute.knn(data$genotypes.na005,4,map,k=8))

impdosage<-geno2
impgeno<-impute.knn(data$genotypes.na005,4,map,30)

ploidy<-4
map<-data$map
K=T #distance matrix
Q=NULL #population effect matrix
Z=NULL
cofactor=NULL
cofactor.type=NULL
dosage=impdosage #dosage matrix
cM=1 #
Qpco=2 #number of axis used for pco decomposition
no_cores=parallel::detectCores()-1
P3D=T
EMMAX=T
permutation = NULL #permutation strategy: "pop" or "fam"
nperm = NULL #number of permutations
alpha = 0.95


if(is.null(K)){
  linear<-T
}else{
  linear<-F
}

markers<-rownames(genotypes)

#We need to standardize the phenotypes
#Put them in a matrix if they come in a vector
#It allows missing values.
#ATN: Shouldn't we put this function outside of mapQTL?
stndrdz.pheno <- function(phenotypes) {
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
phenotypes <- stndrdz.pheno(phenotypes)



# Curation of genotype matrix and dosage matrix
## check whether genotypes contain snp dosages or haplotypes
if (ncol(genotypes)==nrow(phenotypes)) {
  geno.type <- "dosage"
  genotypes <- inputCheck_dos(genotypes, integer=T, ploidy=ploidy)
  cat("SNP dosages have been detected in the genotype matrix\n")

  if (!all(rownames(phenotypes)==colnames(genotypes))) {
    stop("phenotype individuals and genotype individuals have different names")
  }

} else if ((ncol(genotypes) %% ploidy) == 0){
  geno.type <- "haplo"
  cat("Haplotypes have been detected in the genotype matrix\n")

} else {
  stop("If haplotypes have been provided: the number of genotypes columns
       is not a multiple of ploidy.
       If SNP dosages have been provided: the number of individual dosages
       is not equal to the number of individual phenotypes.")
}


## check dosage matrix
if (is.null(dosage)) {
  if (geno.type=="dosage") {
    dosage <- genotypes
  } else {
    stop("If genotypes are not SNP marker dosages, SNP dosages must be provided")
  }
} else {
  if (!all(rownames(phenotypes)==colnames(dosage))) {
    #I see your point with these, but what if phenotypes does not have rownames,
    #or the ind names are a column in a phenotype matrix or something...
    #Also, if both have no rownames/colnames, it will not promt anything
    stop("Phenotype individuals and dosage individuals have different names")
  }
  dosage <- inputCheck_dos(dosage, integer=T, ploidy=ploidy)
}

#bi.geno
#NArate
#MAF
#adjust map accordingly


#gt: Snp dosage matrix and map are needed to calculate K and Q, when
#they are not provided. The map is also used for plotting the association
#results. However, when genotypes are haplotypes (e.g. inferred by
#PolyHaplotyper), the haplotype map might be different from the SNp dosage
#map. For the moment, we can avoid issues providing K and Q, but we need
#to think about a solution.


#DEFINITION OF K
#There are four options
#1) K=NULL, no K is used. We go linear.
#2) K=T, K is calculated using homogeneously distributed markers along a map.
#   using sample.cM and calc.K
#3) K=K matrix, check dimensions, message is printed and K is used.
if(all(K==T)){ #gt: why are you using all?
  #Option 2)
  K<-sample.cM(dosage,map,cM = cM)
  K<-calc.K(t(K))

}else if(nrow(K)==nrow(phenotypes)){
  #Option 3)
  if(ncol(K)!=nrow(phenotypes)) stop("K matrix was specified, but is not square")
  print("Square K matrix has been provided. Using as it is.")
}

#DEFINITION OF Q
#4 options:
#1) No Q is used (Q=NULL), this block is skipped
#2) Q=T, Q is calculated using the K matrix and cmdscale
#3) Q=vector identifying groups, Q is calculated using Q.mat
#4) Q=Q design matrix, message is printed.
if(is.null(Q)) print("No Q matrix will be used")

if(!is.null(Q)){
  if(all(Q==T)){ #gt: why are you using all?
    #atn: because if Q is a vector/matrix Q==T will be a vector of F
    #and that gives an error in the if.

    #Option 2)
    if(is.null(K)){
      if(all(rownames(phenotypes)==colnames(genotypes))){
        dosage<-genotypes
      }else if(is.null(dosage)){
        stop("If genotypes are not SNP marker dosage, SNP dosage must be specified")
      }

      #first we sample homogeneously markers along the genome
      K<-sample.cM(dosage,map,cM = cM)
      #then we calculate the distance
      K<-calc.K(t(K))
    }
    Q <- cmdscale(1-K, k=Qpco, eig = F, add = FALSE, x.ret = FALSE)

  }else if(is.vector(Q)){
    #Option 3)
    Q<-Q.mat(Q)

  }else if(nrow(Q)==nrow(phenotypes)){
    #Option 4)
    print("Q matrix has been provided. Using as it is.")

  }else{
    #Q was not NULL, T, a vector or a cofactor design matrix.
    stop("Wrong Q specification")}
}

##DEFINITION OF C, cofactor matrix.
#Two variables must be provided, cofactor and cofactor.type
#The first contains values, the second specifies whether the cofactor
#should be treated as a numeric regressor or as a categorical (ANOVA-type) parameter
if (!is.null(cofactor)){

  #Make cofactors into a matrix
  if(is.vector(cofactor)) cofactors<-matrix(cofactor,ncol=1)

  #Check cofactor tpyes are specified for cofactors.
  if(is.null(cofactor.type)|length(cofactor.type)!=ncol(cofactor)){
    stop(paste0("All cofactor types must be specified using the cofactor.type parameter.",
                "\nA char vector containing \"numerical\"",
                " or \"categorical\" for each cofactor must be given."))
  }

  #Which are num and which are cat
  c.num<-pmatch(cofactor.type,"numerical")
  c.cat<-pmatch(cofactor.type,"categorical")

  C<-lapply(1:ncol(cofactor),function(i){

    if(!is.na(c.num)){
      result<-as.numeric(cofactor[,i])

    }else if(!is.na(c.cat)){

      #This creates a design matrix with 1s and 0s (representing T and F)
      result<-sapply(unique(cofactor[,i]),function(c){
        cofactor[,i] == c
      })+1-1

      #Last column is taken out to be able to calculate
      result<-result[,-ncol(result)]
    }
  })
  C<-do.call(cbind,C)

}else{C<-NULL}



# Prepare permuted phenotypes
if (!is.null(permutation)) {
  npheno <- ncol(phenotypes) #number of phenotypes
  nind <- nrow(phenotypes)   #number of individuals
  if (permutation=="pop") { #permutation over the whole population
    permid <- cbind(1:nind,
                    replicate(nperm, sample(1:nind)))
  }
  if (permutation=="fam") { #permutation within families
    famid <- lapply(unique(fam), function(x) {
      which(fam==x)
    })
    permid <- cbind(1:nind,
                    replicate(nperm, do.call("c",lapply(famid, sample))))
  }
  phenotypes <- do.call("cbind",lapply(1:ncol(permid), function(i) {
    phenotypes[permid[,i],]
  }))
}

####
# Mixed Models.
if(is.null(Z)){
  Z<-diag(nrow(phenotypes))
}

if(!P3D|!EMMAX){
  Hinv<-NULL
}else{
  Hinv<-calc.Hinv(phenotypes,
                  #X=vector of 1s, to apply P3D/EMMAX algorithm #gt: check#
                  X=matrix(rep(1,nrow(phenotypes))),
                  Z,K)
}

cl<-parallel::makeCluster(no_cores)
if(is.null(Q)){
  print("Mixed model will be used")
}else{
  print("Mixed model with Q correction will be used")
}

#atn: added object ploidy
export<-c("phenotypes","Z","K","Hinv","genotypes","ploidy",
          "mm.solve","dosage.X","Q","C","test.compatibility","comp.vec") #gt

#extra functions need to be exproted
#if we don't want to use P3D approximation
if(!P3D|!EMMAX) export<-c(export,"calc.Hinv")

parallel::clusterExport(cl,export,
                        envir=environment()) #I dunno why, but without this it doesnt work
# #gt: same reason explained before. Here, doesn't even work because Z, K and Hinv exist only
# in this executing environment (and not in the global one).

# #gt version: markers in parallel (parallelization works even with one phenotype)
result<-lapply(1:ncol(phenotypes),function(w){#NOT parallely over each phenotype
  result<-parallel::parSapply(cl,1:nrow(genotypes),FUN=function(k){#parallely over markers, calculate the pvalue for each marker

    #DEFINITION OF X (matrix of fixed effects)
    g<-genotypes[k,]

    #Tweak



    X<-dosage.X(as.matrix(g),ploidy=ploidy,normalize=T)
    # if(any(is.na(genotypes[k,]))) X<-X[,-ncol(X)] #Method 2
    
    if(ncol(X)>2){X<-X[,-2]} #when ancestral/parental model we need to prevent singularity
    nparX<-ncol(X) #total number of genetic parameters
    X<-cbind(Q,C,X) #add population and cofactor parameters
    X<-matrix(as.numeric(X),ncol=ncol(X))
    no.test<-ncol(X)-nparX #number of non genetic parameters

    #If not P3D/EMMAX, then Hinv will be NULL (and cannot be subindexed)
    if(!is.null(Hinv)){
      H<-list(Hinv[[w]][missD,missD]) #gt: Hinv[[w]] will be a matrix, but mm.solve is expecting a list!
      # Hinv[w] will be a list of length 1
    }else{
      H<-NULL
    }
    #gt: phenotypes[,w,drop=F] to avoid conversion to vector
    #gt: added [[1]], after replacement of sapply with lapply in mm.solve
    return(mm.solve(phenotypes[,w,drop=F],X,Z,K,H,no.test = no.test))
  })

  #NEW OUTPUT
  #Create a list of lists (in df form) from all results.
  #It's actually just 6 lists with k elements where k is markers
  res<-as.data.frame(do.call(mapply,c(list,result)))
  #All columns are turned into vectors except the beta column
  for(r in 2:ncol(res)) res[,r]<-do.call(c,res[,r])
  #dimnames[1] is set as markers. Works a bit weird.
  rownames(res)<-markers
  return(res)
})
#Here we obtain a list of df, each df containing all results
parallel::stopCluster(cl)

#The results are given phenotype names
#We obtain the structure we talked about
names(result)<-colnames(phenotypes)

saveRDS(result,"Result_nona.RDS")
saveRDS(result,"Result_na_as_factor.RDS")
saveRDS(result,"Result_na_as_nofactor.RDS")
saveRDS(result,"Result_na_as_omit.RDS")
saveRDS(result,"Result_na_imputed.RDS")

resNA0<-readRDS("Result_nona.RDS")
resNA1<-readRDS("Result_na_as_factor.RDS")
resNA2<-readRDS("Result_na_as_nofactor.RDS")
resNA3<-readRDS("Result_na_as_omit.RDS")
resNA4<-readRDS("Result_na_imputed.RDS")

na0<-resNA0$h08$pval
na1<-resNA1$h08$pval
na2<-resNA2$h08$pval
na3<-resNA3$h08$pval
na4<-resNA4$h08$pval
pvals<-cbind(na0,na1,na2,na3,na4)

plot(-log10(na0))
plot(-log10(na1))
plot(-log10(na2))
plot(-log10(na3))
plot(-log10(na4))

png("Plots/NA_treatment_mapQTL.png",res=150,height = 1000,width=1400)
comp.QQplot(pvals,
            main="NA treatment comparison")
dev.off()

plot(-log10(na0),-log10(na3))
abline(a=0,b=1,col="red")
plot(-log10(na1),-log10(na2))
abline(a=0,b=1,col="red")
plot(-log10(na1),-log10(na3))
abline(a=0,b=1,col="red")
plot(-log10(na2),-log10(na3))
abline(a=0,b=1,col="red")
plot(-log10(na0),-log10(na4))
abline(a=0,b=1,col="red")
