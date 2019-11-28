# QTL mapping functions ---------------------

#' This file includes functions for mapping QTLs using GWAS association model.
#' The association model has been adapted to use in different ploidies, and can
#' use SNP-dosage matrix or haplotype matrices to perform the GWAS analysis.
#' The GWAS model is implemented both as a linear and mixed model, with the
#' approximation known as P3D/EMMAX, which greatly reduces time usage.
#' Additionally, the program can easily handle cofactors and various
#' parameter specifications. Inspirations for this program include:
#'   - Unified mixed model (Yu et al 2006)
#'   - GWASpoly (Rosyara et al 2012)
#'   - mppR (Vincent Garin 2019)
#'
#' Developed by Alejandro Thérèse Navarro & Giorgio Tumino
#' October 2019

#Main wrapper -------------------------


#' QTL mapping of a matrix of phenotypes
#'
#' @param phenotypes A numeric matrix of phenotypes,
#' rows are individuals and columns are different phenotypes.
#' Not updated for vector phenotypes yet. Must follow same
#' order of individuals as genotypes.
#' @param genotypes A matrix of genotypes, rows are markers and
#' columns may be either (1) individual dosages of biallelic markers, with column
#' names coinciding with phenotype individual names or (2) chromosome alleles of each
#' individual (#gt: for multiallelic markers, such as haplotypes).
#' Must follow same order of individuals as phenotypes.
#' @param ploidy a number indicating the ploidy level. All the individuals
#' must have same ploidy.
#' @param K NULL, T or distance matrix. If NULL, no relatedness
#' matrix will be used (i.e., a linear model will be applied). If T
#' a K distance matrix will be calculated. A distance matrix may also
#' be directly specified.
#' @param Q NULL, T or vector identifying populations. If NULL, no Q will
#' be included in the model (i.e., a model without Q correction). If T, a pco
#' decomposition will be used to estimate population differentiation. If a
#' vector specifying population of each individual is passed,
#' it will be used to construct a Q matrix. Vector may contain numerical
#' or character.
#' @param map A table with a genetic map containing at least a "chromosome"
#' and a "position" columns, specifying chromosome number and cM position.
#' Used to sample marker dosages every 1 cM.
#' AN ADRESS of a Pedsim map file: at least with a "chromosome"
#' and "position" columns, specifying chromosome number and cM position
#' of the marker at that chromosome
#' @param cM Numeric. K distance matrix (#gt: and Q matrix) will be calculated using markers every
#' cM centimorgans. Defaults to 1.
#' @param seed An integer to set a seed for random number generation in \code{sample.cM}.
#' @param no_cores Numeric. Number of cores to be used in parallel computing
#' of the p-values. Defaults to number of cores -1.
#' @param Z Identity matrix indicating which individuals correspond to which
#' genotypic effects. For instance, if multiple samples correspond to the same
#' individual, this matrix should indicate so.
#' @param cofactor Possible cofactor matrix (where each column is a cofactor).
#' @param cofactor.type If a cofactor matrix is specified, the type of cofactor
#' ("numerical" or "categorical") to be used for each. It accepts partial strings
#' such as "cat" and "num".
#' @param Qpco Logical value indicating whether a PCo-based population factor should be
#' estimated and included. Is useful for detecting and correcting for population
#' structure.
#' @param approximate Logical value indicating whether P3D/EMMAX approach should be
#' used for computational efficiency (if F, expect much longer waiting times)
#' @param permutation permutation strategy. Can be "pop" (permutation over the
#' whole population) or "fam" (permutation within families). If it is NULL, no
#' permutation will be run.
#' @param alpha Permutation threshold alpha value, defaults to 0.05. Alpha is the
#' chance of observing a false positive experiment-wise.
#' @param impute Logical value indicating whether missing genotypes should be
#' imputed using \code{impute.knn}. Defaults to True, but False is recommended (imputation algorithm
#' needs to be improved and with few missing values it has a small effect
#' on QTL detection)
#' @param nperm number of permutations, if \code{permutation} is not null.
#' @param k Number of neighbours to use in \code{\link{impute.knn}}
#' @param linear logical. If True, linear model (without structure correction)
#' is applied. If False, mixed model (with structure correction) is applied. If
#' not specified, it looks at the value of K, and only applied linear model
#' if K = NULL.
#' @param K_identity logical. If TRUE, an identity matrix is used in the
#' random term of a mixed model. This is to run a mixed model with no kinship
#' correction, even if a kinship matrix is provided or calculated for other
#' purposes (e.g. imputation of missing genotypes, calculation of the Q matrix)
#'
#' @return a pvalue matrix containing the pvalue of each marker with each phenotype passed.
map.QTL<-function(
  phenotypes,
  genotypes, #genotype matrix
  ploidy,
  map, #genetic map table
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
  nperm = 1000, #number of permutations
  alpha = 0.95,
  impute=T,
  k=20,
  linear = NULL,
  K_identity = F
){

  if(is.null(linear)){
    if(is.null(K)){
      linear <- T
    }else{
      linear <- F
    }
  }

  markers <- map$marker

  std_phe <- stndrdz.pheno(phenotypes)
  phenotypes <- std_phe$pheno

  # Curation of genotype matrix and dosage matrix
  ## check whether genotypes contain snp dosages or haplotypes
  if (ncol(genotypes)==nrow(phenotypes)) {
    haplo<-F
    genotypes <- inputCheck_dos(genotypes, integer=T, ploidy=ploidy)
    cat("SNP dosages have been detected in the genotype matrix,\n")

    # if (!all(rownames(phenotypes)==colnames(genotypes))) {
    #   stop("phenotype individuals and genotype individuals have different names")
    # }

  } else if ((ncol(genotypes) %% ploidy) == 0){
    haplo<-T
    cat("Haplotypes have been detected in the genotype matrix.\n")

  } else {
    stop("If haplotypes have been provided: the number of genotypes columns
         is not a multiple of ploidy.
         If SNP dosages have been provided: the number of individual dosages
         is not equal to the number of individual phenotypes.")
  }

  #bi.geno
  #MAF

  #DEFINITION OF K -------------------
  #There are four options
  #1) K=NULL, no K is used. We go linear.
  #2) K=T, K is calculated using homogeneously distributed markers along a map.
  #   using sample.cM and calc.K
  #3) K=K matrix, check dimensions, message is printed and K is used.


  if(!is.null(K)){ #gt: added because when K=NULL the condition below is TRUE
    if(all(K==T)){ #gt: why are you using all?
      #Option 2)
      K<-sample.cM(genotypes,map,cM = cM)
      K<-calc.K(t(K),ploidy = ploidy,haplotypes = haplo)

    }else if(nrow(K)==nrow(phenotypes)){
      #Option 3)
      if(ncol(K)!=nrow(phenotypes)) stop("K matrix was specified, but is not square")
      cat("Square K matrix has been provided. Using as it is.\n")

    }
  }


  ### Na check and impute ----------------
  naperc <- sum(is.na(genotypes))/length(genotypes)
  cat(round(naperc*100,2),"% missing genotypes detected.\n")
  if(naperc == 0) impute<-F
  if(impute){
    if(is.null(K)){
      #first we sample homogeneously markers along the genome
      K <- sample.cM(genotypes,map,cM = cM)
      #then we calculate the distance
      K <- calc.K(t(K),ploidy=ploidy,haplotypes = haplo)
    }
    genotypes<-impute.knn(genotypes,ploidy,map,kneighbors=k,K=K)
    cat("Imputation performed.\n")
  }else{
    cat("No imputation will be performed.\n")
  }

  #DEFINITION OF Q ------------------------
  #4 options:
  #1) No Q is used (Q=NULL), this block is skipped
  #2) Q=T, Q is calculated using the K matrix and cmdscale
  #3) Q=vector identifying groups, Q is calculated using Q.mat
  #4) Q=Q design matrix, message is printed.
  if(is.null(Q)) cat("No Q matrix will be used.\n")

  if(!is.null(Q)){
    if(all(Q==T)){ #gt: why are you using all?
      #atn: because if Q is a vector/matrix Q==T will be a vector of F
      #and that gives an error in the if.

      #Option 2)
      if(is.null(K)){
        #first we sample homogeneously markers along the genome
        K<-sample.cM(genotypes,map,cM = cM)
        #then we calculate the distance
        K<-calc.K(t(K),ploidy=ploidy,haplotypes = haplo)
      }
      Q <- cmdscale(1-K, k=Qpco, eig = F, add = FALSE, x.ret = FALSE)

    }else if(is.vector(Q)){
      #Option 3)
      Q<-calc.Q(Q)

    }else if(nrow(Q)==nrow(phenotypes)){
      #Option 4)
      cat("Q matrix has been provided. Using as it is.\n")

    }else{
      #Q was not NULL, T, a vector or a cofactor design matrix.
      stop("Wrong Q specification.")}
  }

  ##DEFINITION OF C, cofactor matrix. ---------------------------
  #Two variables must be provided, cofactor and cofactor.type
  #The first contains values, the second specifies whether the cofactor
  #should be treated as a numeric regressor or as a categorical (ANOVA-type) parameter
  if (!is.null(cofactor)){

    #Make cofactors into a matrix
    if(is.vector(cofactor)) cofactor <- matrix(cofactor,ncol=1)

    #Check cofactor tpyes are specified for cofactors.
    if(is.null(cofactor.type)|length(cofactor.type)!=ncol(cofactor)){
      stop(paste0("All cofactor types must be specified using the cofactor.type parameter.",
                  "\nA char vector containing \"numerical\"",
                  " or \"categorical\" for each cofactor must be given."))
    }

    #Which are num and which are cat
    c.num <- pmatch(cofactor.type,"numerical")
    c.cat <- pmatch(cofactor.type,"categorical")

    C <- lapply(1:ncol(cofactor),function(i){

      if(!is.na(c.num[i])){
        result <- as.numeric(cofactor[,i])

      }else if(!is.na(c.cat[i])){

        #This creates a design matrix with 1s and 0s (representing T and F)
        vals <- na.omit(unique(cofactor[,i]))
        result <- sapply(vals,function(c){
          cofactor[,i] == c
        })+1-1

        #Last column is taken out to be able to calculate
        result <- result[,-ncol(result)]
      }
    })
    C <- do.call(cbind,C)

  }else{C<-NULL}



  # Prepare permuted phenotypes ------------------------
  npheno <- ncol(phenotypes)
  if (!is.null(permutation)) {
    if (is.null(nperm)) stop("nperm is NULL. Provide the number of permutations.")
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
  # Linear Model ------------
  #with no structure correction or with Q correction
  if(linear){
    #First set up cluster

    cluster <- parallel::makeCluster(no_cores)
    export <- c("phenotypes","genotypes","dosage.X","Q","C","lm_compare","haplo") #gt
    if(is.null(Q)){
      # export<-c("phenotypes","genotypes","dosage.X")
      cat("Linear model will be used.\n")
    }else{
      # export<-c("phenotypes","genotypes","dosage.X","Q")
      cat("Linear model with Q correction will be used.\n")
    }

    parallel::clusterExport(cl=cluster,export,
                            envir=environment()) #needed again? #gt: yes, because
    # we want to use variables defined in this executing environment, not in the
    # global one. Although the function dosage.X is not present in this
    # environment, clusterExport will search it and will find it in the parent environments

    #Testing
    result <- lapply(1:ncol(phenotypes),function(w){
      res <- parallel::parLapply(cluster,1:nrow(genotypes),function(k){
        X <- dosage.X(genotypes[k,],ploidy = 4,haplotype = haplo,normalize = F)
        if(haplo) X <- X[,-ncol(X)]
        #We create the naive and full matrices. The naive
        #model only includes Q, C and an intercept
        N <- cbind(matrix(1,nrow=nrow(X)),Q,C) #naive model
        X <- cbind(1,Q,C,X) #total model
        y <- phenotypes[,w,drop=F]
        solveout <- try(lm_compare(X,N,y),silent = T)

        ## mantain the usual output structure in case of error
        if (class(solveout) == "try-error") {
          write(solveout, file="lm_compare_errors.txt", append = T)

          #ATN the "original structure" requires iterators with as many elements
          #as there are phenotypes. That is, NA matrices/ vectors with length =
          #ncol(phenotypes)
          miss <- rep(NA,ncol(y))
          solveout <- list(
            beta = matrix(miss,nrow=1),
            Fstat = miss,
            pval = miss,
            se = miss
          )
        }

        return(solveout)

      })

      #NEW OUTPUT
      #Create a list of lists from all results.
      #It's actually just 6 lists with k elements where k is markers
      res <- do.call(mapply,c(list,res))
      res <- as.list(as.data.frame(res))
      #All columns are turned into vectors except the beta column

      for(r in 2:length(res)) res[[r]] <- unlist(res[[r]])
      for(r in 1:length(res)) names(res[[r]]) <- markers

      ## for permuted phenotypes store the minimum pvalue only
      if (w > npheno) {
        res <- min(res$pval, na.rm = T)
      }

      return(res)
    })

    parallel::stopCluster(cluster)

    # #Reformatting of results so that they fit standard output
    # result <- lapply(1:ncol(phenotypes),function(w){
    #
    #   #Separates results from each phenotype
    #   res <- lapply(result,function(r) lapply(r,function(e){
    #     if(is.matrix(e)) return(e[,w])
    #     return(e[w])
    #   }))
    #
    #   res <- do.call(mapply,c(list,res))
    #   res <- as.list(as.data.frame(res))
    #   for(r in 2:length(res)) res[[r]] <- unlist(res[[r]]) #gt: convert to vectors Fstat, pval and se (not beta)
    #   for(r in 1:length(res)) names(res[[r]]) <- markers
    #
    #   ## for permuted phenotypes store the minimum pvalue only
    #   if (w > npheno) {
    #     res <- min(res$pval, na.rm = T)
    #   }
    #   return(res)
    # })


  }else{ #Ergo, K must be defined, we must use Mixed Models


    ####
    # Mixed Model --------------
    if (K_identity == T) {
      K <- matrix(1, nrow=nrow(phenotypes), ncol=nrow(phenotypes))
    }

    if(is.null(Z)){
      Z<-diag(nrow(phenotypes))
    }

    if(is.null(Q)){
      cat("Mixed model will be used.\n")
    }else{
      cat("Mixed model with Q correction will be used.\n")
    }

    if(!approximate){
      Hinv<-NULL
    }else{

      X <- rep(1,nrow(phenotypes))
      if(!is.null(C)) X <- cbind(C,Q,X)
      Hinv <- calc.Hinv(phenotypes,
                        #X=vector of 1s, to apply P3D/EMMAX algorithm #gt: check#
                        X = X,
                        Z,K)
    }

    cl<-parallel::makeCluster(no_cores)
    #atn: added object ploidy
    export<-c("phenotypes","Z","K","Hinv","genotypes","ploidy","haplo","npheno",
              "mm.solve","dosage.X","Q","C","test.compatibility","comp.vec") #gt

    #extra functions need to be exproted
    #if we don't want to use P3D approximation
    if(!approximate) export<-c(export,"calc.Hinv")

    parallel::clusterExport(cl,export,
                            envir=environment()) #I dunno why, but without this it doesnt work
    # #gt: same reason explained before. Here, doesn't even work because Z, K and Hinv exist only
    # in this executing environment (and not in the global one).

    # #gt version: markers in parallel (parallelization works even with one phenotype)
    result<-lapply(1:ncol(phenotypes),function(w){#NOT parallely over each phenotype
      result<-parallel::parSapply(cl,1:nrow(genotypes),FUN=function(k){#parallely over markers, calculate the pvalue for each marker

        #DEFINITION OF X (matrix of fixed effects)
        X <- dosage.X(as.matrix(genotypes[k,]),
                    ploidy=ploidy,
                    normalize=T,
                    haplotype = haplo)

        if(ncol(X)>1){X<-X[,-1,drop=F]} #when ancestral/parental model we need to prevent singularity
        # if(any(is.na(genotypes[k,]))) X<-X[,-ncol(X)] #Eliminate NA as factor
        nparX <- ncol(X) #total number of genetic parameters
        X <- cbind(Q,C,X) #add population and cofactor parameters
        X <- matrix(as.numeric(X),ncol=ncol(X))
        no.test<-ncol(X)-nparX #number of non genetic parameters

        #If not P3D/EMMAX, then Hinv will be NULL (and cannot be subindexed)
        if(!is.null(Hinv)){
          H<-Hinv[w] #gt: Hinv[[w]] will be a matrix, but mm.solve is expecting a list!
          # Hinv[w] will be a list of length 1
        }else{
          H<-NULL
        }
        #gt: phenotypes[,w,drop=F] to avoid conversion to vector
        #gt: added [[1]], after replacement of sapply with lapply in mm.solve
        #solveout <- mm.solve(phenotypes[,w,drop=F],X,Z,K,H,no.test = no.test)
        solveout <- try(mm.solve(phenotypes[,w,drop=F],X,Z,K,H,no.test = no.test), silent = T)

        ## mantain the usual output structure in case of error
        if (class(solveout) == "try-error") {
          write(solveout, file="mm.solve_errors.txt", append = T)
          solveout <- list(list(
            beta = matrix(NA),
            Fstat = NA,
            residual = rep(NA,nrow(X)),
            pval = NA,
            se = NA,
            wald = NA,
            real.df = NA
          ))
        }

        return(solveout)
      })

      #NEW OUTPUT
      #Create a list of lists from all results.
      #It's actually just 6 lists with k elements where k is markers
      res<-do.call(mapply,c(list,result))
      res<-as.list(as.data.frame(res))
      #All columns are turned into vectors except the beta column

      for(r in 2:length(res)) res[[r]] <- unlist(res[[r]])
      for(r in 1:length(res)) names(res[[r]]) <- markers
      res$residual <- split(res$residual,rep(1:nrow(genotypes),
                                             each=nrow(phenotypes)))
      names(res$residual) <- markers
      ## for permuted phenotypes store the minimum pvalue only
      if (w > npheno) {
        res <- min(res$pval, na.rm = T)
      }

      return(res)
    })
    #Here we obtain a list of df, each df containing all results
    parallel::stopCluster(cl)
  }

  #The results are given phenotype names
  #We obtain the structure we talked about
  if(!is.null(colnames(phenotypes))){
    names(result) <- make.names(colnames(phenotypes),unique=T)
  }else{
    names(result) <- paste0("pheno",1:ncol(phenotypes))
  }


  cat("Association completed.\n")

  if (!is.null(permutation)) {
    cat("Permutation threshold will now be calculated.\n")
    minpval <- sapply(result[(npheno+1):length(result)], function(x) x)
    minpvalMat <- matrix(minpval, ncol = npheno, byrow = T)
    maxlogpvalMat <- -log10(minpvalMat)
    thr <- sapply(1:ncol(maxlogpvalMat), function(i) {
      quantile(maxlogpvalMat[,i], probs = alpha) # alpha can be a vector
    })
    if (length(alpha)==1) thr <- t(thr)
    colnames(thr) <- colnames(phenotypes)[1:npheno]
    rownames(thr) <- alpha

    #gt: a new element called 'perm.thr' is added to each phenotype in result
    phenonames <- names(result)
    result <- lapply(1:npheno, function(i) {
      result[[i]]$perm.thr <- thr[,i]
      return(result[[i]])
    })
    names(result) <- phenonames[1:npheno]
  }

  #To rescale the effects
  for(i in seq_along(result)){
    beta <- result[[i]]$beta
    if(linear){
      new_beta <- lapply(beta,function(b) {
        b[1] <- b[1] + std_phe$mean[i]
        b[-1] <- b[-1]* std_phe$sd[i]
        return(b)
      })
    }else{
      new_beta <- lapply(beta,function(b) {
        b* std_phe$sd[i]
      })
    }

    result[[i]]$beta <- new_beta
  }

  return(result)
}


# Calculation functions -------------------

#' Calculation of realized distance matrix (K)
#'
#' @description Using dosage scores, a distance matrix is calculated
#' such that the average distance of an individual with itself is 1, and
#' the average with an unrelated individual is 0. Based on the "Realized
#' Relationship" model found in \href{https://dl.sciencesocieties.org/publications/tpg/abstracts/9/2/plantgenome2015.08.0073}
#' {Rosyara et al. 2016}
#'
#' @param matrix Numeric matrix with individuals on rows and markers on columns.
#'
#' @return A numeric matrix nxn where n is the number of rows.
#' @export
#'
#' @examples
#' #Create 10 tetraploid individuals with 200 markers each
#' inds <- lapply(1:10,function(i) round(runif(200,0,4)) )
#'
#' #Put them in a matrix or data.frame with individuals in rows
#' geno <- do.call(rbind,inds)
#'
#' K<-calc.K(geno)
#'
calc.K<-function(
  matrix,
  haplotypes=F,
  ploidy=NULL
){
  #When haplotypes are given we can still perform K!
  if(haplotypes){
    if(is.null(ploidy))
      stop("For haplotype distance calculation ploidy must be defined.")
    #Create an ANOVA type matrix for all haplotypes
    matrix<-lapply(1:ncol(matrix),function(i){
      dosage.X(matrix[,i],haplotype = T,ploidy=ploidy,normalize = F)
    })
    matrix<-do.call(cbind,matrix)
  } else {
    #impute NAs. Let's leave this here, useful for our imputator
    matrix <- imputeNA(matrix)
  }

  #Substract mean dosage for each marker
  M<-apply(matrix,2,function(x) x-mean(x))
  K<-M%*%t(M) #Calculate distance matrix
  K<-K/mean(diag(K)) #average the center
  colnames(K)<-rownames(matrix)
  rownames(K)<-rownames(matrix)
  return(K)
}

#' Calculate a Q design matrix based on a vector
#'
#' @description Using a vector of values as input, it creates a cofactor
#' design matrix (\eqn{Q}) that identifies each individual as belonging to a group.
#'
#' @param pop Vector where each element is a population identifier
#' @param names Optional. A vector of names for the \eqn{Q} matrix
#'
#' @return A matrix with as many columns as population groups
#' [ncol=unique(pop)] identifying each individual belonging to
#' one population.
#' @export
#'
#' @examples
calc.Q<-function(
  pop, #vector identifying a population
  names=NULL #optionally, a vector defining the names
){
  match <- 1*sapply(unique(na.omit(pop)),function(x) x==pop)

  if(!is.null(names)){rownames(match)<-names}

  match<-apply(match,2,function(x) (x-mean(x))/sd(x))

  return(match[,-1]) #we take out the first column, to avoid singularity
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

#' Calculates \eqn{H^-1} matrix
#'
#' @description The core of a Ridge-Regression solution for a mixed model
#' is the calculation of a Hat matrix, or projection matrix \eqn{H}, which contains the
#' variance components of the random part of the model. In the EMMAX/P3D approach,
#' the inverse of such matrix, \eqn{H^-1}, is calculated only once and recycled for the analysis
#' of each marker, which obtains an approximate result but speeds up computational time
#' considerably. For more information consult the articles that simulatenously presented this method:
#' \href{https://www.nature.com/articles/ng.546}{P3D: Zhang et al. 2010}
#' \href{https://www.nature.com/articles/ng.548}{EMMAX: Kang et al. 2010}
#'
#' @param y Response vector or matrix. Each column is taken as a different response.
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix
#' @param K Optional (defaults to NULL).
#' Variance-covariance square matrix (Genetic distance matrix).
#' @param bounds Bounds of the Ridge Regression parameter
#' @param method "REML" or "ML". If "ML", \eqn{H^-1} is calculated using maximum likelihood
#' equations. If REML, Restricted Maximum Likelihood.
#'
#' @return Returns a list \eqn{H^-1},
#' where each element corresponds to the \eqn{H^-1} of each column of y. An
#' attribute "lambda" is included in each \eqn{H^-1} matrix that contains the optimization
#' paramater of ridge regression.
#' @export
#'
#' @examples
calc.Hinv<-function(
  y,
  X,
  Z,
  K=NULL,
  bounds=c(1e-09,1e+09),
  method="REML"
){

  #Obtain the results without NA's at y or X
  #K should probably be recalculated for each y,
  #but we just approximate it. (Might not work with lots of NA's)
  data <- test.compatibility(y,X,Z,K)

  Hinv<-lapply(1:length(data),function(i){
    dat<-data[[i]]
    y <- dat$y
    X <- dat$X
    K <- dat$K
    Z <- dat$Z

    n <- nrow(y) #number of observations
    p <- ncol(X) #number of fixed parameters
    m <- ncol(Z) #number of random variance components

    Xinv <- solve(crossprod(X))
    S <- diag(n) - tcrossprod(X %*% Xinv, X)

    #In case n<(m+p) the following part is done with eigenvalue decomposition.
    #Otherwise Cholesky decomposition is usable. Because of our dataset,
    #only the first method is implemented. Both require positive semidefinite
    #matrices (all eigenvalues are >=0).
    #This part of the code is copied from rrBLUP, so I haven't touched it almost
    offset <- sqrt(n)
    if (is.null(K)) {
      Hb <- tcrossprod(Z) + offset * diag(n)
    } else{
      Hb <- tcrossprod(Z %*% K, Z) + offset * diag(n)
    }

    Hb.system <- eigen(Hb, symmetric = TRUE)
    phi <- Hb.system$values - offset
    if (min(phi) < -1e-06) {
      stop("K not positive semi-definite.")
    }
    U <- Hb.system$vectors
    SHbS <- S %*% Hb %*% S
    SHbS.system <- eigen(SHbS, symmetric = TRUE)
    theta <- SHbS.system$values[1:(n - p)] - offset
    Q <- SHbS.system$vectors[, 1:(n - p)]

    #If we provide multiple columns of y we will calculate one Hinv
    #for each y. Otherwise we just obtain one y.
    #This list set-up allows us to recicle Q calculation, which
    #can decrease considerably the time when dealing
    #with multiple phenotypes.
    Hinv <- lapply(1:ncol(y),function(i){

      suby<-y[,i]
      omega.sq<-crossprod(Q,suby)^2

      if (method == "ML") {
        f.ML <- function(lambda, n, theta, omega.sq, phi) {
          n * log(sum(omega.sq/(theta + lambda))) + sum(log(phi+lambda))
        }
        soln <- optimize(f.ML, interval = bounds, n, theta, omega.sq, phi)
        lambda.opt <- soln$minimum
        df <- n
      } else {
        f.REML <- function(lambda, n.p, theta, omega.sq) {
          n.p * log(sum(omega.sq/(theta + lambda))) + sum(log(theta + lambda))
        }
        soln <- optimize(f.REML, interval = bounds, n - p, theta, omega.sq)
        lambda.opt <- soln$minimum
        df <- n - p
      }

      Hinv <- U %*% (t(U)/(phi + lambda.opt))

      attributes(Hinv)$lambda<-lambda.opt
      return(Hinv)

    })

    return(Hinv)
  })

  #This is to prevent a nested list when calculating Hinv.
  Hinv <- unlist(Hinv,recursive = F)
  Hinv <- Hinv[original_order(y,X)]
  return(Hinv)
}

#' Mixed Model Solver
#' @description Mixed model solver based on the rrBLUP package mixed.solve function
#'
#' @param y Numeric vector of response variable
#' @param X Fixed effect design matrix. Must include an intercept. Generally this
#' matrix corresponds to 1) cofactors and 2) marker(s) genotype(s).
#' @param Z Random effects design matrix. Generally this matrix relates observations with
#' each genotype/individual.
#' @param K Variance-covariance structure. If not specified, \code{\link{calc.K}} is used
#' @param Hinv \eqn{H^-1} matrix. If not provided \code{\link{calc.K}} will be used. Defaults
#' to NULL.
#' @param random Logical value indicating whether random effects should be returned (slows
#' down function significantly).
#' @param no.test Number of parameters in the X matrix that will not be taken into
#' account in the pvalue calculation. In practice, this means that the first \code{no.test}
#' columns of the X matrix will not be used to calculate p-value. For instance, if 3
#' cofactors are added, the X matrix should have 5 columns: one intercept, three cofactors
#' and one genotype column. If we wish to investigate the effect of the genotype
#' only, no.test should take value 4. If we would set \code{no.test} to 1, the p-value
#' would reflect whether the cofactors have any influence on the response \eqn{y}.
#'
#' @return a data.frame where each column corresponds to a column of \eqn{y}. There are six
#' rows on the df: beta (estimates of parameters), Fstat (approximation of
#' F value), residual (vector of residuals), pval (p-value), se (standard error of the model),
#' wald (Wald test values), and real.df (realized degrees of freedom, df approximation).
#' If \code{random} is set to \code{True} an extra row "random" is added, containing the
#' random terms -at the cost of computational efficiency.
#' @export
#'
#' @examples
mm.solve <- function(
  y,
  X,
  Z,
  K,
  Hinv = NULL,
  random = F,
  no.test = 0
){

  #Obtain the results without NA's at y or X
  #K should probably be recalculated for each y,
  #but we just approximate it. (Might not work with lots of NA's)
  data <- test.compatibility(y,X,Z,K,Hinv)


  result<-lapply(1:length(data),function(i){
    dat<-data[[i]]
    y <- dat$y
    X <- dat$X
    K <- dat$K
    Z <- dat$Z
    if(!is.null(Hinv)) Hinv <- dat$Hinv

    n <- nrow(y)
    #number of genetic parameters (not counting intercept, cofactors, etc.)
    p<-ncol(X)-no.test
    if(p==0) stop("no.test is too low, no parameters will be tested")
    df<-n-ncol(X)

    if(no.test==0){
      X2<-Z%*%X
    }else{
      X2<-cbind(X[,1:no.test],Z%*%X[,-1:-no.test,drop=F])#do not multiply cofactors
    }

    #Hat inverse matrix calculation, i.e. rrBLUP estimation.
    if(all(is.null(Hinv)))  Hinv<-calc.Hinv(y=y,X=X,Z=Z,K=K)

    #If multiple phenotypes are provided, Hinv should be
    #a list with one matrix for each phenotype (for each column of y)
    #gt: I commented the check below because we provide mm.solve with one
    #phenotype at time with its Hinv matrix
    if(!length(Hinv)==ncol(y)) stop(
      paste("More than one column for y was provided.",
            "Hinv must be a list of matrices, one for each column of y.",
            "Try setting Hinv to NULL (?)"))


    result<-lapply(1:length(Hinv),function(i){ #gt: sapply -> lapply
      Hi<-Hinv[[i]]
      yi<-y[,i]
      lambda<-attr(Hi,"lambda")

      W <- crossprod(X2, Hi %*% X2)
      Winv <- solve(W) #gt singularity tolerance?
      beta <- Winv %*% crossprod(X2, Hi %*% yi) #estimates
      resid <- yi - X2 %*% beta #residual
      s2 <- as.double(crossprod(resid, Hi %*% resid))/df #variance

      #We only take the Winv col and row of the parameters we want to test
      #(not intercept, so ncol(X)-1)
      Q <- s2*Winv[(no.test+1):ncol(X),(no.test+1):ncol(X)]
      Tt <- solve(Q) #Estimates/se

      #Wald test (??)
      ypred <- X%*%beta
      Vcov <- solve(diag(length(yi))*s2)
      wald <- t(yi)%*%Vcov%*%ypred

      #Realized degrees of freedom
      dj2 <- eigen(t(X)%*%X)$values
      p2 <- sum(dj2/(dj2+lambda))

      #F-statistic of the genetic model?
      V <- beta[(no.test +1):ncol(X),drop=F]
      Fstat <- crossprod(V,Tt %*% V)/(p)

      x <- df/(df + p * Fstat)
      pval <- pbeta(x, df/2, p/2)

      # df<-df^2/(df^2+lambda)
      # pval<-pf(Fstat,p2,df,lower.tail=F) #pvalue of genetic model?

      result<-list(
        beta = beta,
        Fstat = Fstat,
        residual = resid,
        pval = pval,
        se = sqrt(s2),
        wald = wald,
        real.df = p2
      )

      if (random) {
        KZt <- tcrossprod(K, Z)
        u <- crossprod(KZt %*% Hinv, (y - X %*% beta))#random effects
        result$random<-u
      }

      return(result)
    })

    return(result)
  })

  result<-unlist(result,recursive = F)
  if(!is.null(colnames(y))){
    names(result)<-make.unique(colnames(y))
  }
  return(result)
}

#' Linear model comparison
#'
#' @description This function allows to perform a comparison test
#' between two linear models. This is useful to test groups of fixed
#' effects. In our program, this is applied to test the "genetic effects"
#' without testing the cofactors. For instance, if we have a model such as:
#' y = c1 + c2 + g1 + g2 + e, where c1 and c2 are cofactors, and g1 and
#' g2 are genetic effects, we would like to know what is the explanatory
#' value of only g1 and g2. To do so, we compare the model above with
#' the following model: y = c1 + c2 + e. If g1 and g2 add explanatory value,
#' then the first model should have higher explanatory value. In this
#' function, the first model, or full model, is represented by the
#' matrix X, while the second model, or naive model, is represented by
#' the matrix N.
#'
#' @param X matrix of fixed effects in the full model
#' @param N matrix of fixed effects in the naive model
#' @param y matrix or vector of phenotypes
#'
#' @return list with beta (fixed effects), Ftest (F statistic)
#' pval (p-value of the difference between models), se (standard error
#' of fixed effect estimates)
#' @export
#'
#' @examples
lm_compare <- function(X,N,y){

  if(is.vector(y)) y <- matrix(y,ncol=1)

  #Process followed by lm ----
  not_na_y <- rowSums(is.na(y)) == 0
  not_na_X <- rowSums(is.na(X)) == 0
  not_na <- not_na_y & not_na_X
  X <- X[not_na,,drop=F]
  N <- N[not_na,,drop=F]
  y <- y[not_na,,drop=F]
  y_mean <- colMeans(y)

  #Calculations naive ----
  betaN <- solve(t(N)%*%N)%*%t(N)%*%y
  fittedN <- N%*%betaN
  residualsN <- y - fittedN
  SSR_naive <- colSums(residualsN^2)
  df_naive <- nrow(N)-ncol(N)

  #Calculations total ----
  beta <- solve(t(X)%*%X)%*%t(X)%*%y
  fitted <- X%*%beta
  residuals <- y - fitted
  SSR_full <- colSums(residuals^2)
  df_full <- nrow(X) - ncol(X)

  #Last variable calculation ----
  dif_SSR <- SSR_naive - SSR_full
  dif_df <- df_naive - df_full
  MSR_full <- SSR_full/df_full

  Ftest <- (dif_SSR/dif_df)/MSR_full
  pval <- pf(Ftest,dif_df,df_full,lower.tail = F)
  se <- sqrt(SSR_full/(df_full))

  list(beta = beta, Ftest = Ftest,pval = pval,se = se)
}

#Threshold calculations -------------------
## used in LiJi
#' Pooled correlation between two multiallelic markers
#'
#' @description \code{rpool} computes a pooled correlation coefficient between
#' two sets of variables (e.g. two multiallelic markers) using the equation in
#' Zhao et al. 2005, Genet. Res. 86:77-87.
#' Correlation of each pair of alleles is weighted by allele frequency.
#'
#' @param A matrix of haplotype dosages of block A. Haplotypes in rows and
#' individuals in columns.
#' @param B matrix of haplotype dosages of block B.
#' @param ploidy numeric indicating the ploidy level.
#'
#' @return pooled correlation
#' @export
#'
#' @examples
rpool <- function(A,B,ploidy) {

  A <- t(A)
  B <- t(B)

  nhom <- nrow(A)*ploidy  #n homologues
  pA <- colSums(A)/nhom   #haplotype frequency
  pB <- colSums(B)/nhom

  ## for loop semi-vectorised
  k <- ncol(A)  #n haplotypes
  Srsq <- NULL
  for (i in 1:k) {
    Srsq <- sum(Srsq, pA[i] * pB * (cor(A[,i], B))^2)
  }

  ## alternative sapply loop (slower)
  # combid <- expand.grid(1:ncol(A),1:ncol(B))
  # Srsq <- sum(sapply(1:nrow(combid), function(x){
  #   pA[combid[x,1]] * pB[combid[x,2]] * (cor(A[,combid[x,1]], B[,combid[x,2]]))^2
  # }))

  return(sqrt(Srsq))
}

## used in LiJi
#' Pairwise pooled correlation between multiallelic markers
#'
#' @description \code{pairwise_rpool} generates a pairwise correlation matrix
#' between a set of multiallelic markers. It is based on \code{rpool}.
#'
#' @param hapdos a list of matrices with allele dosages. Each matrix contains
#' allele dosages of one multiallelic marker, with alleles in rows and
#' individuals in columns.
#' @param ploidy numeric indicating the ploidy level.
#'
#' @return a correlation matrix between multiallelic markers
#' @export
#'
#' @examples
pairwise_rpool <- function(hapdos,
                           ploidy) {

  nblk <- length(hapdos)
  prws <- combn(nblk,2)

  rtr <- sapply(1:ncol(prws), function(x) {
    rpool(A = hapdos[[prws[1,x]]], B = hapdos[[prws[2,x]]], ploidy=ploidy)
  })

  mcor <- matrix(1, ncol = nblk, nrow = nblk)
  mcor[lower.tri(mcor)] <- rtr
  mcor <- t(mcor)
  mcor[lower.tri(mcor)] <- rtr
  colnames(mcor) <- rownames(mcor) <- names(hapdos)
  return(mcor)
}


#' Li and Ji method for significance level in multiple testing
#'
#' @description This function implements a method to adjust the significance
#' level for multiple testing, as proposed by Li and Ji (Li and Ji, 2005, Heredity).
#' The Bonferroni adjusted threshold is obtained dividing a significance level
#' by the number of tests, assuming tests are independent to each other.
#' The Li and Ji method is more appropriate when the assumption of
#' independence is not met. In this method, the eigenvalues of a correlation
#' matrix are used to estimate the effective number of independent tests.
#'
#' @param m a matrix of SNP dosages or a list of matrices with haplotype dosages,
#' with markers (or haplotypes) in rows and individuals in columns.
#' No missing values are allowed. No monomorphic markers are allowed.
#' @param chrom a vector defining to which chromosome markers belong.
#' @param alpha significance level (use 0.05 for 5% significance).
#' @param ploidy numeric indicating the ploidy level.
#'
#' @return value of adjusted alpha.
#' @export
thr.LiJi <- function(m,
                     chrom,
                     alpha,
                     ploidy) {

  # check input data format
  if (is.matrix(m) || inherits(m, "data.frame")) {
    m <- inputCheck_dos(m, integer=F, ploidy=NULL)
  }

  ### data preparation
  # check missing values
  if(is.matrix(m) || inherits(m, "data.frame")) { #for SNP dosages
    NAsnp <- apply(m, 1, function(x) sum(is.na(x))/length(x))
    if (any(NAsnp > 0.25)) {
      stop(paste("There are markers with more than 25% of missing values",
                 "Remove them before proceeding."))
    } else if (any(NAsnp > 0)) {
      message(paste("The level of missingness per marker is acceptable (<25%).\n",
                    "Before calculating correlations, missing values will be",
                    "imputed using the average dosage per marker."))
      m <- imputeNA(m=m)
    }
  } else if(inherits(m, "list")) { #for haplotype dosages
    NAblock <- sapply(m, function(x) sum(is.na(x[1,]))/ncol(x))
    if (any(NAblock > 0.25)) {
      stop(paste("There are blocks with more than 25% of missing values",
                 "Remove them before proceeding."))
    } else if (any(NAblock > 0)) {
      message(paste("The level of missingness per block is acceptable (<25%).\n",
                    "Before calculating correlations, missing values will be",
                    "imputed using the average dosage per haplotype."))
      m <- lapply(m, imputeNA)
    }
  }

  # check monomorphic markers
  if (inherits(m, "list")) {
    mono <- sapply(m, function(n) {
      apply(n,1, function(x) sd(x, na.rm = T)==0)
    })
    polym <- lapply(mono, function(x) which(!x))
    m <- lapply(1:length(m), function(x) m[[x]][polym[[x]],])
    names(m) <- names(polym)
  } else if (is.matrix(m) || inherits(m, "data.frame")) {
    mono <- apply(m,1, function(x) sd(x, na.rm = T)==0)
    polym <- which(!mono)
    m <- m[polym,]  # discard monomorphic markers
  }
  if (any(unlist(mono))) warning(paste("There are",sum(unlist(mono)),
                                       "monomorphic markers.\n",
                                       "Remove them for subsequent analyses."))


  # define chrom
  # If missing, markers are considered as they were from the same chromosome
  if(missing(chrom)) {
    chrom <- rep(1, length(polym))
  } else if(is.matrix(m) || inherits(m, "data.frame")) {
    chrom <- chrom[polym]
  } # for haplotype lists it is not necessary, because we didn't remove
  # any block, but only monomorphic haplotypes within blocks.


  ####  Loop over chromosomes
  ## If chromosome membership of markers is unknown, a pairwise correlation
  ## matrix is calculated using all the markers, as if they were in the same
  ## chromosome.
  ## If chromosome membership of markers is known, we can assume no linkage
  ## (correlation) between markers of different chromosomes. Therefore, we can
  ## estimate the number of independent tests (Meff) for each chromosome and
  ## sum them to calculate the genome-wide Meff.
  ## TODO: if chromosome membership is missing for a subset of markers only?
  ## For the moment, markers with missing membership could be treated as a
  ## separate chromosome, though this approach could overestimate Meff (those
  ## markers could be in physical linkage with markers with membership).

  chrom[is.na(chrom)] <- "unknown"
  chromnames <- unique(chrom)
  Meff <- NULL
  for(i in chromnames) {
    ## step 1: correlation matrix between tests (markers)
    if(is.matrix(m) || inherits(m, "data.frame")) {
      cm <- cor(t(m[which(chrom == i),]))
    } else if(inherits(m, "list")) {
      cm <- pairwise_rpool(m[which(chrom == i)], ploidy=ploidy)  #pooled correlation between blocks
    }

    ## step 2: estimate the effective number of independent tests (Meff)
    eval <- eigen(cm, symmetric=T, only.values=T)   # eigenvalues
    epos <- eval$values[eval$values>=0]     # select positive only
    epos2 <- as.numeric(epos>=1) + (epos - floor(epos))   # LiJi core method
    tMeff <- sum(epos2)    # number of independent tests

    # add to a vector
    Meff <- c(Meff, tMeff)

    #progress
    cat("Meff of chromosome",i,"=",tMeff,"\n")
  }

  ## step 3: adjust alpha
  alpha.adj <- 1-(1-alpha)^(1/sum(Meff))  # similar to (alpha/Meff)


  return(list(threshold=alpha.adj, Meff=Meff))
}

LD_decay <- function(
  dos,
  map,
  win_size = 0.1,
  max_dist = NULL,
  per_chr = F,
  percentile = c(0.5,0.8,0.9,0.95)
  ){

  #First we split dosages per chromosome and calculate correlations
  #only within chromosomes
  non_segregant <- apply(dos,1,function(d) length(unique(d)) == 1)
  dos <- dos[!non_segregant,]
  map <- map[!non_segregant,]

  dos_per_chr <- split(1:nrow(dos),map$chromosome)
  dos_per_chr <- lapply(dos_per_chr,function(x) dos[x,])

  if(as.character(0) %in% names(dos_per_chr)){
    non <- which(names(dos_per_chr) == as.character(0))
    dos_per_chr <- dos_per_chr[-non]
    map <- map[map$chromosome != 0,]
  }

  LD <- lapply(dos_per_chr,function(g){
    ld <- cor(t(g))
    not_dup <- lower.tri(ld,diag=F)
    return(ld[not_dup]^2)
  })

  #then we calculate the distance between markers
  pos_per_chr <- split(map$position,map$chromosome)
  dis <- lapply(pos_per_chr,function(d){
    res <- sapply(d,function(p) abs(p-d))
    not_dup <- lower.tri(res,diag=F)
    return(res[not_dup])
  })

  if(!per_chr){

    #Per window we calculate the percentile estimates
    LD <- do.call(c,LD)
    dis <- do.call(c,dis)

    if(is.null(max_dist)) max_dist <- max(dis,na.rm=T)

    windows <- seq(0,max_dist,win_size)
    ld_estimates <- t(sapply(seq_along(windows),function(w){
      sel <- dis >= windows[w] & dis < windows[w] + win_size
      return(quantile(LD[sel],percentile,na.rm = T))
    }))

    #We calculate the background correlation between chromosomes
    dos_sample <- lapply(dos_per_chr,function(d) d[sample(1:nrow(d),100),] )
    dos_sample <- do.call(rbind,dos_sample)
    back_ld <- cor(t(dos_sample))^2
    back_ld <- back_ld[lower.tri(back_ld,diag=F)]
    back_ld <- quantile(back_ld,percentile,na.rm = T)

    #We add some extra features
    res <- list(LD = data.frame(ld_estimates,
                                distance = windows),
                background = back_ld)
    attr(res,"max_dist") <- max_dist
    attr(res,"class") <- c("LD","list")

  }else{
    res <- lapply(1:length(dos_per_chr),function(i){
      ld <- LD[[i]]
      d <- dis[[i]]

      windows <- seq(0,max(d,na.rm=T),win_size)
      ld_estimates <- t(sapply(seq_along(windows),function(w){
        sel <- d > windows[w] & d < windows[w] + win_size
        return(quantile(ld[sel],percentile,na.rm=T))
      }))

      back_ld <- quantile(ld[d > max(d,na.rm=T)*0.9],percentile)

      #We add some extra features
      res <- list(LD = data.frame(ld_estimates,distance = windows),
                  background = back_ld)
      attr(res,"max_dist") <- max(d,na.rm=T)
      attr(res,"class") <- c("LD","list")

      return(res)
    })
    names(res) <- names(LD)
  }

  return(res)
}


# Input handling ----------------------

#' Comparison of dimensions between matrices of a mixed model.
#'
#' @description Dimensions of \eqn{y}, \eqn{X}, \eqn{K} and \eqn{Z}
#' matrices are checked to ensure that
#' they are mathematically compatible. That is to say: \eqn{y} is a vector of
#' length \eqn{n}, \eqn{X} is a matrix of \eqn{n*m} where \eqn{m} is the number of parameters,
#' \eqn{K} is a square matrix of \eqn{n*n} and Z is a matrix of
#' \eqn{k*n}.
#'
#' @param y
#' @param X
#' @param K
#' @param Z
#'
#' @return
#' @export
#'
#' @examples
test.compatibility<-function(
  y,
  X,
  Z,
  K,
  Hinv=NULL
){
  # Check y -------
  if(is.vector(y)) y <- matrix(y,ncol=1)
  #Normalize y, i ncase it is not
  if(any(round(colMeans(y,na.rm=T)) != 0)){
    y<-apply(y,2,function(p) (p-mean(p,na.rm=T))/sd(p,na.rm=T))
  }
  dimy <- dim(y)

  # Check X -------
  if(is.vector(X)) X <- matrix(X,ncol=1)
  #Normalize X (in case it is not)
  X <- apply(X,2,function(x){
    if(length(unique(x)) == 1) return(x)
    else return((x-mean(x,na.rm = T))/sd(x,na.rm=T))
  })
  dimX <- dim(X)
  if(!comp.vec(dimy,dimX)[1])
    stop("y and X dimensions do not match. Maybe transpose X?")

  # Check Z ------
  dimZ <- dim(Z)
  yZ <- comp.vec(dimy,dimZ)
  if(!yZ[1]) stop("y and Z dimensions do not match. Maybe transpose Z?")

  # Check K ------
  if (dim(K)[1] != dim(K)[2])  stop("K must be square matrix")
  dimK <- dim(K)
  yK <- comp.vec(dimy,dimK)
  if(!all(yK[c(1,2)])) stop("y and K dimensions do not match.")

  # Clean NAs -----
  clean <- lapply(1:ncol(y),function(p){
    nay <- is.na(y[,p])
    nax <- apply(is.na(X),1,any)
    nas <- nax|nay
    #Changing dimensions for NAs
    y <- y[!nas,p,drop=F]
    X <- X[!nas,,drop=F]
    K <- K[!nas, !nas]
    Z <- Z[!nas, ]
    Z <- Z[,!colSums(Z)==0]

    return(list(y=y,X=X,K=K,Z=Z,na.omit=nas))
  })

  miss<-sapply(clean,function(d){
    paste(which(d$na.omit),collapse=" ")
  })

  #This matrix indicates how many combinations of
  #missing values we have. Columns = final list length
  #rows = number of phenotypes
  miss<-lapply(unique(miss),function(m) m==miss)
  miss<-do.call(cbind,miss)

  result<-lapply(1:ncol(miss),function(i){
    dat<-clean[miss[,i]]
    y <- sapply(dat,function(d) d$y)
    #We take the first one of each because X, K and Z should be
    #the same if y's are equivalent
    X <- dat[[1]]$X
    K <- dat[[1]]$K
    Z <- dat[[1]]$Z
    nas <- dat[[1]]$na.omit
    res<-list(y=y,X=X,K=K,Z=Z,na.omit=nas)

    if(!is.null(Hinv)){
      Hi <- Hinv[miss[,i]]
      res$Hinv <- Hi
    }

    return(res)
  })

  return(result)
}

#' Original order
#'
#' calc.Hinv calculates Hinv for every set of different dimensions
#' that are present in the data to improve speed (certain calculations
#' only depend on the dimensionality of the data). As a result, the Hinv
#' matrices are returned in an order that does not correspond the phenotype
#' order that is provided, causing incompatibilities later. This function
#' returns an index vector that can be used to reorder the Hinv matrix list.
#'
#' @param y
#' @param X
#'
#' @return
#' @export
#'
#' @examples
original_order <- function(y,X){
  #With the lines below we obtain how many types of
  #missing matrices we have in our data
  if(is.vector(y)) y <- matrix(y, ncol = 1)
  if(is.vector(X)) X <- matrix(X, ncol = 1)
  miss <- sapply(1:ncol(y),function(p){
    nay <- is.na(y[,p])
    nax <- apply(is.na(X),1,any)
    nas <- nax|nay
    nas <- paste(which(nas),collapse=" ")
    return(nas)
  })
  miss <- lapply(unique(miss),function(m) m==miss)
  miss <- do.call(cbind,miss)

  #This little piece of code here allows us to recover
  #in which position does each matrix end
  #It's hard to understand, but it works well
  miss[miss] <- 1:sum(miss)
  reorder_vec <- rowSums(miss)

  return(reorder_vec)
}


#' Sample markers every X centimorgans
#' @description Obtains a subset of markers homogeneously distributed across
#' a genetic map based on a cM distance. For instance, if cM=1, will get 1 marker
#' for every 1 cM window if it can. Markers added on chromosome 0 are not returned
#' as their genetic position is not known (and thus cannot be cM sampled)
#'
#' @param genotypes Either a vector where each element is a marker, or a matrix
#' where each row is a marker.
#' @param map table containing a "chromosome" (numerical) and a "position" (in cM) column.
#' Markers with chromosome number == 0 will not be taken into account
#' @param cM Size of bin to sample a marker from.
#' @param seed An integer to set a seed for random number generation.
#'
#' @return A table as genotypes, with as many rows as can be obtained by sampling
#' markers every X centimorgans.
#' @export
#'
#' @examples
sample.cM<-function(
  genotypes,
  map,
  cM=1,
  seed=NULL
){
  if(is.vector(genotypes)) genotypes<-matrix(genotypes,ncol=1)
  if(nrow(genotypes)!=nrow(map)){
    stop("Number of markers in genotypes and in map does not coincide")
  }

  result <- lapply(unique(map$chromosome),function(x){
    pos <- map$position[map$chromosome == x]
    selpos <- seq(0,max(pos),by = cM)

    markers <- sapply(selpos,function(k){
      #what is the closest position to this marker
      dif <- abs(k-pos)

      #the index of the closest marker
      a <- which(dif==min(dif))

      #in case there is more than one marker at that position, we select randomly 1
      if(length(a) > 1){
        set.seed(seed)
        a <- sample(a,1)
      }
      return(a)
    })

    return(as.matrix(genotypes[map$chromosome == x, ][unique(markers),]))
  })

  #In case there are chromosomes labelled as 0, they should not
  #be used for centimorgan sampling

  if(any(map$chromosome == 0)){
    result <- result[-which(unique(map$chromosome) == 0)]
  }

  result <- do.call(rbind,result)
  return(result)
}

stndrdz.pheno <- function(phenotypes) {
  if(is.vector(phenotypes)) {
    phenotypes <- matrix(phenotypes, ncol=1,
                         dimnames = list(names(phenotypes)))
  } else if (is.data.frame(phenotypes)) {
    phenotypes <- as.matrix(phenotypes)
  }

  means <- colMeans(phenotypes,na.rm = T)
  sds <- apply(phenotypes,2,sd,na.rm=T)

  phenotypes <- apply(phenotypes,2,function(p){
    (p-mean(p, na.rm = T))/sd(p, na.rm = T)
  })

  return(list(pheno = phenotypes,
              mean = means,
              sd = sds))
}



## used in map.QTL and thr.LiJi
#' Check format of a numeric input matrix
#'
#' @param x a matrix or a data frame.
#' @param integer Logical value indicating whether to check that the input
#' matrix contains integers only (e.g. allelic dosages). Default if
#' integer = TRUE. If FALSE, non-integers are allowed and the ploidy check is
#' skipped.
#' @param ploidy A number specifying the ploidy level. If a ploidy level is
#' provided, the function will check whether all the data fit with that ploidy.
#'
#' @return a matrix
#' @export
#'
#' @examples
inputCheck_dos <- function(x, integer=TRUE, ploidy=NULL) {

  # conversion to matrix?
  if(!is.matrix(x)){
    cat("x is not a matrix.")
    if(is.data.frame(x)) {
      cat("I will try to convert it into a matrix.")
      # For converting to matrix I am using sapply(x, as.character) instead of
      # as.matrix(x). This is because as.matrix can introduce some extra
      # spaces (by calling format(x, trim=F)). After conversion to numeric
      # (in the next section), numbers will be correctly converted to numeric
      # elements, even if they contain extra spaces. All the other elements
      # will be set to NA. However, I want to show to the user the problematic
      # elements set to NAs. And I need to show the original ones, without
      # extra spaces.
      # Consider performance.
      x <- sapply(x, as.character)
      if(!is.matrix(x)) stop("x cannot be converted into a matrix.")
    } else {
      stop("x cannot be converted into a matrix.")
    }
  }

  # conversion to numeric?
  if(!is.numeric(x)) {
    cat("\nConversion to numeric.")
    u1 <- unique(c(x))
    u2 <- as.numeric(u1)
    d <- is.na(u2)-is.na(u1) # elements coerced to NA
    if(sum(d)>0) {
      toNA <- u1[as.logical(d)]
      warning(paste("The following elements have been set to NA:\n",
                    paste(toNA, collapse = ","),"\n"),
              immediate. = F)
    }
    # x <- apply(x, 2, as.numeric)  #less efficient
    suppressWarnings(storage.mode(x) <- "numeric")  # conversion to numeric
  }

  if (integer) {
    # conversion to integer?
    if(!is.integer(x)) {
      u <- unique(c(x))
      u <- u[!is.na(u)]
      m.int <- u %% 1 == 0 # check if numbers are mathematically integers
      if(all(m.int, na.rm = T)) {
        # x <- apply(x, 2, as.integer)  #less efficient
        storage.mode(x) <- "integer" # conversion to integer
      } else {
        stop(paste("The following non-integers are not allowed:\n",
                   paste(u[!m.int], collapse = ","),
                   "\n"))
      }
    }

    # check the expected range of dosages
    if(!is.null(ploidy)) {
      u <- unique(c(x))
      u <- u[!is.na(u)]
      expdos <- u %in% 0:ploidy
      if(any(!expdos, na.rm = T)) {
        stop(paste("The following dosages do not fit the ploidy level:\n",
                   paste(u[!expdos], collapse = ","),
                   "\n"))
      }
    }
  }
  x
}

# Imputation ----------------------------------

#' NA replacement by mean score per marker
#'
#' @description For some operations (e.g. calculating correlation) the
#' presence of missing values is not allowed. If the level of missingness
#' per marker is not too high (below 25%), missing values can be imputed
#' using the mean score per marker, so that allele frequency is not affected.
#'
#' @param m a numeric matrix, with markers in rows and
#' individuals in columns.
#' @param miss a number used for missing values instead of NA.
#' The function expects missing values are indicated by NA. Alternatively,
#' since the inpute matrix is a numeric matrix, you could indicate
#' missing values by a number.
#'
#' @return a matrix with no missing values.
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

#' KNN genotype imputator
#'
#' @description K-Nearest Neighbour imputation algorithms are useful for
#' imputing missing data. In this function, we use the kinship matrix
#' calculated with \code\link{calc.K} to determine the most similar
#' individuals, and substitute the missing values of one individual with
#' the most common dosage/haplotype of its nearest neighbours.
#'
#' @param geno Haplotype (accepts characters) or Dosage matrix
#' @param ploidy Ploidy of the organism
#' @param map Optional, genetic map used in \code\link{sample.cM}.
#' Improves efficiency
#' @param kneighbors Number of nearest neighbours to use in the imputation. Default is 20.
#' @param K Optionally, provide the similarity matrix (K) directly.
#'
#' @return A Haplotype or Dosage matrix without missing values.
#' @export
#'
#' @examples
impute.knn <- function(
  geno,
  ploidy,
  map = NULL,
  kneighbors = 30,
  K=NULL
){

  #We determine whether we are looking at dosages or haplotypes
  if (all(unique(geno) %in% c(0:ploidy, NA))){
    cat("Imputation will be performed on dosages.\n")
    haplo<-F
  }else if (ncol(geno) %% ploidy == 0){
    cat("Imputation will be performed on haplotypes.\n")
    haplo<-T
  }else stop(paste("Genotypes not recognized as dosage nor haplotypes.",
                   "Dosage values are not between 0 and ploidy or",
                   "ncol of genotypes is not a multiple of ploidy."))


  #Calculate K distance matrix. If possible using fewer markers
  #homogeneously distributed across the genome.
  if(is.null(K)){
    if (!is.null(map)) genoK <- sample.cM(geno, map)
    else genoK <- geno
    K <- calc.K(t(genoK), ploidy = ploidy, haplotypes = haplo)
  }


  #Imputation on dosages
  #Each individual is a column
  if (!haplo) {
    result <- sapply(1:ncol(K), function(d) {
      imputed <- geno[,d]
      if(all(!is.na(imputed))) return(imputed)

      #Find the closest individual to individual d
      best <- order(K[, d][-d], decreasing = T)[1:kneighbors]


      #Get the genotypes of the best individuals
      #at the position where current individual has missing values
      knei <- geno[is.na(imputed), best]
      if(is.vector(knei)) knei<-matrix(knei,nrow=1)
      #Percentage of missing values per marker on the nearest neighbours
      naperc<-apply(knei,1,function(i) sum(is.na(i))/kneighbors)
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
      #Same goes for the individual of interest
      dG <- (d * ploidy - (ploidy - 1)):(d * ploidy)
      imputed <- geno[, dG]
      if(all(!is.na(imputed))) return(imputed)
      #Best individual selection is done the same
      best <- order(K[-d, d], decreasing = T)[1:kneighbors]
      #but we need to transform individual index into
      #the chromosome column index
      best <- as.vector(sapply(best, function(b)
        (b * ploidy - (ploidy - 1)):(b * ploidy)))

      #Missing values are detected in ANY of the haplotypes
      miss <- is.na(imputed)
      miss <- apply(miss, 1, any)

      #We get the genotypes of the neighbours
      knei <- geno[miss, best]
      if(is.vector(knei)) knei<-matrix(knei,nrow=1)

      #Are they expressed numerically? The following function coerces to char
      asnum <- is.numeric(knei)

      nearest <- t(apply(knei, 1, function(alleles) {
        #If there are over 50% missing values, we return NA
        if(sum(is.na(alleles))/length(alleles)>0.5) return(rep(NA,4))

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

      imputed[miss, ] <- nearest

      #Uncertain values are substituted by the mean population allele
      if(any(is.na(imputed))){
        miss <- is.na(imputed)
        miss <- apply(miss, 1, any)

        markers<-names(which(miss))
        error<-paste("Nearest",kneighbors,"neighbours of individual",d,
                     "have over 50% missing values at:",
                     paste(markers,collapse=" "),
                     "\nThey will be imputed using all population genotypes instead")
        warning(error)

        nearest <- t(apply(geno[miss,,drop=F],1,function(alleles){  #gt: drop=F added
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

        imputed[miss,] <- nearest

      }

      return(imputed)
    })
    result <- do.call(cbind, result)
  }

  return(result)

}


# Miscellaneous -------------------

#' Vector comparison function
#' @description Element per element comparison between two vectors.
#' @param vec1 vector of values
#' @param vec2 vector of values
#'
#' @return Matrix of logical values indicating whether the values of vec1
#' and vec2 are equal. Columns correspond to elements
#' of vec1 and rows to elements of vec2.
#' @export
#'
#' @examples
comp.vec <- function(vec1, vec2) {
  mat <- sapply(vec1, function(x)
    x == vec2)
  colnames(mat) <- vec1
  rownames(mat) <- vec2
  return(mat)
}
