
#' QTL mapping of a matrix of phenotypes
#'
#' @description Main wrapper function for QTL analysis in polyploids, using
#' biallelic (continuous or discrete) or multiallelic markers.
#' \code{map.QTL} implements a single-marker regression mixed model, that enables
#' to account for different levels of population structure \href{https://www.nature.com/articles/ng1702}{(Yu et al. 2006)}:
#' \deqn{y = \mu + X\beta + Qv + Cw + Zu + \epsilon}
#' \deqn{var(u) ~ K\sigma^{2}G}
#' \deqn{var(\epsilon) ~ R\sigma^{2}\epsilon}
#' A phenotype y is modelled using an intercept \eqn{\mu}; a genetic matrix \eqn{X\beta} (containing
#' dosages or haplotypes), a family-based correction matrix Qv; an optional set of
#' cofactors, modelled by Cw; a random term that models kinship using K as the
#' variance structure, and a normal error term \eqn{\epsilon}. We can also talk
#' of this parameters referring to the "genetic term" (\eqn{X\beta}), the "genetic
#' structure correction" terms (Qv and \eqn{var(u)}) and the cofactors (Cw).
#'
#' @param phenotypes A numeric matrix of phenotypes,
#' rows are individuals and columns are different phenotypes.
#' Not updated for vector phenotypes yet. Must follow same
#' order of individuals as genotypes.
#' @param genotypes A matrix of SNP genotypes (continuous or discrete) or
#' haplotypes (multi-allelic markers), where rows are markers and
#' columns may be either (1) individuals (for biallelic markers) or (2)
#' individual homologues (for multi-allelic markers, such as haplotypes).
#' Genotypes must follow the same order of individuals as phenotypes.
#' @param ploidy A number indicating the ploidy level. All the individuals
#' must have the same ploidy.
#' @param K It can be 1) NULL, 2) TRUE or 3) a distance matrix. If NULL, no relatedness
#' matrix will be used (i.e. a linear model will be applied). If TRUE,
#' a K distance matrix will be calculated. Otherwise, a pre-calculated
#' distance matrix may also be directly specified.
#' @param Q It can be 1) NULL, 2) TRUE, 3) a vector identifying populations
#' or 4) a Q design matrix. If NULL, no Q will be included in the model (i.e. a
#' model without Q correction). If TRUE, a PCo decomposition will be used to
#' estimate population differentiation. If a vector specifying population of
#' each individual is passed, it will be used to construct a Q matrix. Vector
#' may contain numerical or character.
#' @param map A data frame with 3 columns containing map information. The first
#' column specifies marker names, the second column chromosome names and the
#' third one marker position (any map unit).
#' Marker position is used for plotting or it could be used to sample a subset
#' of evenly spaced markers to be used for kinship estimation.
#' @param binsize Numeric. K distance matrix will be calculated using markers every
#' map unit (for physical maps use Mb). Defaults to 1.
#' @param seed An integer to set a seed for random number generation in \code{sample.marker}.
#' @param no_cores Numeric. Number of cores to be used in parallel computing
#' of the p-values. Defaults to number of cores -1.
#' @param Z Identity matrix indicating which individuals correspond to which
#' genotypic effects. For instance, if multiple samples correspond to the same
#' individual, this matrix should indicate so.
#' @param cofactor Possible cofactor matrix (where each column is a cofactor).
#' @param cofactor.type If a cofactor matrix is specified, a character vector
#' specifying the type ("numerical" or "categorical") of each cofactor.
#' It accepts partial strings such as "cat" and "num".
#' @param Qpco Logical value indicating whether a PCo-based population factor should be
#' estimated and included. It is useful for detecting and correcting for the
#' main axes of genetic variation (i.e. population structure).
#' @param approximate Logical value indicating whether P3D/EMMAX approach should be
#' used for computational efficiency (if FALSE, expect much longer waiting times)
#' @param permutation permutation strategy. Can be "pop" (permutation over the
#' whole population) or "fam" (permutation within families). If it is NULL, no
#' permutation will be run.
#' @param fam If \code{permutation = "fam"}, provide here a vector indicating the family
#' membership (or a population structure membership) of each individual.
#' @param alpha Permutation threshold alpha value, default to 0.05. Alpha is the
#' chance of observing a false positive experiment-wise.
#' @param impute Logical value indicating whether missing genotypes should be
#' imputed using \code{impute.knn}. Defaults to True, but False is recommended (imputation algorithm
#' needs to be improved and with few missing values it has a small effect
#' on QTL detection).
#' @param nperm number of permutations, if \code{permutation} is not null.
#' @param knn Number of neighbours to use in the internal function \code{impute.knn}
#' @param linear logical. If TRUE, linear model (without structure correction)
#' is applied. If FALSE, mixed model (with structure correction) is applied. If
#' not specified, it looks at the value of K, and only applied linear model
#' if K = NULL.
#' @param K_identity logical. If TRUE, an identity matrix is used in the
#' random term of a mixed model. This is to run a mixed model with no kinship
#' correction, even if a kinship matrix is provided or calculated for other
#' purposes (e.g. imputation of missing genotypes, calculation of the Q matrix).
#'
#' @details The arguments \code{linear} and \code{K} can be used to define either a fixed
#' effects or a mixed effects model. By default they are both NULL, that will
#' result in applying a fixed effects linear model.
#' If \code{K} is not NULL (either TRUE or a distance matrix) a mixed model will
#' be used. If a distance matrix is provided for other purposes than kinship
#' correction (e.g. imputation of missing genotypes, calculation of the Q matrix),
#' specify \code{K_identity} = TRUE to use an identity matrix for structuring
#' the variance of the random term. This would be equivalent to a naive model
#' (unless any additional fixed effect parameter is specified).
#' @return A list of length equal to the number of phenotypes. Each element is
#' a list with model coefficients and test results:
#' \itemize{
#'   \item \code{$beta} A list containing the fixed effects of the model (in this
#'   order: intercept, cofactors, structure terms, genetic term). There is a set
#'   of estimates for each marker.
#'   \item \code{$Fstat} A vector containing the F test results only for the genetic
#'   component of each model.
#'   \item \code{$residual}
#'   \item \code{$pval} A vector containing the p-values only of the genetic model
#'   at each marker.
#'   \item \code{$se} A vector containing the standard error of the estimates.
#'   \item \code{$Wald}
#'   \item \code{$real.df}
#' }
#'
#' @export
#' @examples
#' ## Get example genotypes (haplotypes in this case),
#' ## map and phenotypes for a population of tetraploid individuals
#' data("mphapdose")
#' data("mpmap")
#' data("mppheno")
#'
#' ## fixed effects model (y = m + e)
#' results <- map.QTL(phenotypes = mppheno,
#'                    genotypes = mphapdose,
#'                    ploidy = 4,
#'                    map = mpmap,
#'                    no_cores = 2)
#' names(results$phenotype1)
#' skyplot(-log10(results$phenotype1$pval), map = mpmap)
#'
#' ## naive model (y = m + e) - no K correction
#' results <- map.QTL(phenotypes = mppheno,
#'                    genotypes = mphapdose,
#'                    ploidy = 4,
#'                    map = mpmap,
#'                    K = TRUE,
#'                    K_identity = TRUE,
#'                    no_cores = 2)
#' names(results$phenotype1)
#' skyplot(-log10(results$phenotype1$pval), map = mpmap)
#'
#' ## model with Q correction (y = m + Q + e)
#' results <- map.QTL(phenotypes = mppheno,
#'                    genotypes = mphapdose,
#'                    ploidy = 4,
#'                    map = mpmap,
#'                    Q = TRUE,
#'                    Qpco = 2,
#'                    K = TRUE,
#'                    K_identity = TRUE,
#'                    no_cores = 2)
#' names(results$phenotype1)
#' skyplot(-log10(results$phenotype1$pval), map = mpmap)
#'
#' ## model with K correction (y = m + K + e)
#' results <- map.QTL(phenotypes = mppheno,
#'                    genotypes = mphapdose,
#'                    ploidy = 4,
#'                    map = mpmap,
#'                    K = TRUE,
#'                    no_cores = 2)
#' names(results$phenotype1)
#' skyplot(-log10(results$phenotype1$pval), map = mpmap)
#'
#' ## model with Q + K correction (y = m + Q + K + e)
#' results <- map.QTL(phenotypes = mppheno,
#'                    genotypes = mphapdose,
#'                    ploidy = 4,
#'                    map = mpmap,
#'                    Q = TRUE,
#'                    Qpco = 2,
#'                    K = TRUE,
#'                    no_cores = 2)
#' names(results$phenotype1)
#' skyplot(-log10(results$phenotype1$pval), map = mpmap)
#'
map.QTL.notparallel <- function(phenotypes, genotypes, ploidy, map, K=NULL, Q=NULL, Z=NULL,
                    cofactor=NULL, cofactor.type=NULL, binsize=1, seed=NULL, Qpco=2,
                    no_cores=parallel::detectCores()-1, approximate = TRUE,
                    permutation = NULL, fam=NULL, nperm = 100, alpha = 0.05,
                    impute=TRUE, knn=20, linear = NULL, K_identity = FALSE){


  # input check ---------------------------------------------

  ## phenotypes
  std_phe <- stndrdz.pheno(phenotypes)
  phenotypes <- std_phe$pheno

  ## genotypes
  ### Check whether genotypes contain SNPs or haplotypes and
  ### ensure requirements are met.
  ### This part needs to be changed to allow the new list for
  ### probabilistic haplotypes
  if (ncol(genotypes)==nrow(phenotypes)) {
    haplo <- FALSE
    genotypes <- inputCheck_snp(genotypes, integer=TRUE, ploidy=ploidy)
    if (is.integer(genotypes)) {
      cat("SNP dosages have been detected in the genotype matrix,\n")
    } else {
      cat("Continuous SNP genotypes have been detected,\n")
    }
  } else if ((ncol(genotypes) %% ploidy) == 0){
    haplo <- TRUE
    cat("Haplotypes have been detected in the genotype matrix.\n")

  } else {
    stop("If haplotypes have been provided: the number of genotypes columns
         is not a multiple of ploidy.
         If SNP dosages have been provided: the number of individual dosages
         is not equal to the number of individual phenotypes.")
  }

  #bi.geno
  #MAF

  ## map
  colnames(map) <- c("marker","chromosome","position")
  markers <- map$marker

  ## K
  ### check K is either 1) NULL, 2) TRUE or 3) a numeric square matrix.
  input_K <- inputCheck_K(K)
  K <- input_K$K
  input_K$K <- NULL

  ## Q
  ### check Q is either 1) NULL, 2) TRUE, 3) a membership vector or
  ### 4) a design matrix
  input_Q <- inputCheck_Q(Q)
  Q <- input_Q$Q
  input_Q$Q <- NULL

  ## linear
  ### assign linear based on K
  if(is.null(linear)){
    if(is.null(K)){
      linear <- TRUE
    }else{
      linear <- FALSE
    }
  }



  # input order --------------------------------------------------
  # Match and re-order individulas and markers of input matrices
  # ordered_input <- inputOrder(genotypes, pheno=phenotypes, map=map,
  #                             cof=cofactor, Q=Q, K=K)








  #DEFINITION OF K -------------------
  #There are three options
  #1) K=NULL, no K is used. We go linear.
  #2) K=TRUE, K is calculated using homogeneously distributed markers along a map.
  #   using sample.marker and calc.K
  #3) K=K matrix, check dimensions, message is printed and K is used.


  if(!is.null(K)){ #gt: added because when K=NULL the condition below is TRUE
    if(all(K==TRUE)){ #gt: why are you using all?
      #Option 2)
      K<-sample.marker(genotypes,map,binsize = binsize, seed=seed)
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
  if(naperc == 0) impute<-FALSE
  if(impute){
    if(is.null(K)){
      #first we sample homogeneously markers along the genome
      K <- sample.marker(genotypes,map,binsize = binsize, seed=seed)
      #then we calculate the distance
      K <- calc.K(t(K),ploidy=ploidy,haplotypes = haplo)
    }
    genotypes<-impute.knn(genotypes,ploidy,map,kneighbors=knn,K=K)
    cat("Imputation performed.\n")
  }else{
    cat("No imputation will be performed.\n")
  }

  #DEFINITION OF Q ------------------------
  #4 options:
  #1) No Q is used (Q=NULL), this block is skipped
  #2) Q=TRUE, Q is calculated using the K matrix and cmdscale
  #3) Q=vector identifying groups, Q is calculated using Q.mat
  #4) Q=Q design matrix, message is printed.
  if(is.null(Q)) cat("No Q matrix will be used.\n")

  if(!is.null(Q)){
    if(all(Q==TRUE)){ #gt: why are you using all?
      #atn: because if Q is a vector/matrix Q==TRUE will be a vector of FALSE
      #and that gives an error in the if.

      #Option 2)
      if(is.null(K)){
        #first we sample homogeneously markers along the genome
        K<-sample.marker(genotypes,map,binsize = binsize, seed=seed)
        #then we calculate the distance
        K<-calc.K(t(K),ploidy=ploidy,haplotypes = haplo)
      }
      Q <- cmdscale(1-K, k=Qpco, eig = FALSE, add = FALSE, x.ret = FALSE)
      colnames(Q) <- paste0("Q",1:ncol(Q))

    }else if(is.vector(Q)){
      #Option 3)
      Q<-calc.Q(Q)

    }else if(nrow(Q)==nrow(phenotypes)){
      #Option 4)
      cat("Q matrix has been provided. Using as it is.\n")

    }else{
      #Q was not NULL, TRUE, a vector or a cofactor design matrix.
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
    c.num <- pmatch(cofactor.type,"numerical", duplicates.ok = TRUE)
    c.cat <- pmatch(cofactor.type,"categorical", duplicates.ok = TRUE)

    C <- lapply(1:ncol(cofactor),function(i){

      if(!is.na(c.num[i])){
        result <- as.numeric(cofactor[,i])

      }else if(!is.na(c.cat[i])){

        #This creates a design matrix with 1s and 0s (representing TRUE and FALSE)
        vals <- na.omit(unique(cofactor[,i]))
        result <- sapply(vals,function(c){
          cofactor[,i] == c
        })+1-1

        #Last column is taken out to be able to calculate
        #gt: is this correct? Check!
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
    export <- c("phenotypes","genotypes","dosage.X","Q","C","lm_compare","haplo")

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
        X <- dosage.X(genotypes[k,],ploidy = 4,haplotype = haplo,normalize = FALSE)
        if(haplo) X <- X[,-ncol(X)]
        #We create the naive and full matrices. The naive
        #model only includes Q, C and an intercept
        N <- cbind(matrix(1,nrow=nrow(X)),Q,C) #naive model
        X <- cbind(1,Q,C,X) #total model
        y <- phenotypes[,w,drop=FALSE]
        solveout <- try(lm_compare(X,N,y),silent = TRUE)

        ## mantain the usual output structure in case of error
        if (class(solveout) == "try-error") {
          # write(solveout, file="lm_compare_errors.txt", append = TRUE)

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
        res <- min(res$pval, na.rm = TRUE)
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
    #     res <- min(res$pval, na.rm = TRUE)
    #   }
    #   return(res)
    # })


  }else{ #Ergo, K must be defined, we must use Mixed Models


    ####
    # Mixed Model --------------
    if (K_identity == TRUE) {
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

    # cl<-parallel::makeCluster(no_cores)
    # #atn: added object ploidy
    # export<-c("phenotypes","Z","K","Hinv","genotypes","ploidy","haplo","npheno",
    #           "mm.solve","dosage.X","Q","C","test.compatibility","comp.vec","std_phe") #gt
    #
    # #extra functions need to be exproted
    # #if we don't want to use P3D approximation
    # if(!approximate) export<-c(export,"calc.Hinv")
    #
    # parallel::clusterExport(cl,export,
    #                         envir=environment()) #I dunno why, but without this it doesnt work
    # # #gt: same reason explained before. Here, doesn't even work because Z, K and Hinv exist only
    # # in this executing environment (and not in the global one).

    # #gt version: markers in parallel (parallelization works even with one phenotype)
    result<-lapply(1:ncol(phenotypes),function(w){#NOT parallely over each phenotype
      result<-sapply(cl,1:nrow(genotypes),FUN=function(k){#parallely over markers, calculate the pvalue for each marker

        #DEFINITION OF X (matrix of fixed effects)
        X <- dosage.X(as.matrix(genotypes[k,]),
                      ploidy=ploidy,
                      normalize=TRUE,
                      haplotype = haplo)

        if(ncol(X)>1){X<-X[,-1,drop=FALSE]} #when ancestral/parental model we need to prevent singularity
        # if(any(is.na(genotypes[k,]))) X<-X[,-ncol(X)] #Eliminate NA as factor
        nparX <- ncol(X) #total number of genetic parameters
        X <- cbind(Q,C,X) #add population and cofactor parameters
        dimnameX <- dimnames(X)
        X <- matrix(as.numeric(X),ncol=ncol(X))
        dimnames(X) <- dimnameX
        no.test<-ncol(X)-nparX #number of non genetic parameters

        #If not P3D/EMMAX, then Hinv will be NULL (and cannot be subindexed)
        if(!is.null(Hinv)){
          H<-Hinv[w] #gt: Hinv[[w]] will be a matrix, but mm.solve is expecting a list!
          # Hinv[w] will be a list of length 1
        }else{
          H<-NULL
        }
        #gt: phenotypes[,w,drop=FALSE] to avoid conversion to vector
        #gt: added [[1]], after replacement of sapply with lapply in mm.solve
        #solveout <- mm.solve(phenotypes[,w,drop=FALSE],X,Z,K,H,no.test = no.test)
        solveout <- try(mm.solve(phenotypes[,w,drop=FALSE],X,Z,K,H,no.test = no.test), silent = TRUE)

        ## mantain the usual output structure in case of error
        if (class(solveout) == "try-error") {
          # write(solveout, file="mm.solve_errors.txt", append = TRUE)
          solveout <- list(list(
            beta = matrix(NA),
            Fstat = NA,
            residual = rep(NA,nrow(X)),
            pval = NA,
            se = NA,
            wald = NA,
            real.df = NA
          ))
        }else{
          #This turns the standardized effects into un-standardized and adds
          #the mean of the population as a parameter (intercept). Technically this
          #is not very accurate...
          solveout <- lapply(1:length(solveout),function(sol){
            res <- solveout[[sol]]
            par_names <- rownames(res$beta)
            res$beta <- rbind(std_phe$mean[sol], solveout[[sol]]$beta*std_phe$sd[sol])
            rownames(res$beta) <- c("",par_names)
            return(res)
          })
          names(solveout) <- colnames(phenotypes)[w]
        }

        return(solveout)
      })

      #NEW OUTPUT
      #Create a list of lists from all results.
      #It's actually just 6 lists with k elements where k is markers
      res<-do.call(mapply,c(list,result))
      res<-as.list(as.data.frame(res))
      #All columns are turned into vectors except the beta column
      #gt: and except residuals (to fix the bug below)
      for(r in c(2,4:length(res))) res[[r]] <- unlist(res[[r]])
      # for(r in 1:length(res)) names(res[[r]]) <- markers
      for(r in c(1:2,4:length(res))) names(res[[r]]) <- markers
      # names(res[[4]]) <- markers #gt: it should be individuals?
      #gt: here there is a bug. When a phenotype contain NAs, the
      # number of residuals will differ, including only residuals
      # of non-missing phenotypes. The line below will split wrongly!
      # res$residual <- split(res$residual,rep(1:nrow(genotypes),
      #                                        each=nrow(phenotypes)))
      names(res$residual) <- markers
      ## for permuted phenotypes store the minimum pvalue only
      if (w > npheno) {
        res <- min(res$pval, na.rm = TRUE)
      }

      return(res)
    })
    #Here we obtain a list of df, each df containing all results
    # parallel::stopCluster(cl)
  }

  #The results are given phenotype names
  #We obtain the structure we talked about
  if(!is.null(colnames(phenotypes))){
    names(result) <- make.names(colnames(phenotypes),unique=TRUE)
  }else{
    names(result) <- paste0("pheno",1:ncol(phenotypes))
  }


  cat("Association completed.\n")

  if (!is.null(permutation)) {
    cat("Permutation threshold will now be calculated.\n")
    minpval <- sapply(result[(npheno+1):length(result)], function(x) x)
    minpvalMat <- matrix(minpval, ncol = npheno, byrow = TRUE)
    maxlogpvalMat <- -log10(minpvalMat)
    thr <- sapply(1:ncol(maxlogpvalMat), function(i) {
      quantile(maxlogpvalMat[,i], probs = 1-alpha) # alpha can be a vector
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









#' QTL mapping of a matrix of phenotypes
#'
#' @param phenotypes A (#gt: numeric) matrix of phenotypes,
#' rows are individuals and columns are different phenotypes.
#' Not updated for vector phenotypes yet #gt: check#. Must follow same
#' order of individuals as genotypes #gt: add check for same order#.
#' @param genotypes A matrix of genotypes, rows are markers and
#' columns may be either (1) individual dosages (#gt: of biallelic markers), with column
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
#' @param no_cores Numeric. Number of cores to be used in parallel computing
#' of the p-values. Defaults to number of cores -1.
#' @param Z
#' @param cofactor
#' @param cofactor.type
#' @param Qpco
#' @param P3D
#' @param EMMAX
#' @param permutation permutation strategy. Can be "pop" (permutation over the
#' whole population) or "fam" (permutation within families). If it is NULL, no
#' permutation will be run.
#' @param alpha
#' @param impute
#' @param nperm number of permutations, if \code{permutation} is not null.
#'
#' @return a pvalue matrix containing the pvalue of each marker with each phenotype passed.
map.QTL2<-function(
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
  Qpco=2, #number of axis used for pco decomposition
  no_cores=parallel::detectCores()-1,
  P3D=T,
  EMMAX=T,
  permutation = NULL, #permutation strategy: "pop" or "fam"
  nperm = NULL, #number of permutations
  alpha = 0.95,
  impute=T,
  k=20
){

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
    haplo<-F
    genotypes <- inputCheck_dos(genotypes, integer=T, ploidy=ploidy)
    cat("SNP dosages have been detected in the genotype matrix\n")

    # if (!all(rownames(phenotypes)==colnames(genotypes))) {
    #   stop("phenotype individuals and genotype individuals have different names")
    # }

  } else if ((ncol(genotypes) %% ploidy) == 0){
    haplo<-T
    cat("Haplotypes have been detected in the genotype matrix\n")

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

  if(all(K==T)){ #gt: why are you using all?
    #Option 2)
    K<-sample.cM(genotypes,map,cM = cM)
    K<-calc.K(t(K),ploidy = ploidy,haplotypes = haplo)

  }else if(nrow(K)==nrow(phenotypes)){
    #Option 3)
    if(ncol(K)!=nrow(phenotypes)) stop("K matrix was specified, but is not square")
    cat("Square K matrix has been provided. Using as it is.")
  }

  ### Na check and impute ----------------
  naperc<-sum(is.na(genotypes))/length(genotypes)
  cat(round(naperc*100,2),"% missing genotypes detected.\n")
  if(naperc==0) impute<-F
  if(impute){
    if(is.null(K)){
      #first we sample homogeneously markers along the genome
      K<-sample.cM(genotypes,map,cM = cM)
      #then we calculate the distance
      K<-calc.K(t(K),ploidy=ploidy,haplotypes = haplo)
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
      stop("Wrong Q specification")}
  }

  ##DEFINITION OF C, cofactor matrix. ---------------------------
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
  # Linear Model ------------
  #with no structure correction or with Q correction
  if(linear){
    #First set up cluster

    cluster<-parallel::makeCluster(no_cores)
    export<-c("phenotypes","genotypes","dosage.X","Q","C") #gt
    if(is.null(Q)){
      # export<-c("phenotypes","genotypes","dosage.X")
      cat("Linear model will be used")
    }else{
      # export<-c("phenotypes","genotypes","dosage.X","Q")
      cat("Linear model with Q correction will be used.\n")
    }

    parallel::clusterExport(cl=cluster,export,
                            envir=environment()) #needed again? #gt: yes, because
    # we want to use variables defined in this executing environment, not in the
    # global one. Although the function dosage.X is not present in this
    # environment, clusterExport will search it and will find it in the parent environments

    #For each row of genotypes, calculate lm for all phenotypes
    pval<-parallel::parSapply(cl=cluster,1:nrow(genotypes),function(k){
      X<-dosage.X(genotypes[k,])
      nparX<-ncol(X)
      #if Q is NULL, X will not change
      X<-cbind(X,Q,C)

      test1<-lm(phenotypes~X)
      SSR1<-colSums((test1$residuals)^2) #Sum of Squares Residuals Full

      test3<-lm(phenotypes~X[,2:nparX]) #Sum of Squares of Markers
      SSR3<-colSums((test3$residuals)^2)

      SST<-apply(phenotypes,2,function(x){
        sum((x-mean(x))^2)
      })
      SSM<-SST-(SSR3-SSR1)-SSR1 #Calculate the amount of variation explained by the markers
      #in a model with a Q structure

      df2<-apply(test1$coefficients,2,function(x) nparX-sum(is.na(x)))

      Ftest<-(SSM/df2)/(SSR1/test1$df.residual)
      return(pf(Ftest,df2,test1$df.residual,lower.tail = F))

    })
    parallel::stopCluster(cluster)
    pval<-t(pval)

  }else{ #Ergo, K must be defined, we must use Mixed Models


    ####
    # Mixed Model --------------
    if(is.null(Z)){
      Z<-diag(nrow(phenotypes))
    }

    if(is.null(Q)){
      cat("Mixed model will be used")
    }else{
      cat("Mixed model with Q correction will be used")
    }

    if(!P3D|!EMMAX){
      Hinv<-NULL
    }else{
      Hinv<-calc.Hinv(phenotypes,
                      #X=vector of 1s, to apply P3D/EMMAX algorithm #gt: check#
                      X=matrix(rep(1,nrow(phenotypes))),
                      Z,K)
    }

    # cl<-parallel::makeCluster(no_cores)
    # #atn: added object ploidy
    # export<-c("phenotypes","Z","K","Hinv","genotypes","ploidy","haplo",
    #           "mm.solve","dosage.X","Q","C","test.compatibility","comp.vec") #gt
    #
    # #extra functions need to be exproted
    # #if we don't want to use P3D approximation
    # if(!P3D|!EMMAX) export<-c(export,"calc.Hinv")
    #
    # parallel::clusterExport(cl,export,
    #                         envir=environment()) #I dunno why, but without this it doesnt work
    # # #gt: same reason explained before. Here, doesn't even work because Z, K and Hinv exist only
    # # in this executing environment (and not in the global one).

    # #gt version: markers in parallel (parallelization works even with one phenotype)
    result<-lapply(1:ncol(phenotypes),function(w){#NOT parallely over each phenotype
      result<-sapply(1:nrow(genotypes),FUN=function(k){#parallely over markers, calculate the pvalue for each marker

        #DEFINITION OF X (matrix of fixed effects)
        X<-dosage.X(as.matrix(genotypes[k,]),
                    ploidy=ploidy,
                    normalize=T,
                    haplotype = haplo)

        if(ncol(X)>1){X<-X[,-1,drop=F]} #when ancestral/parental model we need to prevent singularity
        # if(any(is.na(genotypes[k,]))) X<-X[,-ncol(X)] #Eliminate NA as factor
        nparX<-ncol(X) #total number of genetic parameters
        X<-cbind(Q,C,X) #add population and cofactor parameters
        X<-matrix(as.numeric(X),ncol=ncol(X))
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
        cat(k,"\n")
        #solveout <- mm.solve(phenotypes[,w,drop=F],X,Z,K,H,no.test = no.test)
        solveout <- try(mm.solve(phenotypes[,w,drop=F],X,Z,K,H,no.test = no.test),
                        silent = T)
        if (class(solveout) == "try-error") {
          solveout <- list(list(
            beta = NA,
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
      for(r in 2:length(res)) res[[r]]<-do.call(c,res[[r]])
      for(r in 1:length(res)) names(res[[r]])<-markers
      return(res)
    })
    # #Here we obtain a list of df, each df containing all results
    # parallel::stopCluster(cl)
  }
  #The results are given phenotype names
  #We obtain the structure we talked about
  names(result)<-make.names(colnames(phenotypes),unique=T)

  if (!is.null(permutation)) {
    pval <- sapply(result, function(x) x$pval) #extract pval for all phenoytpes
    minpval <- apply(pval[,(npheno+1):ncol(pval)],2, min,na.rm=T)
    minpvalMat <- matrix(minpval, ncol = npheno, byrow = T)
    maxlogpvalMat <- -log10(minpvalMat)
    thr <- sapply(1:ncol(maxlogpvalMat), function(i) {
      quantile(maxlogpvalMat[,i], probs = alpha) # alpha can be a vector
    })
    if (length(alpha)==1) thr <- t(thr)
    colnames(thr) <- colnames(phenotypes)[1:npheno]
    rownames(thr) <- alpha

    #gt: a new element called 'perm.thr' is added to each phenotype in result
    result <- lapply(1:npheno, function(i) {
      result[[i]]$perm.thr <- thr[,i]
      return(result[[i]])
    })
    names(result) <- colnames(phenotypes)[1:npheno]
  }

  return(result)
  }
