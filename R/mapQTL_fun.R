# QTL mapping functions ---------------------

# This file includes functions for mapping QTLs using GWAS association model.
# The association model has been adapted to use in different ploidies, and can
# use SNP-dosage matrix or haplotype matrices to perform the GWAS analysis.
# The GWAS model is implemented both as a linear and mixed model, with the
# approximation known as P3D/EMMAX, which greatly reduces time usage.
# Additionally, the program can easily handle cofactors and various
# parameter specifications. Inspirations for this program include:
#   - Unified mixed model (Yu et al 2006)
#   - GWASpoly (Rosyara et al 2012)
#   - mppR (Vincent Garin 2019)
#
# Developed by Alejandro Therese Navarro & Giorgio Tumino
# October 2019

# Main wrapper -------------------------

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
#' estimate population stratification. If a vector specifying population of
#' each individual is passed, it will be used to construct a Q matrix. Vector
#' may contain numerical or character.
#' @param Qpco If Q = TRUE, integer indicating the number of Principal Coordinates to be
#' estimated and included. It is useful for detecting and correcting for the
#' main axes of genetic variation (i.e. population structure).
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
#'   imputed using \code{impute.knn}. Defaults to True, but False is recommended
#'   (imputation algorithm needs to be improved and with few missing values it
#'   has a small effect on QTL detection).
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
#' \dontrun{
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
#' }
#'
map.QTL <- function(phenotypes, genotypes, ploidy, map, K=NULL, Q=NULL, Z=NULL,
                    cofactor=NULL, cofactor.type=NULL, binsize=1, seed=NULL, Qpco=2,
                    no_cores= 1, approximate = TRUE,
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
      phenotypes[permid[,i],,drop=FALSE]
    }))
  }

  ####
  # Linear Model ------------
  #with no structure correction or with Q correction
  if(linear){
    #First set up cluster

    cluster <- parallel::makeCluster(no_cores)
    export <- c("phenotypes","genotypes","dosage.X","Q","C","lm_compare","haplo","ploidy")

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
        X <- dosage.X(genotypes[k,],ploidy = ploidy,haplotype = haplo,normalize = FALSE)
        if(haplo) X <- X[,-ncol(X)]
        #We create the naive and full matrices. The naive
        #model only includes Q, C and an intercept
        N <- cbind(matrix(1,nrow=nrow(X)),Q,C) #naive model
        X <- cbind(1,Q,C,X) #total model
        y <- phenotypes[,w,drop=FALSE]
        solveout <- try(lm_compare(X,N,y),silent = TRUE)

        ## mantain the usual output structure in case of error
        if (inherits(solveout, "try-error")) {
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

    cl<-parallel::makeCluster(no_cores)
    #atn: added object ploidy
    export<-c("phenotypes","Z","K","Hinv","genotypes","ploidy","haplo","npheno",
              "mm.solve","dosage.X","Q","C","test.compatibility","comp.vec","std_phe") #gt

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
        if (inherits(solveout, "try-error")) {
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
    parallel::stopCluster(cl)
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


# Calculation functions -------------------

#' Calculation of realized distance matrix (K)
#'
#' @description Using SNP genotypes (discrete or continuous) or haplotype
#' scores, a distance matrix is calculated
#' such that the average distance of an individual with itself is 1, and
#' the average with an unrelated individual is 0. Based on the "Realized
#' Relationship" model found in Rosyara et al. (2016).
#'
#' @param matrix A matrix of SNP genotypes (continuous or discrete) or
#' haplotypes (multi-allelic markers), where columns are markers and
#' rows may be either (1) individuals (for biallelic markers) or (2)
#' individual homologues (for multi-allelic markers, such as haplotypes).
#' @param haplotypes Logical, whether haplotypes are present in the matrix. If TRUE,
#' a matrix is expected to have p rows per individual, where p is ploidy.
#' @param ploidy integer, the ploidy number p of all individuals. Only used if
#' haplotypes = TRUE.
#'
#' @return A numeric matrix nxn where n is the number of rows.
#' @references
#' Rosyara U.R., De Jong W.S., Douches D.S., Endelman J.B. (2016). Software for
#' Genome-Wide Association Studies in Autopolyploids and Its Application to
#' Potato. *Plant Genome* 9:2. doi: 10.3835/plantgenome2015.08.0073
#' @export
#' @examples
#' ## Create SNP dosages for 10 tetraploid individuals (in rows) and 100 markers
#' snpdose <- matrix(sample(0:4,1000, replace=TRUE), nrow=10)
#'
#' ## Calculate K matrix
#' K <- calc.K(snpdose)
#'
#' ## Create an haplotype matrix for 5 tetraploid individuals (in rows) and 100 markers
#' ## Each individual genotype is represented by four rows (its four homologues)
#' hapdose <- matrix(sample(1:6,2000, replace=TRUE), nrow=20)
#'
#' ## K matrix
#' K <- calc.K(hapdose, haplotype = TRUE, ploidy = 4)
#'
#' ## Haplotypes can be named by characters as well
#' hapdose <- matrix(letters[hapdose], nrow=20)
#' K <- calc.K(hapdose, haplotype = TRUE, ploidy = 4)
#'
calc.K <- function(  matrix,  haplotypes=FALSE,  ploidy=NULL){
  #When haplotypes are given we can still perform K!
  if(haplotypes){
    if(is.null(ploidy))
      stop("For haplotype distance calculation ploidy must be defined.")
    #Create an ANOVA type matrix for all haplotypes
    matrix <-lapply(1:ncol(matrix),function(i){
      dosage.X(matrix[,i],haplotype = TRUE,ploidy=ploidy,normalize = FALSE)
    })
    matrix <- do.call(cbind,matrix)
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
#' @description Using a vector of values as input, it creates a normalized cofactor
#' design matrix (\eqn{Q}) that identifies each individual as belonging to a group.
#'
#' @param pop Vector where each element is a population identifier
#' @param names Optional. A vector of names for the \eqn{Q} matrix
#' @return A normalized matrix with as many columns as population groups
#' [ncol=unique(pop)] identifying each individual belonging to
#' one population.
#' @noRd
calc.Q<-function(pop, names=NULL
){
  match <- 1*sapply(unique(na.omit(pop)),function(x) x==pop)

  if(!is.null(names)){rownames(match)<-names}

  match<-apply(match,2,function(x) (x-mean(x))/sd(x))

  return(match[,-1]) #we take out the first column, to avoid singularity
}

#' X matrix calculator
#'
#' @description Using a vector of dosages per individual, or p haplotypes per
#' individual (where p is ploidy), it creates a design matrix X, with or without
#' normalization.
#'
#' @param genotypes vector of dosages per individual, or p haplotypes per individual
#' (where p is ploidy).
#' @param ploidy integer indicating ploidy.
#' @param haplotype logical, whether the vector contains haplotype or genotype data
#' @param normalize logical, should the X matrix be normalized?
#'
#' @return if haplotypes = FALSE, a matrix is returned whith a single column
#' containing the same values that were provided. If haplotypes = TRUE, a design matrix
#' is created with ncol = number of haplotypes and nrow = number of individuals.
#' @export
#' @keywords internal
#' @examples
#'
#' ## SNP dosages of 10 tetraploid individuals
#' snpdose <- sample(0:4, 10, replace = TRUE)
#' dosage.X(snpdose)
#'
#' ## Haplotypes of 10 tetraploid individuals for a locus with 5 haplotypes
#' hapdose <- sample(1:5,40, replace = TRUE)
#' dosage.X(genotypes = hapdose, haplotype = TRUE, ploidy = 4)
dosage.X <- function(genotypes, haplotype=FALSE, ploidy=NULL, normalize = FALSE ){

  if(!haplotype){
    alcount <- matrix(genotypes,ncol=1)
    if(normalize) alcount <- (alcount-mean(alcount))/sd(alcount)
    rownames(alcount) <- rownames(genotypes)
    colnames(alcount) <- "marker"
  }else{
    #we obtain the different alleles present
    rn <- rownames(genotypes)
    unals <- unique(unlist(genotypes))
    #we obtain a design matrix indicating the allele of each chromosome
    #atn: changed this function to add a NA column
    match <- sapply(unals, function(x) {
      if (!is.na(x)) {
        m <- x == genotypes
        m[is.na(m)] <- FALSE
      } else{
        m <- is.na(genotypes)
      }
      return(m)
    })
    colnames(match) <- unals
    rownames(match) <- rownames(genotypes)

    #we count the number of each allele for each individual
    #For haplotypes we need to combine ploidy columns into one count
    n <- length(genotypes)/ploidy
    alcount <- t(sapply(1:n, function(x){
      colSums(match[1:ploidy + (x - 1) * ploidy, ,drop=FALSE])
    }))

    #this line will not give the correct answer if we have a single individual
    if(nrow(alcount)==1) alcount <- t(alcount)

    if(normalize & ncol(alcount)!= 1) alcount <- apply(alcount,2,function(a) (a-mean(a))/sd(a))
    else if(normalize) alcount <- apply(alcount,2,function(a) (a-mean(a)))

    inds <- unique(substr(rownames(genotypes), 1, nchar(rownames(genotypes)) - 2))
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
#' @keywords internal
calc.Hinv<-function( y, X, Z, K=NULL, bounds=c(1e-09,1e+09), method="REML" ){

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
  Hinv <- unlist(Hinv,recursive = FALSE)
  Hinv <- Hinv[original_order(y,X)]
  return(Hinv)
}

#' Mixed Model Solver
#'
#' @description Mixed model solver. This function is partly based on the
#' \code{mixed.solve} function from the rrBLUP package (Endelman 2011).
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
#' @references
#' Endelman J.B. (2011). Ridge regression and other kernels for genomic
#' selection with r package rrBLUP. *Plant Genome* 4, 250-255.
#' doi: 10.3835/plantgenome2011.08.0024
#' @keywords internal
mm.solve <- function( y,  X, Z, K, Hinv = NULL, random = FALSE, no.test = 0){

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
      X2<-cbind(X[,1:no.test],Z%*%X[,-1:-no.test,drop=FALSE])#do not multiply cofactors
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
      V <- beta[(no.test +1):ncol(X),drop=FALSE]
      Fstat <- crossprod(V,Tt %*% V)/(p)

      x <- df/(df + p * Fstat)
      pval <- pbeta(x, df/2, p/2)

      # df<-df^2/(df^2+lambda)
      # pval<-pf(Fstat,p2,df,lower.tail=FALSE) #pvalue of genetic model?

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

  result<-unlist(result,recursive = FALSE)
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
#' pval (p-value of the difference between models) and se (standard error
#' of fixed effect estimates)
#' @keywords internal
lm_compare <- function(X,N,y){

  if(is.vector(y)) y <- matrix(y,ncol=1)

  #Process followed by lm ----
  not_na_y <- rowSums(is.na(y)) == 0
  not_na_X <- rowSums(is.na(X)) == 0
  not_na <- not_na_y & not_na_X
  X <- X[not_na,,drop=FALSE]
  N <- N[not_na,,drop=FALSE]
  y <- y[not_na,,drop=FALSE]
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
  pval <- pf(Ftest,dif_df,df_full,lower.tail = FALSE)
  se <- sqrt(SSR_full/(df_full))

  list(beta = beta, Ftest = Ftest,pval = pval,se = se)
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
#' This internal function is not intended for general use.
#'
#' @param y matrix of response y
#' @param X matrix of fixed effects X
#' @param K matrix of variance structure K
#' @param Z matrix of random effects Z
#' @return a list of dimension-compatible matrices without missing values.
#' @keywords internal
test.compatibility <- function(
  y,
  X,
  Z,
  K,
  Hinv=NULL
){
  # Check y -------
  if(is.vector(y)) y <- matrix(y,ncol=1)
  #Normalize y, i ncase it is not
  if(any(round(colMeans(y,na.rm=TRUE)) != 0)){
    y<-apply(y,2,function(p) (p-mean(p,na.rm=TRUE))/sd(p,na.rm=TRUE))
  }
  dimy <- dim(y)

  # Check X -------
  if(is.vector(X)) X <- matrix(X,ncol=1)
  #Normalize X (in case it is not)
  X <- apply(X,2,function(x){
    if(length(unique(x)) == 1) return(x)
    else return((x-mean(x,na.rm = TRUE))/sd(x,na.rm=TRUE))
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
    y <- y[!nas,p,drop=FALSE]
    X <- X[!nas,,drop=FALSE]
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
#' @description calc.Hinv calculates Hinv for every set of different dimensions
#' that are present in the data to improve speed (certain calculations
#' only depend on the dimensionality of the data). As a result, the Hinv
#' matrices are returned in an order that does not correspond the phenotype
#' order that is provided, causing incompatibilities later. This function
#' returns an index vector that can be used to reorder the Hinv matrix list.
#'
#' @param y phenotype matrix
#' @param X dosage matrix
#'
#' @return a vector that can be used to reorder the output of test.compatibility()
#' @keywords internal
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


#' Sample markers every bin of map unit (cM or Mb)
#'
#' @description Obtains a subset of markers homogeneously distributed across
#' a genetic map based on a map distance (in any unit). For instance, for if binsize=1,
#' will get 1 marker for every 1 unit window, if it can. Markers added on chromosome 0
#' are not returned as their position is not known (and thus they cannot be sampled per window).
#'
#' @param genotypes Either a character vector where each element is a marker
#' name, or a matrix of genotypes where each row is a marker.
#' @param map A data frame with 3 columns containing map information. The first
#' column specifies marker names, the second column chromosome names and the
#' third one marker position (any map unit).
#' Markers with chromosome number == 0 will not be taken into account
#' @param binsize Size of bin to sample a marker from.
#' @param seed An integer to set a seed for random number generation.
#'
#' @return A table as genotypes, with as many rows as can be obtained by sampling
#' markers every bin.
#' @export
#' @examples
#' ## Create a random map
#' map <- data.frame(marker = paste0("mrk", 1:100),
#'                   chromosome = c(rep(1,50), rep(2,50)),
#'                   position = c(sort(runif(50,0,60)), sort(runif(50,0,80))))
#'
#' ## Create a matrix of random genotypes, with markers in rows
#' snpdose <- matrix(sample(0:4,1000, replace = TRUE), nrow = 100)
#' rownames(snpdose) <- map$marker
#'
#' ## Sample evenly spaced markers
#' sample.marker(genotypes = snpdose, map, binsize=1, seed=3)
#'
#' ## A vector of marker names can be provided instead of genotypes
#' sample.marker(genotypes = map$marker,
#'               map, binsize=1, seed=3)
#'
sample.marker<-function(
  genotypes,
  map,
  binsize=1,
  seed=NULL
){
  if(is.vector(genotypes)) genotypes<-matrix(genotypes,ncol=1)
  if(nrow(genotypes)!=nrow(map)){
    stop("Number of markers in genotypes and in map does not coincide")
  }
  colnames(map) <- c("marker","chromosome","position")

  result <- lapply(unique(map$chromosome),function(x){
    pos <- map$position[map$chromosome == x]
    selpos <- seq(0,max(pos),by = binsize)

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

    return(as.matrix(genotypes[map$chromosome == x,,drop=FALSE][unique(markers),,drop=FALSE]))
  })

  #In case there are chromosomes labelled as 0, they should not
  #be used for centimorgan sampling

  if(any(map$chromosome == 0)){
    result <- result[-which(unique(map$chromosome) == 0)]
  }

  result <- do.call(rbind,result)
  return(result)
}

#' Standaridze phenotypes
#'
#' @description Given a matrix of numeric values, returns a list with
#' the mean and sd per column as well as standaridzed matrix
#'
#' @param phenotypes numeric matrix
#' @keywords internal
#' @noRd
stndrdz.pheno <- function(phenotypes) {
  if(is.vector(phenotypes)) {
    phenotypes <- matrix(phenotypes, ncol=1,
                         dimnames = list(names(phenotypes)))
  } else if (is.data.frame(phenotypes)) {
    phenotypes <- as.matrix(phenotypes)
  }

  means <- colMeans(phenotypes,na.rm = TRUE)
  sds <- apply(phenotypes,2,sd,na.rm=TRUE)

  phenotypes <- apply(phenotypes,2,function(p){
    (p-mean(p, na.rm = TRUE))/sd(p, na.rm = TRUE)
  })

  return(list(pheno = phenotypes,
              mean = means,
              sd = sds))
}



#' Check format of a numeric input matrix
#'
#' Used in map.QTL and thr.LiJi
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
#' @keywords internal
inputCheck_snp <- function(x, integer=TRUE, ploidy=NULL) {

  # conversion to matrix?
  if(!is.matrix(x)){
    cat("x is not a matrix.")
    if(is.data.frame(x)) {
      cat("I will try to convert it into a matrix.")
      # For converting to matrix I am using sapply(x, as.character) instead of
      # as.matrix(x). This is because as.matrix can introduce some extra
      # spaces (by calling format(x, trim=FALSE)). After conversion to numeric
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
              immediate. = FALSE)
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
      if(all(m.int, na.rm = TRUE)) {
        # x <- apply(x, 2, as.integer)  #less efficient
        storage.mode(x) <- "integer" # conversion to integer
      }
      # #gt: the following 'else' has been commented to allow continuous genotypes
      # else {
      #   stop(paste("The following non-integers are not allowed:\n",
      #              paste(u[!m.int], collapse = ","),
      #              "\n"))
      # }
    }

    # check the expected range of dosages
    if(is.integer(x) & !is.null(ploidy)) {
      u <- unique(c(x))
      u <- u[!is.na(u)]
      expdos <- u %in% 0:ploidy
      if(any(!expdos, na.rm = TRUE)) {
        stop(paste("The following dosages do not fit the ploidy level:\n",
                   paste(u[!expdos], collapse = ","),
                   "\n"))
      }
    }
  }
  x
}


#' Check input requirements for the argument K of map.QTL
#'
#' @description Check that the argument K is either 1) NULL, 2) TRUE or
#' 3) a numeric square matrix.
#'
#' @param K the input of the argument K of map.QTL
#' @return When requirements are met, it returns a list with:
#' * $type: number specifying the input type found: 1 for NULL, 2 for TRUE
#'  or 3 for a numeric square matrix;
#' * $K: the K object.
#' @keywords internal
inputCheck_K <- function(K) {

  if (is.null(K)) {
    return(list(type=1,K=K))
  } else {
    if (length(K)==1 && is.logical(K)) {
      if (K==TRUE) {
        return(list(type=2,K=K))
      } else stop("K must be either NULL, TRUE or a numeric square matrix")
    } else {
      d <- dim(K)
      if (length(K)>1 && is.null(d)) {
        stop("K must be either NULL, TRUE or a numeric square matrix")
      } else {
        if (!is.matrix(K)) {K <- as.matrix(K); warning("K coerced to class matrix\n",
                                                       immediate. = TRUE)}
        if (!is.numeric(K)) stop("K matrix is not numeric")
        if (d[1]!=d[2]) stop("K is not a square matrix (i.e. row number differs from column number)")
        return(list(type=3,K=K))
      }
    }
  }
}

#' Check input requirements for the argument Q of map.QTL
#'
#' @description Check that the argument Q is either 1) NULL, 2) TRUE,
#' 3) a vector (num or char) indicating sub-population membership or 4) a user
#' provided design matrix.
#'
#' @param Q the input of the argument Q of map.QTL
#' @return When requirements are met, it returns a list with:
#' * $type: number specifying the input type found: 1 for NULL, 2 for TRUE,
#'  3 for a population membership vector or 4 for a design matrix;
#' * $Q: the Q object.
#' @keywords internal
inputCheck_Q <- function(Q) {

  if (is.null(Q)) {
    return(list(type=1,Q=Q))
  } else {
    if (length(Q)==1 && is.logical(Q)) {
      if (Q==TRUE) {
        return(list(type=2,Q=Q))
      } else stop("Q must be either NULL, TRUE, a membership vector or a design matrix\n")
    } else {
      if (length(Q)>1 && is.vector(Q)) {
        return(list(type=3,Q=Q))
      } else if (length(Q)>1 && length(dim(Q))==2) {
        if (!is.matrix(Q)) {Q <- as.matrix(Q); warning("Q coerced to class matrix\n",
                                                       immediate. = TRUE)}
        if (nrow(Q)==1 | ncol(Q)==1) {
          dim(Q) <- NULL
          warning("Q is a column matrix or a row matrix, coerced to vector\n",
                  immediate. = TRUE)
          return(list(type=3,Q=Q))
        } else {
          if (!is.numeric(Q)) stop("Q design matrix is not numeric\n")
          return(list(type=4,Q=Q))
        }
      }
    }
  }
}


#' Subset and order input matrices based on individual names and marker names
#'
#' @description This function can be used to prepare the input matrices
#' or vectors for \code{map.QTL}. It will match individual names and marker
#' names of the input matrices, selecting the ones in common and putting them
#' in the same order. Individual names and marker names are required.
#'
#' @inheritParams map.QTL
#'
#' @return A list containing the re-ordered input matrices.???
#' @keywords internal
inputOrder <- function(geno, pheno=NULL, map=NULL,
                       cof=NULL, Q=NULL, K=NULL) {

  # geno must be provided, since is required both for individuals and
  # marker alignment

  if(!is.null(geno) & (is.null(rownames(geno)) | is.null(colnames(geno)))) {
    stop("Genotypes must have row names and column names")
  }


  consistent_input <- list(geno=NULL,
                           pheno=NULL,
                           map=NULL,
                           cofactor=NULL,
                           Q=NULL,
                           K=NULL)

  # individuals
  commonind <- NULL #to check wheter we order individuals
  if (!is.null(pheno)) {
    if(is.vector(pheno)) {
      pheno <- matrix(pheno, ncol=1,
                           dimnames = list(names(pheno)))
    } else if (is.data.frame(pheno)) {
      pheno <- as.matrix(pheno)
    }
    ## identify phenotyped and genotyped individuals
    notpheno <- apply(pheno, 1, function(x) all(is.na(x)))
    pheno <- pheno[!notpheno,,drop=FALSE]

    commonind <- intersect(rownames(pheno), colnames(geno))
    ungen_ind <- setdiff(rownames(pheno), colnames(geno))
    unphen_ind <- setdiff( colnames(geno), rownames(pheno))
    cat(length(commonind),"genotyped and phenotyped individuals will be used\n")
    # if (length(ungen_ind)>0) cat("Ungenotyped individuals excluded: ", ungen_ind)
    # if (length(unphen_ind)>0) cat("Unphenotyped individuals excluded: ", unphen_ind)
    ## subset pheno and geno accordingly
    consistent_input$pheno <- pheno[commonind,,drop=FALSE]
    consistent_input$geno <- geno[,commonind,drop=FALSE]
    ## subset cofactors
    if (!is.null(cof)) {
      consistent_input$cofactor <- cof[commonind,,drop=FALSE]
    }
    ## subset K
    if (!is.null(K) & length(K)>1) { #if K is a kinship matrix
      Ksub <- match(commonind, colnames(K))
      consistent_input$K <- K[Ksub,Ksub,drop=FALSE]
    }
    ## subset Q
    if (!is.null(Q) & length(Q)>1) {
      Qsub <- match(commonind, names(Q))
      consistent_input$Q <- Q[Qsub]
    }
  }


  # markers
  if (!is.null(map)) {

    ## Create a new map including all the markers in geno.
    ## For plotting reasons, unmapped markers are assigned to chrom
    ## 0 with fake postions.
    colnames(map) <- c("marker","chromosome","position")

    newmap <- merge(data.frame(marker = rownames(geno)), map,
                    by = "marker", all.x = TRUE, sort = FALSE)
    unmpd <- is.na(newmap$chromosome)
    newmap$chromosome[unmpd] <- 0
    avglen <- mean(tapply(newmap$position, newmap$chromosome, max), na.rm = TRUE)
    newmap$position[unmpd] <- seq(0, avglen, length.out = sum(unmpd))
    ## reorder map and geno
    consistent_input$map <- newmap[order(newmap$chromosome, newmap$position),, drop=FALSE]
    if (!is.null(commonind)) {
      consistent_input$geno <- geno[match(consistent_input$map$marker, rownames(geno)), commonind, drop=FALSE]
    } else {
      consistent_input$geno <- geno[match(consistent_input$map$marker, rownames(geno)), , drop=FALSE]
    }
  }
  return(consistent_input)
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
#' @export
#' @keywords internal
imputeNA <- function(m,
                     miss = NULL) {

  if(!is.null(miss) & is.numeric(miss)) m[m==miss] <- NA

  meanDose <- rowMeans(m, na.rm=TRUE)
  na <- is.na(m)
  means <- na * meanDose
  m[na] <- means[na]

  return(m)
}

#' KNN genotype imputator
#'
#' K-Nearest Neighbour imputation algorithms are useful for
#' imputing missing data. In this function, we use the kinship matrix
#' calculated with \code{\link{calc.K}} to determine the most similar
#' individuals, and substitute the missing values of one individual with
#' the most common dosage/haplotype of its nearest neighbours.
#' For continuous SNP genotypes, the nearest neighbours mean is used instead of
#' the most common dosage/haplotype.
#'
#' @param geno A matrix of SNP genotypes (continuous or discrete) or
#' haplotypes (multiallelic markers), where rows are markers and
#' columns may be either (1) individuals (for biallelic markers) or (2)
#' individual homologues (for multiallelic markers, such as haplotypes).
#' @param ploidy Numeric indicating the ploidy level.
#' @param map Optional. A data frame of three columns (marker, chromosome,
#' position) with map information. To be used in \code{\link{sample.marker}}
#' to sample a subset of evenly distributed markers.
#' @param kneighbors Number of nearest neighbours to use in the imputation. Default is 20.
#' @param K Optionally, provide the similarity matrix (K) directly.
#'
#' @return A complete matrix of haplotypes, SNP dosages or continuous SNP genotypes.
#' @export
#' @examples
#' \dontrun{
#' ## Get simulated genotypes for tetraploid individuals
#' data("mpsnpdose") # SNP dosages
#' data("mphapdose") # haplotypes
#' mpsnpdose[sample(1:length(mpsnpdose),500)] <- NA  # add NAs randomly
#' mphapdose[sample(1:length(mphapdose),500)] <- NA  # add NAs randomly
#'
#' ## Imputation
#' mpsnpdose <- impute.knn(mpsnpdose, ploidy = 4, kneighbors = 50)
#' mphapdose <- impute.knn(mphapdose, ploidy = 4, kneighbors = 50)
#' }
#'
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
    haplo<-FALSE
    contg<-FALSE
  #gt: add a condition for continuous genotypes
  }else if (all((unique(geno) >= 0 & unique(geno) <=1) | is.na(unique(geno)))){
    cat("Imputation will be performed on continuous genotypes.\n")
    haplo<-FALSE
    contg<-TRUE
  }else if (ncol(geno) %% ploidy == 0){
    cat("Imputation will be performed on haplotypes.\n")
    haplo<-TRUE
  }else stop(paste("Genotypes not recognized as dosage nor haplotypes.",
                   "Dosage values are not between 0 and ploidy or",
                   "ncol of genotypes is not a multiple of ploidy."))


  #Calculate K distance matrix. If possible using fewer markers
  #homogeneously distributed across the genome.
  if(is.null(K)){
    if (!is.null(map)) genoK <- sample.marker(geno, map)
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
      best <- order(K[, d][-d], decreasing = TRUE)[1:kneighbors]


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
                     "\nNon confident imputation. Try increasing the number of neighbours? Remove markers?")
        warning(error)
      }

      #Are they expressed numerically? The following function coerces to char
      asnum <- is.numeric(knei)
      #For each marker, select the most common dosage
      if (!contg) {
        nearest <- apply(knei, 1, function(i) {
          res<-names(sort(table(i), decreasing = TRUE)[1])
          if(asnum) res<-as.numeric(res)
          return(res)
        })
      } else {
        #gt: Alternatively, instead of the most common dosage, we could use the
        #mean. This could be more conservative in the case the second most
        #common dosage is not so unfrequent.
        #In any case, at the moment this is necessary for continuous genotypes.
        nearest <- apply(knei, 1, mean, na.rm=TRUE)
      }

      #the nearest neigbhours are put in. Careful, if too many missing markers
      #less neigbhours are being used to impute.
      imputed[is.na(imputed)] <- nearest
      return(imputed)
    })
    if (!is.null(colnames(K))) colnames(result) <- colnames(K) #gt: to keep individual names
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
      best <- order(K[-d, d], decreasing = TRUE)[1:kneighbors]
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
        common <- sort(table(alleles) / length(alleles), decreasing = TRUE)

        #We need to select the four most common alleles
        #This probably needs to be improved
        res <- c()
        for (i in 1:ploidy) {
          res[i] <- names(common)[1]
          common[1] <- common[1] - 1 / ploidy
          common <- sort(common, decreasing = TRUE)
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

        nearest <- t(apply(geno[miss,,drop=FALSE],1,function(alleles){  #gt: drop=FALSE added
          common <- sort(table(alleles) / length(alleles), decreasing = TRUE)
          #We need to select the four most common alleles
          #This probably needs to be improved
          res <- c()
          for (i in 1:ploidy) {
            res[i] <- names(common)[1]
            common[1] <- common[1] - 1 / ploidy
            common <- sort(common, decreasing = TRUE)
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


# Threshold calculations -------------------

# #' Pooled correlation between two multiallelic markers
# #'
# #' @description Used in \code{thr.LiJi}.
# #' \code{rpool} computes a pooled correlation coefficient between
# #' two sets of variables (e.g. two multiallelic markers) using the equation in
# #' Zhao et al. 2005, Genet. Res. 86:77-87.
# #' Correlation of each pair of alleles is weighted by allele frequency.
# #'
# #' @param A matrix of haplotype dosages of block A. Haplotypes in rows and
# #' individuals in columns.
# #' @param B matrix of haplotype dosages of block B.
# #' @param ploidy numeric indicating the ploidy level.
# #'
# #' @return pooled correlation
# #' @noRd
# rpool <- function(A,B,ploidy) {
#
#   A <- t(A)
#   B <- t(B)
#
#   nhom <- nrow(A)*ploidy  #n homologues
#   pA <- colSums(A)/nhom   #haplotype frequency
#   pB <- colSums(B)/nhom
#
#   ## for loop semi-vectorised
#   k <- ncol(A)  #n haplotypes
#   Srsq <- NULL
#   for (i in 1:k) {
#     Srsq <- sum(Srsq, pA[i] * pB * (cor(A[,i], B))^2)
#   }
#
#   ## alternative sapply loop (slower)
#   # combid <- expand.grid(1:ncol(A),1:ncol(B))
#   # Srsq <- sum(sapply(1:nrow(combid), function(x){
#   #   pA[combid[x,1]] * pB[combid[x,2]] * (cor(A[,combid[x,1]], B[,combid[x,2]]))^2
#   # }))
#
#   return(sqrt(Srsq))
# }


# #' Pairwise pooled correlation between multiallelic markers
# #'
# #' @description Used in \code{thr.LiJi}.
# #' \code{pairwise_rpool} generates a pairwise correlation matrix
# #' between a set of multiallelic markers. It is based on \code{rpool}.
# #'
# #' @param hapdos a list of matrices with allele dosages. Each matrix contains
# #' allele dosages of one multiallelic marker, with alleles in rows and
# #' individuals in columns.
# #' @param ploidy numeric indicating the ploidy level.
# #'
# #' @return a correlation matrix between multiallelic markers
# #' @noRd
# pairwise_rpool <- function(hapdos,
#                            ploidy) {
#
#   nblk <- length(hapdos)
#   prws <- combn(nblk,2)
#
#   rtr <- sapply(1:ncol(prws), function(x) {
#     rpool(A = hapdos[[prws[1,x]]], B = hapdos[[prws[2,x]]], ploidy=ploidy)
#   })
#
#   mcor <- matrix(1, ncol = nblk, nrow = nblk)
#   mcor[lower.tri(mcor)] <- rtr
#   mcor <- t(mcor)
#   mcor[lower.tri(mcor)] <- rtr
#   colnames(mcor) <- rownames(mcor) <- names(hapdos)
#   return(mcor)
# }


#' Calculate significance level based on the effective number of independent tests
#'
#' This function implements a method to adjust the significance
#' level for multiple testing based on the effective number of independent
#' tests, as proposed by \href{https://www.nature.com/articles/6800717?report=reader}{Li and Ji (2005, Heredity)}.
#' The Bonferroni adjusted threshold is obtained dividing a significance level
#' by the number of tests, assuming tests are independent to each other.
#' The Li and Ji method is more appropriate when the assumption of
#' independence is not met. In this method, the eigenvalues of a correlation
#' matrix are used to estimate the effective number of independent tests.
#'
#' @param m A matrix of SNP genotypes (discrete or continuous) with markers
#' in rows and individuals in columns.
#' @param chrom A vector defining to which chromosome markers belong.
#' @param alpha Significance level (use 0.05 for 5\% significance).
#' @param ploidy Numeric indicating the ploidy level.
#'
#' @return Value of adjusted alpha.
#' @export
#' @examples
#' ## Get example SNP dosages and a map
#' data("mpsnpdose")
#' data("mpmap")
#'
#' mpsnpdose[1:5,1:5]
#' head(mpmap)
#' tapply(mpmap$marker, mpmap$chromosome, length)
#' identical(rownames(mpsnpdose), mpmap$marker)
#'
#' ## excluding unmapped markers (chromosome 0)
#' mapped <- mpmap$chromosome != 0
#' thr <- thr.LiJi(m = mpsnpdose[mapped,],
#'                 chrom = mpmap$chromosome[mapped])
#' -log10(thr$threshold) #log scale
thr.LiJi <- function(m,
                     chrom,
                     alpha = 0.05,
                     ploidy = NULL) {

  # check input data format
  if (is.matrix(m) || inherits(m, "data.frame")) {
    m <- inputCheck_snp(m, integer=FALSE, ploidy=NULL)
  } else {
    # For the moment I disabled the method for haplotypes. The method used here
    # for calculating correlation between blocks (pooled correlation) needs
    # further checks/developments.
    # In case of haplotypes, m would be a list of matrices with haplotype
    # dosages.
    # To enable the method, 1) source rpool and pairwise_rpool,
    # 2) uncomment all the if(inherits(m, "list")) chunks of code.
    # See an alternative method in the Li and Ji paper, that is based on
    # joint entropy and mutual information.
    stop("m is not a matrix of SNP dosages")
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
  } #else if(inherits(m, "list")) { #for haplotype dosages
  #   NAblock <- sapply(m, function(x) sum(is.na(x[1,]))/ncol(x))
  #   if (any(NAblock > 0.25)) {
  #     stop(paste("There are blocks with more than 25% of missing values",
  #                "Remove them before proceeding."))
  #   } else if (any(NAblock > 0)) {
  #     message(paste("The level of missingness per block is acceptable (<25%).\n",
  #                   "Before calculating correlations, missing values will be",
  #                   "imputed using the average dosage per haplotype."))
  #     m <- lapply(m, imputeNA)
  #   }
  # }

  # check monomorphic markers
  # if (inherits(m, "list")) {
  #   mono <- sapply(m, function(n) {
  #     apply(n,1, function(x) sd(x, na.rm = TRUE)==0)
  #   })
  #   polym <- lapply(mono, function(x) which(!x))
  #   m <- lapply(1:length(m), function(x) m[[x]][polym[[x]],])
  #   names(m) <- names(polym)
  # } else
  if (is.matrix(m) || inherits(m, "data.frame")) {
    mono <- apply(m,1, function(x) sd(x, na.rm = TRUE)==0)
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
    }
    # else if(inherits(m, "list")) {
    #   cm <- pairwise_rpool(m[which(chrom == i)], ploidy=ploidy)  #pooled correlation between blocks
    # }

    ## step 2: estimate the effective number of independent tests (Meff)
    eval <- eigen(cm, symmetric=TRUE, only.values=TRUE)   # eigenvalues
    epos <- eval$values[eval$values>=0]     # select positive eigens only
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





# Miscellaneous -------------------

#' LD decay calculation
#'
#' @description Calculation of linkage disequilibrium per windows across
#' the genome can help paint a picture of the distribution of recombinations
#' in a population. This function, provided a set of dosages and a genetic
#' map, is able to obtain LD estimates for any percentile, by default
#' at 50%, 80% 90% and 95%. It also calculates the background LD (between unlinked
#' markers) and can be performed over the whole genome, or per chromosome.
#'
#' @param dos integer matrix with markers on rows and individuals in columns.
#' @param map data.frame genetic map with columns "chromosome" "marker" and "position"
#' @param win_size numeric indicating the window size to calculate LD estimates.
#' Must be in the same unit as the "position" column of map.
#' @param max_dist numeric, maximum distance to consider between markers. Usually, beyond
#' 50 cM the estimates are not informative.
#' @param per_chr logical, whether the estimations should be performed per chromosome.
#' @param percentile numeric vector, which percentiles to calculate?  By default
#' 0.5, 0.8, 0.9 and 0.95. Any value between 0 and 1 is allowed.
#'
#' @return If per_chr = FALSE, an LD object (list with LD estimates dataframe and
#' background LD for each percentile); if per_chr = TRUE a list of LD objects, one
#' per chromosome. Background LD represents the LD between unlinked markers. If per_chr = FALSE,
#' background LD is the LD between markers on different chromosomes. If per_chr = TRUE,
#' background LD is the LD between markers on the opposite ends of a chromosome, and thus
#' is different per each chromosome.
#' @export
#' @examples
#'
#' ## Get example SNP dosages and a map
#' data("mpsnpdose")
#' data("mpmap")
#'
#' mpsnpdose[1:5,1:5]
#' head(mpmap)
#' tapply(mpmap$marker, mpmap$chromosome, length)
#' identical(rownames(mpsnpdose), mpmap$marker)
#'
#' ## LD decay calculation
#' LD <- LD_decay(mpsnpdose, mpmap, win_size = 2,
#'                max_dist = 50, percentile = c(0.9,0.95))
#' LD <- LD_decay(mpsnpdose, mpmap, win_size = 2,
#'                per_chr = TRUE,
#'                max_dist = 50, percentile = c(0.9,0.95))
LD_decay <- function(
  dos,
  map,
  win_size = 0.1,
  max_dist = NULL,
  per_chr = FALSE,
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
    not_dup <- lower.tri(ld,diag=FALSE)
    return(ld[not_dup]^2)
  })

  #then we calculate the distance between markers
  pos_per_chr <- split(map$position,map$chromosome)
  dis <- lapply(pos_per_chr,function(d){
    res <- sapply(d,function(p) abs(p-d))
    not_dup <- lower.tri(res,diag=FALSE)
    return(res[not_dup])
  })

  if(!per_chr){

    #Per window we calculate the percentile estimates
    LD <- do.call(c,LD)
    dis <- do.call(c,dis)

    if(is.null(max_dist)) max_dist <- max(dis,na.rm=TRUE)

    windows <- seq(0,max_dist,win_size)
    ld_estimates <- t(sapply(seq_along(windows),function(w){
      sel <- dis >= windows[w] & dis < windows[w] + win_size
      return(quantile(LD[sel],percentile,na.rm = TRUE))
    }))

    #We calculate the background correlation between chromosomes
    #In this case, background LD represents LD between markers on different chromosomes
    back_ld <- lapply(1:length(dos_per_chr),function(chr){
      cors <- lapply(dos_per_chr[-chr],function(d) cor(t(d),t(dos_per_chr[[chr]]))^2 )
      do.call(c,cors)
    })
    back_ld <- do.call(c,back_ld)
    back_ld <- quantile(back_ld,percentile,na.rm = TRUE)

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

      windows <- seq(0,max(d,na.rm=TRUE),win_size)
      ld_estimates <- t(sapply(seq_along(windows),function(w){
        sel <- d > windows[w] & d < windows[w] + win_size
        return(quantile(ld[sel],percentile,na.rm=TRUE))
      }))

      #Background LD represents the LD between the ends of the chromosomes
      #More specifically, between markers at a distance of >90% of max distance
      back_ld <- quantile(ld[d > max(d,na.rm=TRUE)*0.9],percentile,na.rm=TRUE)

      #We add some extra features
      res <- list(LD = data.frame(ld_estimates,distance = windows),
                  background = back_ld)
      attr(res,"max_dist") <- max(d,na.rm=TRUE)
      attr(res,"class") <- c("LD","list")

      return(res)
    })
    names(res) <- names(LD)
  }

  return(res)
}



#' Vector comparison function
#'
#' @description Element per element comparison between two vectors.
#'
#' @param vec1 vector of values
#' @param vec2 vector of values
#'
#' @return Matrix of logical values indicating whether the values of vec1
#' and vec2 are equal. Columns correspond to elements
#' of vec1 and rows to elements of vec2.
#' @noRd
comp.vec <- function(vec1, vec2) {
  mat <- sapply(vec1, function(x)
    x == vec2)
  colnames(mat) <- vec1
  rownames(mat) <- vec2
  return(mat)
}
