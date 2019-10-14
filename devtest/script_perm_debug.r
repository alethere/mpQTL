
source("../mpQTL/R/mpQTL_fun_clean.R")



# data ----------------------------------------
set.seed(3)
phenotypes <- matrix(c(sample(c(5:10,NA), 200, replace = T),
                  sample(c(5:10), 200, replace = T),
                  sample(c(15:100), 200, replace = T)), ncol = 3)
phenotypes <- cbind(phenotypes, phenotypes, phenotypes)


set.seed(5)
genotypes <- matrix(sample(0:4, 5000*200, replace = T), ncol = 200)


map <- data.frame(marker=paste0("snp",1:5000),
                  chromosome=c(rep(1,3000),
                               rep(2,2000)),
                  position=c(seq(0,100, length.out = 3000),
                             seq(0,100, length.out = 2000)))

cof <- sample(c(1,2,3),50,replace=T,prob=c(0.30,0.30,0.4))
cof <- cbind(cof,runif(50))
cof <- cbind(cof,cof,cof)

cof
plot(-log10(res[,1]),ylim=c(0,10),xlim=c(0,500),type="n")
for(i in 1:ncol(res)){
  points(-log10(res[,i]),pch=i)
}



# permutation --------------------
res <- map.QTL(
  phenotypes = phenotypes, 
  genotypes = genotypes,  
  ploidy = 4,
  map = map, 
  K = T,
  cM = 2,
  no_cores = 6)

y = phenotypes
X = matrix(1,nrow=nrow(phenotypes),ncol=1)
Z = diag(nrow(phenotypes))
K <- sample.cM(genotypes,map,cM = cM)
K <- calc.K(t(K),ploidy = ploidy,haplotypes = haplo)

sapply(data,function(d) dim(d$y))

Hinv <- calc.Hinv(y,X,Z,K)
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

  Hinv <- lapply(1:length(data),function(i){
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

data <- test.compatibility(y,X,Z,K)

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
    res <- list(y=y,X=X,K=K,Z=Z,na.omit=nas)

    return(res)
  })
  
  #This indicates the NA combinations of each dataset
  #Eventually it will indicate the final order
  miss <- sapply(clean,function(d){
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
  # Clean NAs -----
  miss <- sapply(1:ncol(y),function(p){
    nay <- is.na(y[,p])
    nax <- apply(is.na(X),1,any)
    nas <- nax|nay
    nas <- paste(which(nas),collapse=" ")
    return(nas)
  })
  
  #This indicates the NA combinations of each dataset
  #Eventually it will indicate the final order
  
  #This matrix indicates how many combinations of
  #missing values we have. Columns = final list length
  #rows = number of phenotypes
  miss <- lapply(unique(miss),function(m) m==miss)
  miss <- do.call(cbind,miss)
  miss[miss] <- 1:sum(miss)
  reorder_vec <- rowSums(miss)
  return(reorder_vec)
}

Hinv <- Hinv[original_order(y,X)]
sapply(Hinv2,dim)
