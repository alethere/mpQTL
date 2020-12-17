


# Consistency of input matrices: phenotypes, genotypes, map, cofactor, K, Q.
# For the wrapper and for the functions that can be used outside it.

# Currently, mpQTL accepts input data without rownames an colnames for marker
# names and individual names. Those matrices are expected to be in the right
# order:
# - same individuals order in pheno, geno, cofactor, K, Q;
# - same marker order in geno and map.


# inputCheck_dos is just checking the right data format for the genotypes





# example data ------------------------
## pheno
pheno <- matrix(rnorm(20,0,1), ncol = 2)
rownames(pheno) <- paste0("ind",1:10)
colnames(pheno) <- c("pheno1","pheno2")
pheno
pheno[1] <- NA

## geno
geno <- matrix(sample(0:4, 20*12, replace = T), ncol=12)
rownames(geno) <- paste0("snp",1:20)
colnames(geno) <- paste0("ind",1:12)
geno

## map
map <- data.frame(markers = paste0("snp",c(1:10,16:35)),
                  chrom = rep(1:3, each=10),
                  position = rep(seq(0,60, length.out=10), times=3))
map

## cof (if provided)
cof <- matrix(rnorm(20,0,1), ncol = 2)
rownames(cof) <- paste0("ind",1:10)
colnames(cof) <- c("cof1","cof2")
cof

## K (if provided)
K <- matrix(1, nrow = ncol(geno), ncol = ncol(geno))
K
colnames(K) <- colnames(geno)

## Q
Q <- c(rep(1:5, each=2))
names(Q) <- rownames(pheno)
Q

## Z
# diag(nrow(pheno))




# map.QTL input ----------------------

#There are four options
#1) K=NULL, no K is used. We go linear.
#2) K=T, K is calculated using homogeneously distributed markers along a map.
#   using sample.cM and calc.K
#3) K=K matrix, check dimensions, message is printed and K is used.

# create some possible vaues for argument K
Klist <- list(K1 = NULL,
              K2 = T,
              K3 = F,
              K4 = c(1,1,1,1),
              K5 = matrix(1, nrow = ncol(geno), ncol = ncol(geno)), #identity matrix,
              K6 = matrix(rnorm(ncol(geno)^2,0.3,1),
                          nrow = ncol(geno), ncol = ncol(geno)) #a random kinship matrix
)

inputCheck_K_test <- function(K) {

  # check that the argument K is either
  # 1) NULL or
  # 2) TRUE or
  # 3) a numeric square matrix

  # for testing I substituted stop with cat and return

  if (is.null(K)) {
    return(list(type=1,K=K))
  } else {
    if (length(K)==1 && is.logical(K)) {
      if (K==TRUE) {
        return(list(type=2,K=K))
      } else cat("K must be either NULL, TRUE or a numeric square matrix\n"); return()
            # stop("K must be either NULL or TRUE or a numeric square matrix")
    } else {
      d <- dim(K)
      if (length(K)>1 && is.null(d)) {
        cat("K must be either NULL, TRUE or a numeric square matrix\n"); return()
        # stop("K must be either NULL or TRUE or a numeric square matrix")
      } else {
        if (!is.matrix(K)) {K <- as.matrix(K); warning("K coerced to class matrix\n",
                                                       immediate. = T)}
        if (!is.numeric(K)) {cat("K matrix is not numeric\n"); return()}  #stop("K matrix is not numeric")
        if (d[1]!=d[2]) {cat("K is not a square matrix (i.e. nrow(K)!=ncol(K))\n"); return()}   #stop("K is not a square matrix (i.e. nrow(K)!=ncol(K))")
        return(list(type=3,K=K))
      }
    }
  }
}

# test
for(i in 1:length(Klist)) {
  cat(paste0(names(Klist)[i],"\n"))
  print(inputCheck_K_test(K=Klist[[i]]))
}



inputCheck_Q_test <- function(Q) {

  # check that the argument Q is either
  # 1) NULL or
  # 2) TRUE or
  # 3) a vector
  # 4) a numeric matrix

  # for testing I substituted stop with cat and return

  if (is.null(Q)) {
    return(list(type=1,Q=Q))
  } else {
    if (length(Q)==1 && is.logical(Q)) {
      if (Q==TRUE) {
        return(list(type=2,Q=Q))
      } else {cat("Q must be either NULL, TRUE, a membership vector or a design matrix\n"); return()} #stop("Q must be either NULL, TRUE, a membership vector or a design matrix\n")
    } else {
      if (length(Q)>1 && is.vector(Q)) {
        return(list(type=3,Q=Q))
      } else if (length(Q)>1 && length(dim(Q))==2) {
        if (!is.matrix(Q)) {Q <- as.matrix(Q); warning("Q coerced to class matrix\n",
                                                       immediate. = T)}
        if (nrow(Q)==1 | ncol(Q)==1) {
          dim(Q) <- NULL
          warning("Q is a column matrix or a row matrix, coerced to vector\n",
                  immediate. = T)
          return(list(type=3,Q=Q))
        } else {
          if (!is.numeric(Q)) {cat("Q design matrix is not numeric"); return()} #stop("Q design matrix is not numeric")
          return(list(type=4,Q=Q))
        }
      }
    }
  }
}


# create some possible vaues for argument Q
Qlist <- list(Q1 = NULL,
              Q2 = T,
              Q3 = F,
              Q4 = c(1,1,1,1), #num vector,
              Q5 = letters[1:5], #char vector
              Q6 = matrix(c(1,1,1,1), ncol = 1), #num one-column matrix,
              Q7 = matrix(rnorm(20), ncol = 4) #a numeric matrix
)

# test
for(i in 1:length(Qlist)) {
  cat(paste0(names(Qlist)[i],"\n"))
  print(inputCheck_Q_test(Q=Qlist[[i]]))
}

length(dim(Qlist$Q6))
h <- data.frame(1,2,3,4)
dim(as.matrix(h))
dim(h) <- NULL


# Order input matrices -------------------------------

# test
newinput <- inputOrder(geno, pheno=pheno, map=map,
                       cof=cof, Q=Q, K=K)
newinput

newinput2 <- inputOrder(geno, map=map,
                       cof=cof, Q=Q, K=K)
newinput2


