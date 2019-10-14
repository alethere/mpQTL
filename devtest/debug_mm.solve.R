

myf <- function(phenotypes,
                genotypes,
                dosage,
                map,
                ploidy=4,
                C=NULL,
                Q=NULL,
                w) {
  
  K<-sample.cM(dosage,map,cM = 1)
  K<-calc.K(t(K))
  Z<-diag(nrow(phenotypes)) 
  Hinv<-calc.Hinv(phenotypes,
                  #X=vector of 1s, to apply P3D/EMMAX algorithm #gt: check#
                  X=matrix(rep(1,nrow(phenotypes))),
                  Z,K)
  markers<-rownames(genotypes)
  # k=1
  pval <- sapply(1:nrow(genotypes), function(k) {
    #DEFINITION OF X (matrix of fixed effects)
    X<-dosage.X(as.matrix(genotypes[k,]),ploidy=ploidy,normalize=T) 
    if(ncol(X)>2){X<-X[,-2]} #when ancestral/parental model we need to prevent singularity
    nparX<-ncol(X) #total number of genetic parameters 
    X<-cbind(Q,C,X) #add population and cofactor parameters
    X<-matrix(as.numeric(X),ncol=ncol(X))
    no.test<-ncol(X)-nparX #number of non genetic parameters
    
    #If not P3D/EMMAX, then Hinv will be NULL (and cannot be subindexed)
    if(!is.null(Hinv)){
      H<-Hinv[w] #gt: H will be a matrix, but mm.solve is expecting a list!
      # Hinv[w] will be a list of length 1
    }else{
      H<-NULL
    }
    #g: phenotypes[,w,drop=F] to avoid conversion to vector
    b <- mm.solve(phenotypes[,w,drop=F],X,Z,K,H,no.test = no.test)
    b[[1]]$pval
    return(mm.solve(phenotypes[,w,drop=F],X,Z,K,H,no.test = no.test)[[1]]$pval)
  })
  # pval <- matrix(pval, ncol = 1)
  # colnames(pval)<-colnames(phenotypes)
  # rownames(pval)<-markers
  return(pval)
}
