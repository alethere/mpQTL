#### QTL map dev ####

# Script to test and evaluate the different points of QTL calculation
source("mpQTL_fun.R")

### 1 Phenotyping ###
cfiles<-paste0("../mpQTL/PedigreeSIM/NAM_crosses/",c(1,3,7,10),
               "_ancestral/cross",c(0,2,6,9),"00_founderalleles.dat")
dfiles<-paste0("../mpQTL/PedigreeSIM/NAM_crosses/",c(1,3,7,10),
               "_ancestral/cross",c(0,2,6,9),"00_alleledose.dat")


data<-lapply(1:length(cfiles),function(w){
  #Read the founder genotypes
  genotypes<-link_NAM(crossfile=cfiles[w],
                      totallele="../mpQTL/PedigreeSIM/Parents/Total_pop.txt")
  dosage<-data.table::fread(dfiles[w])[,-1]
  rownames(dosage)<-rownames(genotypes)
  
  QTL<-genotypes[c(1348,4023,892),]
  #Obtain a set of QTL effects
  effects.QTL <- multi.effects2(
    gen = QTL,
    anc_alleles = matrix(0:399, nrow = 40),
    seed = 7,
    anc_sd = 100
  )
  
  set.seed(7)
  polygenic<-sample(1:nrow(genotypes),50)
  polygen<-genotypes[polygenic,]
  
  effects.polyg<-apply(polygen,1,function(gen){
    e<-rnorm(length(unique(gen)),0,10) #mu should be zero, sd depends on main effects (anc_range).
    names(e)<-unique(gen)
    return(e)
  })
  
  phenotypes<-pheno2(genotypes = QTL,
                     effects = effects.QTL,
                     polygen=polygen,
                     effects.polygen=effects.polyg,
                     ploidy=4,
                     herit = 0.8,
                     mu = 100,
                     return.effects = T)
  
  return(list(geno=genotypes, #ancestral allele genotype
              dosage=dosage, #SNP marker genotypes
              pheno=phenotypes[[1]], #phenotypes
              pos=list(QTL=c(1348,4023,892), #position of genetic effects
                       polyg=polygenic), 
              effects=list(QTL=phenotypes[[2]], #size of genetic effects
                           polyg=effects.polyg)
              ))
  
})
names(data)<-paste0("A",c(1,3,7,10))
saveRDS(data,"data_QTL_3.RDS")
# data<-readRDS("data_QTL_dev.RDS")

pops<-c(0:9,rep(1:9,each=50))
for(i in 1:length(data)){
  boxplot(as.numeric(data[[i]]$pheno[-1:-10])~pops[-1:-10],
          ylim=c(min(as.numeric(data[[i]]$pheno)),max(as.numeric(data[[i]]$pheno))),
          xlab="Cross index",ylab="phenotypes",main=paste("Cross",names(data)[i]))
  points(pops[2:10],data[[i]]$pheno[2:10],cex=1.5,pch=19)
  abline(h=data[[i]]$pheno[1],lty=2)
  text(9.5,data[[i]]$pheno[1],"central parent",xpd=T,cex=0.7)  
}

#### 2 QTL calculation
#In here, the different processes contained within map.QTL function are 
#disentangled to evaluate them one by one

#Arguments

map=data.table::fread("../mpQTL/PedigreeSIM/Potato/Potato.map")
K=T #distance matrix
Q=c(0:9,rep(1:9,each=50)) #population effect matrix
Q=c(rep(0,10),rep(1,100),rep(2,150),rep(3,200)) #ancestral matrix effect
Z=NULL
cofactor=NULL
cofactor.type=NULL
cM=1 #centimorgans to sample the K matrix
k=3 #number of axis used for pco decomposition
no_cores=parallel::detectCores()-1
P3D=T
EMMAX=T

genotypes<-as.matrix(data.table::fread(crossdosage[i]))
markers<-result[,1]; 
genotypes<-result[,-1] #take out marker column
rownames(genotypes)<-markers
genotypes<-bi.geno(dosage)

for(z in 1:length(data)){
  phenotypes<-data[[z]]$pheno
  phenotypes<-(phenotypes-mean(phenotypes))/sd(phenotypes)
  genotypes<-data[[z]]$geno
  
  dosage<-data[[z]]$dosage
  genotypes<-bi.geno(dosage)
  K=T #distance matrix
  Q=c(rep(0,10),rep(1,100),rep(2,150),rep(3,200))
  Z=NULL
  
  # ##### TEST SUBSETTING ####
  # sel<-c(1,4,5,6,111:260)
  # selc<-sapply(sel,function(i) (i-1)*4+1:4 )
  # selc<-as.vector(selc) 
  # phenotypes<-phenotypes[sel]
  # genotypes<-genotypes[,selc] #because the genotype matrix is per chromosome
  # dosage<-dosage[,..sel]
  # Q<-Q[sel]
  # ##########################
  
  if(is.vector(phenotypes)) {
    inds<-names(phenotypes)
    phenotypes<-matrix(phenotypes,ncol=1)
    rownames(phenotypes)<-inds
  }
  
  #DEFINITION OF K (kinship matrix)
  if(all(K==T)){ #if K is T, we calculate K distance matrix
    
    if(is.null(dosage)){
      if(all(rownames(phenotypes)==colnames(genotypes))){
        dosage<-genotypes
      }else{
        stop("If genotypes are not marker dosages, dosages must be specified") 
      }
    }
    
    #first we sample homogeneously markers along the genome
    K<-sample.cM(dosage,map,cM = cM)
    #then we calculate the distance
    K<-calc.K(t(K))
    
  }else if(!all(ncol(K)==length(phenotypes),nrow(K)==length(phenotypes))){
    #in case K is already defined, we check dimensions
    stop("K matrix does not have adequate dimensions")
  }
  
  #DEFINITION OF Q (population structure matrix)
  if(!is.null(Q)){
    if(all(Q==T)){
      if(is.null(K)){
        if(all(rownames(phenotypes)==colnames(genotypes))){
          dosage<-genotypes
        }else if(is.null(dosage)){
          stop("If genotypes are not marker dosage, dosage must be specified")
        }
        
        #first we sample homogeneously markers along the genome
        K<-sample.cM(dosage,map,cM = cM)
        #then we calculate the distance
        K<-calc.K(t(K))
      }
      Q <- cmdscale(1-K, k=k, eig = F, add = FALSE, x.ret = FALSE)
      
    }else if(is.vector(Q)){#if Q is a vector identifying populations, a Q matrix is designed
      Q<-Q.mat(Q)
    }else{ stop("Wrong Q specification")}
  }
  
  #DEFINITION OF C (cofactor matrix)
  if (!is.null(cofactor)){
    if(is.vector(cofactor)) cofactors<-matrix(cofactor,ncol=1)
    
    if(is.null(cofactor.type)|length(cofactor.type)!=ncol(cofactor)){
      stop("All cofactor types must be specified (by column if more than one) 
           either as \"numerical\" or \"categorical\"")
    }
    
    C<-lapply(1:ncol(cofactor),function(i){
      if(!is.na(pmatch(cofactor.type[i],"numerical"))){
        result<-as.numeric(cofactor[,i])
      }else if(!is.na(pmatch(cofactor.type[i],"categorical"))){
        result<-sapply(unique(cofactor[,i]), function(c)
          #+1 -1 transforms logical values into numerical values (more useful than as.numeric)
          cofactor[,i] == c) + 1 - 1
        result<-result[,-ncol(result)]
      }
    })
    C<-do.call(cbind,C)
    C<-apply(C,2,function(x) (x-mean(x))/sd(x))
    
  }else{C<-NULL}
  
  
  #Calculation of mixed models
  if(is.null(Z))  Z<-diag(nrow(phenotypes)) 
  
  if(!P3D|!EMMAX){
    #if we don't want to use P3D approximation
    Hinv<-NULL
  }else{
    Hinv<-calc.Hinv(phenotypes,
                    #X=vector of 1s, to apply P3D/EMMAX algorithm #gt: check#
                    X=X,
                    Z,K,lambda=T) 
    lambda<-Hinv$rr.parameter
    Hinv<-Hinv[1]
  }

  cl<-parallel::makeCluster(no_cores)
  export<-c("phenotypes","Z","K","Hinv","genotypes",
            "mm.solve","dosage.X","Q","C","lambda") #gt
  
  #extra functions need to be exproted 
  #if we don't want to use P3D approximation
  if(!P3D|!EMMAX) export<-c(export,"calc.Hinv","comp.vec")
  
  parallel::clusterExport(cl,export,
                          envir=environment()) #I dunno why, but without this it doesnt work
  # gt: same reason explained before. Here, doesn't even work because Z, K and Hinv exist only
  # in this executing environment (and not in the global one).
  # gt version: markers in parallel (parallelization works even with one phenotype)
  
  result<-lapply(1:ncol(phenotypes),function(w){#NOT parallely over each phenotype
    
    #parallely over markers, calculate the pvalue for each marker
    parallel::parLapply(cl,1:nrow(genotypes),function(k){      
      
      #DEFINITION OF X (matrix of fixed effects)
      X<-dosage.X(as.matrix(genotypes[k,]),ploidy=4,normalize=F) 
      if(ncol(X)>2){X<-X[,-2]} #when ancestral/parental model we need to prevent singularity
      nparX<-ncol(X) #total number of genetic parameters 
      X<-cbind(Q,C,X) #add population and cofactor parameters
      X<-matrix(as.numeric(X),ncol=ncol(X))
      no.test<-ncol(X)-nparX #number of non genetic parameters
      
      #If not P3D/EMMAX, then Hinv will be NULL (and cannot be subindexed)
      if(!is.null(Hinv)){
        H<-Hinv[[w]]
      }else{
        H<-NULL
      }
      
      res<-mm.solve(phenotypes[,w],X,Z,K,H,no.test = no.test,lambda = lambda)
    })
  })
  
  parallel::stopCluster(cl)
  names(result[[1]])<-rownames(genotypes)
  
  data[[z]]$result<-result[[1]]
}


plot(0, type = "n", 
     xlim = c(0, ncol(X) * space),
     ylim = c(min(y), max(y)),
     main=paste(names(data)[z],"QTL=",QTL[i]))

col <- viridis::viridis(11)
col2<- viridis::magma(12)
for (i in 1:ncol(X)) {
  x <- X[, i] / 4 + (i - 1) * space
  no<-X[,i]!=0
  points(x[no],y[no],
         pch = 19,
         cex = 0.7,
         col = col[an[i]+1])
  # rect(min(x),min(y),
  #      min(x)+space,
  #      min(y)+0.05*(max(y)-min(y)),
  #      col=col2[an[i]+1])

}

pval<-sapply(1:nrow(genotypes),function(k){
  X<-dosage.X(genotypes[k,],ploidy=4,F)
  test<-lm(phenotypes~X)

  R<-sum(test$residuals^2)/test$df.residual
  M<-sum(test$fitted.values^2)/(ncol(X)-1)
  p<-pf(M/R,ncol(X)-1,test$df.residual,lower.tail = F)
  return(p)
})
plot(-log10(pval))


#### Phenotype Variations ####
# In this section we will generate phenotypes using different
# parameters and then evaluate its effect on QTL detection
A3<-readRDS("data_A3.RDS")

#Positions of QTLs and polygenic effects
QTL<-c(1348,4023,892)
set.seed(7)
polyg<-sample(1:nrow(A3$geno),50)

#Genotypes (per chromosome) of the QTLs and polygenic effects
#the genotypes at the QTL positions
Qg<-A3$geno[QTL,] 
pg<-A3$geno[polyg,]
#Genetic map
map<-data.table::fread("../mpQTL/PedigreeSIM/Potato/Potato.map")

#mu should be zero, sd depends on main effects (less than 1 is better).
effects.polyg<-apply(pg,1,function(gen){
  e<-rnorm(length(unique(gen)),0,0.01) 
  names(e)<-unique(gen)
  return(e)
})

params<-c(0.01,0.1,1)
phenos<-lapply(1:length(params),function(w){
  #Obtain a set of QTL effects
  effects.QTL <- multi.effects2(
    gen = Qg,
    anc_alleles = matrix(0:399, nrow = 40),
    seed = 7,
    anc_sd = params[w]
  )
  
  phenotypes<-pheno2(genotypes = Qg,
                     effects = effects.QTL,
                     polygen=pg,
                     polygen.effects = effects.polyg,
                     ploidy=4,
                     herit = 0.8,
                     mu = 100,
                     return.effects = T)
  return(phenotypes)

  
})
names(phenos)<-params

phsum(Qg,phenos$`100`$effects,ploidy = 4,table=T)
phsum(pg,effects.polyg,ploidy = 4,table=T)

res<-lapply(1:length(phenos),function(i){
  p<-phenos[[i]]$pheno
  r<-map.QTL(phenotypes = p,
          genotypes = A3$geno,
          map = map,
          dosage = A3$dosage,
          K = T,Q = T )
  return(r)
})

cols<-colorspace::sequential_hcl(5,h=240,alpha = 0.7)
plot(0,type="n",
     xlim=c(0,length(res[[1]])),
     ylim=c(0,max(-log10(unlist(res)))),
     main="Manplot",xlab="markers",ylab="-log10(pval)"
     )
abline(v=QTL)
abline(v=polyg,col="green",lty=2)
for( i in 1:3){
  points(-log10(res[[i]]),pch=19,cex=0.5,main=paste("Manplot",3),col=cols[i])
}

colQQplot(do.call(cbind,res))
r
plot(-log10(r))
sapply(effects.polyg,max)*4
sapply(phenos$`1`$effects,max)

for(j in 1:3){
  r<-res[1:3+(j-1)*3]
  par(mfrow=c(1,2))
  plot(0,type="n",
       xlim=c(0,length(r[[1]])),
       ylim=c(0,max(-log10(unlist(r)))),
       main="Manplot",xlab="markers",ylab="-log10(pval)"
  )
  abline(v=QTL)
  abline(v=polyg,col="grey",lty=2)
  for(i in 1:3){
    points(-log10(r[[i]]),pch=19,cex=0.5,main=paste("Manplot",3),col=cols[i])
  }
  comp.QQplot(do.call(cbind,r),coltype = "sequential",legend=paste("Div",1:3))
  dev.off()
}


