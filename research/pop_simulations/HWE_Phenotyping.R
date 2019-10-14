
# Phenotype generation with realistic gene frequencies ---------------
# Alejandro Thérèse NAvarro
# July 2019

#' In this file we will generate phenotypes using the new cross files
#' that are based on HWE-realistic stuff
source("mpQTL_fun.R")
source("mpQTL_fun_clean.R")

cfiles<-paste0("PedigreeSim/HWE_NAM/",c(1,2,3),
               "_ancestral/cross0",c(0,1,2),"0_founderalleles.dat")

phe.poly<-list()
phe.non<-list()
for(w in cfiles){
    genotypes<-link_NAM(crossfile=w,
                        totallele="PedigreeSim/HWE_Parents/HWE_Totallele.txt")
    #main QTLs
    QTL <- genotypes[c(438,1348,4023),]
    effects.QTL <- multi.effects2(gen=QTL,
                                anc_alleles=matrix(0:119,nrow=12),
                                seed=7,
                                anc_sd=1)
    
    set.seed(7)
    polygenic<-sample(1:nrow(genotypes),50)
    polygen<-genotypes[polygenic,]
    
    effects.polyg <- apply(polygen,1,function(gen){
      e <- rnorm(length(unique(gen)),0,0.1) #mu should be zero, sd depends on main effects (anc_range).
      names(e)<-unique(gen)
      return(e)
    })
    
    phenotypes<-pheno2(genotypes = rbind(QTL,polygen),
                       effects = c(effects.QTL,effects.polyg),
                       ploidy=4,
                       herit = 0.8,
                       mu = 100)
    phe.poly[[w]]<-phenotypes
    phe.non[[w]]<-pheno2(genotypes = rbind(QTL),
                         effects = c(effects.QTL),
                         ploidy=4,
                         herit = 0.8,
                         mu = 100)

}


pops<-c(0:9,rep(1:9,each=50))
for(i in 1:3){
  boxplot(as.numeric(phe.poly[[i]][-1:-10])~pops[-1:-10],
          ylim=c(min(as.numeric(phe.poly[[i]])),max(as.numeric(phe.poly[[i]]))),
          xlab="Cross index",ylab="phenotypes",main="Cross phenotypes")
  points(pops[2:10],phe.poly[[i]][2:10],cex=1.5,pch=19)
  abline(h=phe.poly[[i]][1],lty=2)
  text(9.5,phe.poly[[i]][1],"central parent",xpd=T,cex=0.7)  
}

dfiles<-paste0("PedigreeSim/HWE_NAM/",c(1,2,3),
               "_ancestral/cross0",c(0,1,2),"0_alleledose.dat")
mapfile<-"PedigreeSim/Potato.map"
map<-data.table::fread(mapfile)
source("mpQTL_fun_clean.R")
pval<-lapply(1:3,function(i){
  genotypes<-link_NAM(cfiles[i],totallele="PedigreeSim/HWE_Parents/HWE_Totallele.txt")
  genotypes <- as.matrix(read.table(dfiles[i], header=T)[,-1])
  pval<-map.QTL(phenotypes = phe.poly[[i]],
                genotypes = genotypes,
                K=T,
                Q=T,
                map=map,
                ploidy= ploidy)
  
  pval2<-map.QTL(phenotypes = phe.non[[i]],
                 genotypes = genotypes,
                 K=T,
                 Q=T,k=5,
                 map=map,
                 ploidy = ploidy)
  return(list(polygen=pval,main=pval2))
})

