#Creating a sample dataset for mpQTL package


map<-data.table::fread("PedigreeSim/Potato.map")
less.markers<-map$chromosome==1|map$chromosome==2

dosage<-data.table::fread("../mpQTL/PedigreeSIM/NAM_crosses/3_ancestral/cross200_alleledose.dat")
dosage<-as.matrix(dosage[less.markers,-1])
"../mpQTL/PedigreeSIM/NAM_crosses/3_ancestral/"
geno<-link_NAM(
  crossfile="../mpQTL/PedigreeSIM/NAM_crosses/3_ancestral/cross200_founderalleles.dat",
                    totallele="../mpQTL/PedigreeSIM/Parents/Total_pop.txt")
geno<-geno[less.markers,]

map.less<-map[less.markers,]


QTL<-matrix(geno[489,],nrow=1)
effects.QTL<-multi.effects2(gen=QTL,
                            anc_alleles=matrix(0:399,nrow=40),
                            seed=7,
                            anc_sd=1)

set.seed(7)
polygenic<-sample(1:nrow(geno),50)
polygen<-geno[polygenic,]

effects.polyg<-apply(polygen,1,function(gen){
  e<-rnorm(length(unique(gen)),0,0.1) #mu should be zero, sd depends on main effects (anc_range).
  names(e)<-unique(gen)
  return(e)
})

phenotypes<-pheno2(genotypes = rbind(QTL,polygen),
                   effects = c(effects.QTL,effects.polyg),
                   ploidy=4,
                   herit = 0.8,
                   mu = 100)
phenotypes2<-pheno2(genotypes = rbind(QTL,polygen),
                   effects = c(effects.QTL,effects.polyg),
                   ploidy=4,
                   herit = 0.5,
                   mu = 100)
phenotypes3<-pheno2(genotypes = rbind(QTL,polygen),
                    effects = c(effects.QTL,effects.polyg),
                    ploidy=4,
                    herit = 0.2,
                    mu = 100)

cofactors<-cbind(
  soil=c(rep("optimal",10),
         rep(rep(c("optimal","lowP","lowK","lowC","highNa"),each=10),9)
         ),
  pop=c("Europe","America","Asia")[c(c(1,1,1,2,2,2,3,3,3,3),
        rep(1,100),rep(2,150),rep(3,200))]
)

#5% of missing markers
genotype<-geno
na.matrix<-function(dosage,genotype,miss=0.05,ploidy=4){
  
  #How many markers should be missing?
  reps<-round(length(dosage)*miss)
  
  d<-dim(dosage)
  miss.geno<-genotype
  miss.dosage<-dosage
  missAB<-c()
  
  for(i in 1:reps){
    
    repeat{
      A<-sample(1:d[1],1)
      B<-sample(1:d[2],1)
      AB<-paste0(A,"-",B)
      if(!AB%in%missAB) break
    }
    missAB[i]<-AB
    
    miss.dosage[A,B]<-NA
    anc.B<-(B*ploidy-(ploidy-1)):(B*ploidy)
    miss.geno[A,anc.B]<-NA
  }
  
  return(list(na.geno=miss.geno,
              na.dosage=miss.dosage))
}


na005<-na.matrix(dosage,geno,0.05)
na0001<-na.matrix(dosage,geno,0.001)


data<-list(
  genotypes=geno,
  dosage=dosage,
  phenotype=cbind(h08=phenotypes,
                  h05=phenotypes2,
                  h02=phenotypes3),
  cofactors=cofactors,
  map=map.less,
  genotypes.na005=na005$na.geno,
  dosage.na005=na005$na.dosage,
  genotypes.na001=na0001$na.geno,
  dosage.na001=na0001$na.dosage,
  QTL.pos=489,
  polgen.pos=polygenic
)
saveRDS(data,"sample_data.RDS")
data<-readRDS("sample_data.RDS")

map<-data$map

result<-map.QTL(phenotypes = data$phenotype,genotypes = data$dosage,ploidy = 4,K=T,map=data$map)
resultNA01<-map.QTL(data$phenotype,data$genotypes.na001,ploidy=4,K=T,map=data$map)
resultNA01d<-map.QTL(data$phenotype,data$dosage.na001,ploidy = 4,K=T,map=data$map)
skyplot(resultNA01d$h08$pval,data$map)
pvals<-cbind(Ancestral=result$h08$pval,
             Ance_NA001=resultNA01$h08$pval,
             Dosage_NA001=resultNA01d$h08$pval)
comp.skyplot(pvals,data$map)
comp.QQ(pvals)

impgeno<-impute.knn(data$genotypes.na005,ploidy = 4,map=data$map)
miss<-as.vector(is.na(data$genotypes.na005))
imp<-impgeno[miss]
true<-data$genotypes[miss]

similarity<-sapply(1:(length(imp)/4),function(g){
  i<-(g*4-(4-1)):(g*4)
  sum(imp[i]%in%true[i])/4
})
hist(similarity)
