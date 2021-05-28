### Figures of MPP article -----------------
QTLpos<-c(453,1348,4023) #QTL positions
set.seed(50)
polyg<-sample(1:6910,50) #polygenic effect positions
herit <- c(0.2,0.5,0.8)
source("mpQTL_fun.R")

files<-paste0("../mpQTL v0/PedigreeSIM/NAM_crosses/",c(1,3,7,10),
               "_ancestral/cross",c(0,2,6,9),"00_founderalleles.dat")

phe<-list()
for(w in files){
  genotypes<-link_NAM(crossfile=w,
                      totallele="../mpQTL v0//PedigreeSIM/Parents/Total_pop.txt")
  #main QTLs
  QTL<-genotypes[QTLpos,]
  effects.QTL<-multi.effects2(gen=QTL,
                              anc_alleles=matrix(0:399,nrow=40),
                              seed=7,
                              anc_sd=1)

  polygen<-genotypes[polyg,]
  set.seed(7)
  effects.polyg<-apply(polygen,1,function(gen){
    e<-rnorm(length(unique(gen)),0,0.1) #mu should be zero, sd depends on main effects (anc_range).
    names(e)<-unique(gen)
    return(e)
  })

  phenotypes<-lapply(herit,function(h2){
      pheno2(genotypes = rbind(QTL,polygen),
             effects = c(effects.QTL,effects.polyg),
             ploidy=4,
             herit = h2,
             mu = 100)
    })
  names(phenotypes)<-herit
  phe[[w]]<-phenotypes
}
names(phe)<-paste0("cross",c(0,2,6,9),"00")
saveRDS(phe,"phenotypes.RDS")

phe<-readRDS("phenotypes.RDS")
# Results 1. Pop Simulation ---------------

parents<-as.data.frame(data.table::fread("../mpQTL v0/PedigreeSIM/Parents/Total_pop.txt"))
sapply(paste0("A",0:9),function(name){
  AG<-parents[,grep(name,colnames(parents))]
  alleles<-apply(AG,1,function(x) length(unique(x)))
  return(summary(alleles))
})

NAM3<-as.data.frame(data.table::fread("../mpQTL v0/PedigreeSIM/NAM_crosses/3_ancestral/cross200_alleledose.dat"))
K<-calc.K(t(NAM3[,-1]))
pc<-prcomp(K)

#Define colours for parents/offsprings
id<-sapply(c(paste0("A",0:9),paste0("_",1:9)),grepl,colnames(NAM3)[-1])
id<-apply(id,2,function(i){
  if(!any(i)) return(NULL)
  else return(i)
})
id<-do.call(cbind,id)

AGs<-apply(id[1:10,],1,which)
# col<-colorspace::divergex_hcl(7,palette="Tropic")
col<-colorspace::sequential_hcl(length(unique(AGs[-1])),h=240)
col<-col[c(unique(AGs[-1]),AGs[-1])]
cols<-col[id%*%1:ncol(id)]



var.pc<-paste0("PC",1:length(pc$sdev),"  ",
               round(pc$sdev^2/sum(pc$sdev^2)*100,2),"% variance")
comp<-c(1,2)

png("Plots/NAM3 genetic structure.png",width = 800,height=800,res=130)
plot(pc$rotation[1:10,comp],
     pch=c("X",LETTERS[AGs[-1]]),
     col=cols[1:10],cex=2,
     main="NAM3 Genetic Structure",
     xlab=var.pc[comp[1]],ylab=var.pc[comp[2]])
points(pc$rotation[11:nrow(pc$rotation),comp],
       pch=1,col=cols[-1:-10],cex=3)
legend("topright",pch=LETTERS[unique(AGs)],
       legend=paste("AG",unique(AGs)),
       col=unique(cols[1:10]),
       cex=1.2,bty="n")
dev.off()


iK<-apply(K,1,function(k) (k-mean(k))/sd(k))
iK<-iK[,ncol(iK):1]
offspring <- 50
nparents <- 10
png("Plots/NAM3 heatmap.png",width = 800, height = 800, res=130)
par(mar=c(4,3,5,2),oma=c(0,3,0,0))
image(iK,axes=F,main="Heatmap of NAM3 K matrix",xlab="Crosses",
      col = colorspace::sequential_hcl(200,h=240,c=80))
mtext(side = 2,xpd=T,text = "Ancestral Groups",outer=T,line=1.8)
#Y axis
#Code to work out which offsprings are bound together
lab<-sapply(unique(AGs[-1]),function(A){
  mid<- sum(A==AGs[-1])*offspring/2
  before <- sum(A>AGs[-1])*offspring
  return(mid+before+nparents)
})
at1<-(ncol(K)-c(nparents/2,lab))/ncol(K)
labs <- c("P",paste0("A1xA",unique(AGs[-1])))
axis(2,at=at1,labels=labs,line = 0.3,col="white",las=2)

lab<-sapply(unique(AGs[-1]),function(A){
  mid<- sum(A==AGs[-1])*offspring
  before <- sum(A>AGs[-1])*offspring
  return(mid+before+nparents)
})
at2 <- (ncol(K)-c(1,nparents,lab))/ncol(K)
axis(2,at=at2,line=0.3,labels=NA)

#X axis
#Code to work out which corresponds to which cross
at1<-c(nparents/2,seq(nparents+offspring/2,ncol(K),offspring))/ncol(K)
axis(1,at=at1,line=0.3,
     labels=c("P",paste0("C",1:(nparents-1))),col="white")
at2<-c(0,nparents,seq(nparents,ncol(K),offspring))/ncol(K)
axis(1,at=at2,labels=NA,line=0.3)
dev.off()


# Results 2. Phenotype simulation -------------------------
phe <- readRDS("phenotypes.RDS")

AGs
crossid<-id[,-1:-3]%*%AGs[-1]
source("mpQTL_fun_clean.R")
boxplot(phe$cross200$`0.5`~crossid)
cor(phe$cross200$`0.8`,phe$cross200$`0.8`)
dev.off()

# Appendix 1. Phenotype variations
div<-readRDS("results_diversity_test.RDS")
str(div)

dev.off()
comp.QQ(pval.anc,coltype = "sequential",
        h=240,main="QQplot comparison")
comp.skyplot(pval.anc,map=map,
             coltype="qualitative",h=240,main="Manhattan plot comparison")
dev.off()

#Poster Figure. All models.
haplo <- readRDS("output_inferredHaps/haplo_res.rds")
pheno <- phe$`../mpQTL v0/PedigreeSIM/NAM_crosses/3_ancestral/cross200_founderalleles.dat`$`0.5`

hblock <- grepl("chr",rownames(haplo$cross200$genotypes))
geno <- haplo$cross200$genotypes[hblock,]
maph <- haplo$cross200$map[hblock,]
resh <- map.QTL(pheno,geno,4,map,K=T,permutation = "pop",nperm = 10,impute=F)
geno <- link_NAM(files[2],totallele="../mpQTL v0//PedigreeSIM/Parents/Total_pop.txt")
map <- data.table::fread("PedigreeSim/Potato.map")
resa <- map.QTL(pheno,geno,4,map,K=T,permutation = "pop",nperm = 10,impute=F)
geno <- data.table::fread("../mpQTL v0/PedigreeSIM/NAM_crosses/3_ancestral/cross200_alleledose.dat")
geno <- as.matrix(geno[,-1])
resb <- map.QTL(pheno,geno,4,map,K=T,permutation = "pop",nperm = 10,impute=F)

pvals <-cbind(resa[[1]]$pval,resb[[1]]$pval)
colnames(pvals) <- c("True IBD","SNP")

qp<-sapply(QTLpos,function(q){
  map$chromosome == map$chromosome[q] & map$position == map$position[q]
})
qtl.pval <- qp[,1] | qp[,2] | qp[,3]

colors = colorspace::hcl_palettes("YlGnBu",3)
colors = colorspace::sequential_hcl(palette="YlGnBu",9)[c(1,5,8)]


#First we calculate the positions based on cumulative cM
chroms<-unique(map$chromosome)
ch.length<-sapply(chroms,function(i) max(map$position[map$chromosome==i]))
space<-sum(ch.length)*0.1/(length(chroms)-1)
spaces<-sapply(1:length(ch.length),function(i) space*(i-1))
#absolute chromosome length (cumulative)
abschlen<-sapply(1:length(ch.length),function(i) sum(ch.length[1:i]))

abschlen<-c(0,abschlen)
ch.start<-abschlen[-length(abschlen)]+spaces
ch.end<-abschlen[-1]+spaces
abspos <- map$position+ch.start[map$chromosome]
absposh<-maph$position+ch.start[maph$chromosome]


jpeg("Plots/skyplot_model_comparison_bigger.jpg",
    width = 1980*2 ,height = 800*2 ,res=280)
skyplot(-log10(pvals[!qtl.pval,2]),
        map[!qtl.pval,],col=colors[1],ylim=c(0,12))
abline(v=abspos[QTLpos],lty=2)
abline(h=c(resb[[1]]$perm.thr,resa[[1]]$perm.thr,resh[[1]]$perm.thr),lty=2,
       lwd=2,col=colors)
skyplot(-log10(pvals[!qtl.pval,1]),map[!qtl.pval,],add=T,col=colors[2])
points(absposh,-log10(resh[[1]]$pval),cex=1.5,pch=21,bg=colors[3])
legend("topright",legend=c("True IBD","SNP","Haplo"),cex = 1.5,
       col=c(colors[1:2],"black"),pt.bg=colors[5],pch=c(19,19,21),bty="n")
dev.off()

QTLb <- (-log10(pvals[,2]) > resb[[1]]$perm.thr) & !qtl.pval
QTLa <- (-log10(pvals[,1]) > resa[[1]]$perm.thr) & !qtl.pval
QTLh <- which(-log10(resh[[1]]$pval) > resh[[1]]$perm.thr)

map[QTLpos,]
View(map[QTLa,])
View(map[QTLb,])
View(maph[QTLh,])

pos.calc <- function(map,space=0.1,ch.length=NULL){
  #First we calculate the positions based on cumulative cM
  if(is.null(ch.length)){
    chroms<-unique(map$chromosome)
    ch.length<-sapply(chroms,function(i) max(map$position[map$chromosome==i]))
  }
  space<-sum(ch.length)*space/(length(chroms)-1)
  spaces<-sapply(1:length(ch.length),function(i) space*(i-1))
  #absolute chromosome length (cumulative)
  abschlen<-sapply(1:length(ch.length),function(i) sum(ch.length[1:i]))

  abschlen<-c(0,abschlen)
  ch.start<-abschlen[-length(abschlen)]+spaces
  ch.end<-abschlen[-1]+spaces

  abspos<-map$position+ch.start[map$chromosome]
}
