source("mpQTL_fun.R")
source("mpQTL_fun_clean.R")
phe <- readRDS("phenotypes.RDS")
gen <- data.table::fread("../mpQTL v0/PedigreeSIM/NAM_crosses/3_ancestral/cross200_genotypes.dat",
                  header = T)
dos <- data.table::fread("../mpQTL v0/PedigreeSIM/NAM_crosses/3_ancestral/cross200_alleledose.dat",
                         header = T)
dos <- dos[,-1]
gen <- gen[,-1]
anc <- link_NAM("../mpQTL v0/PedigreeSIM/NAM_crosses/3_ancestral/cross200_founderalleles.dat",
                totallele = "../mpQTL v0/PedigreeSIM/Parents/Total_pop.txt",ploidy = 4)
map <- read.table("PedigreeSim/Potato.map",header=T)
QTLpos <- c(453,1348,4023)

small_g <- sample.cM(gen,map,cM=0.5,seed = 8)
small_d <- sample.cM(dos,map,cM=0.5, seed = 8)
small_map <- data.frame(sample.cM(map,map,cM=0.5, seed = 8),
                           stringsAsFactors = F)
small_anc <- sample.cM(anc,map,cM = 0.5,seed = 8)

small_map$position <- as.numeric(small_map$position)
small_map$chromosome <- as.numeric(small_map$chromosome)
# Create haplotypes -------------------------------------------------------

#First we create the indeces of the SNPs we must concatenate
hap_indexer <- function(gen,l){
  t(sapply(1:floor(nrow(gen)/l),function(n) 1:l+(n-1)*l ))
}


#Then we calculate the position as the average of the position of each marker
hap_pos <- function(hap_index,map){
  position <- apply(hap_index,1,function(n) mean(map$position[n]))
  chromosome <- apply(hap_index,1,function(n) mean(map$chromosome[n]))
  marker <- paste0("hap",sprintf("%04.f",1:length(position)))
  return(data.frame(marker,position,chromosome))
}
#Now we have an updated haplotype map 

#Now we must create a set of haplotypes
happer <- function(gen,hap_index){
  t(apply(hap_index,1,function(h_i){
    Reduce(paste0,data.frame(t(gen[h_i,])))
  }))
  
}

#' Haplotype corrector
#' 
#' Sometimes haplotypes might be built accross chromosomes. In order
#' to avoid this (which results in some haplotypes belonging to 
#' chromosome 7.4), we can remove those haplotypes that bridge multiple
#' chromosomes. This will, at most, be as many haplotypes as chromosomes.
#'
#' @param hap_index 
#' @param hap_map 
#'
#' @return
#' @export
#'
#' @examples
hap_correct <- function(hap_index,hap_map){
  wrong <- hap_map$chromosome%%1 != 0
  return(list(hap_index = hap_index[!wrong,],
              map = hap_map[!wrong,]))
}

# Wrapper function that creates haplotypes based on a map, length and genotype
haplotyper <- function(gen,map,l){
  hap_i <- hap_indexer(gen,l)
  new_map <- hap_pos(hap_i,map)
  corrected <- hap_correct(hap_i,new_map)
  haps <- happer(gen,corrected$hap_index)
  return(list(haplotypes = haps,
              map = corrected$map))
}

gen <- small_g
l <- 4
map <- small_map

haplos4 <- haplotyper(gen,map,4)
haplos6 <- haplotyper(gen,map,6)
haplos8 <- haplotyper(gen,map,8)

small_hap <- haplotyper(small_g,small_map,4)

haplo_res <- map.QTL(phenotypes = data.frame(phe$cross200),
          genotypes = small_hap$haplotypes,
          map = small_hap$map,
          ploidy = 4,K = T)


dos_QTL <- map.QTL(phenotypes = data.frame(phe$cross200),
                     genotypes = small_d,
                     map = small_map,
                     ploidy = 4,K = T)

anc_res <- map.QTL(phenotypes = data.frame(phe$cross200),
                  genotypes = small_anc,
                  map = small_map,
                  ploidy = 4,K = T)

qp<-sapply(QTLpos,function(q){
  map$chromosome == map$chromosome[q] & map$position == map$position[q]
})
qtl.pval <- qp[,1] | qp[,2] | qp[,3]

colors = colorspace::sequential_hcl(palette="YlGnBu",5)[c(1,3,5)]

#First we calculate the positions based on cumulative cM 
map <- small_map
chroms<-unique(map$chromosome)
ch.length<-sapply(chroms,function(i) max(map$position[map$chromosome==i]))
space<-sum(ch.length)*0.1/(length(chroms)-1)
spaces<-sapply(1:length(ch.length),function(i) space*(i-1))
#absolute chromosome length (cumulative)
abschlen<-sapply(1:length(ch.length),function(i) sum(ch.length[1:i]))

abschlen<-c(0,abschlen)
ch.start<-abschlen[-length(abschlen)]+spaces
ch.end<-abschlen[-1]+spaces

absposh<-small_hap$map$position+ch.start[small_hap$map$chromosome]
abspos <- map$position+ch.start[map$chromosome]
thresh <- c(anc_res$X0.5$perm.thr,dos_QTL$X0.5$perm.thr,haplo_res$X0.5$perm.thr)
thresh <- thresh[c(2,1,3)]
for(pic in 1:3){
  jpeg(paste0("Plots/True_haplo",pic,".jpeg"),width = 2800,height = 1400,res=220)
  if(pic>0){
    skyplot(-log10(dos_QTL$X0.5$pval),
            map,col=colors[2],ylim=c(0,15))
     }
  if(pic>1){
    skyplot(-log10(anc_res$X0.5$pval),add=T,
            map,col=colors[1],ylim=c(0,15))
      }
  if(pic>2){
    points(absposh,-log10(haplo_res$X0.5$pval)
           ,cex=1.5,pch=21,bg=colors[3])
  }
  abline(v=abspos[QTLpos],lty=2)
  legend("topright",legend=c("True IBD","SNP","Haplo"),cex = 1.5,
         col=c(colors[1:2],"black"),pt.bg=colors[3],pch=c(19,19,21),bty="n")
  # abline(h=thresh[1:pic],col=c(colors[2:1],"orange"),lty=2,lwd=2)
  dev.off()
}


QTLb <- (-log10(dos_QTL$X0.8$pval) > dos_QTL$X0.8$perm.thr) & !qtl.pval
QTLa <- (-log10(anc_res$X0.8$pval) > anc_res$X0.8$perm.thr) & !qtl.pval
QTLh <- (-log10(haplo_res$X0.8$pval) > haplo_res$X0.8$perm.thr) 
QTLh[is.na(QTLh)] <- F

map[QTLpos,]
for(m in list(map[QTLa,],map[QTLb,],haplos4$map[QTLh,])){
  a <- sapply(unique(m$chromosome),function(x){
    range(m$position[m$chromosome == x],na.rm = T)
  })
  
  #The range of significant positions
  print(a)
}

max_pval <- function(pval,m){
  data <- cbind(m,pval)
  max_mark <- sapply(unique(m$chromosome),function(x){
    data[data$chromosome == x,][which.max(data$pval[data$chromosome == x]),]
  })
  return(t(max_mark))
}

max_pval(-log10(dos_QTL$X0.8$pval[QTLb]),map[QTLb,])
max_pval(-log10(anc_res$X0.8$pval[QTLa]),map[QTLa,])
best_pos <- max_pval(-log10(haplo_res$X0.8$pval)[QTLh],haplos4$map[QTLh,])

good_markers <- as.character(unlist(best_pos[,1]))



saveRDS(list(haplo4= haplo_res,
             snp = dos_QTL,
             anc = anc_res),"true_haplotypes_result.RDS")
data <- readRDS("true_haplotypes_result.RDS")

haplo_res <- data$haplo4
dos_QTL <- data$snp
anc_res <- data$anc

h_i <- hap_indexer(gen = anc,4)
which(apply(h_i,1,function(i) 1348%in%i))

# Scenario simulation --> Genetic diversity ----------------
anc <- link_NAM("../mpQTL v0/PedigreeSIM/NAM_crosses/3_ancestral/cross200_founderalleles.dat",
                totallele = "../mpQTL v0/PedigreeSIM/Parents/Total_pop.txt",ploidy = 4)
map <- read.table("PedigreeSim/Potato.map",header=T)
QTLpos <- c(453,1348,4023)

#Create sets of effects with different variations
effects <- sapply(c(0.01,0.1,1),function(sd){
  multi.effects2(anc[QTLpos[1],,drop=F],
                 anc_alleles = matrix(0:399,ncol=10),
                 anc_sd = sd,seed = 7)
})
eff <- do.call(cbind,effects)

#All the alleles present at this QTL
alleles <- unique(anc[QTLpos[1],])

#Calculate the x values for each effect
anc_alleles <- mat2list(matrix(0:399,ncol=10))
anc_matrix <- sapply(anc_alleles,function(c) alleles%in%c)
ancestors <- apply(anc_matrix,1,which)
ancestors <- sapply(ancestors,function(i) which(i==unique(ancestors)))
xvals <- sapply(1:length(unique(ancestors)),function(i){
  r <- ancestors/max(ancestors)+(i-1)*1.5
})

#Transform the values on the same scale
eff <- apply(eff,2,function(e) (e-min(e))/max(e-min(e)))

colors = colorspace::sequential_hcl(palette="YlGnBu",5)[c(1,3,5)]

jpeg("Plots/Effect_comparison.jpeg",width = 2400,height = 1000,res=220)
plot(xvals,eff,pch=19,col=rep(color,each=nrow(eff)),
     axes=F,xlab="",
     ylab="Relative allelic effect size",
     main="19 alleles from 3 AG",cex=1.5)
#We use the xvals to create axis
at <- apply(xvals,2,function(i) mean(unique(i)))
mtext(side=1,text = paste("Scenario",1:3),at = at,line=3)
for(i in 1:ncol(eff)){
  axis(1,at=unique(xvals[,i]),
       labels = paste0("AG",1:3),cex.axis=0.85,padj = -1,
       line = 1)
}
dev.off()

#We use a reduced dataset to generate phenotypes,
#and QTL mapping. 
red_phe <- sapply(1:3,function(i){
  pheno2(anc[QTLpos[1],,drop=F],effects = effects[[i]],
         ploidy = 4,herit = 0.8)
})
red_ge <- anc[1:1000,]
red_map <- map[1:1000,]
red_res <- map.QTL(red_phe,genotypes = red_ge,4, red_map,K = T)
pvals <- sapply(red_res, function(x) -log10(x$pval))

#Calculate where to place marker values on the x axis
abspos <- sapply(unique(red_map$chromosome),function(x) max(map$position[map$chromosome == x]))
abspos<-c(0,abspos)
pos <- red_map$position + abspos[red_map$chromosome]

#Plotting
for(pic in 1:3){
  jpeg(paste0("Plots/Effect_comparison_Manhattan_",pic,".jpeg"),
       width = 3200,height = 1800,res=250)
  plot(0,type="n",xlim=c(0,max(pos)),ylim=c(0,max(pvals)),
       ylab="-log10(pval)",xlab="cM",main="Association with IBD alleles")
  abline(v=pos[QTLpos[1]],lty=2,lwd=1.5)
  for(i in 1:pic){
    points(pos,-log10(red_res[[i]]$pval),col=color[i],pch=19,cex=1.3)
  }
  legend("topright",legend=paste("Scenario",1:3),pch=19,col=colors,bty="n")
  dev.off()
  
}
jpeg("Plots/Effect_comparison_Manhattan.jpeg",width = 3200,height = 1800,res=250)
plot(0,type="n",xlim=c(0,max(pos)),ylim=c(0,max(pvals)),
     ylab="-log10(pval)",xlab="cM",main="Association with IBD alleles")
abline(v=pos[QTLpos[1]],lty=2,lwd=1.5)
for(i in 1:pic){
  points(pos,-log10(red_res[[i]]$pval),col=color[i],pch=19,cex=1.3)
}
legend("topright",legend=paste("Scenario",1:3),pch=19,col=colors,bty="n")
dev.off()

