#Alejandro Thérèse Navarro
#March 2021

#This file allows to utilize the mpQTL package to obtain the results
#shown in the article "Multiallelic models for QTL mapping in diverse polyploid populations"
#by Therese Navarro et al.

#0 Loading of mpQTL package
source("mapQTL_fun.R")

#1 Genetic distances of NAM------------
#In the article the genetic distance of a single NAM3 population is shown
#In the form of a heatmap and in a pca

NAM3 <- data.table::fread("Genotypes/NAM3_1_biallelic.txt")

#The matrix K contains the genetic distances
K <- calc.K(t(NAM3[,-1]),ploidy = 4)
heatmap(K,Colv = NA,Rowv = NA)

#The principal components was made using the matrix K
pc <- prcomp(K)
var.pc <- paste0("PC",1:length(pc$sdev),"  ",
               round(pc$sdev^2/sum(pc$sdev^2)*100,2),"% variance")

#We plot first and second components
comp <- 1:2
plot(pc$x[,comp],xlab=var.pc[comp[1]],ylab=var.pc[comp[2]],cex = 2)

#2 QTL analyses --------------------
#2.0 Data loading ------------------
#Please note these will take a while

#We first load the maps
map <- read.table("Potato.map",header = T)
hapmap <- read.table("haplotype.map",header = T)

#We read the file names for each type of marker
snp_files <- list.files("Genotypes",pattern = "biallelic",full.names = T)
ibd_files <- list.files("Genotypes",pattern = "IBD",full.names = T)
hap_files <- list.files("Genotypes",pattern = "haplotype",full.names = T)

#Then we load all the genotype files into a geno list
#each element will represent a cross, and will contain a list
#of three genotype files (snp, ibd and hap)
geno <- lapply(1:length(snp_files),function(i){
  snp <- read.table(snp_files[i])
  ibd <- read.table(ibd_files[i])
  hap <- read.table(hap_files[i])
  return(list(snp = snp,ibd = ibd, hap = hap))
})

#We also load the phenotypes
#Each column corresponds to one population
#Each row corresponds to each individual in that population
pheno <- read.table("Phenotypes.txt",header = T)

#2.1 P-value calculation ------------------
#For each population we run the three methods
#Please note this is a long step, will take some time
result <- lapply(1:length(geno),function(pop){
  phe <- pheno[,pop]
  snp <- geno[[pop]]$snp
  ibd <- geno[[pop]]$ibd
  hap <- geno[[pop]]$hap

  #Note that the argument no_cores controls the number of cores
  #assigned to the job (for parallel computation). Do not set this number
  #higher than the number of cores in your computer (can be checked using parallel::detectCores())
  bi <- map.QTL(phe,snp,ploidy = 4,map = map, K = T,no_cores = 10)
  anc <- map.QTL(phe,ibd,ploidy = 4,map = map, K = T,no_cores = 10)
  haplo <- map.QTL(phe,hap,ploidy = 4,map = hapmap, K = T,no_cores = 10)
  cat("

      Population",pop,"done ----------

      ")
  return(list(snp = bi, ibd = anc, hap = haplo))
})
saveRDS(result,"qtl_analyses.RDS",compress = T)

#2.2 Threshold calculation -----------
#This is a permutation-based approach. Can be very lengthy for
#all datasets. This loop performs it once per model type.

result <- lapply(1,function(pop){
  phe <- pheno[,pop]
  snp <- geno[[pop]]$snp
  ibd <- geno[[pop]]$ibd
  hap <- geno[[pop]]$hap

  bi <- map.QTL(pe,snp,ploidy = 4,nc,map = map, K = T,
                permutation = "pop",nperm = 100)
  anc <- map.QTL(pe,ibd,ploidy = 4,map = map, K = T,
                 permutation = "pop",nperm = 100)
  haplo <- map.QTL(pe,hap,ploidy = 4,map = hapmap,K = T,
                   permutation = "pop",nperm = 100)
  cat("

      Population",pop,"done ----------

      ")
  return(list(snp = bi, ibd = anc, hap = haplo))
})
thresh <- sapply(result[[1]],function(mod) mod$pheno1$perm.thr)
names(thresh) <- c("snp","ibd","hap")
saveRDS(thresh,"threshold.RDS")

#3 Power studies ------------
source("research/data_availability/power_fun.R")
qtl_res <- readRDS("research/data_availability/qtl_analyses.RDS")
thresh <- readRDS("research/data_availability/threshold.RDS")
map <- read.table("research/data_availability/Potato.map",header = T)
hapmap <- read.table("research/data_availability/haplotype.map",header = T)


pval <- lapply(c("snp","ibd","hap"),function(model){
  res <- sapply(qtl_res,function(qr){
    qr[[model]]$pheno1$pval
  })
  #Uncomment this line if populations names are missing in the result
  #colnames(res) <- colnames(pheno)
  return(res)
})
names(pval) <- c("snp","ibd","hap")

#These are the true QTL positions in all populations
QTLpos <- c(430,1348,2800)
trueQTL <- map[QTLpos,]

linkings <- seq(0,10,length.out = 21)
linkings <- 3
#For each linking distance
power_study <- lapply(linkings,function(space){
  #And for each model type
  power_data <- lapply(1:length(pval),function(model_type){
    #We obtain the threshold and map corresponding
    #to the model type
    if(model_type == 1){
      th <- 10^-thresh[1]
      m <- map
    }else if(model_type == 2){
      th <- 10^-thresh[2]
      m <- map
    }else{
      th <- 10^-thresh[3]
      m <- hapmap
    }

    #And then return the QTL analyses
    apply(pval[[model_type]],2,function(pv){
      QTL <- read_QTL(pv,m,threshold = th,space = space)

      q <- summary(QTL)
      pvm <- cbind(m,pval = pv[m$marker])

      fp_rate <- lapply(1:nrow(trueQTL),function(i){
        true_q <- trueQTL[i,]
        q_range <- data.frame(min = true_q$position - 3, max = true_q$position + 3)
        # in_chr <- q$chr == true_q$chromosome
        # in_interval <- (q$min <= true_q$position) & (q$max >= true_q$position)
        # q_range <- q[in_interval & in_chr,c("min","max")]
        negs <- !((pvm$position >= q_range$min) & (pvm$position <= q_range$max) &
                    (pvm$chromosome == true_q$chromosome))
        false_pos <- negs & pvm$pval <= th
        poses <- !negs
        true_pos <- poses & pvm$pval <= th

        return(list(fp = false_pos,
                    neg = negs,
                    tp =true_pos ,
                    pos = poses))
      })
      fp <- apply(sapply(fp_rate,'[[',"fp"),1,any)
      neg <- apply(sapply(fp_rate,'[[',"neg"),1,all)
      tp <- apply(sapply(fp_rate,'[[',"tp"),1,any)
      pos <- apply(sapply(fp_rate,'[[',"pos"),1,any)
      fprate <- round(sum(fp,na.rm = T)/sum(neg,na.rm = T),2)
      tprate <- round(sum(tp,na.rm = T)/sum(pos,na.rm = T),2)

      res <- power_QTL(QTL,trueQTL)
      mean_pv <- mean(res$data[ ,"mean_pv"],na.rm=T)
      return(c(res$power,
               fp_rate = fprate,
               tp_rate = tprate,
               mean_pv = mean_pv,
               n_QTL = length(QTL) ))
    })
  })

  names(power_data) <- c("snp","ibd","hap")
  extract <- function(code) sapply(power_data,'[',code,)

  power_measures <- rownames(power_data$snp)
  power_data <- lapply(power_measures,extract)
  names(power_data) <- power_measures

  return(power_data)
})
names(power_study) <- linkings

#This loop helps summarizing the results of the power study
summaries <- lapply(power_study,function(pw){
  means <- lapply(pw,function(p){
    res <- sapply(1:4,function(nam){
      colMeans(p[1:11 + 11*(nam-1),])
    })
    colnames(res) <- paste0("NAM",c(1,3,7,10))
    return(t(res))
  })
  return(means)
})

#Table 1 in the article is showing the following results:
#QTL sensitivity is detection power
#mean_dist is accuracy in cM from true position
lapply(summaries$`3`[c("QTL_sensitivity","QTL_precision","mean_dist")],round,2)

#The percentiles of figure 3 can be obtained as such:
power_quantiles <- lapply(1:3,function(model_type){
  pq <- lapply(names(power_study[[1]]),function(pow_stat){
    pol_edges <- sapply(seq_along(linkings),function(i){
      #These are all the measures of a particular pow_stat for a particular model_type
      pow <- as.vector(power_study[[i]][[pow_stat]][,model_type])
      quantile(pow,c(0.2,0.8),na.rm = T,type=1)
    })
    colnames(pol_edges) <- linkings
    return(pol_edges)
  })
  names(pq) <- names(power_study[[1]])
  return(pq)
})
names(power_quantiles) <- c("snp","ibd","hap")
power_quantiles

#4 Example QTL detection -------------
#Figure 5 example was taken from chromosome 2 of population NAM7_4
fig5 <- lapply(seq_along(pval),function(model){
  if(model == 3){
    #Then it's the haplotype model
    m <- hapmap
  }else{
    m <- map
  }
  sel <- m$chromosome == 2 & m$position > 50 & m$position < 70
  cbind(m[sel,],pval = pval[[model]][sel,"NAM7_4"])
})

par(mfrow = c(1,3))
for(data in fig5){
  plot(data$position,-log10(data$pval))
}

