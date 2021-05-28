#Defining cross names
per_nam <- 0:10
pop_nums <- as.vector(sapply(c(0,200,600,900),function(i) i+ per_nam))
pop_nums <- sprintf("%03d",pop_nums)
pop_files <- paste0("research/PedigreeSim/","cross",pop_nums)
found_files <- paste0(pop_files,"_founderalleles.dat")

#Getting stored phenotypes
pheno <- readRDS("research/phenotyping/manuscript_pheno.RDS")
names(pheno) <- paste0("cross",pop_nums)

#Defining the QTL positions and polygenic positions
QTLpos <- c(430,1348,2800)
geno <- data.table::fread(found_files[1])
polyg <- which(geno$marker %in% names(pheno$cross000$effects))


herits <- lapply(1:length(found_files),function(i){
  IBD <- link_NAM(found_files[i],totallele = "research/PedigreeSim/Total_pop.txt",
                  marker = c(QTLpos,polyg))
  names(pheno[[i]]$effects)[1:3] <- rownames(IBD)[1:3]

  IBD <- IBD[names(pheno[[i]]$effects),]
  phe <- phsum(IBD,effects = pheno[[i]]$effects,ploidy = 4,table = T)
  phe <- cbind(phe[,1:3],rowSums(phe[,-1:-3]))
  colnames(phe) <- c(paste0("QTL",1:3),"polygen")
  tot_var <- var(pheno[[i]]$pheno)

  #Calculate all variances and co-variance
  varcov <- var(phe)
  variances <- diag(varcov)

  names(variances) <- c(paste0("QTL",1:3),"polygen")


  narrow_sense <- variances/tot_var
  broad_sense <- colSums(varcov)/tot_var
  return(data.frame(narrow = narrow_sense,broad = broad_sense))
})

names(herits) <- names(pheno)

lapply(c("narrow","broad"),function(type){
  h <- t(sapply(herits,'[',1:3,type))
  perc <- t(apply(h,1,function(n) n/sum(n)))

  list(h = h,
       mean_perc )
})
narrow <- t(sapply(herits,'[',1:3,"narrow"))
broad <- t(sapply(herits,'[',1:3,"broad"))

plot_narrow_broad <- function(narrow,broad,...){
  col <- hcl.colors(3,"Zissou 1")
  plot(narrow,broad,
       xlab = "option A, without covariance",
       ylab = "option B, with covariance",
       pch = 19,
       col = rep(col,each = 44),cex = 2.2,...)
  legend("topright",
         bty = "n",
         legend = paste0("QTL",1:3),
         col = col,
         pch = 19)

  mbroad <- colMeans(broad)
  mnarrow <- colMeans(narrow)
  abline(h = mbroad,
         v = mnarrow,
         col = col,lty = 2)

  labs <- round(c(mbroad,mnarrow),3)
  text(x = c(0,0,0,mnarrow+0.01),y = c(mbroad+0.01,0,0,0),
       labels = labs, col = col,cex = 0.7,adj = 0)
}


plot_narrow_broad(narrow,broad,xlim = c(0,0.5),ylim = c(0,0.5),main = "QTL heritabilities")
plot_narrow_broad(perc_n,perc_b,xlim = c(0,0.7),ylim = c(0,0.7),main = "Percentage of QTL variance explained by each QTL")

#Detection power per QTL
d_test <- seq(0,5,by = 0.5)
power_study <- lapply(d_test,function(space){
  power_data <- lapply(1:length(pval),function(model_type){
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


    lapply(1:ncol(pval[[model_type]]),function(pvi){
      pv <- pval[[model_type]][,pvi]
      QTL <- read_QTL(pv,m,threshold = th,space = space)

      res <- sapply(1:3,function(i){
        chrom <- c(1,2,4)
        qi <- grepl(paste0(chrom[i],"\\."),names(QTL))
        if(!any(qi)) return(c(QTL_precision = 0,
                              QTL_sensitivity = 0,
                              mark_precision = 0,
                              mean_dist = NA))
        q <- QTL[which(qi)]
        class(q) <- "QTL_list"
        power_QTL(q,trueQTL[i,])$power
      })
      colnames(res) <- paste0("QTL",1:3)
      return(res)
    })
  })

  names(power_data) <- c("snp","ibd","hap")
  extract <- function(code) sapply(power_data,'[',code,)

  power_data <- lapply(1:3,function(model){
    res <- sapply(1:3,function(qn){
      qtl_power <- data.frame(t(sapply(power_data[[model]],'[',,qn)))
      pops <- rep(sprintf("NAM%02d",c(1,3,7,10)),each = 11)
      qtl_ls <- split(qtl_power,pops)
      sapply(qtl_ls,colMeans,na.rm = T)["QTL_sensitivity",]
    })
    colnames(res) <- c("QTL1","QTL2","QTL3")
    return(res)
  })
 names(power_data)<- c("snp","ibd","hap")

  return(power_data)
})
names(power_study) <- d_test
lapply(power_study$`3`,round,2)

#Calculating false negative rate--------------

