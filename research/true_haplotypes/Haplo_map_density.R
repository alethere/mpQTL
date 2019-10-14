
# Functions ---------------------------------------------------------------
source("mapQTL_fun.R")
source("haplo_fun.R")
source("viz_fun.R")

# Data loading ------------------------------------------------------------
phe <- readRDS("phenotypes.RDS")
gen <- data.table::fread(
  "../mpQTL v0/PedigreeSIM/NAM_crosses/3_ancestral/cross200_genotypes.dat",
  header = T)[,-1]
dos <- data.table::fread(
  "../mpQTL v0/PedigreeSIM/NAM_crosses/3_ancestral/cross200_alleledose.dat",
  header = T)[,-1]
anc <- link_NAM(
  "../mpQTL v0/PedigreeSIM/NAM_crosses/3_ancestral/cross200_founderalleles.dat",
  totallele = "../mpQTL v0/PedigreeSIM/Parents/Total_pop.txt",
  ploidy = 4)
map <- read.table("PedigreeSim/Potato.map",header=T)

data <- list(gen = as.data.frame(gen),
             dos = as.data.frame(dos),
             anc = as.data.frame(anc),
             map = map)
# Create haplotype sets ---------------------------------------------------

## Clean maps (wo QTL positions) --------
QTLpos <- c(453,1348,4023)
true_markers <- map[QTLpos,]
QTLmarkers <- sapply(1:3,function(m){
  m <- unlist(true_markers[m,])
  which(m[2] == map$chromosome & m[3] == map$position)
})
QTLmarkers <- do.call(c,QTLmarkers)
#ilter out QTL markers
data <- lapply(data,'[',-QTLmarkers, )

## Create new haplotype/map sets
#Start with 6876 markers, we'll go down by steps
marker_num <- c(6000, 4000, 2000, 1000)

cl <- parallel::makeCluster(11)
parallel::clusterExport(cl,c(".hap_indexer",".hap_pos",".happer",
                             ".hap_correct","haplotyper","data"))
hap_data <- parallel::parLapply(cl,marker_num,function(m){
  less_markers <- round(seq(1,nrow(data$map),length.out = m))
  new_data <- lapply(data,'[',less_markers, )
  hap4 <- haplotyper(new_data$gen,new_data$map,l = 4,method = "sliding")
  hap6 <- haplotyper(new_data$gen,new_data$map,l = 6,method = "sliding")
  hap8 <- haplotyper(new_data$gen,new_data$map,l = 8,method = "sliding")
  
  return(list(data=new_data,
              hap4=hap4,
              hap6=hap6,
              hap8=hap8))
})

parallel::stopCluster(cl)
saveRDS(hap_data,"Haplo_map_density.RDS")

hap_data <- readRDS("Haplo_map_density.RDS")
pheno <- phe$cross200$`0.5`
hap_res <- lapply(marker_num,function(m){
  i <- which(m == marker_num)
  
  dat <- hap_data[[i]]
  
  haplo_4 <- map.QTL(phenotypes = pheno,
                       genotypes = dat$hap4$haplotypes,
                       map = dat$hap4$map,
                       ploidy = 4,K = T)
  haplo_6 <- map.QTL(phenotypes = pheno,
                     genotypes = dat$hap6$haplotypes,
                     map = dat$hap6$map,
                     ploidy = 4,K = T)
  haplo_8 <- map.QTL(phenotypes = pheno,
                     genotypes = dat$hap8$haplotypes,
                     map = dat$hap8$map,
                     ploidy = 4,K = T)
  
  res <- list(maps = list(hap4 = dat$hap4$map,
                          hap6 = dat$hap6$map,
                          hap8 = dat$hap8$map),
              res = list(hap4 = haplo_4,
                         hap6 = haplo_6,
                         hap8 = haplo_8))
  print(paste(m,"Completed"))
  return(res)
})
names(hap_res) <- marker_num
saveRDS(hap_res,"Results_haplotypes_densities.RDS")
dim(hap_res$`6000`$maps$hap8)

res_dos <- map.QTL(phenotypes = pheno,
        genotypes = data$dos,
        map = data$map,
        ploidy = 4,K = T)
res_dos_hap <- lapply(marker_num,function(m){
  less_markers <- round(seq(1,nrow(data$map),length.out = m))
  pv <- res_dos$pheno1$pval[less_markers]
  mp <- data$map[less_markers,]
  return(list(pvals = pv, map= mp))
})
names(res_dos_hap) <- paste0(marker_num,".dos")

res_anc <- map.QTL(phenotypes = pheno,
                   genotypes = data$anc,
                   map = data$map,
                   ploidy = 4,K = T)

# QTLs detected ----------------------


pvals <- lapply(hap_res,function(h){
  sapply(h$res,function(nm){
    nm[[1]]$pval
  })
})
pvals <- unlist(pvals,recursive = F)
pvals <- c(lapply(res_dos_hap,'[[',"pvals"),pvals)


hap_maps <- lapply(hap_res,function(h){
  h$maps
})
hap_maps <- unlist(hap_maps, recursive = F)
hap_maps <- c(lapply(res_dos_hap,'[[',"map"),hap_maps)

saveRDS(hap_maps,"Haplotype_maps_large_list.RDS")

#Function for extracting significant QTLs -----------------
#threshold = 3.6, distance between markers = 1cM
QTLs <- lapply(1:length(pvals),function(i){
  pv <- pvals[[i]]
  hm <- hap_maps[[i]]
  sig_m <- which(-log10(pv) > 3.6)
  
  tit <- names(pvals)[i]
  png(paste0("Plots/Hap_density_",tit,".png"),width = 2000,height = 800,res=200)
  skyplot(-log10(pv),hm,main=tit,threshold = 3.6,ylim=c(0,15))
  dev.off()
  
  #With these two lines we identify what is the distance between
  #significant markers (expressed as a vector of  "distance 
  #with the previous marker"). The first value of each vector is Inf
  sig_p <- split(hm$position[sig_m],hm$chromosome[sig_m])
  difs <- do.call(c,lapply(sig_p,function(p){
    c(Inf,p[-1]-p[-length(p)])
  }))
  
  # Loop for generating QTL positions ------
  sig <- cbind(hm[sig_m,],difs,pval=pv[sig_m])
  QTLs <- list(); li_count <- 0; current_ch <- sig[1,]$chromosome
  for(i in 1:nrow(sig)){

    if(sig[i,]$chromosome != current_ch){
      li_count <- li_count + 1
      current_ch <- sig[i,]$chromosome
      QTLs[[li_count]] <- sig[i,]
    }else if(round(sig[i,]$difs) > 1){
      li_count <- li_count + 1
      QTLs[[li_count]] <- sig[i,]
    }else{
      QTLs[[li_count]] <- rbind(QTLs[[li_count]],sig[i,])
    }
  }
  
  summary_QTL <- t(sapply(QTLs,function(q){
    return(data.frame(
      chr = mean(q$chromosome),
      min = min(q$position),
      max = max(q$position),
      mean = mean(q$position),
      median = median(q$position),
      peak = q$position[which.min(q$pval)],
      peak_pval = min(q$pval),
      mean_pval = mean(q$pval),
      median_pval = median(q$pval),
      support = nrow(q)
    ))
  }))
  rownames(summary_QTL) <- paste0("QTL",1:length(QTLs))
  
  return(summary_QTL)
})
names(QTLs) <- names(pvals)
sink("QTLs_map_density.txt")
print(QTLs)
sink()

saveRDS(h,"Results_haplotype_QTL_detection.RDS")

# Number of missing values ------------------------------------------------
hap_res[[4]]$res$hap8$pheno1$pval
X <- dosage.X(hap_data[[4]]$hap8$haplotypes[911,],haplotype=T,ploidy = 4)
K <- calc.K(t(hap_data[[4]]$hap8$haplotypes),haplotypes = T,ploidy = 4)
heatmap(K,Colv = NA,Rowv = NA)
sdX <- apply(hap_data[[4]]$hap8$haplotypes,1,function(g){
  X <- dosage.X(g,haplotype = T,ploidy = 4)
  sd(colSums(X))
})

hist(sdX)

# LD decay study ------------


LD_decay <- function(dos,map,win_size = 0.1,max_dist = NULL){
  
  #LD calculation
  LD <- cor(t(dos))
  not_dup <- lower.tri(LD,diag=F)
  LD <- LD[not_dup]
  
  #centimorgan distance between markers
  cm_d <- sapply(map$position,function(p) abs(p-map$position))
  cm_d <- cm_d[not_dup]
  
  #but distances between markers at different chromosomes are meaningless
  same_c <- sapply(map$chromosome, function(p) p == map$chromosome)
  same_c <- same_c[not_dup]
  
  cm_d <- cm_d[same_c]
  LD <- LD[same_c]
  
  if(is.null(max_dist)) max_dist <- max(cm_d)
  
  #This is inefficient, first, markers at max_dist should 
  #be selected and then LD calculated in order to minimize the number
  #of operations. Optimization would be marginal, though
  windows <- seq(0,max_dist,win_size)
  ld_estimates <- t(sapply(seq_along(windows),function(w){
    sel <- cm_d > windows[w] & cm_d < windows[w] + win_size
    return(quantile(LD[sel]^2,c(0.5,0.8,0.9,0.95)))
  }))
  
  res <- data.frame(ld_estimates, distance = windows)
  attr(res,"class") <- c("LD","data.frame")
  attr(res,"max_dist") <- max_dist
  
  return(res)
}

LD_decay <- function(dos,map,win_size = 0.1,max_dist = NULL,per_chr = F){
  
  #First we split dosages per chromosome and calculate correlations
  #only within chromosomes
  dos_per_chr <- split(dos,map$chromosome)
  LD <- lapply(dos_per_chr,function(g){
    ld <- cor(t(g))
    not_dup <- lower.tri(ld,diag=F)
    return(ld[not_dup]^2)
  })
  
  #then we calculate the distance between markers
  pos_per_chr <- split(map$position,map$chromosome)
  dis <- lapply(pos_per_chr,function(d){
    res <- sapply(d,function(p) abs(p-d))
    not_dup <- lower.tri(res,diag=F)
    return(res[not_dup])
  })
  
  if(!per_chr){
    
    #Per window we calculate the percentile estimates
    LD <- do.call(c,LD)
    dis <- do.call(c,dis)
    
    if(is.null(max_dist)) max_dist <- max(dis)
    
    windows <- seq(0,max_dist,win_size)
    ld_estimates <- t(sapply(seq_along(windows),function(w){
      sel <- dis >= windows[w] & dis < windows[w] + win_size
      return(quantile(LD[sel],c(0.5,0.8,0.9,0.95)))
    }))
    
    #We calculate the background correlation between chromosomes
    dos_sample <- lapply(dos_per_chr,function(d) d[sample(1:nrow(d),100),] )
    dos_sample <- do.call(rbind,dos_sample)
    back_ld <- cor(t(dos_sample))^2
    back_ld <- back_ld[lower.tri(back_ld,diag=F)]
    back_ld <- quantile(back_ld,c(0.5,0.8,0.9,0.95))
    
    #We add some extra features
    res <- list(LD = data.frame(ld_estimates,
                                distance = windows), 
                background = back_ld)
    attr(res,"max_dist") <- max(dis,na.rm=T)
    attr(res,"class") <- c("LD","list")
    
  }else{
    res <- lapply(1:length(unique(map$chromosome)),function(i){
      ld <- LD[[i]]
      d <- dis[[i]]
      
      windows <- seq(0,max(d,na.rm=T),win_size)
      ld_estimates <- t(sapply(seq_along(windows),function(w){
        sel <- d > windows[w] & d < windows[w] + win_size
        return(quantile(ld[sel],c(0.5,0.8,0.9,0.95)))
      }))
      
      back_ld <- quantile(ld[d > max(d,na.rm=T)*0.9],c(0.5,0.8,0.9,0.95))
      
      #We add some extra features
      res <- list(LD = data.frame(ld_estimates,distance = windows), 
                  background = back_ld)
      attr(res,"max_dist") <- max(d,na.rm=T)
      attr(res,"class") <- c("LD","list")
      return(res)
    })
    
  }
  
  return(res)
}

LD <- LD_decay(dos,map)

max_dist = 50 #for very large values of max_dist the spline doesn't work
#because there are not enough observations to calculate the ld for large
#distances
LD <- LD_decay(dos,map,win_size = 0.5)

plot.LD <- function(LD,max_dist = NULL,main=NULL){
  if(is.null(max_dist)) max_dist <- attr(LD,"max_dist")
  if(max_dist > attr(LD,"max_dist")) 
    stop("Cannot choose a larger max_dist than ",attr(LD,"max_dist"))
  if(!"LD" %in% class(LD)){
    stop("Object provided is not an LD list")
  }
  
  linkage <- LD$LD
  
  plot(0,type="n",
       xlab = "cM",ylab = "r2",main = main,
       xlim = c(0,max_dist),ylim=c(0,0.5))
  
  cols <- colorspace::qualitative_hcl(4)
  halflife <- c()
  for(i in 1:4){
    ld <- linkage[,i][linkage$distance <= max_dist]
    dis <- linkage$distance[linkage$distance <= max_dist]
    
    #We take out missing percentile values because spline cannot handle it
    if(any(is.na(ld))){
      ld <- na.omit(ld)
      dis <- dis[-attr(ld,"na.action")]
    }
    
    points(dis,ld,
           pch=19,cex=0.85,col=cols[i])
    
    sp <- smooth.spline(dis,ld,spar = 0.6)
    lines(sp,lwd=2,col= lighten(cols[i]))
    
    top <- max(ld,na.rm=T) 
    halfpoint <- (top - LD$background[i])/2 + LD$background[i]
    abline(v = sp$x[which.min(abs(sp$y - halfpoint))],
           col=cols[i],
           lty=2,lwd=1.5)
    halflife[i] <- round(sp$x[which.min(abs(sp$y - halfpoint))],1)
    
  }
  legend("topright",pch=19,
         legend=paste(c("50th","80th","90th","95th")," ld1/2=",halflife)
                      ,col=cols,bty="n")
  
}

plot(LD,main="NAM3 LD decay",max_dist = 50)
LD$background



