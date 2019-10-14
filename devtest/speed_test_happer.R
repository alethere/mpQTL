gen <- data.table::fread("../mpQTL v0/PedigreeSIM/NAM_crosses/3_ancestral/cross200_genotypes.dat",
                         header = T)
gen <- gen[,-1]
map <- read.table("PedigreeSim/Potato.map",header=T)


# Create haplotypes -------------------------------------------------------

#First we create the indeces of the SNPs we must concatenate
hap_indexer <- function(gen,l){
  t(sapply(1:floor(nrow(gen)/l),function(n) 1:l+(n-1)*l ))
}
hap_index <- hap_indexer(gen,5)


happer <- function(gen,hap_index){
  t(apply(hap_index,1,function(h_i){
    Reduce(paste0,data.frame(t(gen[h_i,])))
  }))
  
}

happer2 <- function(gen,hap_index){
  res <- apply(hap_index,1,function(h_i){
    gen[,Map(paste,gen[h_i,],MoreArgs = list(collapse=""))] 
  })
  do.call(rbind,res)
  
}

happer(gen,hap_index)
happer2(gen,hap_index)