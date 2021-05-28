
geno <- readRDS("research/phenotyping/manuscript_geno.RDS")
snpmap <- read.table("research/Potato.map",header = T)

geno$snp <- lapply(geno$snp,function(gs){
  rownames(gs) <- subset(snpmap, chromosome %in% 1:5)$marker
  return(gs)
})

hapmap <- geno$haplo[[1]]$map
write.table(hapmap,file = "research/data_availability/haplotype.map",quote = F,row.names = F,col.names = T)
geno$haplo <- lapply(1:length(geno$haplo),function(i){
  gh <- geno$haplo[[i]]$haplotypes
  gi <- geno$ibd[[i]]
  rownames(gh) <- hapmap$marker
  colnames(gh) <- colnames(gi)
  return(gh)
})

geno$haplo[[1]][1:10,1:10]

nam_names <- rep(paste("NAM",c(1,3,7,10),sep = ""),each = 11)
nam_names <- paste(nam_names,rep(1:11,4),sep = "_")

for(n in seq_along(nam_names)){
  write.table(geno$snp[[n]],file = paste0("research/data_availability/",nam_names[n],"_biallelic.txt"),
              quote = F, col.names = T, row.names = T)
  write.table(geno$ibd[[n]],file = paste0("research/data_availability/",nam_names[n],"_IBD.txt"),
              quote = F, col.names = T, row.names = T)
  write.table(geno$haplo[[n]],file = paste0("research/data_availability/",nam_names[n],"_haplotype.txt"),
              quote = F, col.names = T, row.names = T)
}

pheno <- readRDS("research/phenotyping/manuscript_pheno.RDS")
phens <- sapply(pheno,'[[',"pheno")
rownames(phens) <- NULL
colnames(phens) <- nam_names
write.table(phens,file = "research/data_availability/Phenotypes.txt",
            quote = F, col.names = T, row.names = F)
