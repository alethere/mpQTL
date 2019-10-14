


source("mpQTL_fun_clean.R")
source("D:/OneDrive - WageningenUR/polyploids/software/haploutils/NAMhaplo_utils.R")
source("D:/OneDrive - WageningenUR/polyploids/software/gtfun/pkg_rasso/rasso/R/GWASplot.R")
source("D:/OneDrive - WageningenUR/polyploids/software/PolyHaplotyper/PolyHaplotyper_20190128.R")


outf <- "output_inferredHaps"
dir.create(outf)


# # data for testing ------------
# 
# sample_data <- readRDS("sample_data.rds")
# names(sample_data)
# sample_data$phenotype[1:5,]
# sample_data$genotypes[1:5,1:5]
# dim(sample_data$genotypes)
# dim(sample_data$dosage)
# sample_data$dosage[1:5,1:5]
# rownames(sample_data$phenotype)[1:20]
# colnames(sample_data$dosage)[1:20]
# 
# ## individual names in phenotypes are wrong
# rownames(sample_data$phenotype) <- colnames(sample_data$dosage)

# ### snpdose file
# snpdose <- sample_data$dosage
# dim(snpdose)
# snpdose[1:5,1:5]
# rownames(snpdose) <- map$marker
# # snpdose <- as.matrix((snpdose))
# # check map markers are identical to snpdose markers
# if(!all(map$marker %in% rownames(snpdose))) stop("error 1")
# 
# ### ped file
# ped <- read.table("cross200.ped", header = T)
# 
# ### parents
# colnames(ped) <- c("genotype","mother","father")
# parents <- data.frame(mother=unique(ped$mother[!is.na(ped$mother)]),
#                       father=unique(ped$father[!is.na(ped$father)]))
# parents <- as.matrix(parents)
# 
# F1 <- list()
# for (p in seq_len(nrow(parents))) {
#   F1[[p]] <- ped$genotype[!is.na(ped$mother) & ped$mother==parents[p, 1] &
#                             !is.na(ped$father) & ped$father==parents[p, 2]]
# }
# 
# # # noPedigree
# # if(noPedigree) {
# #   parents <- NULL
# #   F1 <- NULL
# # }
# 
# ### subset snpdose according to hb_list
# snpdose <- snpdose[match(unlist(hb_list),rownames(snpdose)),]
# dim(snpdose)
# 
# 
# ## haplotype inference
# load(paste0(outf,"/hb_list.RData"))
# 
# options(PolyHaplotyper_ahcdir=paste0("D:/polyploids/software/PolyHaplotyper/",
#                                      "ahccompletelist_8snps"))
# ### infer
# write(paste0("Start time: ", Sys.time()),
#       file=paste0(outf,"/log.txt"), append = T)
# 
# results <- inferHaplotypes(mrkDosage=snpdose, ploidy=4,
#                            haploblock=hb_list[1:10],
#                            parents=parents, F1=F1,
#                            maxmrk=6,
#                            maxparcombs=150000)
# 
# write(paste0("Finish time: ", Sys.time()),
#       file=paste0(outf,"/log.txt"), append = T)






# make haploblocks ------------------------------

f <- "D:/OneDrive - WageningenUR/polyploids/workpackage4/data/PedigreeSim/NAMcrosses"
map <- read.table(paste0(f,"/Potato.map"),
                  header = T)
# map <- sample_data$map
# map[1:5,]
# dim(map)


### binning by genetic distance
### sliding window
hb_list1 <- map2blocks(map=map,
                       winsize=0.01,
                       sldpace=0.01)

length(hb_list1)
print(table(sapply(hb_list1, length)))
names(hb_list1)[1:50]

hb_list <- refineBlocks(hb_list1,  # haploblock list
                         nmrk=2:8, # numeric vector with block length
                         mrkDosage=NULL,
                         method="random")

names(hb_list)[1:100]
print(table(sapply(hb_list, length)))

save(hb_list, file = paste0(outf,"/hb_list.RData"))


### complete block list, including blocks with one SNP
# set.seed(3)
hb_list2 <- refineBlocks(hb_list1,  # haploblock list
                         nmrk=1:8, # numeric vector with block length
                         mrkDosage=NULL,
                         method="random")

hb_list2_1snp <- hb_list2[sapply(hb_list2, length) == 1]
hb_list_0.01 <- c(hb_list2_1snp, hb_list)
save(hb_list_0.01, file = paste0(outf,"/hb_list_0.01.RData"))


### sliding window 0.1 cM
hb_list_0.1 <- map2blocks(map=map,
                          winsize=0.1,
                          sldpace=0.1)

length(hb_list_0.1)
print(table(sapply(hb_list_0.1, length)))
names(hb_list_0.1)[1:50]

set.seed(3)
hb_list_0.1a <- refineBlocks(hb_list_0.1,  # haploblock list
                             nmrk=1:8, # numeric vector with block length
                             mrkDosage=NULL,
                             method="random")

save(hb_list_0.1a, file = paste0(outf,"/hb_list_0.1.RData"))
rm(hb_list_0.1a)

### sliding window 0.5 cM
hb_list_0.5 <- map2blocks(map=map,
                          winsize=0.5,
                          sldpace=0.5)

length(hb_list_0.5)
print(table(sapply(hb_list_0.5, length)))
names(hb_list_0.5)[1:50]

set.seed(3)
hb_list_0.5a <- refineBlocks(hb_list_0.5,  # haploblock list
                             nmrk=1:8, # numeric vector with block length
                             mrkDosage=NULL,
                             method="random")

save(hb_list_0.5a, file = paste0(outf,"/hb_list_0.5.RData"))
rm(hb_list_0.5a)






# haplotyping ----------------------------

f2 <- paste0(f,"/group")


load(paste0(outf,"/hb_list_0.1.RData"))
print(table(sapply(hb_list_0.1a, length)))
hb1snp <- sapply(hb_list_0.1a, length) == 1


options(PolyHaplotyper_ahcdir=paste0("D:/polyploids/software/PolyHaplotyper/",
                                     "ahccompletelist_8snps"))


write(paste0("Start time: ", Sys.time()),
      file=paste0(outf,"/log.txt"), append = T)

NAM_PolyHaplotyper(folder=f2,
                   map=map,
                   hb_list=hb_list_0.1a[!hb1snp],
                   maxmrk=8,
                   nhb=NULL,
                   outfld = outf)

write(paste0("Finish time: ", Sys.time()),
      file=paste0(outf,"/log.txt"), append = T)


rm(hb_list_0.1a, ahccompletelist)






# preparation to QTL mapping ----------------------------------

crossfolder <- list(cross000 = "/group/1ancestral",
                    cross200 = "/group/3ancestral",
                    cross600 = "/group/7ancestral",
                    cross900 = "/group/10ancestral")


truehapres <- list()
hapres <- list()

# hb_list
load(paste0(outf,"/hb_list_0.1.RData"))

for (crs in c("cross000","cross200","cross600","cross900")) { #crosses
  
  # load data
  f2 <- paste0(f,crossfolder[[crs]])
  load(paste0(outf,"/",crs,".results_0.1.RData"))
  
  ## SNP alleles
  snpdose <- read.delim(paste0(f2,"/",crs,"_alleledose.dat"),
                        row.names = 1, check.names = F, stringsAsFactors = F)
  snpdose <- as.matrix((snpdose))
  
  ## phased genotypes
  geno <- read.delim(paste0(f2,"/",crs,"_genotypes.dat"),
                     row.names = 1, check.names = F, stringsAsFactors = F)
  
  
  
  
  # inferred haplotypes
  hapres[[crs]] <- HapCurate(results, #results of PolyHaplotyper
                            hb_list = hb_list_0.1a, #containing also 1 SNP blocks
                            snpdose = snpdose,
                            map = map,
                            ploidy = 4,
                            na.rate = 0.2)
  
  
  # true haplotypes
  truehapres[[crs]] <- TrueGeno2TrueHap(geno, hb_list_0.1a)
  
  # remove
  rm(f2,results,snpdose,geno)
}

# check
namrk <- apply(hapres[["cross200"]]$genotypes, 1, function(x) sum(is.na(x))/length(x))
plot(sort(namrk))
mono <- apply(hapres[["cross200"]]$genotypes, 1, function(x) length(unique(x)) )
plot(sort(mono))
sort(mono)[1:10]
namrk <- apply(truehapres[["cross600"]], 1, function(x) sum(is.na(x))/length(x))
plot(sort(namrk))
mono <- apply(truehapres[["cross900"]], 1, function(x) length(unique(x)) )
sort(mono)[1:10]

#save
saveRDS(truehapres, paste0(outf,"/truehaplo_res.rds"))
saveRDS(hapres, paste0(outf,"/haplo_res.rds"))


# maps for true haplotypes
truehaplomap <- data.frame(snp = sapply(hb_list_0.1a, function(x) x[1]),
                           marker = names(hb_list_0.1a))
head(truehaplomap)
truehaplomap <- merge(truehaplomap, map, by.x = "snp", by.y = "marker")[,-1]
truehaplomap <- truehaplomap[order(truehaplomap$chromosome, truehaplomap$position),]

## blocks with more than one snp
hb_list_0.1a[[4]]
blks_len <- sapply(hb_list_0.1a, length)
blks_no1snp <- names(hb_list_0.1a)[blks_len > 1] 

##
rm(hb_list_0.1a)








# QTL mapping -------------------------------

phenotypes2 <- readRDS("Phenotypes.rds")


## SNP alleles
crs <- "cross000"
f2 <- paste0(f,crossfolder[[crs]])
snpdose <- read.delim(paste0(f2,"/",crs,"_alleledose.dat"),
                      row.names = 1, check.names = F, stringsAsFactors = F)
snpdose <- as.matrix((snpdose))

## adjust phenotypes names
phenotypes2 <- lapply(phenotypes2, function(cross) {
  lapply(cross, function(h) {
    names(h) <- colnames(snpdose)
    return(h)
  })
})




## QTL positions
map[1:5,]
QTLpos <- c(453,1348,4023)
mapQTLs <- map[QTLpos,]
mapQTLs[1,3]
lqtl1 <- map[map$chromosome == mapQTLs[1,2] & 
      map$position > (mapQTLs[1,3] - 3) & 
      map$position < (mapQTLs[1,3] + 3),]
lqtl2 <- map[map$chromosome == mapQTLs[2,2] & 
      map$position > (mapQTLs[2,3] - 3) & 
      map$position < (mapQTLs[2,3] + 3),]
lqtl3 <- map[map$chromosome == mapQTLs[3,2] & 
      map$position > (mapQTLs[3,3] - 3) & 
      map$position < (mapQTLs[3,3] + 3),]
qtl1 <- map[map$chromosome == mapQTLs[1,2] & 
              map$position == mapQTLs[1,3],]
qtl2 <- map[map$chromosome == mapQTLs[2,2] & 
              map$position == mapQTLs[2,3],]
qtl3 <- map[map$chromosome == mapQTLs[3,2] & 
              map$position == mapQTLs[3,3],]


qtl1_snps <- as.character(qtl1$marker)
qtl2_snps <- as.character(qtl2$marker)
qtl3_snps <- as.character(qtl3$marker)



load(paste0(outf,"/hb_list_0.1.RData"))
hb_list_0.1a[[4]]

qtl1_blks <- names(hb_list_0.1a[sapply(hb_list_0.1a, function(x) {
  any(qtl1[,1] %in% x)
})])
qtl2_blks <- names(hb_list_0.1a[sapply(hb_list_0.1a, function(x) {
  any(qtl2[,1] %in% x)
})])
qtl3_blks <- names(hb_list_0.1a[sapply(hb_list_0.1a, function(x) {
  any(qtl3[,1] %in% x)
})])


mapQTLs
hb_list_0.1a[sapply(hb_list_0.1a, function(x) {
  any(lqtl1[,1] %in% x)
})]


load(paste0(outf,"/hb_list_0.5.RData"))
hb_list_0.5a_qtls <- hb_list_0.5a[sapply(hb_list_0.5a, function(x) {
  any(x %in% c(as.character(lqtl1[,1]),as.character(lqtl2[,1]),as.character(lqtl3[,1])))
})]

# sample(hb_list_0.5a_qtls$chr1_85,6)
hb_list_0.5a_nmrk6_qtls <- lapply(hb_list_0.5a_qtls, function(x) {
  if (length(x) > 6) {
    sample(x,6)
  } else x
})

##
rm(hb_list_0.1a, hb_list_0.5a)



# ## to check singularity in the X matrix
# check.sing <- function(geno) {
#   dets <- sapply(1:nrow(geno), function(i) {
#     #DEFINITION OF X (matrix of fixed effects)
#     cat(i,"\n")
#     if (length(unique(geno[i,])) == 1) {
#       NA
#     } else {
#       X<-dosage.X(as.matrix(geno[i,]),
#                   ploidy=4,
#                   normalize=T,
#                   haplotype = T)
#       if(ncol(X)>1){X<-X[,-1,drop=F]} #when ancestral/parental model we need to prevent singularity
#       det(t(X) %*% X)
#     }
#   }) 
# }



# d1 <- check.sing(g2)
# sort(d1)[1:50]
# plot(sort(d1[d1 < 1e-30]))
# which(d1==0)
# 
# monoa <- apply(truehapres$cross000, 1, function(x) length(unique(x)))
# d000 <- check.sing(truehapres$cross000[monoa>1,])
# d000b <- check.sing(truehapres$cross000)
# which(is.na(d000b))
# monoa <- apply(truehapres$cross200, 1, function(x) length(unique(x)))
# d200 <- check.sing(truehapres$cross200[monoa>1,])
# d200b <- check.sing(truehapres$cross200)
# 
# monoa <- apply(truehapres$cross600, 1, function(x) length(unique(x)))
# d600 <- check.sing(truehapres$cross600[monoa>1,])
# monoa <- apply(truehapres$cross900, 1, function(x) length(unique(x)))
# d900 <- check.sing(truehapres$cross900[monoa>1,])
# 
# sort(d000)[1:10]
# sort(d000b)[1:10]
# sort(d200)[1:10]
# sort(d200b)[1:10]
# sort(d600)[1:10]
# sort(d900)[1:10]
# 
# which(d000 < 1e-30)
# which(d200b < 1e-30)





## QTL mapping
crossfolder <- list(cross000 = "/group/1ancestral",
                    cross200 = "/group/3ancestral",
                    cross600 = "/group/7ancestral",
                    cross900 = "/group/10ancestral")



res.dose <- list()
res.truehap <- list()
res.fall <- list()
res.hap <- list()

for (crs in c("cross000","cross200","cross600","cross900")) { #crosses
  
  # load data
  f2 <- paste0(f,crossfolder[[crs]])
  # load(paste0(outf,"/",crs,".results_0.1.RData"))
  
  ## SNP alleles
  snpdose <- read.delim(paste0(f2,"/",crs,"_alleledose.dat"),
                        row.names = 1, check.names = F, stringsAsFactors = F)
  snpdose <- as.matrix((snpdose))
  
  ## founder alleles
  crossfile <- paste0(f2,"/",crs,"_founderalleles.dat")
  totallele <- "D:/OneDrive - WageningenUR/polyploids/workpackage4/data/PedigreeSim/Total_pop.txt"
  
  fall <- link_NAM(#function to connect ancestral founder alleles and NAM alleles
    crossfile=crossfile, #founderallele file of the cross
    marker=NULL, #index of marker or markers to get
    totallele=totallele, #pedigreeSim-type table with chromosomes and individuals, specifying alleles for each parental chromosome
    parental=F, #F:return the whole cross NAM, T:returns only the parental alleles
    parents=10,
    ploidy=4)
  
  
  # QTL mapping
  for (h in "0.5") { # c("0.2","0.5","0.8") #heritability
  
    ### QTL mapping SNP dosages
    
    res.dose[[crs]][[h]] <- map.QTL(
      phenotypes = phenotypes2[[crs]][[h]],
      genotypes = snpdose, #[monod>1,] #genotype matrix
      ploidy = 4,
      map = map, #[monod>1,] #genetic map table
      K=T, #distance matrix
      Q=NULL, #population effect matrix
      Z=NULL,
      cofactor=NULL,
      cofactor.type=NULL,
      cM=1, #
      Qpco=2, #number of axis used for pco decomposition
      no_cores=6,
      P3D=T,
      EMMAX=T,
      permutation = NULL, #permutation strategy: "pop" or "fam"
      nperm = NULL, #number of permutations
      impute=T,
      k=20
    )
    
    
    
    ## QTL mapping haplotypes
    
    ### True haplotypes
    
    res.truehap[[crs]][[h]] <- map.QTL(
      phenotypes = phenotypes2[[crs]][[h]],
      genotypes = truehapres[[crs]], # [!is.na(monoa) & monoa > 1e-30,] #genotype matrix
      ploidy = 4,
      map = truehaplomap, #[!is.na(monoa) & monoa > 1e-30,] #genetic map table
      K=T, #distance matrix
      Q=NULL, #population effect matrix
      Z=NULL,
      cofactor=NULL,
      cofactor.type=NULL,
      cM=1, #
      Qpco=2, #number of axis used for pco decomposition
      no_cores=6,
      P3D=T,
      EMMAX=T,
      permutation = NULL, #permutation strategy: "pop" or "fam"
      nperm = NULL, #number of permutations
      impute=T,
      k=20
    )
    
    
    ### Founder alleles
    
    res.fall[[crs]][[h]] <- map.QTL(
      phenotypes = phenotypes2[[crs]][[h]],
      genotypes = fall, #[monoa>1,] #genotype matrix
      ploidy = 4,
      map = map, #[monoa>1,] #genetic map table
      K=T, #distance matrix
      Q=NULL, #population effect matrix
      Z=NULL,
      cofactor=NULL,
      cofactor.type=NULL,
      cM=1, #
      Qpco=2, #number of axis used for pco decomposition
      no_cores=6,
      P3D=T,
      EMMAX=T,
      permutation = NULL, #permutation strategy: "pop" or "fam"
      nperm = NULL, #number of permutations
      impute=T,
      k=20
    )
    
    
    ### Inferred haplotypes
    #### K not provided
    ##### knn imputation
    
    res.hap[[crs]][[h]] <- map.QTL(
      phenotypes = phenotypes2[[crs]][[h]],
      genotypes = hapres[[crs]]$genotypes, #[monoa>1,] #genotype matrix
      ploidy = 4,
      map = hapres[[crs]]$map, #[monoa>1,] #genetic map table
      K=T, #distance matrix
      Q=NULL, #population effect matrix
      Z=NULL,
      cofactor=NULL,
      cofactor.type=NULL,
      cM=1, #
      Qpco=2, #number of axis used for pco decomposition
      no_cores=6,
      P3D=T,
      EMMAX=T,
      permutation = NULL, #permutation strategy: "pop" or "fam"
      nperm = NULL, #number of permutations
      impute=T,
      k=20
    )
    
    
    # #### no knn imputation
    # res.gen3.cross000 <- map.QTL(
    #   phenotypes = phenotypes2$cross000[[h]],
    #   genotypes = cross000_0.1$genotypes[monoa>1,], #genotype matrix
    #   ploidy = 4,
    #   map = cross000_0.1$map[monoa>1,], #genetic map table
    #   K=T, #distance matrix
    #   Q=NULL, #population effect matrix
    #   Z=NULL,
    #   cofactor=NULL,
    #   cofactor.type=NULL,
    #   cM=1, #
    #   Qpco=2, #number of axis used for pco decomposition
    #   no_cores=6,
    #   P3D=T,
    #   EMMAX=T,
    #   permutation = NULL, #permutation strategy: "pop" or "fam"
    #   nperm = NULL, #number of permutations
    #   impute=F,
    #   k=20
    # )
    # 
    # png(paste0(outf,"/cross000_hap_noKnn_",h,".png"))
    # skyplot(-log10(res.gen3.cross000[[1]]$pval), cross000_0.1$map[monoa>1,])
    # dev.off()
    
  }
  rm(f2,snpdose,crossfile,fall)
}





# plots

## my plotting function
# crs="cross000"
# h="0.5"

for (crs in c("cross000","cross200","cross600","cross900")) {
  for (h in c("0.5")) {
    
    ### dosages
    exclQTLs <- names(res.dose[[crs]][[h]][[1]]$pval) %in% c(qtl1_snps,qtl2_snps,qtl3_snps)
    res.dose[[crs]][[h]][[1]]$pval[exclQTLs]
    
    png(paste0(outf,"/",crs,"_snp_",h,"_bis.png"))
    skyplot(-log10(res.dose[[crs]][[h]][[1]]$pval[!exclQTLs]), map[!exclQTLs,],
            ylim = c(0,10))
    dev.off()
    
    Manhplot(scores= -log10(res.dose[[crs]][[h]][[1]]$pval[!exclQTLs]),     # a 1-column matrix of gwas p-values
             map = map[!exclQTLs,],            # custom map with 3 columns for marker, chrom, position
             Chr=NULL,            # chromosome name
             ChrCol = c("dodgerblue4","dodgerblue"),
             threshold=NULL,
             yuplim=10,
             marks=mapQTLs[,2:3],
             MapUnit="Mbp",
             filename=paste0(outf,"/",crs,"_snp_",h))
    
    
    ### founder alleles
    exclQTLs <- names(res.fall[[crs]][[h]][[1]]$pval) %in% c(qtl1_snps,qtl2_snps,qtl3_snps)
    res.fall[[crs]][[h]][[1]]$pval[exclQTLs]
    
    png(paste0(outf,"/",crs,"_fall_",h,"_bis.png"))
    skyplot(-log10(res.fall[[crs]][[h]][[1]]$pval[!exclQTLs]), map[!exclQTLs,],
            ylim = c(0,10))
    dev.off()
    
    Manhplot(scores= -log10(res.fall[[crs]][[h]][[1]]$pval[!exclQTLs]),     # a 1-column matrix of gwas p-values
             map = map[!exclQTLs,],            # custom map with 3 columns for marker, chrom, position
             Chr=NULL,            # chromosome name
             ChrCol = c("dodgerblue4","dodgerblue"),
             threshold=NULL,
             yuplim=10,
             marks=mapQTLs[,2:3],
             MapUnit="Mbp",
             filename=paste0(outf,"/",crs,"_fall_",h))
    
    # Manhplot(scores= -log10(res.fall[[crs]][[h]][[1]]$pval),     # a 1-column matrix of gwas p-values
    #          map = map,            # custom map with 3 columns for marker, chrom, position
    #          Chr=NULL,            # chromosome name
    #          ChrCol = c("dodgerblue4","dodgerblue"),
    #          threshold=NULL,
    #          yuplim=10,
    #          marks=mapQTLs[,2:3],
    #          MapUnit="Mbp",
    #          filename=paste0(outf,"/",crs,"_fall_",h,"withQTLs"))
    
    
    ### true haplotypes
    exclQTLs <- names(res.truehap[[crs]][[h]][[1]]$pval) %in% c(qtl1_blks,qtl2_blks,qtl3_blks)
    res.truehap[[crs]][[h]][[1]]$pval[exclQTLs]
    
    png(paste0(outf,"/",crs,"_truehap_",h,"_bis.png"))
    skyplot(-log10(res.truehap[[crs]][[h]][[1]]$pval[!exclQTLs]), truehaplomap[!exclQTLs,],
            ylim = c(0,10))
    dev.off()
    
    Manhplot(scores= -log10(res.truehap[[crs]][[h]][[1]]$pval[!exclQTLs]),     # a 1-column matrix of gwas p-values
             map = truehaplomap[!exclQTLs,],            # custom map with 3 columns for marker, chrom, position
             Chr=NULL,            # chromosome name
             ChrCol = c("dodgerblue4","dodgerblue"),
             threshold=NULL,
             yuplim=10,
             marks=mapQTLs[,2:3],
             MapUnit="Mbp",
             filename=paste0(outf,"/",crs,"_truehap_",h))
    
    #### excluding 1-snp blocks
    # truehaplomap[1:5,]
    blks <- truehaplomap$marker %in% blks_no1snp
    exclQTLs <- names(res.truehap[[crs]][[h]][[1]]$pval[blks]) %in% c(qtl1_blks,qtl2_blks,qtl3_blks)
    res.truehap[[crs]][[h]][[1]]$pval[blks][exclQTLs]
    
    png(paste0(outf,"/",crs,"_truehaponly_",h,"_bis.png"))
    skyplot(-log10(res.truehap[[crs]][[h]][[1]]$pval[blks][!exclQTLs]), truehaplomap[blks,][!exclQTLs,],
            ylim = c(0,10))
    dev.off()
    
    Manhplot(scores= -log10(res.truehap[[crs]][[h]][[1]]$pval[blks][!exclQTLs]),     # a 1-column matrix of gwas p-values
             map = truehaplomap[blks,][!exclQTLs,],            # custom map with 3 columns for marker, chrom, position
             Chr=NULL,            # chromosome name
             ChrCol = c("dodgerblue4","dodgerblue"),
             threshold=NULL,
             yuplim=10,
             marks=mapQTLs[,2:3],
             MapUnit="Mbp",
             filename=paste0(outf,"/",crs,"_truehaponly_",h))
    
    ### inferred haplotypes
    exclQTLs <- names(res.hap[[crs]][[h]][[1]]$pval) %in% c(qtl1_blks,qtl2_blks,qtl3_blks)
    res.hap[[crs]][[h]][[1]]$pval[exclQTLs]
    
    png(paste0(outf,"/",crs,"_hap_",h,"_bis.png"))
    skyplot(-log10(res.hap[[crs]][[h]][[1]]$pval[!exclQTLs]), hapres[[crs]]$map[!exclQTLs,],
            ylim = c(0,10))
    dev.off()
    
    Manhplot(scores= -log10(res.hap[[crs]][[h]][[1]]$pval[!exclQTLs]),     # a 1-column matrix of gwas p-values
             map = hapres[[crs]]$map[!exclQTLs,],            # custom map with 3 columns for marker, chrom, position
             Chr=NULL,            # chromosome name
             ChrCol = c("dodgerblue4","dodgerblue"),
             threshold=NULL,
             yuplim=10,
             marks=mapQTLs[,2:3],
             MapUnit="Mbp",
             filename=paste0(outf,"/",crs,"_hap_",h))
    
    #### excluding 1-snp blocks
    # hapres[[crs]]$map[1:5,]
    blks <- grep("^chr", hapres[[crs]]$map$marker)
    exclQTLs <- names(res.hap[[crs]][[h]][[1]]$pval[blks]) %in% c(qtl1_blks,qtl2_blks,qtl3_blks)
    res.hap[[crs]][[h]][[1]]$pval[blks][exclQTLs]
    
    png(paste0(outf,"/",crs,"_haponly_",h,"_bis.png"))
    skyplot(-log10(res.hap[[crs]][[h]][[1]]$pval[blks][!exclQTLs]), hapres[[crs]]$map[blks,][!exclQTLs,],
            ylim = c(0,10))
    dev.off()
    
    Manhplot(scores= -log10(res.hap[[crs]][[h]][[1]]$pval[blks][!exclQTLs]),     # a 1-column matrix of gwas p-values
             map = hapres[[crs]]$map[blks,][!exclQTLs,],            # custom map with 3 columns for marker, chrom, position
             Chr=NULL,            # chromosome name
             ChrCol = c("dodgerblue4","dodgerblue"),
             threshold=NULL,
             yuplim=10,
             marks=mapQTLs[,2:3],
             MapUnit="Mbp",
             filename=paste0(outf,"/",crs,"_haponly_",h))
    
  }
}



# focus on QTL regions -------------------------------------------
crs="cross000"
h="0.5"


res.dose.logpval <- -log10(res.dose[[crs]][[h]][[1]]$pval)
res.fall.logpval <- -log10(res.fall[[crs]][[h]][[1]]$pval)

all(map$marker == names(res.dose.logpval))
all(map$marker == names(res.fall.logpval))
pvalmap <- cbind(map, res.dose.logpval, res.fall.logpval)
pvalmap[1:10,]





qtl1_snps
pvalmap[pvalmap$marker %in% qtl1_snps,]
reg1 <- pvalmap[pvalmap$chromosome == 1 & (pvalmap$position > 74 & pvalmap$position < 78),]

qtl2_snps
pvalmap[pvalmap$marker %in% qtl2_snps,]
reg2 <- pvalmap[pvalmap$chromosome == 2 & (pvalmap$position > 55 & pvalmap$position < 65),]

qtl3_snps
pvalmap[pvalmap$marker %in% qtl3_snps,]
reg3 <- pvalmap[pvalmap$chromosome == 7 & (pvalmap$position > 0 & pvalmap$position < 10),]







res.fall[[crs]][[h]][[1]]$pval["PotVar0006158"]
res.fall[[crs]][[h]][[1]]$beta["PotVar0006158"]
res.fall[[crs]][[h]][[1]]$Fstat["PotVar0006158"]
res.fall[[crs]][[h]][[1]]$se["PotVar0006158"]
res.fall[[crs]][[h]][[1]]$wald["PotVar0006158"]
res.fall[[crs]][[h]][[1]]$real.df["PotVar0006158"]


# load data
f2 <- paste0(f,crossfolder[[crs]])
# load(paste0(outf,"/",crs,".results_0.1.RData"))

## SNP alleles
snpdose <- read.delim(paste0(f2,"/",crs,"_alleledose.dat"),
                      row.names = 1, check.names = F, stringsAsFactors = F)
snpdose <- as.matrix((snpdose))

## founder alleles
crossfile <- paste0(f2,"/",crs,"_founderalleles.dat")
totallele <- "D:/OneDrive - WageningenUR/polyploids/workpackage4/data/PedigreeSim/Total_pop.txt"

fall <- link_NAM(#function to connect ancestral founder alleles and NAM alleles
  crossfile=crossfile, #founderallele file of the cross
  marker=NULL, #index of marker or markers to get
  totallele=totallele, #pedigreeSim-type table with chromosomes and individuals, specifying alleles for each parental chromosome
  parental=F, #F:return the whole cross NAM, T:returns only the parental alleles
  parents=10,
  ploidy=4)

## phased genotypes
geno <- read.delim(paste0(f2,"/",crs,"_genotypes.dat"),
                   row.names = 1, check.names = F, stringsAsFactors = F)





cross000.par.snpdose <- snpdose[,1:10]
cross000.par.snpdose[rownames(cross000.par.snpdose) %in% reg1$marker,]

cross000.par.fall <- fall[,1:40]
cross000.par.fall[rownames(cross000.par.fall) %in% reg1$marker,]

cross000.par.phasedsnpdose <- geno[,1:40]
cross000.par.phasedsnpdose[rownames(cross000.par.phasedsnpdose) %in% reg1$marker,]

cross000.par.fall[rownames(cross000.par.fall) %in% "PotVar0006158",]
cross000.par.phasedsnpdose[rownames(cross000.par.phasedsnpdose) %in% "PotVar0006158",]
table(cross000.par.fall[rownames(cross000.par.fall) %in% "PotVar0006158",],
      cross000.par.phasedsnpdose[rownames(cross000.par.phasedsnpdose) %in% "PotVar0006158",])

reg1tab <- lapply(reg1$marker, function(x) {
  table(cross000.par.fall[rownames(cross000.par.fall) == x,],
        cross000.par.phasedsnpdose[rownames(cross000.par.phasedsnpdose) == x,])
})
names(reg1tab) <- reg1$marker

phenotypes2$cross000$`0.5`[1:10]



####################
crs="cross200"
h="0.5"


res.dose.logpval <- -log10(res.dose[[crs]][[h]][[1]]$pval)
res.fall.logpval <- -log10(res.fall[[crs]][[h]][[1]]$pval)

all(map$marker == names(res.dose.logpval))
all(map$marker == names(res.fall.logpval))
pvalmap <- cbind(map, res.dose.logpval, res.fall.logpval)
pvalmap[1:10,]





qtl1_snps
pvalmap[pvalmap$marker %in% qtl1_snps,]
reg1 <- pvalmap[pvalmap$chromosome == 1 & (pvalmap$position > 72 & pvalmap$position < 90),]
# pvalmap[pvalmap$chromosome == 1 & (pvalmap$position > 86 & pvalmap$position < 90),]

qtl2_snps
pvalmap[pvalmap$marker %in% qtl2_snps,]
reg2 <- pvalmap[pvalmap$chromosome == 2 & (pvalmap$position > 55 & pvalmap$position < 65),]

qtl3_snps
pvalmap[pvalmap$marker %in% qtl3_snps,]
reg3 <- pvalmap[pvalmap$chromosome == 7 & (pvalmap$position > 0 & pvalmap$position < 10),]







res.fall[[crs]][[h]][[1]]$pval["PotVar0006158"]
res.fall[[crs]][[h]][[1]]$beta["PotVar0006158"]
res.fall[[crs]][[h]][[1]]$Fstat["PotVar0006158"]
res.fall[[crs]][[h]][[1]]$se["PotVar0006158"]
res.fall[[crs]][[h]][[1]]$wald["PotVar0006158"]
res.fall[[crs]][[h]][[1]]$real.df["PotVar0006158"]


# load data
f2 <- paste0(f,crossfolder[[crs]])
# load(paste0(outf,"/",crs,".results_0.1.RData"))

## SNP alleles
snpdose <- read.delim(paste0(f2,"/",crs,"_alleledose.dat"),
                      row.names = 1, check.names = F, stringsAsFactors = F)
snpdose <- as.matrix((snpdose))

## founder alleles
crossfile <- paste0(f2,"/",crs,"_founderalleles.dat")
totallele <- "D:/OneDrive - WageningenUR/polyploids/workpackage4/data/PedigreeSim/Total_pop.txt"

fall <- link_NAM(#function to connect ancestral founder alleles and NAM alleles
  crossfile=crossfile, #founderallele file of the cross
  marker=NULL, #index of marker or markers to get
  totallele=totallele, #pedigreeSim-type table with chromosomes and individuals, specifying alleles for each parental chromosome
  parental=F, #F:return the whole cross NAM, T:returns only the parental alleles
  parents=10,
  ploidy=4)

## phased genotypes
geno <- read.delim(paste0(f2,"/",crs,"_genotypes.dat"),
                   row.names = 1, check.names = F, stringsAsFactors = F)





cross200.par.snpdose <- snpdose[,1:10]
cross200.par.snpdose[rownames(cross200.par.snpdose) %in% reg1$marker,]

cross200.par.fall <- fall[,1:40]
cross200.par.fall[rownames(cross200.par.fall) %in% reg1$marker,]

cross200.par.phasedsnpdose <- geno[,1:40]
cross200.par.phasedsnpdose[rownames(cross200.par.phasedsnpdose) %in% reg1$marker,]

# cross200.par.fall[rownames(cross200.par.fall) %in% "PotVar0006158",]
# cross200.par.phasedsnpdose[rownames(cross200.par.phasedsnpdose) %in% "PotVar0006158",]
# table(cross200.par.fall[rownames(cross200.par.fall) %in% "PotVar0006158",],
#       cross200.par.phasedsnpdose[rownames(cross200.par.phasedsnpdose) %in% "PotVar0006158",])

reg1tab <- lapply(reg1$marker, function(x) {
  table(cross200.par.fall[rownames(cross200.par.fall) == x,],
        cross200.par.phasedsnpdose[rownames(cross200.par.phasedsnpdose) == x,])
})
names(reg1tab) <- reg1$marker

c("PotVar0006158","PotVar0006356","PotVar0029143")
phenotypes2$cross200$`0.5`[1:10]

