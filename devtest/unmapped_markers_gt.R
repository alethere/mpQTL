


# prepare map and genotypes -----------------------

## map is a data.frame with 3 columns ("marker","chromosome","position")
## I would suggest to include in this map also a chromosome 0 containing unmapped markers
## Add also a fake position for unmapped markers using e.g. seq(0,100, length.out = length(unmapped markers))

## subset the map and genotype table, so that they include only mapped markers
mapped <- map[map$chromosome != 0,]
genotypes_mapdmrkr <- genotypes[match(mapped$marker, rownames(genotypes)),]
dim(genotypes_mapdmrkr)



# calculate kinship -------------------------------------
K <- sample.cM(genotypes_mapdmrkr,mapped,cM = 1) #sample markers every 1 cM
K <- calc.K(t(K),ploidy = 4,haplotypes = F) #calculate a kinship matrix



# QTL mapping --------------------------------------
res <- map.QTL(
  phenotypes = phenotype, #phenotype matrix
  genotypes = genotypes,  #genotype containg mapped and unmapped markers
  ploidy = 4,
  map = NULL, 
  K = K,  #here provide your kinship 
  P3D = T #approximation, to reduce computation time
)


# plot -------------------------------------------
skyplot(pval = -log10(res$yourPheno$pval),
        map = map) #here use the complete map, with mapped and unmapped markers