

# random data -----------------------------
pheno <- rnorm(10)
names(pheno) <- paste0("ind",1:10)
pheno

snpdose <- matrix(sample(0:4, 10*3, replace = TRUE), nrow = 3)
rownames(snpdose) <- paste0("snp",1:3)
colnames(snpdose) <- names(pheno)
snpdose

hapdose <- matrix(sample(letters[1:4], 40*3, replace = TRUE), nrow = 3)
rownames(hapdose) <- paste0("snp",1:3)
colnames(hapdose) <- sapply(names(pheno), function(x) paste(x, 1:4, sep = "_"))
hapdose

hapdose_2alleles <- matrix(sample(letters[1:2], 40*3, replace = TRUE), nrow = 3)
rownames(hapdose_2alleles) <- paste0("snp",1:3)
colnames(hapdose_2alleles) <- sapply(names(pheno), function(x) paste(x, 1:4, sep = "_"))
hapdose_2alleles

map <- data.frame(marker = rownames(hapdose),
                  chromosome = rep(1,3),
                  position = 1:3)
map



##
dosage.X(genotypes = snpdose[1,], ploidy = 4)
dosage.X(genotypes = hapdose[1,], haplotype = TRUE, ploidy = 4)
dosage.X(genotypes = hapdose_2alleles[1,], haplotype = TRUE, ploidy = 4)



## snpdose linear
results <- map.QTL(phenotypes = pheno,
                   genotypes = snpdose,
                   ploidy = 4,
                   map = map)

results$pheno1$beta

## snpdose mixed
results <- map.QTL(phenotypes = pheno,
                   genotypes = snpdose,
                   ploidy = 4,
                   map = map,
                   K = TRUE,
                   K_identity = TRUE)

results$pheno1$beta

## hapdose linear
results <- map.QTL(phenotypes = pheno,
                   genotypes = hapdose,
                   ploidy = 4,
                   map = map)

results$pheno1$beta

## hapdose mixed
results <- map.QTL(phenotypes = pheno,
                   genotypes = hapdose,
                   ploidy = 4,
                   map = map,
                   K = TRUE,
                   K_identity = TRUE)

results$pheno1$beta




# mpqtl example data -----------------------------------
data("mpsnpdose")
data("mphapdose")
data("mppheno")
data("mpmap")

# source("../backup_R_folder_20210510/R/mapQTL_fun.R")



results <- map.QTL(phenotypes = mppheno,
                   genotypes = mpsnpdose,
                   ploidy = 4,
                   map = mpmap)

results$phenotype1$beta[1:3]

results <- map.QTL(phenotypes = mppheno,
                   genotypes = mpsnpdose,
                   ploidy = 4,
                   map = mpmap,
                   K = TRUE,
                   K_identity = TRUE)

results$phenotype1$beta[1:3]

results <- map.QTL(phenotypes = mppheno,
                   genotypes = mphapdose,
                   ploidy = 4,
                   map = mpmap)

results$phenotype1$beta[1:3]

results <- map.QTL(phenotypes = mppheno,
                   genotypes = mphapdose,
                   ploidy = 4,
                   map = mpmap,
                   K = TRUE,
                   K_identity = TRUE)

results$phenotype1$beta[1:3]
