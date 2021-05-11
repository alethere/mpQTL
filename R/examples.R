
# map.QTL ---------------------------------
# ## Get example genotypes (haplotypes in this case),
# ## map and phenotypes for a population of tetraploid individuals
# data("mphapdose")
# data("mpmap")
# data("mppheno")
#
# ## fixed effects model (y = m + e) LINEAR
# results <- map.QTL(phenotypes = mppheno,
#                    genotypes = mphapdose,
#                    ploidy = 4,
#                    map = mpmap)
# names(results$phenotype1)
# skyplot(-log10(results$phenotype1$pval), map = mpmap)
#
# ## naive model (y = m + e) - no K correction
# results <- map.QTL(phenotypes = mppheno,
#                    genotypes = mphapdose,
#                    ploidy = 4,
#                    map = mpmap,
#                    K = TRUE,
#                    K_identity = TRUE)
# names(results$phenotype1)
# skyplot(-log10(results$phenotype1$pval), map = mpmap)
#
# ## model with Q correction (y = m + Q + e) LINEAR
# results <- map.QTL(phenotypes = mppheno,
#                    genotypes = mphapdose,
#                    ploidy = 4,
#                    map = mpmap,
#                    Q = TRUE,
#                    Qpco = 2)
# names(results$phenotype1)
# results$phenotype1$beta[1]
# skyplot(-log10(results$phenotype1$pval), map = mpmap)
#
# ## model with Q correction (y = m + Q + e)
# results <- map.QTL(phenotypes = mppheno,
#                    genotypes = mphapdose,
#                    ploidy = 4,
#                    map = mpmap,
#                    Q = TRUE,
#                    Qpco = 2,
#                    K = TRUE,
#                    K_identity = TRUE)
# names(results$phenotype1)
# results$phenotype1$beta[1]
# skyplot(-log10(results$phenotype1$pval), map = mpmap)
#
# ## model with K correction (y = m + K + e)
# results <- map.QTL(phenotypes = mppheno,
#                    genotypes = mphapdose,
#                    ploidy = 4,
#                    map = mpmap,
#                    K = TRUE)
# names(results$phenotype1)
# skyplot(-log10(results$phenotype1$pval), map = mpmap)
#
# ## model with Q + K correction (y = m + Q + K + e)
# results <- map.QTL(phenotypes = mppheno,
#                    genotypes = mphapdose,
#                    ploidy = 4,
#                    map = mpmap,
#                    Q = TRUE,
#                    Qpco = 2,
#                    K = TRUE)
# names(results$phenotype1)
# skyplot(-log10(results$phenotype1$pval), map = mpmap)



# calc.K ------------------------------------
# ## Create SNP dosages for 10 tetraploid individuals (in rows) and 100 markers
# snpdose <- matrix(sample(0:4,1000, replace=TRUE), nrow=10)
#
# ## Calculate K matrix
# K <- calc.K(snpdose)
#
# ## Create an haplotype matrix for 5 tetraploid individuals (in rows) and 100 markers
# ## Each individual genotype is represented by four rows (its four homologues)
# hapdose <- matrix(sample(1:6,2000, replace=TRUE), nrow=20)
#
# ## K matrix
# K <- calc.K(hapdose, haplotype = TRUE, ploidy = 4)
#
# ## Haplotypes can be named by characters as well
# hapdose <- matrix(letters[hapdose], nrow=20)
# K <- calc.K(hapdose, haplotype = TRUE, ploidy = 4)



# impute.knn -------------------------------
# ## Get simulated genotypes for tetraploid individuals
# data("mpsnpdose") # SNP dosages
# mpsnpdose[sample(1:length(mpsnpdose),500)] <- NA  # add NAs randomly
# data("mphapdose") # haplotypes
# mphapdose[sample(1:length(mphapdose),500)] <- NA  # add NAs randomly
#
# ## Imputation
# mpsnpdose <- impute.knn(mpsnpdose, ploidy = 4, kneighbors = 50)
# mphapdose <- impute.knn(mphapdose, ploidy = 4, kneighbors = 50)



# sample.marker --------------------------------
# ## Create a random map
# map <- data.frame(marker = paste0("mrk", 1:100),
#                   chromosome = c(rep(1,50), rep(2,50)),
#                   position = c(sort(runif(50,0,60)), sort(runif(50,0,80))))
#
# ## Create a matrix of random genotypes, with markers in rows
# snpdose <- matrix(sample(0:4,1000, replace = TRUE), nrow = 100)
# rownames(snpdose) <- map$marker
#
# ## Sample evenly spaced markers
# sample.marker(genotypes = snpdose, map, binsize=1, seed=3)
#
# ## A vector of marker names can be provided instead of genotypes
# sample.marker(genotypes = map$marker,
#               map, binsize=1, seed=3)



# dosage.X ---------------------------------
# ## SNP dosages of 10 tetraploid individuals
# snpdose <- sample(0:4, 10, replace = TRUE)
# dosage.X(snpdose)
#
# ## Haplotypes of 10 tetraploid individuals for a locus with 5 haplotypes
# hapdose <- sample(1:5,40, replace = TRUE)
# dosage.X(genotypes = hapdose, haplotype = TRUE, ploidy = 4)


# pcoa.plot --------------------------------
# ## Get example SNP dosages
# data("mpsnpdose")
# ## Calculate distances between individuals
# K <- as.matrix(dist(t(mpsnpdose)))
# ## Plot Principal Coordinates
# pcoa.plot(K)
#
# ## Colors according to a grouping variable
# pop <- substr(colnames(K),1,2)
# pcoa.plot(K, col = pop)
# pcoa.plot(K, col = pop, h = c(0,360))


# thr.LiJi -------------------------------
# ## Get example SNP dosages and a map
# data("mpsnpdose")
# data("mpmap")
#
# mpsnpdose[1:5,1:5]
# head(mpmap)
# tapply(mpmap$marker, mpmap$chromosome, length)
# identical(rownames(mpsnpdose), mpmap$marker)
#
# ## excluding unmapped markers (chromosome 0)
# mapped <- mpmap$chromosome != 0
# thr <- thr.LiJi(m = mpsnpdose[mapped,],
#                 chrom = mpmap$chromosome[mapped])
# -log10(thr$threshold) #log scale


# LD_decay ----------------------------
# ## Get example SNP dosages and a map
# data("mpsnpdose")
# data("mpmap")
#
# mpsnpdose[1:5,1:5]
# head(mpmap)
# tapply(mpmap$marker, mpmap$chromosome, length)
# identical(rownames(mpsnpdose), mpmap$marker)
#
# ## LD decay calculation
# LD <- LD_decay(mpsnpdose, mpmap, win_size = 2,
#                max_dist = 50, percentile = c(0.9,0.95))
# LD <- LD_decay(mpsnpdose, mpmap, win_size = 2,
#                per_chr = TRUE,
#                max_dist = 50, percentile = c(0.9,0.95))


# plot_LD ------------------------------
# # ## Get example SNP dosages and a map
# data("mpsnpdose")
# data("mpmap")
#
# mpsnpdose[1:5,1:5]
# head(mpmap)
# tapply(mpmap$marker, mpmap$chromosome, length)
# identical(rownames(mpsnpdose), mpmap$marker)
#
# ## LD decay calculation and plotting
# LD <- LD_decay(mpsnpdose, mpmap, win_size = 1,
#                max_dist = 100, percentile = c(0.9,0.95))
#
# plot_LD(LD, max_dist = NULL, main=NULL)
#
# LD <- LD_decay(mpsnpdose, mpmap, win_size = 1,
#                max_dist = 100, percentile = c(0.9,0.95), per_chr = TRUE)
# plot_LD(LD[[1]])



# HapdoseToHapname -----------------------
# ## Create random haplotype dosage data for ten tetraploid
# ## individuals and three loci.
# hapdose <- lapply(1:3, function(l) {
#   sapply(1:10, function(x) {
#     tabulate(sample(1:3, 4, replace = TRUE), nbins = 3)
#   })
# })
# names(hapdose) <- paste0("locus",1:3)
#
# for(i in 1:3){ #add names to markers and individuals
#   d <- dim(hapdose[[i]])
#   rownames(hapdose[[i]]) <- paste0("hap",1:d[1])
#   colnames(hapdose[[i]]) <- paste0("ind",1:d[2])
# }
#
# colSums(hapdose$locus1) #each column sums to ploidy
#
# ## Convert to Hapname format
# hapname <- HapdoseToHapname(hapdose, ploidy = 4)
#
# ## and viceversa
# HapnameToHapdose(hapname)



# dosP2Bfreq --------------------------
# ## Create an example matrix of dosage probabilities for a tetraploid individual
# dosP <- matrix(c(1,0,0,0,0,
#                  0,1,0,0,0,
#                  0,0,0,0,1,
#                  0,0,0.5,0.5,0,
#                  0,0,0.2,0.8,0,
#                  0,0.95,0.05,0,0), ncol = 5, byrow = TRUE)
# colnames(dosP) <- paste0("D",0:4)
#
# ## Convert dosage probabilities into B allele probabilities
# dosP2Bfreq(dosP)


# QQ.plot ------------------
# ## One set of random p-values
# pvals <- runif(100)
# QQ.plot(pvals)
#
# ## two sets of p-values
# pvals <- matrix(runif(200), ncol = 2) # matrix
# colnames(pvals) <- paste0("pval",1:2)
# QQ.plot(pvals)
#
# pvals <- list(pval1 = runif(100), # list
#               pval2 = runif(100))
# QQ.plot(pvals)


# skyplot -----------------------
# ## Create a map and set of random p-values
# map <- data.frame(marker = paste0("mrk", 1:100),
#                   chromosome = c(rep(1,50), rep(2,50)),
#                   position = c(sort(runif(50,0,60)), sort(runif(50,0,80))))
# pvals <- runif(100)
#
# ## Plot -log10(pvals)
# skyplot(pval = -log10(pvals), map)
# skyplot(pval = -log10(pvals), map, small = FALSE)
# skyplot(pval = -log10(pvals), map, small = FALSE, threshold = 3)
# skyplot(pval = -log10(pvals), map, small = FALSE, threshold = 3, ylim = c(0,4))



# comp.skyplot -----------------------
# ## Create a map and two sets of random p-values
# map <- data.frame(marker = paste0("mrk", 1:100),
#                   chromosome = c(rep(1,50), rep(2,50)),
#                   position = c(sort(runif(50,0,60)), sort(runif(50,0,80))))
#
# pvals <- matrix(runif(200), ncol = 2) # matrix
# colnames(pvals) <- paste0("pval",1:2)
#
# ## Plot -log10(pvals)
# comp.skyplot(-log10(pvals), map)
# comp.skyplot(-log10(pvals), map, small = FALSE)
# comp.skyplot(pval = -log10(pvals), map, small = FALSE, threshold = 3, ylim = c(0,4))



# pheno_box ---------------------
# ## Get example genotypes and phenotypes of tetraploid individuals
# data("mpsnpdose")
# mpsnpdose[1:5,1:5]
# data("mphapdose")
# mphapdose[1:5,1:9]
# data("mppheno")
# head(mppheno)
#
# ## pheno_box of the first SNP
# pheno_box(phe = mppheno, gen = mpsnpdose[1,])
# pheno_box(phe = mppheno, gen = mphapdose[1,],
#           haplotype = TRUE, ploidy = 4)


# hue_wheel --------------------------
# hue_wheel() #luminance = 60 by default
# hue_wheel(90) #luminance = 90
