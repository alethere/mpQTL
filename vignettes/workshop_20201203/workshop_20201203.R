

# Load mpQTL ----------------------------------------------
library(mpQTL)



# Load example data ------------------------
load("workshop_data.RData")

# it contains
head(mrksel)
head(indsel)
head(map)
head(phenotypes)


# Chunkwise data import by 'ff' -----------------------------
library("ff")
library("ffbase")

## preview
readLines("scores.dat", 10)

scores_preview <- read.delim.ffdf(file="scores.dat",
                                  nrows = 10)
scores_preview


## import data and subset based on marker names (or sample names)
names_ffdf <- read.delim.ffdf(file="scores.dat",
                                     transFUN = function(x, n = mrksel, m = indsel) {
                                       idx <- x[,2] %in% n & x[,3] %in% m
                                       x[idx, 1:3]
                                     })
names_ffdf
rownames(names_ffdf) <- NULL

data_ffdf <- read.delim.ffdf(file="scores.dat",
                             transFUN = function(x, n = mrksel, m = indsel){
                               idx <- x[,2] %in% n & x[,3] %in% m
                               x[idx, c(4:9,12)]
                             })
data_ffdf
rownames(data_ffdf) <- NULL

data_ffmat <- ffbase:::as.ff_matrix.ffdf(data_ffdf)
data_ffmat



## get the RAM object
ramattribs(data_ffmat[1:10,])
data <- data_ffmat[]
data[1:5,]
names <- names_ffdf[]
head(names)

# ## save and load
# ffsave(data_ffmat, file=""vignettes/workshop_20201203/ff_files/data_ffmat")
# ffload("vignettes/workshop_20201203/ff_files/data_ffmat")

## delete ff files
delete(names_ffdf, data_ffdf, data_ffmat, scores_preview)
rm(names_ffdf, data_ffdf, data_ffmat, scores_preview)
gc()
detach("package:ffbase", unload=TRUE)
detach("package:ff", unload=TRUE)
detach("package:bit", unload=TRUE)


# Transform and reshape data ---------------------------
# + reshape ratios -------------------
ratios <- matrix(data[,1], nrow = length(unique(names$MarkerName)), byrow = T)
dim(ratios)
ratios <- round(ratios,4)

rownames(ratios) <- unique(names$MarkerName)
colnames(ratios) <- unique(names$SampleName)

ratios[1:5,1:5]

# + dosage probabilities ------------------
Bfreq <- dosP2Bfreq(data[,2:6])
dim(Bfreq)
Bfreq[1:9]

attributes(Bfreq)
class(Bfreq)
typeof(Bfreq)
Bfreq <- round(Bfreq,4)

Bfreq <- matrix(Bfreq, nrow = length(unique(names$MarkerName)), byrow = T)
dim(Bfreq)
Bfreq[1:5,1:5]

rownames(Bfreq) <- unique(names$MarkerName)
colnames(Bfreq) <- unique(names$SampleName)


Bfreq[1:5,1:5]
range(c(Bfreq))

# Simulate filtration steps ----------------------------
set.seed(3)
naidx <- sample(1:length(ratios), 100)
ratios[naidx] <- NA
Bfreq[naidx] <- NA

# Kinship calculation and genotype imputation ------------------
## + Kinship matrix -----------
## Used for knn imputation or structure correction.
## To calculate a kinship using a subset of markers see sample.cM
Kr <- calc.K(t(ratios))
Kr[1:5,1:5]

Kp <- calc.K(t(Bfreq))
Kp[1:5,1:5]


## visualize kinship
# pcoa.plot(Kr)

## + geno imputation --------------------
## Since we provide K, map is not needed
sum(is.na(ratios))
ratios2 <- impute.knn(ratios,
                      ploidy = 4,
                      kneighbors = 20,
                      K=Kr)
ratios2[1:5,1:5]

sum(is.na(Bfreq))
Bfreq2 <- impute.knn(Bfreq,
                     ploidy = 4,
                     kneighbors = 20,
                     K=Kp)


# Run association models --------------------------------
# + Input matrices order -------------------
ordInp_ratios <- inputOrder(ratios2, pheno = phenotypes, map = map, K = Kr)
str(ordInp_ratios)

ordInp_Bfreq <- inputOrder(Bfreq2, pheno = phenotypes, map = map, K = Kp)


# + Li & Ji threshold --------------------
ratios_thr <- thr.LiJi(ordInp_ratios$geno,
                       chrom = ordInp_ratios$map$chromosome,
                       ploidy = 4)
ratios_thr




Bfreq_thr <- thr.LiJi(ordInp_Bfreq$geno,
                       chrom = ordInp_Bfreq$map$chromosome,
                       ploidy = 4)

# + Naive model and K model -------------------------
## using ratios
ratios_naive <- map.QTL(phenotypes = ordInp_ratios$pheno,
                        genotypes = ordInp_ratios$geno,
                        ploidy = 4,
                        K = ordInp_ratios$K,
                        K_identity = T, #to force naive
                        # impute = F, #not needed, since we provided imputed genotypes
                        map = ordInp_ratios$map,
                        no_cores = 6)

ratios_K <- map.QTL(phenotypes = ordInp_ratios$pheno,
                    genotypes = ordInp_ratios$geno,
                    ploidy = 4,
                    K = ordInp_ratios$K,
                    map = ordInp_ratios$map,
                    no_cores = 6)




## using dosage probabilities
Bfreq_naive <- map.QTL(phenotypes = ordInp_Bfreq$pheno,
                        genotypes = ordInp_Bfreq$geno,
                        ploidy = 4,
                        K = ordInp_Bfreq$K,
                        K_identity = T, #to force naive
                        # impute = F, #not needed, since we provided imputed genotypes
                        map = ordInp_Bfreq$map,
                        no_cores = 6)

Bfreq_K <- map.QTL(phenotypes = ordInp_Bfreq$pheno,
                   genotypes = ordInp_Bfreq$geno,
                   ploidy = 4,
                   K = ordInp_Bfreq$K,
                   map = ordInp_Bfreq$map,
                   no_cores = 6)

# + Manhattan plots -------------------------
## plots for ratios
skyplot(-log10(ratios_naive$pheno01$pval),
        map = ordInp_ratios$map,
        chromspace = 0,
        threshold = -log10(ratios_thr$threshold),
        main="naive model using ratios")

skyplot(-log10(ratios_K$pheno01$pval),
        map = ordInp_ratios$map,
        chromspace = 0,
        threshold = -log10(ratios_thr$threshold),
        main="K model using ratios")


## plots for dosage prob.
skyplot(-log10(Bfreq_naive$pheno01$pval),
        map = ordInp_Bfreq$map,
        chromspace = 0,
        threshold = -log10(Bfreq_thr$threshold),
        main="naive model using dosage prob.")

skyplot(-log10(Bfreq_K$pheno01$pval),
        map = ordInp_Bfreq$map,
        chromspace = 0,
        threshold = -log10(Bfreq_thr$threshold),
        main="K model using dosage prob.")

# + Genotype-phenotype plots -----------------------------
## peak marker for ratios
ratios_K$pheno01$pval[which.max(-log10(ratios_K$pheno01$pval))]
plot(ordInp_ratios$geno["mrk050949",],
     ordInp_ratios$pheno[,1],
     xlim = c(0,1),
     xlab = "ratios",
     ylab = "phenotype")

## peak marker for dosage prob.
Bfreq_K$pheno01$pval[which.max(-log10(Bfreq_K$pheno01$pval))]
plot(ordInp_Bfreq$geno["mrk050949",],
     ordInp_Bfreq$pheno[,1],
     xlim = c(0,1),
     xlab = "dosage Prob.",
     ylab = "phenotype")


# B allele frequency plots ---------------------------
library(scales)

ordMapRatio <- inputOrder(ratios[rownames(ratios) %in% map$marker[map$chromosome != 0],], map = map)
# colnames(ordMapRatio$geno)
# str(ordMapRatio)

skyplot(ordMapRatio$geno[,"ind010"],
        map = ordMapRatio$map, #[physmap$physical_chr!=0,]
        ylim = c(0,1.1),
        pch = 20,
        cex = 0.9,
        ylab = "ratios",
        col = alpha(c("dodgerblue4","dodgerblue"),0.2),
        chromspace = 0,
        main="ind010")

skyplot(ordMapRatio$geno[,"ind011"],
        map = ordMapRatio$map, #[physmap$physical_chr!=0,]
        ylim = c(0,1.1),
        pch = 20,
        cex = 0.9,
        ylab = "ratios",
        col = alpha(c("dodgerblue4","dodgerblue"),0.2),
        chromspace = 0,
        main="ind011")



hist(ordMapRatio$geno[,"ind010"],
     breaks = 500,
     xlim = c(0,1),
     main="ind010",
     xlab = "ratios")
hist(ordMapRatio$geno[,"ind011"],
     breaks = 500,
     xlim = c(0,1),
     main="ind011",
     xlab = "ratios")

ordMapBfreq <- inputOrder(Bfreq[rownames(Bfreq) %in% map$marker[map$chromosome != 0],], map = map)

skyplot(ordMapBfreq$geno[,"ind011"],
        map = ordMapBfreq$map, #[physmap$physical_chr!=0,]
        ylim = c(0,1.1),
        pch = 20,
        cex = 0.9,
        ylab = "Dosage Prob.",
        col = alpha(c("dodgerblue4","dodgerblue"),0.2),
        chromspace = 0,
        main="ind011")

