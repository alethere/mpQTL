
#### prepare workshop

devtools::load_all()



# write tab sep data PhenoGeno -------------------------------
load("C:/Users/tumin001/OneDrive - WageningenUR/polyploids/rose/GWASdata_PhenoGeno/geno/data/20161216-175242_scores.RData")
head(scores)

## anonymize data
sample_names <- as.character(unique(scores$SampleName))
indnames <- sprintf("ind%03d", 1:length(sample_names))
indnamestab <- data.frame(orig = sample_names,
                          anonym = indnames)
head(indnamestab)

rm(sample_names)

marker_names <- as.character(unique(scores$MarkerName))
length(marker_names)
mrknames <- sprintf("mrk%06d", 1:length(marker_names))
mrknamestab <- data.frame(orig = marker_names,
                          anonym = mrknames)
head(mrknamestab)
rm(marker_names)


head(scores)
scores$MarkerName <- mrknamestab$anonym[match(scores$MarkerName, mrknamestab$orig)]
scores$SampleName <- indnamestab$anonym[match(scores$SampleName, indnamestab$orig)]






## prepare map
map <- readRDS("C:/Users/tumin001/OneDrive - WageningenUR/polyploids/rose/Rose_Korea/Map/physicalSNPposition/map_gen_phys.rds")
head(map)
map <- map[order(map$physical_chr, map$physical_pos),]
map[1:5,]
unique(map$physical_chr)

physmap <- map[,c("MarkerName","physical_chr","physical_pos")]
physmap$physical_pos <- physmap$physical_pos/1e6
physmap[1:10,]
sum(is.na(physmap$physical_chr))
sum(is.na(physmap$physical_pos))
physmap$physical_chr[is.na(physmap$physical_chr)] <- 0
physmap$physical_pos[physmap$physical_chr==0] <- seq(1,50,length.out = sum(physmap$physical_chr==0))
sum(is.na(physmap$physical_chr))
sum(is.na(physmap$physical_pos))

colnames(physmap) <- c("marker","chromosome","position")
head(physmap)

### anonimyze marker names of map
physmap$marker <- mrknamestab$anonym[match(physmap$marker, mrknamestab$orig)]
head(physmap)

map <- physmap
rm(physmap)

# saveRDS(map, "vignettes/workshop_20201203/map.rds")

## sample names and marker names for subsetting
set.seed(3)
mrksel <- sample(mrknames, 10000)
indsel <- indnames[c(1:19,399:478)]

## phenotypes
phenotypes <- read.csv("C:/Users/tumin001/OneDrive - WageningenUR/polyploids/rose/GWASdata_PhenoGeno/pheno/petalnum_accessionMeans_for_workshop.csv",
                  stringsAsFactors = F)
head(phenotypes)

phenotypes <- phenotypes[, c("genoaccs","mean","sd")]
indnamestab[indnamestab$anonym %in% indsel,]

phenotypes$anonym <- indnamestab$anonym[match(phenotypes$genoaccs, indnamestab$orig)]
phenotypes <- phenotypes[, c(4,2,3)]
colnames(phenotypes) <- c("SampleName","pheno01","pheno02")

phenotypes <- phenotypes[phenotypes$SampleName %in% indnamestab$anonym,]
phenotypes[,2:3] <- round(phenotypes[,2:3],1)
rownames(phenotypes) <- phenotypes$SampleName
phenotypes <- phenotypes[,2,drop=F]

# saveRDS(phenotypes, "vignettes/workshop_20201203/phenotypes.rds")


mrksel <- sample(mrknames, 10000)
indsel <- indnames[c(1:19,399:478)]

mrksel_dat <- c(mrksel, sample(mrknames[!mrknames %in% mrksel], 5000))
sum(duplicated(mrksel_dat))

indsel_dat <- c(indsel, sample(indnames[!indnames %in% indsel], 20))
sum(duplicated(indsel_dat))



head(scores)
scores <- scores[scores$MarkerName %in% mrksel_dat,]
scores <- scores[scores$SampleName %in% indsel_dat,]


write.table(scores, file = "vignettes/workshop_20201203/scores.dat",
            sep = "\t", row.names = F, quote = F) # write a tab separated file (to test chunkwise access)
rm(scores)



save(mrksel,indsel,map,phenotypes,
     file="vignettes/workshop_20201203/workshop_data.RData")



# ff approch: subset during import -----------------------------
## Not all the data might be necessary (e.g. one could need to focus
## on a subset of individuals or to use only well performing markers
## or mapped markers).

library("ff")
library("ffbase")

## preview
readLines("vignettes/workshop_20201203/scores.dat", 10)

scores_preview <- read.delim.ffdf(file="vignettes/workshop_20201203/scores.dat",
                                  nrows = 10)
scores_preview


# subset based on marker names (or sample names)
names_ffdf <- read.delim.ffdf(file="vignettes/workshop_20201203/scores.dat",
                                     transFUN = function(x, n = mrksel, m = indsel) {
                                       idx <- x[,2] %in% n & x[,3] %in% m
                                       x[idx, 1:3]
                                     })
names_ffdf
rownames(names_ffdf) <- NULL

# names_ffdf <- read.delim.ffdf(file="vignettes/workshop_20201203/20180712-234650_scores.dat",
#                                # nrows = 10000,
#                                transFUN = function(x){
#                                  x[,1:3]
#                                })
#
#

data_ffdf <- read.delim.ffdf(file="vignettes/workshop_20201203/scores.dat",
                             transFUN = function(x, n = mrksel, m = indsel){
                               idx <- x[,2] %in% n & x[,3] %in% m
                               x[idx, c(4:9,12)]
                             })
data_ffdf
rownames(data_ffdf) <- NULL

# data_ffdf <- read.delim.ffdf(file="vignettes/workshop_20201203/20180712-234650_scores.dat",
#                                # nrows = 10000,
#                                transFUN = function(x){
#                                  x[,c(4:9,12)]
#                                })

data_ffmat <- ffbase:::as.ff_matrix.ffdf(data_ffdf)
data_ffmat



# ## check each marker has the same number of individuals and
# ## each individual has the same number of markers
# mrk_len <- table(names_ffdf$MarkerName)
# unique(mrk_len)
# ind_len <- table(names_ffdf$SampleName)
# unique(ind_len)
#
# ## get alphabetical order based on Markername and Samplename
# newidx <- ffdforder(names_ffdf[c("MarkerName","SampleName")])
#
# names_ffdf <- names_ffdf[newidx,]
# data_ffmat <- data_ffmat[newidx[],]


## If no subscript selection is needed,
## use data_ffmat[] to get the ram object.
## Make sure that row names are not inherited by the ff object
## using ramattribs(data_ffmat[1:10,2:6])
ramattribs(data_ffmat[1:10,])
data <- data_ffmat[]
data[1:5,]
names <- names_ffdf[]
head(names)
# saveRDS(data, "vignettes/workshop_20201203/data.rds")


## remove ff objects and clear reference to flat files
# ## save and load
# ffsave(data_ffmat, file=""vignettes/workshop_20201203/ff_files/data_ffmat")
# ffload("vignettes/workshop_20201203/ff_files/data_ffmat")
delete(names_ffdf, data_ffdf, data_ffmat, scores_preview)
rm(names_ffdf, data_ffdf, data_ffmat, scores_preview)
gc()
detach("package:ffbase", unload=TRUE)
detach("package:ff", unload=TRUE)
detach("package:bit", unload=TRUE)


## transform dosage probabilities
# Bfreq <- dosP2Bfreq(data_ffmat[,2:6])
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


# Bfreq <- Bfreq[match(integrmap_dos$marker, rownames(Bfreq)),]


## reshape ratios (snp array intensities)
ratios <- matrix(data[,1], nrow = length(unique(names$MarkerName)), byrow = T)
dim(ratios)
ratios[1:9]
ratios <- round(ratios,4)

rownames(ratios) <- unique(names$MarkerName)
colnames(ratios) <- unique(names$SampleName)
ratios[1:5,1:5]


## Filtration steps (not included in mpQTL) might add NAs
set.seed(3)
naidx <- sample(1:length(ratios), 100)
ratios[naidx] <- NA
Bfreq[naidx] <- NA


## save
# saveRDS(ratios, "vignettes/workshop_20201203/ratios.rds")
# saveRDS(Bfreq, "vignettes/workshop_20201203/Bfreq.rds")
# ratios0 <- readRDS("vignettes/workshop_20201203/ratios.rds")




# Association study ----------


## + Kinship matrix -----------
## (for knn imputation or structure correction)
## using all the markers (otherwise use sample.cM)
Kr <- calc.K(t(ratios))
Kr[1:5,1:5]

Kp <- calc.K(t(Bfreq))
Kp[1:5,1:5]


# ### visualize kinship
# pcoa.plot(Kr)

# saveRDS(Kr, "vignettes/workshop_20201203/Kr.rds")
# saveRDS(Kp, "vignettes/workshop_20201203/Kp.rds")


## + geno imputation --------------------
## with continuos genotypes impute.knn will use the mean of kneighbors.
## Since we provide K, there is non need to provide a map
sum(is.na(ratios))
ratios2 <- impute.knn(ratios,
                      ploidy = 4,
                      kneighbors = 20,
                      K=Kr)
ratios2[1:5,1:5]


# ratios22 <- imputeNA(ratios)
# ratios0[naidx]
# plot(ratios0[naidx], ratios2[naidx])
# plot(ratios2[naidx], ratios22[naidx])
sum(is.na(ratios))


sum(is.na(Bfreq))
Bfreq2 <- impute.knn(Bfreq,
                     ploidy = 4,
                     kneighbors = 20,
                     K=Kp)
Bfreq2[1:5,1:5]


# + run models -----------------------
# ++ ratios ----------------

## input order
ordInp_ratios <- inputOrder(ratios2, pheno = phenotypes, map = map, K = Kr)
str(ordInp_ratios)
sapply(ordInp_ratios, dim)

## Li & Ji threshold
ratios_thr <- thr.LiJi(ordInp_ratios$geno,
                       chrom = ordInp_ratios$map$chromosome,
                       ploidy = 4)
ratios_thr

## naive model
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


## plots
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


ratios_K$pheno01$pval[which.max(-log10(ratios_K$pheno01$pval))]
plot(ordInp_ratios$geno["mrk050949",],
     ordInp_ratios$pheno[,1],
     xlim = c(0,1),
     xlab = "ratios",
     ylab = "phenotype")



# ++ dosage Prob. (Bfreq) ----------------
## input order
ordInp_Bfreq <- inputOrder(Bfreq2, pheno = phenotypes, map = map, K = Kp)
str(ordInp_Bfreq)

## Li & Ji threshold
# mono <- apply(ordInp_Bfreq$geno,1, function(x) sd(x, na.rm = T)==0)
# ordInp_Bfreq$geno[mono, 1:9]
Bfreq_thr <- thr.LiJi(ordInp_Bfreq$geno,
                       chrom = ordInp_Bfreq$map$chromosome,
                       ploidy = 4)
Bfreq_thr

## naive model
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

## plots
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

Bfreq_K$pheno01$pval[which.max(-log10(Bfreq_K$pheno01$pval))]
plot(ordInp_Bfreq$geno["mrk050949",],
     ordInp_Bfreq$pheno[,1],
     xlim = c(0,1),
     xlab = "dosage Prob.",
     ylab = "phenotype")



# plot ratios -----------------
# source("file:///C:/Users/tumin001/OneDrive - WageningenUR/polyploids/software/gtfun/pkg_rasso/rasso/R/GWASplot.R")
library(scales)

dir.create("vignettes/workshop_20201203/Bfreq_plots")


ordMapRatio <- inputOrder(ratios[rownames(ratios) %in% map$marker[map$chromosome != 0],], map = map)
colnames(ordMapRatio$geno)
str(ordMapRatio)

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
     xlim = c(0,1))
hist(ordMapRatio$geno[,"ind011"],
     breaks = 500,
     xlim = c(0,1))


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

# skyplot(ordMapBfreq$geno[,"ind201"],
#         map = ordMapBfreq$map, #[physmap$physical_chr!=0,]
#         ylim = c(0,1.1),
#         pch = 20,
#         cex = 1,
#         ylab = "ratios",
#         col = alpha("dodgerblue4",0.2),
#         chromspace = 0,
#         main="ind201")



# Manhplot(scores=ratios[,1],         # a vector of gwas p-values
#          map=physmap[physmap$chromosome!=0,],       # custom map with 3 columns for marker, chrom, position
#          width=1200,
#          ChrCol = alpha(c("dodgerblue4","dodgerblue"),0.2),
#          yuplim=1.1,
#          ylab="ratios",
#          pch=20,
#          cex=0.7,
#          filepath=paste0("vignettes/workshop_20201203/Bfreq_plots/","ratios_01"))
#
#
#
#
# Manhplot(scores=Bfreq[,1],         # a vector of gwas p-values
#          map=physmap[physmap$chromosome!=0,],       # custom map with 3 columns for marker, chrom, position
#          width=1200,
#          ChrCol = alpha(c("dodgerblue4","dodgerblue"),0.2),
#          yuplim=1.1,
#          ylab="BfreqP",
#          pch=20,
#          cex=0.7,
#          filepath=paste0("vignettes/workshop_20201203/Bfreq_plots/","BfreqP_01"))
#
#
#
# for(i in colnames(Bfreq)) {
#   Manhplot(scores=Bfreq[,i],         # a vector of gwas p-values
#            map=physmap[physmap$chromosome!=0,],       # custom map with 3 columns for marker, chrom, position
#            width=1200,
#            ChrCol = alpha(c("dodgerblue4","dodgerblue"),0.2),
#            yuplim=1.1,
#            ylab="Bfreq",
#            pch=20,
#            cex=0.7,
#            filepath=paste0("vignettes/workshop_20201203/Bfreq_plots/","Bfreq_",i))
# }
#aneup: 14-225, 14-247, 14-325
#trip: 32-61















