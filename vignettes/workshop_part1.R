


# Installing and loading the package
install.packages("../mpQTL_0.2.0.tar.gz",repos=NULL)
library("mpQTL")

# the object "data" contains an example dataset
data <- readRDS("new_workshop_data.RDS")
names(data)


# Population design ----------------------------
## pedigree file
ped <- read.table("cross_workshop.ped", header = T, stringsAsFactors = F)
head(ped,15)


# Haplotyping ----------------------------------
## parents used in each cross
parents <- as.matrix(unique(ped[!is.na(ped[,2]) & !is.na(ped[,3]),2:3]))
parents

## merge reciprocal crosses
for(i in 1:nrow(parents)) {
  parents[i,] <- parents[i, order(parents[i,])]
}
parents <- as.matrix(unique(parents[,1:2]))
parents

F1 <- list()
for (p in seq_len(nrow(parents))) {
  F1[[p]] <- ped[,1][(!is.na(ped[,2]) & !is.na(ped[,3]) & ped[,2]==parents[p, 1] & ped[,3]==parents[p, 2]) |
                       (!is.na(ped[,2]) & !is.na(ped[,3]) & ped[,2]==parents[p, 2] & ped[,3]==parents[p, 1])]
}
str(F1)


# Make haploblocks
## genetic map
map <- data$map
head(map)
table(map$chromosome) # number of markers per chromosome
tapply(map$position, map$chromosome, range) # chromosome lenght

## create blocks by a sliding window
hb_list <- map2blocks(map=map[map$chromosome != 0,], #excluding chrom 0
                      winsize=1, # window size
                      sldpace=1) # window step (sliding pace)

length(hb_list)
table(sapply(hb_list, length))
names(hb_list)[1:5]

## refine blocks with more than 8 SNPs
set.seed(3)
hb_list <- refineBlocks(hb_list,   # haploblock list
                        nmrk=2:7,  # numeric vector indicating allowed block length
                        method="random",
                        nrand = 2) # number of random selections per block
names(hb_list)[1:5]

# Haplotype inference (DO NOT RUN, we get true haplotypes from our simulation)
# ## infer haplotypes by PolyHaplotyper (DO NOT RUN)
# results <- inferHaplotypes(mrkDosage=snp, ploidy=4,
#                            haploblock=hb_list,
#                            parents=parents, F1=F1,
#                            maxmrk=7)

## get true haplotypes from the simulation
PolyHap_res <- data$PolyHap_res
str(PolyHap_res[[1]])


# Haplotype data cuaration and format conversion -------------------------------------
PolyHap_res[[1]]$hapdos[,1:3] #first block, 21 haplotypes, first 3 individuals
hapres <- HapCurate(PolyHap_res,  # results of PolyHaplotyper
                    hb_list = hb_list, #may contain also 1-SNP blocks
                    map = map,
                    ploidy = 4,
                    na.rate = 1,  # allowed missingness per marker
                    use.SNPs = F, # whether you want to include SNPs of discarded blocks
                    snpdose = snp)

hap <- hapres$genotypes  # haplotypes
hapmap <- hapres$map     # haplotype map


# Types of genetic data accepted by `mpQTL` --------------------------
## SNP dosages
snp <- as.matrix(data$snp)
snp[1:5,1:4]

## Haplotype names
hap[1:5,1:8]

## from (a list of) haplotype dosages to haplotype names
eg1 <- HapdoseToHapname(PolyHap_res[1:3], ploidy = 4)
eg1[,1:8]

## from SNP dosages to haplotype names
eg2 <- SNPdoseToHapname(snp[1:5,1:4], ploidy = 4)
eg2

## from haplotype names to haplotype dosages
eg3 <- HapnameToHapdose(eg1)
str(eg3)
eg3[[1]][,1:3]




