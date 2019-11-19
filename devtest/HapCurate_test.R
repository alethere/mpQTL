
source("D:/OneDrive - WageningenUR/polyploids/software/haploutils/NAMhaplo_utils.R")



outf <- "research/inferred_haplotypes/output_inferredHaps"

f <- "D:/OneDrive - WageningenUR/polyploids/workpackage4/data/PedigreeSim/NAMcrosses"
map <- read.table(paste0(f,"/Potato.map"),
                  header = T)


crs <- "cross200"
h <- "0.5"
ws <- "1"
nmrk <- "6"


# load inferred haplotypes
load(paste0(outf,"/QTL1loop/",crs,".results.",ws,"_",nmrk,"_QTL1peak.RData"))
results[[2]]$hapdos[,1:9]

# extract hb_list
hb_list <- lapply(results, function(x) {
  x$markers
})
names(hb_list)[1:3]





# test HapdoseToHapname ------------------
inherits(results, "list")
inherits(results[[1]], "list")

haplo.temp <- lapply(results, function(x) x[[1]])
haplo.test <- HapdoseToHapname(results, 4)
haplo.test2 <- HapdoseToHapname(haplo.temp, 4)
identical(haplo.test,haplo.test2)
haplo.test[1:5,1:9]
haplo.test2[1:5,1:9]


# test HapCurate -------------------

##
hap1 <- HapCurate(haplo = results, #results of PolyHaplotyper
                  # hb_list, #containing also 1 SNP blocks
                  # map,  #snp map
                  ploidy = 4,
                  na.rate = 1, #from 0 to 1. 1 means no filtration
                  use.SNPs = F
                  # snpdose
)

hap1$genotypes[1:5,1:9]
hap1$map

identical(haplo.test, hap1$genotypes)

##
hap2 <- HapCurate(haplo = haplo.test, #results of PolyHaplotyper
                  # hb_list, #containing also 1 SNP blocks
                  # map,  #snp map
                  ploidy = 4,
                  na.rate = 1, #from 0 to 1. 1 means no filtration
                  use.SNPs = F
                  # snpdose
)

identical(haplo.test, hap2$genotypes)
hap2$genotypes[1:15,1:9]

##
hap3 <- HapCurate(haplo = results, #results of PolyHaplotyper
                  hb_list = hb_list, #containing also 1 SNP blocks
                  map = map,  #snp map
                  ploidy = 4,
                  na.rate = 1, #from 0 to 1. 1 means no filtration
                  use.SNPs = F
                  # snpdose
)


identical(haplo.test[order(rownames(haplo.test)),],
          hap3$genotypes[order(rownames(hap3$genotypes)),])


##
hap4 <- HapCurate(haplo = results, #results of PolyHaplotyper
                  hb_list = hb_list, #containing also 1 SNP blocks
                  map = map,  #snp map
                  ploidy = 4,
                  na.rate = 0.4, #from 0 to 1. 1 means no filtration
                  use.SNPs = F
                  # snpdose
)

nrow(haplo.test) == nrow(hap4$genotypes)

identical(haplo.test[match(rownames(hap4$genotypes), rownames(haplo.test)),],
          hap4$genotypes)

hap4$genotypes[1:15,1:9]
hap4$map[1:20,]


##
hap5 <- HapCurate(haplo = results, #results of PolyHaplotyper
                  hb_list = hb_list, #containing also 1 SNP blocks
                  map = map,  #snp map
                  ploidy = 4,
                  na.rate = 0.4, #from 0 to 1. 1 means no filtration
                  use.SNPs = T,
                  snpdose = snpdose
)

nrow(haplo.test) == nrow(hap5$genotypes)

identical(haplo.test[intersect(rownames(hap5$genotypes), rownames(haplo.test)),],
          hap5$genotypes[intersect(rownames(hap5$genotypes), rownames(haplo.test)),])
all.equal(haplo.test[intersect(rownames(hap5$genotypes), rownames(haplo.test)),],
          hap5$genotypes[intersect(rownames(hap5$genotypes), rownames(haplo.test)),])

hap5$genotypes[1:15,1:9]
hap5$map[1:20,]
