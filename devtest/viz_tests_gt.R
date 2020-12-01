
devtools::load_all()


# example dataset --------------
load("data/data.RData")
names(data)

data$map[1:6,]
data$result$phenotype1$pval[1:6]
dim(data$map)
length(data$result$phenotype1$pval)
identical(data$map$marker, names(data$result$phenotype1$pval))




skyplot(pval = -log10(data$result$phenotype1$pval),
        map = data$map)


# map and pval order ---------

# change order of markers in map and pval
set.seed(3)
neworder <- sample(1:nrow(data$map))
map <- data$map[neworder,]
map[1:6,]
pval <- data$result$phenotype1$pval[neworder]

skyplot(pval = -log10(pval),
        map = map)

skyplot(pval = -log10(data$result$phenotype1$pval),
        map = data$map,
        pch=21)



# performance -----------------------------------
library(microbenchmark)

pval <- matrix(runif(100000), ncol = 1)
rownames(pval) <- paste0("snp",1:nrow(pval))
map <- data.frame(marker = rownames(pval),
                  chromosome = rep(1:10, length.out=nrow(pval)),
                  position = seq(1,100,length.out=nrow(pval)))
map <- map[order(map$chromosome, map$position),]
pval <- pval[match(map$marker, rownames(pval)),]


skyplot(pval = -log10(pval),
        map = map)


Manhplot(scores=-log10(pval),         # a vector of gwas p-values
         map = map,                 # map with 3 columns for marker, chrom, position
         Chr=NULL,            # chromosome name
         ChrCol = c("dodgerblue4","dodgerblue"))




microbenchmark(skyplot(pval = -log10(pval),
                       map = map),
               Manhplot(scores=-log10(pval),         # a vector of gwas p-values
                        map = map,                 # map with 3 columns for marker, chrom, position
                        Chr=NULL,            # chromosome name
                        ChrCol = c("dodgerblue4","dodgerblue")), times=10)

microbenchmark(skyplot(pval = -log10(data$result$phenotype1$pval),
                       map = data$map),
               Manhplot(scores=-log10(data$result$phenotype1$pval),         # a vector of gwas p-values
                        map = data$map,                 # map with 3 columns for marker, chrom, position
                        Chr=NULL,            # chromosome name
                        ChrCol = c("dodgerblue4","dodgerblue")))
