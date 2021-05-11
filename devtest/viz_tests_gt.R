
devtools::load_all()



# example dataset --------------
load("data-raw/rawdata.RData")
names(data)

data$map[1:6,]
data$result$phenotype1$pval[1:6]
dim(data$map)
length(data$result$phenotype1$pval)
identical(data$map$marker, names(data$result$phenotype1$pval))




# skyplot -------------------------
skyplot(pval = -log10(data$result$phenotype1$pval),
        map = data$map)


# + map and pval order ---------

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



# + performance -----------------------------------
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




# + skyplot > draw_chrom_axis: wrong labelling of small axis ----------
# test
testmap <- data.frame(marker=paste0("m",0:10),
                      chromosome=1,
                      position=0:10,
                      pval=rep(1,11))

skyplot(testmap$pval,
        map = testmap,
        chromspace = 0,
        cex.small = 0.7,
        ylim=c(0,2),
        # threshold = -log10(ratios_thr$threshold),
        main="test")


plot(testmap$position,testmap$pval)

##
testmap2 <- data.frame(marker=paste0("g",0:10),
                       chromosome=2,
                       position=0:10,
                       pval=rep(1,11))

testmap2 <- rbind(testmap, testmap2)

skyplot(testmap2$pval,
        map = testmap2,
        chromspace = 0.1,
        cex.small = 0.7,
        ylim=c(0,2),
        # threshold = -log10(ratios_thr$threshold),
        main="test")


plot(testmap2$position,testmap2$pval)


##
testmap3 <- data.frame(marker=paste0("t",0:10),
                       chromosome=3,
                       position=seq(5,10, by=0.5),
                       pval=rep(1,11))



skyplot(testmap3$pval,
        map = testmap3,
        chromspace = 0,
        cex.small = 0.7,
        ylim=c(0,2),
        # threshold = -log10(ratios_thr$threshold),
        main="test")

##
testmap4 <- rbind(testmap2, testmap3)

skyplot(testmap4$pval,
        map = testmap4,
        chromspace = 0.1,
        cex.small = 0.7,
        ylim=c(0,2),
        # threshold = -log10(ratios_thr$threshold),
        main="test")

##
testmap5 <- data.frame(marker=paste0("m",0:50),
                       chromosome=1,
                       position=0:50,
                       pval=rep(1,51))


skyplot(testmap5$pval,
        map = testmap5,
        chromspace = 0,
        cex.small = 0.7,
        ylim=c(0,2),
        # threshold = -log10(ratios_thr$threshold),
        main="test")

##
testmap6 <- data.frame(marker=paste0("m",0:104),
                       chromosome=c(rep(1,53), rep(2,52)),
                       position=c(0:51,52.4,1:52),
                       pval=c(rep(1,53), rep(0.1,52)))


skyplot(testmap6$pval[1:53],
        map = testmap6[1:53,],
        chromspace = 0.1,
        cex.small = 0.7,
        ylim=c(0,2),
        # threshold = -log10(ratios_thr$threshold),
        main="test")

skyplot(testmap6$pval[54:105],
        map = testmap6[54:105,],
        chromspace = 0.1,
        cex.small = 0.7,
        ylim=c(0,2),
        # threshold = -log10(ratios_thr$threshold),
        main="test")

skyplot(testmap6$pval,
        map = testmap6,
        chromspace = 0.1,
        cex.small = 0.7,
        ylim=c(0,2),
        # threshold = -log10(ratios_thr$threshold),
        main="test")


plot(testmap6$position[1:53],testmap6$pval[1:53])

# map = testmap6
# ch_edges = axis_res$chr_edges
# i=2



# pheno_box -------------------
data("mppheno")
data("mphapdose")

pheno_box(phe = mppheno, gen = mphapdose[1,], haplotype = TRUE, ploidy = 4)
pheno_box(phe = mppheno, gen = mphapdose[1,], haplotype = TRUE, h=90)
pheno_box(phe = mppheno, gen = mphapdose[1,], haplotype = TRUE, ploidy=4, h=200)
pheno_box(phe = mppheno, gen = mphapdose[1,], haplotype = TRUE, ploidy=4, hap.select=c("28","9"), h=10)




