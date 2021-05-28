## code to prepare the example dataset


load("data-raw/rawdata.RData")

names(data)
sapply(data, object.size)
nrow(data$map)
tapply(data$map$marker, data$map$chromosome, length)

# sample markers -----------
data$result$phenotype1$pval
skyplot(-log10(data$result$phenotype1$pval), data$map, chrom = 4)
nrow(data$map)

qtls <- which(-log10(data$result$phenotype1$pval)>4)
mrkid <- sample(1:nrow(data$map), 1000)
mrkid <- union(mrkid, qtls)
mrkid <- sort(mrkid)
mrkid <- mrkid[data$map$chromosome[mrkid] %in% c(0,1,2,3,4)]

# mrkid <- which(data$map$chromosome %in% c(2,3,4))

length(mrkid)
mrknames <- data$map$marker[mrkid]

## map
map <- data$map[mrkid,]
object.size(map)
head(map)
nrow(map)

## snp dosages
data$snp[1:5,1:5]
snpdose <- data$snp[mrkid,]
rownames(snpdose) <- mrknames
snpdose[1:5,1:5]
object.size(snpdose)

class(snpdose)
snpdose <- as.matrix(snpdose)
typeof(snpdose)
object.size(snpdose)


## ancestral alleles
data$founder[1:5,1:9]

hapdose <- data$founder[mrkid,]
hapdose[1:5,1:5]
class(hapdose)
typeof(hapdose)
object.size(hapdose)



# rename parents --------------
head(data$pedigree)
ped <- as.matrix(data$pedigree[-1,])
pheno <- data$pheno[,1, drop=F]
head(pheno)

ped[1:15,]

renamevect <- ped[1:9,1]
names(renamevect) <- c("A1_P1","A1_P2","A1_P3","A2_P4","A2_P5","A2_P6","A3_P7","A3_P8","A3_P9")

## ped
for(i in 1:length(renamevect)){
  ped[ped == renamevect[i]] <- names(renamevect)[i]
}
ped[1:15,]

## pheno
pheno[1:5,]
rownames(pheno)
for(i in 1:length(renamevect)){
  rownames(pheno)[rownames(pheno) == renamevect[i]] <- names(renamevect)[i]
}

## snpdose
snpdose[1:5,1:5]
for(i in 1:length(renamevect)){
  colnames(snpdose)[colnames(snpdose) == renamevect[i]] <- names(renamevect)[i]
}

## hapdose
hapdose[1:5,1:9]
for(i in 1:length(renamevect)){
  # colnames(snpdose)[colnames(snpdose) == renamevect[i]] <- names(renamevect)[i]
  colnames(hapdose) <- gsub(renamevect[i],names(renamevect)[i],colnames(hapdose))
}
colnames(hapdose)[1:40]



# save --------------------

sum(object.size(pheno),
    object.size(ped),
    object.size(map),
    object.size(snpdose),
    object.size(hapdose))

# mpdata <- list(pheno=pheno,
#                ped=ped,
#                map=map,
#                snpdose=snpdose,
#                hapdose=hapdose)
# object.size(mpdata)
#
#
# save(mpdata, file="data/mpdata.RData")


mppheno <- pheno
mpped <- ped
mpmap <- map
mpsnpdose <- snpdose
mphapdose <- hapdose

save(mppheno, file="data/mppheno.RData")
save(mpped, file="data/mpped.RData")
save(mpmap, file="data/mpmap.RData", compress = "bzip2")
save(mpsnpdose, file="data/mpsnpdose.RData", compress = "xz")
save(mphapdose, file="data/mphapdose.RData", compress = "xz")
