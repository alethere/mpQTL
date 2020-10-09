
# example dataset --------------
names(data)

data$map[1:6,]
data$result$phenotype1$pval[1:6]
dim(data$map)
length(data$result$phenotype1$pval)
identical(data$map$marker, names(data$result$phenotype1$pval))

source("R/viz_fun.R")


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

