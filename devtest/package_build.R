#Package building
library(devtools)

document(); build_manual()

#Include data
data <- readRDS("vignettes/new_workshop_data.RDS")
save(data,file="data/data.RData")


setwd("D:/Science/PhD/Code/mpQTL/mpQTL/")
build(vignettes = F)
