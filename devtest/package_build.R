#Package building
library(devtools)

document(); build_manual()

data <- readRDS("vignette/workshop_data.RDS")
save(data,file="data/data.RData")

build()
ñ
