#Package building
# library(devtools)

# function documentation
devtools::document()
devtools::build_manual() #to add a pdf manual

# data
# data to be included must be in the folder 'data' with extension .RData
# data <- readRDS("vignettes/new_workshop_data.RDS")
# save(data,file="data/data.RData")

# vignettes
# to create a template  use:
# usethis::use_vignette("mpQTL_newvignette")

devtools::build()


devtools::check_built(paste0("../mpQTL_0.3.0.tar.gz"), cran=TRUE)



