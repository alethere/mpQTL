devtools::load_all()

# Add citation -------------------
usethis::use_citation()






# Package building ------------------------------------
# library(devtools)

## function documentation --------------
devtools::document()
devtools::check_man() #calls document() and runs the checks for manual included in check_built()
# devtools::build_manual() #to add a pdf manual

# ## data ----------------
# data to be included must be in the folder 'data' with extension .RData
# data <- readRDS("vignettes/new_workshop_data.RDS")
# save(data,file="data/data.RData")

# ## vignettes ------------
# to create a template  use:
# usethis::use_vignette("mpQTL_newvignette")

## check before building ------------
devtools::check()
devtools::check_win_devel()

## build -----------------
# devtools::build()
devtools::build(vignettes = F) #skip building vignettes

## check after building -------------
devtools::check_built(paste0("../mpQTL_0.6.2.tar.gz"), cran=TRUE)


# setwd("D:/Science/PhD/Code/mpQTL/mpQTL/")
# build(vignettes = F)

# #To install it
# install.packages("../mpQTL_0.6.2.tar.gz", repos = NULL, type = "source")
