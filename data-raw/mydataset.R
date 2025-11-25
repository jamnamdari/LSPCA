## code to prepare `mydataset` dataset goes here
HC <- load("data-raw/HC.rda")
FEP <- load("data-raw/FEP.rda")
D <- load("data-raw/simulated_data.rda")
usethis::use_data(HC, FEP, D, overwrite = TRUE)
