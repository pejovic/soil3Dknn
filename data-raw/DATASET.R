## code to prepare `DATASET` dataset goes here

library(Surveyor)

load(here::here("Bor_data", "edgeroi.rda"))
load(here::here("Bor_data", "edgeroi.covmaps.rda"))

edgeroi_data <- edgeroi
edgeroi_maps <- cov.maps


usethis::use_data(edgeroi_data, overwrite = TRUE)
usethis::use_data(edgeroi_maps, overwrite = TRUE)


load("C:/R_projects/soil3Dknn/edgeroi_results.rda")
load("C:/R_projects/soil3Dknn/bor_results.rda")


usethis::use_data(edgeroi_results, overwrite = TRUE)
usethis::use_data(bor_results, overwrite = TRUE)
