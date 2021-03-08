library(usethis)
library(devtools)
library(sinew)

sinew::makeOxygen(hs_knn_pred)
sinew::makeOxygen(hs3D_knn)

usethis::use_description()

usethis::use_build_ignore(here::here("Bor_data"))
usethis::use_build_ignore(here::here("R", "deprecated"))
