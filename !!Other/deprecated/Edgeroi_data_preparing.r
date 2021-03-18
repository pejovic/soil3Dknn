
library(rgdal)
library(sp)
library(aqp)
library(tidyverse)
library(tidymodels)
library(here)
library(raster)
library(gower)
library(caret)
library(sf)

library(soil3Dknn)
source(here::here("R", "deprecated", "preprocessing.r"))


# Edgeroi data
load(here::here("Bor_data", "edgeroi.rda"))
load(here::here("Bor_data", "edgeroi.covmaps.rda"))


edgeroi$horizons %<>% dplyr::rename("Top" = "UHDICM", "Bottom" = "LHDICM", "ID" = "SOURCEID") %>%
  dplyr::mutate(ORCperc = (ORCDRC)/10)

edgeroi$sites <-  edgeroi$sites %>% dplyr::rename("ID" = "SOURCEID")
sites <- edgeroi$sites
coordinates(sites) <- ~ LONGDA94 + LATGDA94
proj4string(sites) <- CRS("+proj=longlat +ellps=GRS80 +datum=WGS84 +no_defs")
sites <- spTransform(sites, proj4string(cov.maps))
edgeroi$sites <- data.frame(sites)
edgeroi$sites <- edgeroi$sites %>% dplyr::rename(c("x"="LONGDA94", "y"="LATGDA94"))

edgeroi_sf <- dplyr::inner_join(edgeroi$horizons, edgeroi$sites) %>% dplyr::select(ID, ORCDRC, ORCperc, PHIHO5, Top, Bottom, x, y)

edgeroi_sf <- edgeroi_sf %>% sf::st_as_sf(., coords = c("x", "y"), remove = FALSE) %>% 
  dplyr::group_by(ID) %>% 
  dplyr::mutate(hz_ID = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(!is.na(ORCperc)) %>% 
  dplyr::arrange(ID) %>% 
  dplyr::mutate(ID = as.character(ID), x = as.double(x), y = as.double(y), depth = (Top + Bottom)/2) %>%
  dplyr::select(ID, hz_ID, Top, Bottom, ORCDRC, ORCperc, PHIHO5, x, y, depth, geometry) %>%
  dplyr::filter(!is.na(ORCperc)) %>%
  sf::st_set_crs(proj4string(cov.maps))

edgeroi_sf <- edgeroi_sf %>% dplyr::mutate(x = sf::st_coordinates(.)[,1], y = sf::st_coordinates(.)[,2]) %>% sf::st_drop_geometry() %>% dplyr::select(ID, Top, Bottom, ORCDRC, ORCperc, PHIHO5, x, y, depth)

edgeroi_sp <- edgeroi_sf %>% dplyr::select(-Top, -Bottom) %>% dplyr::distinct(x, y, .keep_all = TRUE) 
coordinates(edgeroi_sp) <- ~x + y 
proj4string(edgeroi_sp) <- proj4string(cov.maps)
edgeroi_ov <- sp::over(edgeroi_sp, cov.maps) %>% dplyr::mutate(ID = edgeroi_sp$ID) %>% dplyr::select(ID, everything())

edgeroi_df <- dplyr::inner_join(edgeroi_sf, edgeroi_ov) %>% dplyr::select(ID, Top, Bottom, ORCperc, MVBSRT6:EV3MOD5, x, y, depth)

#============================== Creating folds ==============================

folds <- stratpart(target.name = c("ORCperc", "x", "y"), sp.data = edgeroi_sp, num.folds = 5, princ.comp = TRUE, kmean.vars = NA, num.means = 3, spatial.cluster = TRUE ,nested = TRUE, seed = 45, cum.prop = 0.9)


edgeroi_list <- list(data = edgeroi_df, folds = folds)

save(edgeroi_list, file = here::here("Bor_data", "edgeroi_data.rda"))



