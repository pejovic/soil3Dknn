
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

load(here::here("Bor_data", "BorData.rda"))
load(here::here("Bor_data", "covmaps.rda"))

source(here::here("R", "deprecated", "preprocessing.r"))




#================== Spatial references ===================================================================
gk_7 <- "+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m"
utm <- "+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m"

#================== Names and abbrevations of covariates =================================================

CovNames <- c("DEM","Aspect","Slope","TWI","ConvInd","CrSectCurv","LongCurv","ChNetBLevel","VDistChNet","NegOp", "PosOp","WEeast","WEnw","clc","SoilType")



#=================== DATA ================================================================================
bor <- bor %>% dplyr::select(ID, Top, Bottom, As, pH, Co, SOM, Soil.Type, x, y ) %>%
  dplyr::mutate(OC = SOM*10/1.724,
                depth = (Top + Bottom)/2)

bor.profs <- bor
depths(bor.profs) <- ID ~ Top + Bottom
site(bor.profs) <- ~ Soil.Type + x + y
coordinates(bor.profs) <- ~ x+y
proj4string(bor.profs) <- CRS(utm)
#=====================================================================================

#proj4string(gridmaps.sm2D) <- gk_7
#==========================================================



#==================================================================================================================

bor_sf <- bor %>% sf::st_as_sf(., coords = c("x", "y"), remove = FALSE) %>% 
  dplyr::group_by(ID) %>% 
  dplyr::mutate(hz_ID = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(!is.na(OC)) %>% 
  dplyr::arrange(ID) %>% 
  dplyr::mutate(ID = as.character(ID), x = as.double(x), y = as.double(y)) %>%
  dplyr::select(ID, hz_ID, Top, Bottom, As:SOM, OC, Soil.Type, x, y, depth, geometry) %>%
  dplyr::filter(!is.na(OC)) %>%
  sf::st_set_crs(utm)

bor_sf <- bor_sf %>% sf::st_transform(gk_7) %>% dplyr::mutate(x = sf::st_coordinates(.)[,1], y = sf::st_coordinates(.)[,2]) %>% sf::st_drop_geometry() %>% dplyr::select(ID, Top, Bottom, As, pH, Co, SOM, OC, x, y, depth)

bor_sp <- bor_sf %>% dplyr::select(-Top, -Bottom) %>% dplyr::distinct(x, y, .keep_all = TRUE) 
coordinates(bor_sp) <- ~x + y 
proj4string(bor_sp) <- gk_7
bor_ov <- sp::over(bor_sp, gridmaps.sm2D) %>% dplyr::mutate(ID = bor_sp$ID) %>% dplyr::select(ID, everything())

bor_df <- dplyr::inner_join(bor_sf, bor_ov) %>% dplyr::select(ID, Top, Bottom, OC, AHils:DEM, LongCurv:VDistChNet, clc, SoilType, x, y, depth)

#============================== Creating folds ==============================

bor.folds <- stratpart(target.name = c("OC", "x", "y"), sp.data = bor_sp, num.folds = 5, princ.comp = TRUE, kmean.vars = NA, num.means = 3, spatial.cluster = TRUE ,nested = TRUE, seed = 45, cum.prop = 0.9)

bor_list <- list(data = bor_df, folds = bor.folds)

save(bor_list, file = here::here("Bor_data", "bor_data.rda"))

