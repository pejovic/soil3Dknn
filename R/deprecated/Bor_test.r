
library(rgdal)
library(aqp)
library(psych)
library(mda)
library(classInt)
library(tidyverse)
library(tidymodels)
library(MASS)
library(magrittr)
library(here)
library(raster)
library(gower)
library(caret)
library(tidymodels)
library(sf)

load(here::here("data", "Bor", "BorData.rda"))
load(here::here("data", "Bor", "covmaps.rda"))

source(here::here("R","preprocessing.r"))


devtools::install_github("pejovic/soil3Dknn")



#================== Spatial references ===================================================================
gk_7 <- "+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +datum=WGS84 +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m"
utm <- "+proj=utm +zone=34 +ellps=GRS80 +datum=WGS84 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m"

#================== Names and abbrevations of covariates =================================================

CovNames <- c("DEM","Aspect","Slope","TWI","ConvInd","CrSectCurv","LongCurv","ChNetBLevel","VDistChNet","NegOp", "PosOp","WEeast","WEnw","clc","SoilType")



#=================== DATA ================================================================================
bor.profs <- bor[,c("ID","x","y","Soil.Type","Top","Bottom","As","pH","Co","SOM")]
bor.profs$logSOM <- log(bor.profs$SOM)
bor.profs$logAs <- log(bor.profs$As)
bor.profs$ORCDRC <- bor.profs$SOM*10/1.724
bor.profs$logORCDRC <- log1p(bor.profs$ORCDRC)

bor.profs$mid <- (bor.profs$Top+bor.profs$Bottom)/2
depths(bor.profs) <- ID ~ Top + Bottom
site(bor.profs) <- ~ Soil.Type + x + y
coordinates(bor.profs) <- ~ x+y
proj4string(bor.profs) <- CRS(utm)
#=====================================================================================
str(gridmaps.sm2D)
proj4string(gridmaps.sm2D) <- gk_7
#==========================================================



#==================================================================================================================

spc <- bor.profs

spc_sf <- sf::st_as_sfc(spc@sp) %>% sf::st_sf(spc@site, geometry = .) %>% dplyr::inner_join(spc@horizons, .) %>% dplyr::group_by(ID) %>% dplyr::mutate(hz_ID = 1:n()) %>% dplyr::ungroup() %>% dplyr::filter(!is.na(ORCDRC))
  
spc_sf <- sf::st_as_sf(spc_sf) %>% sf::st_transform(gk_7) %>% dplyr::mutate(x = sf::st_coordinates(.)[,1], y = sf::st_coordinates(.)[,2]) %>% sf::st_drop_geometry() %>% dplyr::select(ID, Top, Bottom, As, pH, Co, SOM, ORCDRC, x, y)

spc_sp <- spc_sf %>% dplyr::select(-Top, -Bottom) %>% dplyr::distinct(x, y, .keep_all = TRUE) 
coordinates(spc_sp) <- ~x + y 
proj4string(spc_sp) <- gk_7
spc_ov <- sp::over(spc_sp, gridmaps.sm2D) %>% dplyr::mutate(ID = spc_sp$ID) %>% dplyr::select(ID, everything())

spc_df <- dplyr::inner_join(spc_sf, spc_ov) %>% dplyr::select(ID, Top, Bottom, ORCDRC, AHils:DEM, LongCurv:VDistChNet, clc, SoilType, x, y)


# Formulas
ORCDRC.fun <- as.formula(paste("ORCDRC ~", paste(c(names(spc_df)[-c(1:4)]), collapse="+")))

#============================== Creating folds ==============================

folds <- stratpart(target.name = c("ORCDRC", "x", "y"), sp.data = spc_sp, num.folds = 5, princ.comp = TRUE, kmean.vars = NA, num.means = 3, spatial.cluster = TRUE ,nested = TRUE, seed = 45, cum.prop = 0.9)


#============================== Homosoil knn ================================
spc_df <- spc_df %>% dplyr::mutate(depth = (Top+Bottom)/2)
"%ni%" <- Negate("%in%")
train_df_ID <- spc_df %>% dplyr::distinct(x, y, .keep_all = TRUE) %>% .[!is.na(.$ORCDRC), ] %>% .[folds$outer.fold != 3, "ID"]
train_df <- spc_df %>% dplyr::filter(ID %in% train_df_ID$ID)
pred_df_all <- spc_df %>% dplyr::filter(ID %ni% train_df_ID$ID)
#%>% dplyr::distinct(x, y, .keep_all = TRUE)

#spc_sf <- gowers_list2[[1]]
#depth <- as.list(pred_df$depth)[[1]]
#gdist <- unique(gowers_list2[[1]]$gower_dist)
#soil.fun <- ORCDRC.fun

spc_sf = gdist_list5[[1]]; pred_depth = pred_list[[1]]$depth; d_inc = d_inc; n = n

# Uvesti distance kao tezinu

hs_knn_pred <- function(soil.fun, spc_sf, pred_depth, d_inc, n, p){
  target.name <- all.vars(soil.fun)[1]
  gdist = unique(spc_sf$gw)[1:n]
  n_ids <- unique(spc_sf$ID)[1:n] %>% data.frame(ID = ., stringsAsFactors = FALSE)
  spc_sf %>% dplyr::mutate(pred_top = pmax(pred_depth-d_inc, 0), pred_bot = pred_depth+d_inc) %>% 
    dplyr::filter(Top < pred_bot, Bottom > pred_top) %>% 
    mutate(dh = case_when(pred_top >= Top & pred_bot <= Bottom ~ pred_bot - pred_top,
                          pred_top >= Top & pred_bot >= Bottom ~ Bottom - pred_top,
                          pred_top <= Top & pred_bot >= Bottom ~ Bottom - Top,
                          pred_top <= Top & pred_bot <= Bottom ~ pred_bot - Top)) %>%
    dplyr::filter(., ID %in% n_ids$ID) %>% dplyr::left_join(n_ids, ., by = "ID") %>%
    dplyr::group_by(ID) %>% dplyr::summarise(obs_mean = weighted.mean(eval(parse(text = target.name)), w = dh, na.rm = TRUE)) %>%
    dplyr::ungroup() %>% dplyr::mutate(gower_distance = gdist) %>% 
    dplyr::summarise(pred_mean = sum(obs_mean/gower_distance^p)/sum(1/gower_distance^p))
    #dplyr::summarise(pred_mean = weighted.mean(obs_mean, w = gower_distance, na.rm = TRUE)) 
}



hs_knn_prep <- function(soil.fun, spc_sf, pred_depth, d_inc, n){
  target.name <- all.vars(soil.fun)[1]
  gdist = unique(spc_sf$gw)[1:n]
  n_ids <- unique(spc_sf$ID)[1:n] %>% data.frame(ID = ., stringsAsFactors = FALSE)
  spc_sf %>% dplyr::mutate(pred_top = pmax(pred_depth-d_inc, 0), pred_bot = pred_depth+d_inc) %>% 
    dplyr::filter(Top < pred_bot, Bottom > pred_top) %>% 
    mutate(dh = case_when(pred_top >= Top & pred_bot <= Bottom ~ pred_bot - pred_top,
                          pred_top >= Top & pred_bot >= Bottom ~ Bottom - pred_top,
                          pred_top <= Top & pred_bot >= Bottom ~ Bottom - Top,
                          pred_top <= Top & pred_bot <= Bottom ~ pred_bot - Top)) %>%
    dplyr::filter(., ID %in% n_ids$ID) %>% dplyr::left_join(n_ids, ., by = "ID") %>%
    dplyr::group_by(ID) %>% dplyr::summarise(obs_mean = weighted.mean(eval(parse(text = target.name)), w = dh, na.rm = TRUE)) %>%
    dplyr::ungroup() %>% dplyr::mutate(gower_dist = gdist)
    #dplyr::summarise(pred_mean = weighted.mean(obs_mean, w = gdist, na.rm = TRUE)) 
}


soil.fun <- as.formula(paste("ORCDRC ~", paste(c("x", "y", "DEM", "depth"), collapse="+")))
pred_df = pred_df_all 
train_df = train_df
d_inc = 5
n = 5
output <- "prediction"

pred.data <- pred_df
obs.data <- train_df

max(st_distance(obs_sf))

# soil.fun = soil.fun; pred.data = pred_df; obs.data = train_df; radius = 5000; d_inc = 15; n = 10; p = 2; min_obs = 12; output = "prediction"



hs3D <- function(soil.fun, pred.data, obs.data, d_inc, n, p = 1, min_obs = 4, radius, output = as.list("prediction", "preparation")){
   output <- output[1]
   obs_sf <- sf::st_as_sf(obs.data, coords = c("x", "y"), remove = FALSE)
   pred_sf <- sf::st_as_sf(pred.data, coords = c("x", "y"), remove = FALSE)
   
   if(min_obs < n) stop("Each prediction location must have more available obs.data than n")
   
   target.name <- all.vars(soil.fun)[1]
   coord_names <- c("x", "y", "depth")
   cov.names <- all.vars(soil.fun)[-1]
   cov.names <- cov.names[-which(cov.names == "depth")]
   
   pred_depth <- as.list(pred.data$depth) %>% purrr::map(.x = ., .f = ~data.frame(depth = .x))
   
   obs_inds <- nngeo::st_nn(x = pred_sf, y = obs_sf %>% dplyr::distinct(x, y, .keep_all = TRUE), sparse = TRUE, k = length(unique(obs.data$ID)), maxdist = radius)
   
   obs_list <- purrr::map(.x = obs_inds, .f = ~obs.data %>% dplyr::distinct(x, y, .keep_all = TRUE) %>% .[.x, ] %>% .$ID) %>% 
      purrr::map(.x = ., .f = ~dplyr::filter(obs.data[, c("ID", "Top", "Bottom", all.vars(soil.fun))], ID %in% .x))
   
   gdist_list <- purrr::map2(.x = pred_depth, .y = obs_list, .f = ~gower_dist(x = .x, y = .y[, "depth"])) %>%
      purrr::map2(.x = ., .y = obs_list, .f = ~.x < d_inc/diff(range(.y$depth)))
   
   #gdist_list <- purrr::map(.x = pred_depth, .f = ~gower_dist(x = .x, y = data.frame(depth = obs.data$depth))) %>%
   #purrr::map(.x = ., .f = ~.x < d_inc/diff(range(obs.data$depth)))
   
   #gdist_list2 <- purrr::map(.x = gdist_list, .f = ~obs.data[.x, ]) # Ovde se može dodati selekcija u okviru nekog buffer-a
   
   obs_list2 <- purrr::map2(.x = gdist_list, .y = obs_list, .f = ~.y[.x, ]) # Ovde se može dodati selekcija u okviru nekog buffer-a
   
   ex_obs <- unlist(map(.x = obs_list2, .f = ~dim(.x)[1])) < min_obs # TRUE for obs.data with no obs.data
   
   if(any(ex_obs)) warning(paste("There is no enough available obs.data for prediction at location of : ", paste(which(ex_obs), collapse = " ")))
   
   if(sum(!ex_obs) == 0) stop("There is no available obs.data on target depth")
   
   pred.data <- pred.data[!ex_obs, ]
   obs_list2 <- obs_list2[!ex_obs]
   
   pred_list <- pred.data %>% tibble::rowid_to_column() %>% dplyr::group_split(rowid, .keep = FALSE)
   
   #purrr::map2(.x = pred_list, .y = gdist_list2, .f = ~nn2(data = .y[, coord_names[1:2]], query = .x[, coord_names[1:2]], searchtype = "radius", radius = 1000))
   
   #nn2(data = gdist_list2[[1]][, c("x", "y")], query = pred_list[[1]][, c("x", "y")], radius = 1000)
   
   
   gdist_list1 <- purrr::map2(.x = pred_list, .y = obs_list2, .f = ~gower_topn(x = .x[, cov.names], y = .y[, cov.names], n = dim(.y)[1]))
   
   obs_list3 <- purrr::map2(.x = obs_list2, .y = gdist_list1, .f = ~.x[.y$index, ]) %>%
      purrr::map2(.x = ., .y = gdist_list1, .f = ~dplyr::mutate(.x, gw = .y$distance[,1]))
   
   #obs_list4 <- purrr::map(.x = obs_list3, .f = ~obs.data[obs.data$ID %in% .x$ID, ] %>% dplyr::left_join(.x, ., by = "ID"))
   
   if(output == "prediction"){
      prediction <- purrr::map2(.x = obs_list3, .y = pred_list, .f = ~hs_knn_pred(soil.fun = soil.fun, spc_sf = .x, pred_depth = .y$depth, d_inc = d_inc, n = n, p = p)) %>%
         bind_rows()
      output <- rep(NA, length(pred_depth))
      output[!ex_obs] <- prediction$pred_mean
   }else{
      preparation <- purrr::map2(.x = obs_list3, .y = pred_list, .f = ~hs_knn_prep(soil.fun = soil.fun, spc_sf = .x, pred_depth = .y$depth, d_inc = d_inc, n = n))
      output <- as.list(rep(NA, length(pred_depth)))
      output[!ex_obs] <- preparation
   }
   return(output)
}

ORCDRC.fun <- as.formula(paste("ORCDRC ~", paste(c(names(spc_df)[-c(1:4)]), collapse="+")))
ORCDRC.fun <- as.formula(paste("ORCDRC ~", paste(c("x", "y", "DEM", "depth"), collapse="+")))

system.time(aa <- hs3D(soil.fun = ORCDRC.fun, pred.data = pred_df_all, obs.data = train_df, radius = 4000, d_inc = 5, n =3, p = 2, min_obs = 6, output = "prediction"))

aa

ggplot(data = spc_df) + geom_point(aes(x = x, y = y, size = ORCDRC)) + 
  geom_point(data = spc_df %>% dplyr::filter(ID %in% aa[[35]]$ID), aes(x = x, y = y, size = ORCDRC), color = "red") + 
  geom_point(data = spc_df %>% dplyr::filter(ID %ni% train_df_ID$ID) %>%.[35, ], aes(x = x, y = y, size = ORCDRC), color = "blue")



spc_df <- spc_df %>% dplyr::mutate(depth = (Top+Bottom)/2)
"%ni%" <- Negate("%in%")
train_df_ID <- spc_df %>% dplyr::distinct(x, y, .keep_all = TRUE) %>% .[!is.na(.$ORCDRC), ] %>% .[folds$outer.fold != 3, "ID"]
train_df <- spc_df %>% dplyr::filter(ID %in% train_df_ID$ID)
pred_df_all <- spc_df %>% dplyr::filter(ID %ni% train_df_ID$ID)


#[1] "ID"               "Top"              "Bottom"           "ORCDRC"           "AHils"            "Aspect"           "CatchArea"       
#[8] "ChNetBLevel"      "ConvInd"          "CrSectCurv"       "DEM"              "LongCurv"         "LSFactor"         "NegOp"           
#[15] "PosOp"            "RelSlopePosition" "Slope"            "TWI"              "VelleyDepth"      "VDistChNet"       "clc"             
#[22] "SoilType"         "x"                "y"                "depth"     

ORCDRC.fun <- as.formula(paste("ORCDRC ~", paste(c("x", "y", "DEM", "depth"), collapse="+")))

folds <- stratpart(target.name = c("ORCDRC", "x", "y"), sp.data = spc_sp, num.folds = 5, princ.comp = TRUE, kmean.vars = NA, num.means = 3, spatial.cluster = TRUE ,nested = TRUE, seed = 111, cum.prop = 0.9)


prediction <- data.frame()
for(i in 1:5){
  train_df_ID <- spc_df %>% dplyr::distinct(x, y, .keep_all = TRUE) %>% .[!is.na(.$ORCDRC), ] %>% .[folds$outer.fold != i, "ID"]
  train_df <- spc_df %>% dplyr::filter(ID %in% train_df_ID$ID)
  pred_df_all <- spc_df %>% dplyr::filter(ID %ni% train_df_ID$ID)
  pred_df_all$prediction <-   hs3D(soil.fun = ORCDRC.fun, pred.data = pred_df_all, obs.data = train_df, d_inc = 15, n = 5, radius = 7000, p = 1, min_obs = 10, output = "prediction")
  prediction <- rbind(prediction, pred_df_all[, c("ID", "Top", "Bottom", "ORCDRC", "prediction")])
}
 

rsq(prediction, truth = ORCDRC, estimate = prediction)



