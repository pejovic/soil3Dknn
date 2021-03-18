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
#=============================== soil3D_kidw ==================================================

soil3D_kidw <- function(soil.fun, obs.data, pred.data, n.obs, depth_th, p = 1, output = "prediction"){
  target.name <- all.vars(soil.fun)[1]
  coord_names <- c("x", "y")
  cov.names <- all.vars(soil.fun)[-1]
  cov.names <- cov.names[-which(cov.names == "depth")]
  
  in.obs.data <- obs.data %>% dplyr::filter(depth <= pred.data$depth + depth_th & depth >= pmax(pred.data$depth - depth_th))
  in.obs.data.ID <- unique(in.obs.data$ID)
  if(length(in.obs.data.ID) < n.obs){
    while(length(in.obs.data.ID) < n.obs) {
      depth_th <- depth_th + 5
      in.obs.data <- obs.data %>% dplyr::filter(depth <= pred.data$depth + depth_th & depth >= pmax(pred.data$depth - depth_th))
      in.obs.data.ID <- unique(in.obs.data$ID)
    }
  }

  gower_search <- gower::gower_topn(x = pred.data[, cov.names], y = in.obs.data %>% dplyr::distinct(x, y, .keep_all = TRUE) %>% .[, cov.names], n = n.obs)
  
  in.obs.data.ID <- in.obs.data %>% dplyr::distinct(x, y, .keep_all = TRUE) %>% .[as.numeric(gower_search$index), "ID"]
  in.obs.data <- suppressMessages(in.obs.data %>% dplyr::left_join(data.frame(ID = in.obs.data.ID, gower_distance = gower_search$distance), .))

  in.obs.data <- in.obs.data %>% dplyr::mutate(pred_top = pmax(pred.data$depth-depth_th, 0), pred_bot = pred.data$depth+depth_th) %>% 
    dplyr::filter(Top < pred_bot, Bottom > pred_top) %>% 
    mutate(dh = dplyr::case_when(pred_top >= Top & pred_bot <= Bottom ~ pred_bot - pred_top,
                                 pred_top >= Top & pred_bot >= Bottom ~ Bottom - pred_top,
                                 pred_top <= Top & pred_bot >= Bottom ~ Bottom - Top,
                                 pred_top <= Top & pred_bot <= Bottom ~ pred_bot - Top)) %>%
    dplyr::group_by(ID) %>% dplyr::summarise(obs_mean = weighted.mean(eval(parse(text = target.name)), w = dh/(2*depth_th), na.rm = TRUE),
                                             gower_distance = mean(gower_distance)) %>%
    dplyr::ungroup() %>% dplyr::arrange(gower_distance)
  
  if(output == "prediction"){
    results <- in.obs.data %>% dplyr::summarise(pred_mean = sum(obs_mean/gower_distance^p)/sum(1/gower_distance^p)) %>% .$pred_mean
  }else{
    results <- in.obs.data
  }
  return(results)
}


#===================== Tune =============================================================

tune_soil3D_kidw <- function(soil.fun, obs.data, params, folds){
  params <- params %>% dplyr::mutate(rsq = as.numeric(NA))
  for(k in 1:dim(params)[1]){
    prediction <- data.frame()
    for(i in 1:length(unique(folds))){ 
      train_df_ID <- obs.data %>% dplyr::distinct(x, y, .keep_all = TRUE) %>% .[folds != i, "ID"] # OC mora genericki
      train_df <- obs.data %>% dplyr::filter(ID %in% train_df_ID$ID)
      pred_df <- obs.data %>% dplyr::filter(ID %ni% train_df_ID$ID)
      pred_list <- pred_df %>% tibble::rowid_to_column() %>% dplyr::group_split(rowid, .keep = FALSE)
      pred_df$prediction <-   purrr::map(.x = pred_list, .f = ~soil3D_kidw(soil.fun = soil.fun, obs.data = train_df, pred.data = .x, depth_th = params$depth_th[k], n.obs = params$n.obs[k], p = params$p[k])) %>% unlist()
      prediction <- rbind(prediction, pred_df[, c("ID", "Top", "Bottom", all.vars(soil.fun)[1], "prediction")])
    }
    params[k, "rsq"] <- rsq(prediction, truth = eval(parse(text = all.vars(soil.fun)[1])), estimate = prediction)[3]$.estimate
  }
  return(params)
}
#=========================================================================================


#=========================== ncv =========================================================

#soil.fun = m3.fun; obs.data = bor_df; folds = bor_folds; params = bor_params

ncv_soil3D_kidw <- function(soil.fun, obs.data, params, folds){
  outer.results <- as.list(rep(NA, length(unique(folds$outer.fold))))
  for(i in 1:length(unique(folds$outer.fold))){
    train_ID <- obs.data %>% dplyr::distinct(x, y, .keep_all = TRUE) %>% .[folds$outer.fold != i, "ID"]
    train_df <- obs.data %>% dplyr::filter(ID %in% train_ID$ID)
    pred_df <- obs.data %>% dplyr::filter(ID %ni% train_ID$ID)
    inner.folds <- folds[folds$outer.fold != i, ] %>% .$inner.fold
    tune.results <- tune_soil3D_kidw(soil.fun = soil.fun, obs.data = train_df, folds = inner.folds, params = params)
    best.tune <- tune.results[which.max(tune.results$rsq), 1:3]
    pred_list <- pred_df %>% tibble::rowid_to_column() %>% dplyr::group_split(rowid, .keep = FALSE)
    pred_df$prediction <-   purrr::map(.x = pred_list, .f = ~soil3D_kidw(soil.fun = soil.fun, obs.data = train_df, pred.data = .x, depth_th = best.tune$depth_th, n.obs = best.tune$n.obs, p = best.tune$p)) %>% unlist()
    outer.results[[i]] <- pred_df[, c("ID", "Top", "Bottom", all.vars(soil.fun)[1], "prediction")]
  }
  outer.rsq <- rsq(outer.results %>% bind_rows(), truth = eval(parse(text = all.vars(soil.fun)[1])), estimate = prediction)[3]$.estimate
  results <- list(obs.pred = outer.results,  rsq = outer.rsq)
  return(results)
}


rsq(outer.prediction, truth = OC, estimate = prediction)[3]$.estimate

rsq(outer.prediction, truth = eval(parse(text = all.vars(soil.fun)[1])), estimate = prediction)[3]$.estimate


#===========================================================================================


#================== Test for Bor data ======================================================
load(here::here("Bor_data", "bor_data.rda"))
"%ni%" <- Negate("%in%")

bor_df <- bor_list[[1]]
bor_folds <- bor_list[[2]]
bor_params <- expand.grid(depth_th = seq(3, 15, 3), n.obs = seq(4, 8, 1), p = c(1:2))

m1.fun <- as.formula(paste("OC ~", paste(c("x", "y", "depth"), collapse="+")))
ncv.m1.bor <- ncv_soil3D_kidw(soil.fun = m1.fun, obs.data = bor_df, folds = bor_folds, params = bor_params)

m2.fun <- as.formula(paste("OC ~", paste(names(bor_df[-c(1:4)]), collapse="+")))
ncv.m2.bor <- ncv_soil3D_kidw(soil.fun = m2.fun, obs.data = bor_df, folds = bor_folds, params = bor_params)

m3.fun <- as.formula(paste("OC ~", paste(c("x", "y", "DEM", "TWI", "Slope", "depth"), collapse="+")))
ncv.m3.bor <- ncv_soil3D_kidw(soil.fun = m3.fun, obs.data = bor_df, folds = bor_folds, params = bor_params)


library(foreach)
bor_results <- foreach(f = c(m1.fun, m2.fun, m3.fun)) %dopar% ncv_soil3D_kidw(soil.fun = f, obs.data = bor_df, folds = bor_folds, params = bor_params)

#===========================================================================================


#================== Test for Edgeroi data ==================================================
load(here::here("Bor_data", "edgeroi_data.rda"))
"%ni%" <- Negate("%in%")

edgeroi_df <- edgeroi_list[[1]] %>% dplyr::rename(OC = ORCperc)
edgeroi_folds <- edgeroi_list[[2]]
edgeroi_params <- expand.grid(depth_th = seq(3, 15, 3), n.obs = seq(4, 8, 1), p = c(1:2))

m1.fun <- as.formula(paste("OC ~", paste(c("x", "y", "depth"), collapse="+")))
#ncv.m1.edgeroi <- ncv_soil3D_kidw(soil.fun = m1.fun, obs.data = edgeroi_df, folds = edgeroi_folds, params = edgeroi_params)

m2.fun <- as.formula(paste("OC ~", paste(names(edgeroi_df[-c(1:4)]), collapse="+")))
#ncv.m2.edgeroi <- ncv_soil3D_kidw(soil.fun = m2.fun, obs.data = edgeroi_df, folds = edgeroi_folds, params = edgeroi_params)

m3.fun <- as.formula(paste("OC ~", paste(c("x", "y", "DEMSRT5", "TWISRT5", "PMTGEO5", "depth"), collapse="+")))
#ncv.m3.edgeroi <- ncv_soil3D_kidw(soil.fun = m3.fun, obs.data = edgeroi_df, folds = edgeroi_folds, params = edgeroi_params)

library(foreach)
library(doParallel)
registerDoParallel(cores=6)

edgeroi_results <- foreach(f = c(m1.fun, m2.fun, m3.fun)) %dopar% ncv_soil3D_kidw(soil.fun = f, obs.data = edgeroi_df, folds = edgeroi_folds, params = edgeroi_params)

#save(edgeroi_results, file = "edgeroi_results.rda")



plot.obs <- function(soil.fun, obs.data, n.pred, n.obs){
  id <- obs.data[n.pred, ]$ID
  n.obs.data <- filter(obs.data, ID != id)
  n.pred.data <-  obs.data[n.pred, ]
  search.obs <- soil3D_kidw(soil.fun = soil.fun, obs.data = n.obs.data, pred.data = n.pred.data, depth = n.pred.data$depth, n.obs = n.obs, output = "preparation")
  search.obs <- suppressMessages(dplyr::left_join(search.obs, n.obs.data)) %>% dplyr::group_by(ID) %>% 
    dplyr::summarise(obs_mean = mean(obs_mean),gower_dist = mean(gower_distance), x = mean(x), y = mean(y)) %>%
    dplyr::ungroup()
  
  obs.plot <- ggplot(data = n.obs.data %>% dplyr::group_by(ID) %>% dplyr::summarise(OC.mean = mean(OC), x = mean(x), y = mean(y))) + 
    geom_point(aes(x = x, y = y), shape = 1, color = "black") + 
    geom_point(data = n.pred.data, aes(x = x, y = y, size = depth), shape = 17, color = "blue") + 
    geom_point(data = search.obs, aes(x = x, y = y, size = obs_mean), color = "red") +
    theme_bw()
  return(obs.plot)
  
}


plot.obs(soil.fun = m2.fun, obs.data = edgeroi_df, n.pred = 285, n.obs = 6)


#============ plot folds ===================

plotfolds <- function(sp.data, target.name, folds){
  data <- as.data.frame(sp.data)
  data$ID <- c(1:dim(data)[1])
  data <- data %>% dplyr::rename(x = names(.)[which(names(data) %in% colnames(coordinates(sp.data)))[1]], y = names(.)[which(names(data) %in% colnames(coordinates(sp.data)))[2]], target = names(.)[which(names(data) == target.name)[1]] )
  data <- cbind(data, folds)
  data$inner.fold <- factor(data$inner.fold)
  q <- ggplot(data, aes(x = x, y = y, fill = inner.fold))
  r <- q + geom_point(aes(size = sqrt(target/pi)), pch = 21, show.legend = FALSE) + scale_size_continuous(range=c(1,7))
  r <- r + facet_wrap(~ outer.fold) + theme_bw()
  plot(r)
}


edgeroi_sp <- edgeroi_df %>% dplyr::group_by(ID) %>% 
  dplyr::summarise(OC = mean(OC), x = mean(x), y = mean(y), max.depth = max(depth), n.obs = n()) %>%
  dplyr::ungroup()



coordinates(edgeroi_sp) <- ~x + y

edgeroi_folds_plot <- plotfolds(sp.data = edgeroi_sp, target.name = "OC", folds = edgeroi_folds)

ggsave("edgeroi_folds_plot.jpg", plot = edgeroi_folds_plot)

#===========================================================================================



