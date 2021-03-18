

#soil.fun = soil.fun; obs.data = n.obs.data; pred.data = n.pred.data; depth_th = 10; n.obs = n.obs; output = "preparation"

#' @title soil3D_knn
#' @description Perform interpolation of soil variables in 3D by usin knn approach.
#' @param soil.fun formula object that defines target and predictor variables which will be used in searching for neighboring profiles. (eg. `OC~x + y + DEM + TWI + depth`. In this case, `x`, `y`, `DEM` and `TWI` are included in the searching algorithm. `depth` is mandatory for 3D soil mapping.
#' @param obs.data  Observation data. Tidy `data.frame` with the following columns: `ID`, `Top`, `Bottom`, `Target variable` column, `Covariates` columns, `x`, `y`, `depth`.
#' @param pred.data Prediction data. Tidy `data.frame` with the following columns: `ID`, `Covariates` columns, `x`, `y`, `depth`.
#' @param n.obs number of neighboring profiles.Tuning parameter.
#' @param depth_th depth treshold. Serves to define prediction depth interval: `pdi = pred.data$depth +/- depth_th`, Tuning parameter.
#' @param p power of distance. Tuning parameter. Default: 1
#' @param output prediction or preparation. If output = "prediction" prediction will be applied, Default: 'prediction'
#' @return prediction based on IDW and gower distance.
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname soil3D_knn
#' @export 
#' @importFrom dplyr filter distinct left_join mutate case_when group_by summarise ungroup arrange
#' @importFrom gower gower_topn
#' 
soil3D_knn <- function(soil.fun, obs.data, pred.data, n.obs, depth_th, p = 1, output = "prediction"){
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
    dplyr::mutate(dh = dplyr::case_when(pred_top >= Top & pred_bot <= Bottom ~ pred_bot - pred_top,
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



#' @title tune_soil3D_knn
#' @description performs the cross-validation based tuning of meta-parameters.
#' @param soil.fun formula object that defines target and predictor variables which will be used in searching for neighboring profiles. (eg. `OC~x + y + DEM + TWI + depth`. In this case, `x`, `y`, `DEM` and `TWI` are included in the searching algorithm. `depth` is mandatory for 3D soil mapping.
#' @param obs.data  Observation data. Tidy `data.frame` with the following columns: `ID`, `Top`, `Bottom`, `Target variable` column, `Covariates` columns, `x`, `y`, `depth`.
#' @param n.obs number of neighboring profiles.Tuning parameter.
#' @param params `data.frame` of possible values of parameters (`n.obs`, `depth_th`, `p`).
#' @param folds vector indicating which profiles belong to each fold.
#' @return params data.frame with the corresponding accuracy measures
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname tune_soil3D_knn
#' @export 
#' @importFrom dplyr mutate distinct filter group_split
#' @importFrom tibble rowid_to_column
#' @importFrom purrr map
tune_soil3D_knn <- function(soil.fun, obs.data, params, folds){
  params <- params %>% dplyr::mutate(rsq = as.numeric(NA))
  for(k in 1:dim(params)[1]){
    prediction <- data.frame()
    for(i in 1:length(unique(folds))){ 
      train_df_ID <- obs.data %>% dplyr::distinct(x, y, .keep_all = TRUE) %>% .[folds != i, "ID"] # OC mora genericki
      train_df <- obs.data %>% dplyr::filter(ID %in% train_df_ID$ID)
      pred_df <- obs.data %>% dplyr::filter(ID %ni% train_df_ID$ID)
      pred_list <- pred_df %>% tibble::rowid_to_column() %>% dplyr::group_split(rowid, .keep = FALSE)
      pred_df$prediction <-   purrr::map(.x = pred_list, .f = ~soil3D_knn(soil.fun = soil.fun, obs.data = train_df, pred.data = .x, depth_th = params$depth_th[k], n.obs = params$n.obs[k], p = params$p[k])) %>% unlist()
      prediction <- rbind(prediction, pred_df[, c("ID", "Top", "Bottom", all.vars(soil.fun)[1], "prediction")])
    }
    params[k, "rsq"] <- rsq(prediction, truth = eval(parse(text = all.vars(soil.fun)[1])), estimate = prediction)[3]$.estimate
  }
  return(params)
}



#' @title ncv_soil3D_knn
#' @description performs the nested cross-validation procedure for assessing prediction accuracy.
#' @param soil.fun formula object that defines target and predictor variables which will be used in searching for neighboring profiles. (eg. `OC~x + y + DEM + TWI + depth`. In this case, `x`, `y`, `DEM` and `TWI` are included in the searching algorithm. `depth` is mandatory for 3D soil mapping.
#' @param obs.data  Observation data. Tidy `data.frame` with the following columns: `ID`, `Top`, `Bottom`, `Target variable` column, `Covariates` columns, `x`, `y`, `depth`.
#' @param n.obs number of neighboring profiles.Tuning parameter.
#' @param params `data.frame` of possible values of parameters (`n.obs`, `depth_th`, `p`).
#' @param folds data.frame with two columns indicating which profiles belong to each fold in outer and inner loop of the nested cross-validation.
#' @return Accuracy measures
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ncv_soil3D_knn
#' @export 
#' @importFrom dplyr distinct filter group_split
#' @importFrom tibble rowid_to_column
#' @importFrom purrr map
#' 
ncv_soil3D_knn <- function(soil.fun, obs.data, params, folds){
  outer.results <- as.list(rep(NA, length(unique(folds$outer.fold))))
  for(i in 1:length(unique(folds$outer.fold))){
    train_ID <- obs.data %>% dplyr::distinct(x, y, .keep_all = TRUE) %>% .[folds$outer.fold != i, "ID"]
    train_df <- obs.data %>% dplyr::filter(ID %in% train_ID$ID)
    pred_df <- obs.data %>% dplyr::filter(ID %ni% train_ID$ID)
    inner.folds <- folds[folds$outer.fold != i, ] %>% .$inner.fold
    tune.results <- tune_soil3D_knn(soil.fun = soil.fun, obs.data = train_df, folds = inner.folds, params = params)
    best.tune <- tune.results[which.max(tune.results$rsq), 1:3]
    pred_list <- pred_df %>% tibble::rowid_to_column() %>% dplyr::group_split(rowid, .keep = FALSE)
    pred_df$prediction <-   purrr::map(.x = pred_list, .f = ~soil3D_knn(soil.fun = soil.fun, obs.data = train_df, pred.data = .x, depth_th = best.tune$depth_th, n.obs = best.tune$n.obs, p = best.tune$p)) %>% unlist()
    outer.results[[i]] <- pred_df[, c("ID", "Top", "Bottom", all.vars(soil.fun)[1], "prediction")]
  }
  outer.rsq <- rsq(outer.results %>% bind_rows(), truth = eval(parse(text = all.vars(soil.fun)[1])), estimate = prediction)[3]$.estimate
  results <- list(obs.pred = outer.results,  rsq = outer.rsq)
  return(results)
}




#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param partitioning.data PARAM_DESCRIPTION
#' @param num.folds PARAM_DESCRIPTION
#' @param princ.comp PARAM_DESCRIPTION
#' @param num.means PARAM_DESCRIPTION
#' @param spatial.cluster PARAM_DESCRIPTION, Default: TRUE
#' @param seed PARAM_DESCRIPTION, Default: 46
#' @param cum.prop PARAM_DESCRIPTION, Default: 0.9
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname partitioning

partitioning <- function(partitioning.data, num.folds, princ.comp, num.means, spatial.cluster = TRUE, seed = 46, cum.prop = 0.9){
  
  na_count <- sapply(partitioning.data, function(y) sum(length(which(is.na(y)))))
  if(sum(na_count) > 0){stop("there is NA in the data for clustering")}
  names(partitioning.data) <- c("target", names(partitioning.data)[-1])
  
  if(princ.comp){
    if((dim(partitioning.data)[2]-1) <= 2){stop("Principal components analysis need more than two variables")
    }
    
    prc <- prcomp(as.formula(paste(paste(" ~ "), paste(names(partitioning.data[, -1]), collapse="+"))), data = partitioning.data[, -1])
    temp <- c()
    
    if(sum(summary(prc)$importance[3,] >= cum.prop) == 0){
      stop("None of principal components is greater than cum.prop")
    }
    
    for(i in 1:sum(summary(prc)$importance[3,] >= cum.prop)){
      temp <- cbind(temp, as.matrix(partitioning.data[, -1]) %*% as.matrix(prc[[2]][,i]))
    }
    partitioning.data <- data.frame(target = partitioning.data[,1], temp)
  }
  
  km.quality <- c()
  seed.seq <- seed + seq(0, 1000, 1)
  for (i in seed.seq){
    set.seed(i)
    km.temp <- kmeans(partitioning.data[, -1],  centers = num.means)$tot.withinss
    km.quality <- c(km.quality, km.temp)
  }
  seed.list <- seed.seq[which(km.quality == min(km.quality))]
  km <- as.list(rep(NA, length(seed.list)))
  for(i in 1:length(seed.list)){
    set.seed(seed.list[i])
    km[[i]] <- kmeans(partitioning.data[, -1],  centers = num.means)
  }
  min.size.clusters <- sapply(km, function(x) min(x$size))
  max.min.size.cluster <- which(min.size.clusters == max((sapply(km, function(x) min(x$size)))))[1]
  
  partitioning.data$cluster <- as.factor(km[[max.min.size.cluster]]$cluster)
  partitioning.data$ID <- c(1:dim(partitioning.data)[1])
  
  
  # Creating empty list to contain stratified folds
  cluster.list <- as.list(rep(NA,length(unique(partitioning.data$cluster))))
  names(cluster.list) <- paste("cluster",c(1:length(cluster.list)), sep="")
  
  # Each cluster is split to folds, stratified according to profile depth and weighted mean of observed target variable in the profile
  
  for(i in 1:length(cluster.list)){
    set.seed(seed)
    cluster.list[[i]] <- caret::createFolds(partitioning.data[which(partitioning.data$cluster == levels(partitioning.data$cluster)[i]),"target"], k = num.folds)
    if(length(cluster.list[[i]]) < num.folds){stop(paste("There is no enough data in cluster", i, sep = " "))}
    for(j in 1:num.folds){
      cluster.list[[i]][[j]] <- partitioning.data[which(partitioning.data$cluster == levels(partitioning.data$cluster)[i]),"ID"][cluster.list[[i]][[j]]]
    }
  }
  
  # List containing a vector of observation IDs that constitute each fold
  fold.list <- as.list(rep(NA, num.folds))
  names(fold.list) <- paste("fold", c(1:num.folds), sep = "")
  for(i in 1:num.folds){
    fold.list[[i]] <- do.call(c,lapply(cluster.list,function(x) x[[i]]))
    names(fold.list[[i]]) <- NULL
  }
  
  # Each row is augmented by fold number
  data.with.folds <- data.frame()
  for(i in 1:length(fold.list)){
    tmp <- partitioning.data[fold.list[[i]], ]
    tmp$fold <- i
    data.with.folds <- rbind(tmp, data.with.folds) %>% arrange(ID)
  }
  
  return(data.with.folds)
  
}


#' @title stratpart
#' @description perform n-fold stratified partitioning of data 
#' @param target.name names of target variables
#' @param sp.data "SpatialPointsDataFrame" object of data
#' @param num.folds Number of folds, Default: 5
#' @param princ.comp logical, do the principal compoment shoud be performed, Default: TRUE
#' @param kmean.vars Names of variables used for k-means clustering, Default: NA
#' @param num.means Number of clusters, Default: 3
#' @param spatial.cluster logical, if spatial clustering must be done, Default: TRUE
#' @param nested logical, if two step partitioning (outer and inner) must be done, Default: TRUE
#' @param seed random number generator seed, Default: 45
#' @param cum.prop proportion treshold for principal component analysis, Default: 0.9
#' @return vector of length sp.data, with the numbers indicating which point belongs to each fold. In case of nested = TRUE, the output is data.frame with two vectors, one for outer and one for inner partitioning.
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname stratpart
#' @export 
#' @importFrom dplyr rename select
#' 
#' 
stratpart <- function(target.name, sp.data, num.folds = 5, princ.comp = TRUE, kmean.vars = NA, num.means = 3, spatial.cluster = TRUE ,nested = TRUE, seed = 45, cum.prop = 0.9){
  
  if(!is.na(kmean.vars)){
    if(spatial.cluster){
      kmean.vars <- c(colnames(coordinates(sp.data)), kmean.vars)
      temp.data <- as.data.frame(sp.data)[, c(target.name, kmean.vars)]
    }
  }else{
    kmean.vars <- colnames(coordinates(sp.data))
    temp.data <- as.data.frame(sp.data)[, c(target.name, kmean.vars)]
  }
  
  part.data <- partitioning(partitioning.data = temp.data, num.folds = num.folds, princ.comp = princ.comp, num.means = num.means, seed = seed, cum.prop = cum.prop)    
  
  if(nested){
    part.data <- dplyr::rename(part.data, outer.fold = fold)
    part.data$inner.fold <- rep(NA, dim(part.data)[1])
    for(i in 1:num.folds){
      inner.temp <- partitioning(partitioning.data = dplyr::filter(part.data, outer.fold == i) %>% dplyr::select(-c(cluster:inner.fold)), num.folds = num.folds-1, princ.comp = princ.comp, num.means = num.means, seed = seed, cum.prop = cum.prop)
      inner.temp[which(inner.temp$fold == i), "fold"] <- num.folds
      part.data[which(part.data$outer.fold == i),"inner.fold"] <- inner.temp$fold
    }
    return(part.data[,c("outer.fold","inner.fold")])
  }else{
    return(part.data[,"outer.fold"])
  }
  
}
