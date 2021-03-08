
hs_knn_pred <- function(soil.fun, spc_sf, pred_depth, depth_th, n, p, output = list("prediction", "preparation")){
  output <- output[[1]]
  target.name <- all.vars(soil.fun)[1]
  gdist = unique(spc_sf$gw)[1:n]
  n_ids <- unique(spc_sf$ID)[1:n] %>% data.frame(ID = ., stringsAsFactors = FALSE)
  spc_sf %>% dplyr::mutate(pred_top = pmax(pred_depth-depth_th, 0), pred_bot = pred_depth+depth_th) %>% 
    dplyr::filter(Top < pred_bot, Bottom > pred_top) %>% 
    mutate(dh = case_when(pred_top >= Top & pred_bot <= Bottom ~ pred_bot - pred_top,
                          pred_top >= Top & pred_bot >= Bottom ~ Bottom - pred_top,
                          pred_top <= Top & pred_bot >= Bottom ~ Bottom - Top,
                          pred_top <= Top & pred_bot <= Bottom ~ pred_bot - Top)) %>%
    dplyr::filter(., ID %in% n_ids$ID) %>% dplyr::left_join(n_ids, ., by = "ID") %>%
    dplyr::group_by(ID) %>% dplyr::summarise(obs_mean = weighted.mean(eval(parse(text = target.name)), w = dh, na.rm = TRUE)) %>%
    dplyr::ungroup() %>% dplyr::mutate(gower_distance = gdist)
  if(output == "prediction"){
    results <- spc_sf %>% dplyr::summarise(pred_mean = sum(obs_mean/gower_distance^p)/sum(1/gower_distance^p))
  }else{
    results <- spc_sf
  }
    return(results)
}




hs3D_knn <- function(soil.fun, pred.data, obs.data, depth_th, n, p = 2, min_obs = 10, radius = NA, output = list("prediction", "preparation")){
  output <- output[[1]]
  obs_sf <- sf::st_as_sf(obs.data, coords = c("x", "y"), remove = FALSE)
  pred_sf <- sf::st_as_sf(pred.data, coords = c("x", "y"), remove = FALSE)
  if(is.na(radius)){radius <- max(st_distance(obs_sf))/3}
  
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
    purrr::map2(.x = ., .y = obs_list, .f = ~.x < depth_th/diff(range(.y$depth)))

  obs_list2 <- purrr::map2(.x = gdist_list, .y = obs_list, .f = ~.y[.x, ]) # Ovde se može dodati selekcija u okviru nekog buffer-a
  
  ex_obs <- unlist(map(.x = obs_list2, .f = ~dim(.x)[1])) < min_obs # TRUE for obs.data with no obs.data
  
  if(any(ex_obs)) warning(paste("There is no enough available obs.data for prediction at location of : ", paste(which(ex_obs), collapse = " ")))
  
  if(sum(!ex_obs) == 0) stop("There is no available obs.data on target depth")
  
  pred.data <- pred.data[!ex_obs, ]
  obs_list2 <- obs_list2[!ex_obs]
  
  pred_list <- pred.data %>% tibble::rowid_to_column() %>% dplyr::group_split(rowid, .keep = FALSE)

  gdist_list1 <- purrr::map2(.x = pred_list, .y = obs_list2, .f = ~gower_topn(x = .x[, cov.names], y = .y[, cov.names], n = dim(.y)[1]))
  
  obs_list3 <- purrr::map2(.x = obs_list2, .y = gdist_list1, .f = ~.x[.y$index, ]) %>%
    purrr::map2(.x = ., .y = gdist_list1, .f = ~dplyr::mutate(.x, gw = .y$distance[,1]))

  if(output == "prediction"){
    prediction <- purrr::map2(.x = obs_list3, .y = pred_list, .f = ~hs_knn_pred(soil.fun = soil.fun, spc_sf = .x, pred_depth = .y$depth, depth_th = depth_th, n = n, p = p)) %>%
      bind_rows()
    output <- rep(NA, length(pred_depth))
    output[!ex_obs] <- prediction$pred_mean
  }else{
    preparation <- purrr::map2(.x = obs_list3, .y = pred_list, .f = ~hs_knn_prep(soil.fun = soil.fun, spc_sf = .x, pred_depth = .y$depth, depth_th = depth_th, n = n, output = output))
    output <- as.list(rep(NA, length(pred_depth)))
    output[!ex_obs] <- preparation
  }
  return(output)
}
