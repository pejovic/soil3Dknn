# Preprocessing functions
#base.model = lead.fun; sp.data = meuse; coord.trend = TRUE; use.interactions = TRUE; poly.vars = NA; poly.deg = 1; interaction.var = NA; standardize = TRUE
preproc <- function(base.model, sp.data, use.interactions = FALSE, interaction.var = NA, poly.vars = NA, poly.deg = 2, coord.trend = FALSE, coord.names = NA, standardize = TRUE, num.folds = 5, num.means = NA, seed = 46){
  
  "%ni%" <- Negate("%in%")
  if(coord.trend){
    coord.names <- colnames(coordinates(sp.data))
  }
  
  
  # check if the number of variables for polynomial expansion is in line with the total number of explanatory variables.
  if(!is.na(poly.vars[1])){
    if(length(poly.vars) > length(all.vars(base.model)[-1])){
      stop("Number of variables for polynomial expansion is greater than total number of explanatory variables")
    }
  }
  
  # check if length of poly.vars and poly.deg is equal
  if(sum(is.na(poly.vars)) == 0){
    if(length(poly.vars) != length(poly.deg)){
      stop("length of poly.vars and poly.deg must be equal")
    }
  }
  
  target.name <- all.vars(base.model)[1]
  cov.names <- all.vars(base.model)[-1]
  
  # Extracting the names of categorical variagles and removing empty classes
  factor.names <- sp.data@data %>% subset(., select=which(sapply(., is.factor))) %>% names()
  
  if(length(factor.names) > 0){
    for(i in factor.names){
      sp.data@data[,i] <- factor(sp.data@data[,i])
    }
  }
  
  # check if variables for polynomial expansion are in factor.names:
  if(sum(poly.vars %in% factor.names) !=0){
    stop("poly.vars contain factor variables")
  }
  
  if(coord.trend){base.model <- as.formula(paste(target.name,"~", paste(c(cov.names, coord.names), collapse="+"))); cov.names <- c(cov.names, coord.names)}
  
  # check if variables for polynomial expansion are in base model:
  if(!is.na(poly.vars[1])){
    if(sum(poly.vars %in% all.vars(base.model)[-1]) == 0){
      stop("None of the variables for polynomial expansion matches the names in the base model")
    }
  }

  
  sp.data.df <- as.data.frame(sp.data) %>% dplyr::select(all.vars(base.model))
  
  # Adding polynomial depth terms in input data matrix, only if poly.deg > 1
  if(sum(!is.na(poly.vars)) != 0){
    for(i in 1:length(poly.vars)){
      sp.data.df <- cbind(sp.data.df, poly(sp.data.df[,poly.vars[i]], poly.deg[i], raw = TRUE, simple = TRUE)[,-1])
      names(sp.data.df) <- c(names(sp.data.df)[1:(length(names(sp.data.df))-(poly.deg[i]-1))], (paste(poly.vars[i], c(2:poly.deg[i]), sep="")))
      base.model <- as.formula(paste(target.name,"~", paste(c(all.vars(base.model)[-1], paste(poly.vars[i], c(2:poly.deg[i]),sep="")), collapse="+")))
    }
  }
  
  # Dummy coding
  if(length(factor.names) > 0){
    dummy.par <- dummyVars(as.formula(paste("~", paste(c(all.vars(base.model))[-1], collapse="+"))), sp.data.df, levelsOnly = FALSE)
    cov.data <- as.data.frame(predict(dummy.par, newdata = sp.data.df))
  }else{
    cov.data <- sp.data.df[,-1]
  }
  
  # Names
  names(cov.data) <- gsub( "\\_|/|\\-|\"|\\s" , "." , names(cov.data) )
  
  base.model <- as.formula(paste(target.name,"~", paste(names(cov.data), collapse="+")))
  cov.names <- all.vars(base.model)[-1]
  
  if(use.interactions){
    f <- as.formula(~ .^2)
    cov.data <- as.data.frame(model.matrix(f, cov.data)[,-1])
    
    if(!is.na(interaction.var[1])){
      interaction.ind <- c()
      for(i in 1:length(interaction.var)){
        temp<- unique(c(which(startsWith(names(cov.data), interaction.var[i])), grep(paste(":", interaction.var[i], sep = ""), names(cov.data))))
        interaction.ind <- unique(c(interaction.ind, temp))
      }
      interaction.names <- names(cov.data)[interaction.ind]
      interaction.names <- interaction.names[grep(":", interaction.names)]
      
      base.model <- as.formula(paste(target.name,"~", paste(c(all.vars(base.model)[-1], interaction.names), collapse="+")))
      cov.data <- cov.data %>% dplyr::select(c(all.vars(base.model)[-1], interaction.names))
    }
  }
  
  # Standardization of input data
  if(standardize) {
    cnt.par <- as.data.frame(cov.data) %>% preProcess(.,method=c("center", "scale"))
    cov.data <- predict(cnt.par, newdata = cov.data)
  }
  
  data.df <- data.frame(sp.data.df[,1], cov.data)
  names(data.df) <- c(target.name, names(cov.data))
  
  return(data.df)
}


partitioning <- function(partitioning.data, num.folds, princ.comp, num.means, spatial.cluster = TRUE, seed = 46, cum.prop = 0.9){
  
  na_count <-sapply(partitioning.data, function(y) sum(length(which(is.na(y)))))
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
    cluster.list[[i]] <- createFolds(partitioning.data[which(partitioning.data$cluster == levels(partitioning.data$cluster)[i]),"target"], k = num.folds)
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
      inner.temp <- partitioning(partitioning.data = filter(part.data, outer.fold == i) %>% dplyr::select(-c(cluster:inner.fold)), num.folds = num.folds-1, princ.comp = princ.comp, num.means = num.means, seed = seed, cum.prop = cum.prop)
      inner.temp[which(inner.temp$fold == i), "fold"] <- num.folds
      part.data[which(part.data$outer.fold == i),"inner.fold"] <- inner.temp$fold
    }
    return(part.data[,c("outer.fold","inner.fold")])
  }else{
    return(part.data[,"outer.fold"])
  }
  
}


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
