extract.year<-function(character.string){
  out<-as.numeric(strsplit(character.string, ".", fixed=TRUE)[[1]][2])
  return(out)
}

extract.month<-function(character.string){
  out<-strsplit(character.string, ".", fixed=TRUE)[[1]][1]
  return(out)
}

match_year_to_corine <- function(year){
  if (year >= 2003 & year <= 2008) out <- 2006
  if (year >= 2009 & year <= 2014) out <- 2012
  if (year >= 2015 & year <= 2020) out <- 2018
  return(out)
}

identify_clc_index <- function(X_var_corine){
  temp <- as.numeric(strsplit(strsplit(X_var_corine, "_", fixed=TRUE)[[1]][1],"X",fixed=TRUE)[[1]][2])
  return(temp)
}

zero_test_fn <- function(numeric_vector){
  out <- sum(abs(numeric_vector))
  return(out)
}

quick_glm <- function(corine_to_use, month_string){
  index_of_explanatory_variables <- grep(paste0("_",corine_to_use), names(combined.df))
  N_exp_vars <- length(index_of_explanatory_variables)
  index_month_corine_to_use <- which(which_corine_lookup$corine_to_use==corine_to_use
                                     & which_corine_lookup$month == month_string)
  index_of_responses <- which(names(combined.df) %in% 
                                which_corine_lookup$date[index_month_corine_to_use])
  curlew_prop <- apply(combined.df[,index_of_responses], 1, sum)/length(index_of_responses)
  temp <- combined.df[,index_of_explanatory_variables[1:(N_exp_vars-1)]]
  x <- as.matrix(temp)
  N<-dim(x)[1]
  weights <- rep(length(index_of_responses), N)
  temp <- cbind(temp, curlew_prop)
  out.glm <- glm(curlew_prop ~ ., data=temp, family=binomial, weights=weights)
  return(out.glm)
}


breaks_prob_fn <- function(probability, threshold_vector){
  threshold_vector[1] <- -1
  N <- length(threshold_vector)
  for (i in 1:(N-1)){
    if (probability > threshold_vector[i] & probability <= threshold_vector[i+1]) out <- i
  }
  return(out)
}


colour_codes_prob_vector_fn <- function(prob_vector, probs=seq(0.25, 1, 0.25)){
  min_prob <- min(prob_vector)
  non_min_probs <- prob_vector[which(prob_vector>min_prob)]
  prob_quantiles <- quantile(non_min_probs, probs)
  threshold_vector <- c(min_prob, prob_quantiles)
  out <- sapply(prob_vector, breaks_prob_fn, threshold_vector=threshold_vector)
  return(out)
}

ggplot_of_predicted_prob_surface <- function(rasterstack_for_preds, glm_object, title=NULL){
  temp <-predict(rasterstack_for_preds, glm_object, type="response")
  temp2 <- projectRaster(temp, crs=geo.prj)
  temp2_masked <- mask(temp2, polygons_spdf)
  temp2_pts <- rasterToPoints(temp2_masked, spatial=TRUE)
  temp2_df <- data.frame(temp2_pts)
  probs <- c(0, 0.25, 0.5, 0.75, 1)
  temp2_df$col <- colour_codes_prob_vector_fn(temp2_df$layer, probs=probs)
  
  out <-  ggplot() +
    geom_polygon(data=UKmap, aes(x=long, y=lat, group=group),
                 fill='gray90',
                 color='black') +
    geom_tile(data = filter(temp2_df, col>1) , aes(x = x, y = y, fill = as.factor(col))) + 
    ggtitle(title) +
    labs(x = "long", y="lat", fill="Model probability class") +
    scale_fill_manual(values=safe_colorblind_palette[1:(length(probs)-1)],
                      labels=c("Low", "Moderate", "High", "Very high")) +
    theme_void()
  return(out)
}

raster_of_predicted_prob_surface_long_lat <- function(rasterstack_for_preds, 
                                                      glm_object){
  temp <-predict(rasterstack_for_preds, glm_object, type="response")
  temp2 <- projectRaster(temp, crs=geo.prj)
  temp2_masked <- mask(temp2, polygons_spdf)
  temp2_pts <- rasterToPoints(temp2_masked, spatial=TRUE)
  temp2_df <- data.frame(temp2_pts)
  probs <- c(0, 0.25, 0.5, 0.75, 1)
  temp2_df$col <- colour_codes_prob_vector_fn(temp2_df$layer, probs=probs)
  temp2_filtered <- filter(temp2_df, col>1)
  out <- rasterize(temp2_filtered[,c("x","y")], temp2_masked, 
                   field=temp2_filtered$col)
  return(out)
}


quick_glm <- function(corine_to_use, month_string){
  index_of_explanatory_variables <- grep(paste0("_",corine_to_use), names(combined.df))
  N_exp_vars <- length(index_of_explanatory_variables)
  index_month_corine_to_use <- which(which_corine_lookup$corine_to_use==corine_to_use
                                     & which_corine_lookup$month == month_string)
  index_of_responses <- which(names(combined.df) %in% 
                                which_corine_lookup$date[index_month_corine_to_use])
  curlew_prop <- apply(combined.df[,index_of_responses], 1, sum)/length(index_of_responses)
  temp <- combined.df[,index_of_explanatory_variables[1:(N_exp_vars-1)]]
  x <- as.matrix(temp)
  N<-dim(x)[1]
  weights <- rep(length(index_of_responses), N)
  temp <- cbind(temp, curlew_prop)
  out.glm <- glm(curlew_prop ~ ., data=temp, family=binomial, weights=weights)
  return(out.glm)
}

glmnet_particular_month_aggregated <- function(corine_to_use, month_string, use_cv.glmnet=FALSE){
  out_cv.glmnet<-NULL
  index_of_explanatory_variables <- grep(paste0("_",corine_to_use), names(combined.df))
  index_of_clc_classes <- sapply(names(combined.df)[index_of_explanatory_variables],
                                 identify_clc_index)
  N_exp_vars <- length(index_of_explanatory_variables)
  index_month_corine_to_use <- which(which_corine_lookup$corine_to_use==corine_to_use
                                     & which_corine_lookup$month == month_string)
  index_of_responses <- which(names(combined.df) %in% 
                                which_corine_lookup$date[index_month_corine_to_use])
  curlew_prop <- apply(combined.df[,index_of_responses], 1, sum)/length(index_of_responses)
  temp <- combined.df[,index_of_explanatory_variables[1:(N_exp_vars-1)]]
  x <- as.matrix(temp)
  y <- cbind(curlew_prop, 1-curlew_prop)
  N<-dim(x)[1]
  weights <- rep(length(index_of_responses), N)
  out.glmnet <- glmnet(x, y, family=binomial, weights=weights)
  if (use_cv.glmnet) out_cv.glmnet <- cv.glmnet(x, y, family=binomial, weights=weights)
  out <- list (out.glmnet = out.glmnet, out_cv.glmnet = out_cv.glmnet, x=x, y=y)
  return(out)
}

return_coefs_from_set_of_cv.glmnet_objects <- function(list_of_cv.glmnet_objects,
                                                       vector_month_strings,
                                                       vector_corine_to_use){
  N_objects <- length(list_of_cv.glmnet_objects)
  if (length(vector_month_strings)!=N_objects |
      length(vector_corine_to_use)!=N_objects){ 
    stop("vector of months or vector of corine years incorrect length")
  }
  coefs_matrix <- matrix(0, length(clc_labels), N_objects)
  rownames(coefs_matrix)[2:dim(coefs_matrix)[1]] <- as.character(clc_labels[1:(length(clc_labels)-1)])
  rownames(coefs_matrix)[1] <- "Intercept"
  colnames(coefs_matrix) <- paste0(vector_month_strings,"_",vector_corine_to_use)
  for (i in 1:N_objects){
    print(i)
    corine_to_use <- vector_corine_to_use[i]
    cv.glmnet_object <- list_of_cv.glmnet_objects[[i]]
    coefs <- as.matrix(coef(cv.glmnet_object$out_cv.glmnet,  s="lambda.1se"))
    #    index_of_explanatory_variables <- grep(paste0("_",corine_to_use), names(combined.df))
    #    index_of_clc_classes <- sapply(names(combined.df)[index_of_explanatory_variables],
    #                                   identify_clc_index)
    index_of_explanatory_variables <- grep(paste0("_",corine_to_use), vector_of_names)
    index_of_clc_classes <- sapply(vector_of_names[index_of_explanatory_variables],
                                   identify_clc_index)
    N_exp_vars <- length(index_of_explanatory_variables)
    coefs_matrix[c(1,index_of_clc_classes[1:(N_exp_vars-1)]+1),i] <- coefs[,1]
  }
  coefs.df <- as.data.frame(coefs_matrix)
  labels <- c("Intercept", as.character(clc_labels[1:(length(clc_labels)-1)]))
  coefs.df$label <- factor(labels, levels=labels, ordered=TRUE)
  out <- coefs.df
  return(out)
}

return_coefs_from_set_of_glmnet_objects <- function(list_of_glmnet_objects,
                                                       vector_month_strings,
                                                       vector_corine_to_use){
  N_objects <- length(list_of_glmnet_objects)
  if (length(vector_month_strings)!=N_objects |
      length(vector_corine_to_use)!=N_objects){ 
    stop("vector of months or vector of corine years incorrect length")
  }
  coefs_matrix <- matrix(0, length(clc_labels), N_objects)
  rownames(coefs_matrix)[2:dim(coefs_matrix)[1]] <- as.character(clc_labels[1:(length(clc_labels)-1)])
  rownames(coefs_matrix)[1] <- "Intercept"
  colnames(coefs_matrix) <- paste0(vector_month_strings,"_",vector_corine_to_use)
  for (i in 1:N_objects){
    # print(i)
    corine_to_use <- vector_corine_to_use[i]
    glmnet_object <- list_of_glmnet_objects[[i]]
    coefs <- as.matrix(coef(glmnet_object$out.glmnet,  s=0))
    #    index_of_explanatory_variables <- grep(paste0("_",corine_to_use), names(combined.df))
    #    index_of_clc_classes <- sapply(names(combined.df)[index_of_explanatory_variables],
    #                                   identify_clc_index)
    index_of_explanatory_variables <- grep(paste0("_",corine_to_use), vector_of_names)
    index_of_clc_classes <- sapply(vector_of_names[index_of_explanatory_variables],
                                   identify_clc_index)
    N_exp_vars <- length(index_of_explanatory_variables)
    coefs_matrix[c(1,index_of_clc_classes[1:(N_exp_vars-1)]+1),i] <- coefs[,1]
  }
  coefs.df <- as.data.frame(coefs_matrix)
  labels <- c("Intercept", as.character(clc_labels[1:(length(clc_labels)-1)]))
  coefs.df$label <- factor(labels, levels=labels, ordered=TRUE)
  out <- coefs.df
  return(out)
}



cv.glmnet_month_corine <- function(month_string, corine_to_use){
  assign(paste0(month_string,"_",corine_to_use,".cv.glmnet"), 
         glmnet_particular_month_aggregated(corine_to_use, 
                                            month_string, 
                                            use_cv.glmnet = TRUE))
  cv.glmnet_object <- get(paste0(month_string,"_",corine_to_use,".cv.glmnet"))
  filename <- paste0("data/",month_string,"_",corine_to_use,".cv.glmnet.rds")
  saveRDS(cv.glmnet_object, filename)
}


plot_coefs <- function(cv.glmnet_object, month_string=NULL, corine_to_use=NULL){
  if (is.null(month_string) & is.null(corine_to_use)){
    name_of_object <- deparse(substitute(cv.glmnet_object))
    temp <- strsplit(name_of_object, "_",fixed=TRUE)
    month_string <- temp[[1]][1]
    corine_to_use <- strsplit(temp[[1]][2], ".", fixed=TRUE)[[1]][1]
  }
  coefs <- as.matrix(coef(cv.glmnet_object$out_cv.glmnet,  s="lambda.1se"))
  index_of_explanatory_variables <- grep(paste0("_",corine_to_use), names(combined.df))
  index_of_clc_classes <- sapply(names(combined.df)[index_of_explanatory_variables],
                                 identify_clc_index)
  N_exp_vars <- length(index_of_explanatory_variables)
  rownames(coefs)[2:dim(coefs)[1]] <- 
    clc_labels[index_of_clc_classes[1:(N_exp_vars-1)]]
  colnames(coefs) <- "coefficient"
  coefs <- as.data.frame(coefs)
  coefs$label <- rownames(coefs)
  index <- which(coefs$coefficient!=0)
  coefs <- coefs[index,]
  p <- ggplot(coefs, aes(x=label, y=coefficient)) +
    geom_bar(stat="identity") +
    coord_flip() +
    ggtitle(paste(month_string,corine_to_use,"coefficients from CV lasso GLM (Curlew)", sep=" ")) +
    labs(x="Corine land cover class", y="Coefficient")
  return(p)
}

return_coefs_from_set_of_cv.glmnet_objects <- function(list_of_cv.glmnet_objects,
                                                       vector_month_strings,
                                                       vector_corine_to_use){
  N_objects <- length(list_of_cv.glmnet_objects)
  if (length(vector_month_strings)!=N_objects |
      length(vector_corine_to_use)!=N_objects){ 
    stop("vector of months or vector of corine years incorrect length")
  }
  coefs_matrix <- matrix(0, length(clc_labels), N_objects)
  rownames(coefs_matrix)[2:dim(coefs_matrix)[1]] <- as.character(clc_labels[1:(length(clc_labels)-1)])
  rownames(coefs_matrix)[1] <- "Intercept"
  colnames(coefs_matrix) <- paste0(vector_month_strings,"_",vector_corine_to_use)
  for (i in 1:N_objects){
    print(i)
    corine_to_use <- vector_corine_to_use[i]
    cv.glmnet_object <- list_of_cv.glmnet_objects[[i]]
    coefs <- as.matrix(coef(cv.glmnet_object$out_cv.glmnet,  s="lambda.1se"))
    #    index_of_explanatory_variables <- grep(paste0("_",corine_to_use), names(combined.df))
    #    index_of_clc_classes <- sapply(names(combined.df)[index_of_explanatory_variables],
    #                                   identify_clc_index)
    index_of_explanatory_variables <- grep(paste0("_",corine_to_use), vector_of_names)
    index_of_clc_classes <- sapply(vector_of_names[index_of_explanatory_variables],
                                   identify_clc_index)
    N_exp_vars <- length(index_of_explanatory_variables)
    coefs_matrix[c(1,index_of_clc_classes[1:(N_exp_vars-1)]+1),i] <- coefs[,1]
  }
  coefs.df <- as.data.frame(coefs_matrix)
  labels <- c("Intercept", as.character(clc_labels[1:(length(clc_labels)-1)]))
  coefs.df$label <- factor(labels, levels=labels, ordered=TRUE)
  out <- coefs.df
  return(out)
}

# create matrix of non-zero coefficients for individual CORINE years
remove_zero_rows <- function(coef_matrix, cols_to_keep){
  coef_matrix <- coef_matrix[,cols_to_keep]
  coef_matrix$zeroflag <- apply(coef_matrix, 1, zero_test_fn)==0
  coef_matrix_nozeros <- filter(coef_matrix, zeroflag==FALSE)
  zeroflag_index <- which(colnames(coef_matrix_nozeros)=="zeroflag")
  out <- coef_matrix_nozeros[,-zeroflag_index]
  return(out)
}

## Creating unconstrained GLM based on lasso output

unconstrained_glm_from_cv.glmnet_coefs <- function(cv.glmnet_coefs_matrix, corine_to_use, month_string){
  col_of_interest <- which(colnames(cv.glmnet_coefs_matrix) == 
                             paste0(month_string, "_", corine_to_use))
  index_of_exp_vars <- (which(cv.glmnet_coefs_matrix[,col_of_interest] !=0)-1)[-1]
  # The -1 here is because the first row contains the intercept for the glm
  # We then need to drop the initial 0 corresponding to the intercept row
  xnam <- paste0("X", index_of_exp_vars, "_", corine_to_use)
  fmla <- as.formula(paste("curlew_prop ~ ", paste(xnam, collapse= "+")))
  index_month_corine_to_use <- which(which_corine_lookup$corine_to_use==corine_to_use
                                     & which_corine_lookup$month == month_string)
  index_of_responses <- which(names(combined.df) %in% 
                                which_corine_lookup$date[index_month_corine_to_use])
  curlew_prop <- apply(combined.df[,index_of_responses], 1, sum)/length(index_of_responses)
  index_of_exp_vars_in_combined <- which(names(combined.df) %in% xnam)
  index_of_clc_classes <- sapply(names(combined.df)[index_of_exp_vars_in_combined],
                                 identify_clc_index)
  temp <- combined.df[,index_of_exp_vars_in_combined]
  temp <- set_label(temp, label = clc_labels[index_of_clc_classes])
  temp$curlew_prop <- curlew_prop
  temp$weights <- length(index_of_responses)
  out <- list(model=glm(fmla, data=temp, weights=weights, family=binomial),
              var_names = as.character(clc_labels[index_of_clc_classes]))
  return(out)
}

unconstrained_glm_from_glmnet_coefs <- function(glmnet_coefs_matrix, corine_to_use, month_string){
  col_of_interest <- which(colnames(glmnet_coefs_matrix) == 
                             paste0(month_string, "_", corine_to_use))
  index_of_exp_vars <- (which(glmnet_coefs_matrix[,col_of_interest] !=0)-1)[-1]
  # The -1 here is because the first row contains the intercept for the glm
  # We then need to drop the initial 0 corresponding to the intercept row
  xnam <- paste0("X", index_of_exp_vars, "_", corine_to_use)
  fmla <- as.formula(paste("curlew_prop ~ ", paste(xnam, collapse= "+")))
  index_month_corine_to_use <- which(which_corine_lookup$corine_to_use==corine_to_use
                                     & which_corine_lookup$month == month_string)
  index_of_responses <- which(names(combined.df) %in% 
                                which_corine_lookup$date[index_month_corine_to_use])
  curlew_prop <- apply(combined.df[,index_of_responses], 1, sum)/length(index_of_responses)
  index_of_exp_vars_in_combined <- which(names(combined.df) %in% xnam)
  index_of_clc_classes <- sapply(names(combined.df)[index_of_exp_vars_in_combined],
                                 identify_clc_index)
  temp <- combined.df[,index_of_exp_vars_in_combined]
  temp <- set_label(temp, label = clc_labels[index_of_clc_classes])
  temp$curlew_prop <- curlew_prop
  temp$weights <- length(index_of_responses)
  out <- list(model=glm(fmla, data=temp, weights=weights, family=binomial),
              var_names = as.character(clc_labels[index_of_clc_classes]))
  return(out)
}


unconstrained_glms_from_average_sgn_coefs_matrix <- function(average_sgn_coefs_matrix, 
                                                      corine_years=c(2006, 2012, 2018), 
                                                             Months_to_use){
  coefs <- average_sgn_coefs_matrix
  xnam <- paste0(rownames(coefs), "_", corine_to_use)
  fmla <- as.formula(paste("curlew_prop ~ ", paste(xnam, collapse= "+")),
                     "+X101+X102+X103+X104+X105+X201+X202")
  index_month_corine_to_use <- which(which_corine_lookup$corine_to_use==corine_to_use
                                     & which_corine_lookup$month == month_string)
  index_of_responses <- which(names(combined.df) %in% 
                                which_corine_lookup$date[index_month_corine_to_use])
  curlew_prop <- apply(combined.df[,index_of_responses], 1, sum)/length(index_of_responses)
  index_of_exp_vars_in_combined <- which(names(combined.df) %in% xnam)
  index_of_clc_classes <- sapply(names(combined.df)[index_of_exp_vars_in_combined],
                                 identify_clc_index)
  temp <- combined.df[,index_of_exp_vars_in_combined]
  temp <- set_label(temp, label = clc_labels[index_of_clc_classes])
  temp$curlew_prop <- curlew_prop
  temp$weights <- length(index_of_responses)
  out <- list(model=glm(fmla, data=temp, weights=weights, family=binomial),
              var_names = as.character(clc_labels[index_of_clc_classes]))
  return(out)
}





plot_coefs_unconstrained_glm <- function(cv.glmnet_coefs_matrix, corine_to_use, month_string){
  final_glm_object <- unconstrained_glm_from_cv.glmnet_coefs(cv.glmnet_coefs_matrix, corine_to_use, month_string)
  model <- final_glm_object$model
  var_names <- final_glm_object$var_names
  coefs <- as.matrix(model$coefficients)
  rownames(coefs)[2:dim(coefs)[1]] <- var_names
  colnames(coefs) <- "coefficient"
  coefs <- as.data.frame(coefs)
  coefs$label <- rownames(coefs)
  index <- which(coefs$coefficient!=0)
  coefs <- coefs[index,]
  p <- ggplot(coefs, aes(x=label, y=coefficient)) +
    geom_bar(stat="identity") +
    coord_flip() +
    ggtitle(paste(month_string,corine_to_use,"coefficients final GLM (Curlew)", sep=" ")) +
    labs(x="Corine land cover class", y="Coefficient")
  return(p)
}



# functions for Boyce indices

maxpred <-function(cv.glmnet_object){
  out <- max(cv.glmnet_object$pred)
  return(out)
}

find_max_pred_value <- function(list_of_cv.glmnet_objects){
  out <- max(unlist(lapply(list_of_cv.glmnet_objects, maxpred)))
  return(out)
}

data_frame_for_boyce_plots <- function(list_of_cv.glmnet_objects, xlim=c(0,1), N_windows=10, titles=titles){
  N_objects <- length(list_of_cv.glmnet_objects)
  xmax <- find_max_pred_value(list_of_cv.glmnet_objects)
  xlim <- c(0,xmax)
  window_width <- (xlim[2]-xlim[1])/N_windows
  xvalues <- seq(xlim[1]+window_width/2, xlim[2]-window_width/2, by=window_width)
  predicted_frequencies <- matrix(0, length(xvalues), N_objects)
  expected_frequencies <- matrix(0, length(xvalues), N_objects)
  predicted_to_expected <- matrix(0, length(xvalues), N_objects)
  means <- matrix(0, length(xvalues), N_objects)
  boyce <- numeric(N_objects)
  
  for (i in 1:length(list_of_cv.glmnet_objects)){  
    print(i)
    temp <- list_of_cv.glmnet_objects[[i]]
    temp.df <- data.frame(presence = temp$y[,1]!=0, pred=temp$pred[,1])
    total_presences <- sum(temp.df$presence==TRUE)
    N_preds <- dim(temp.df)[1]
    
    for (j in 1:length(xvalues)){
      window <- c(xvalues[j]-window_width/2, xvalues[j]+window_width/2)
      temp2.df <- filter(temp.df, pred>=window[1] & pred<window[2])
      actual_presences <- sum(temp2.df$presence==TRUE)
      N_temp2 <- dim(temp2.df)[1]
      predicted_frequencies[j,i] <- actual_presences/total_presences
      expected_frequencies[j,i] <- N_temp2/dim(temp.df)[1]
      predicted_to_expected[j,i] <- predicted_frequencies[j,i]/expected_frequencies[j,i]
      means[j,i] <- mean(window)
    }
    boyce[i] <- cor(means[,i], predicted_to_expected[,i], use="complete.obs")
  }
  
  
  boyce.df <- data.frame(mean = as.numeric(means), 
                         predicted_to_expected=as.numeric(predicted_to_expected), 
                         label = rep(titles, rep(length(xvalues), length(titles))))
  boyce.df$label <- factor(boyce.df$label, labels=month_string_vector, ordered=TRUE)
  out <- list(boyce.df = boyce.df, boyce = boyce)
  return(out)
}

boyce_calc <- function(corr_probs, obs, xlim=c(0,1), N_windows=10){
  xmax <- max(corr_probs)
  xlim <- c(0,xmax)
  window_width <- (xlim[2]-xlim[1])/N_windows
  xvalues <- seq(xlim[1]+window_width/2, xlim[2]-window_width/2, by=window_width)
  predicted_frequencies <- numeric(length(xvalues))
  expected_frequencies <- numeric(length(xvalues))
  predicted_to_expected <- numeric(length(xvalues))
  means <- numeric(length(xvalues))
  total_presences <- sum(obs)
  temp.df <- data.frame(obs=obs, pred=corr_probs)

      for (i in 1:length(xvalues)){
      window <- c(xvalues[i]-window_width/2, xvalues[i]+window_width/2)
      temp2.df <- filter(temp.df, pred>=window[1] & pred<window[2])
      actual_presences <- sum(temp2.df$obs==1)
      N_temp2 <- dim(temp2.df)[1]
      predicted_frequencies[i] <- actual_presences/total_presences
      expected_frequencies[i] <- N_temp2/dim(temp.df)[1]
      predicted_to_expected[i] <- predicted_frequencies[i]/expected_frequencies[i]
      means[i] <- mean(window)
    }
    boyce <- cor(means, predicted_to_expected, use="complete.obs")
    out <- boyce
    return(out)
}



boyce_plots <- function(boyce_object, year){
  boyce.df <- boyce_object$boyce.df
  boyce <- boyce_object$boyce
  text_values <- paste("BI = ", round(boyce,2), sep="")
  month_labels <- month_string_vector
  facet_labs <- paste0(month_labels, " (", text_values, ")")
  names(facet_labs) <- month_string_vector
  ggplot(boyce.df, aes(x=mean, y=predicted_to_expected)) +
    geom_line() +
    labs(x="mean predicted value", y="Predicted to expected") +
    facet_wrap(~label, labeller=labeller(label=facet_labs))+
    ggtitle(paste("Boyce index plots for ", year, sep=""))
}

# Empirical cdfs

create_df_for_ecdf <- function(list_of_cv.glmnet_objects, 
                               month_corine_vec = paste0(month_string_rep, "_", corine_to_use_rep)){
  list_of_dfs <- vector("list", 36)
  for (i in 1:36){
    temp <- list_of_cv.glmnet_objects[[i]]
    list_of_dfs[[i]] <- data.frame(presence = temp$y[,1]!=0, pred=temp$pred[,1])
    list_of_dfs[[i]]$month_corine <- month_corine_vec[i]
  }
  out <- bind_rows(list_of_dfs)
  out$month_corine <- factor(out$month_corine, levels=month_corine_vec, ordered=TRUE)
  return(out)
}


convert_curlew_obs_to_newCRS <- function(curlew_obs_df, newCRS){
  points <- st_as_sf(curlew_obs_df, coords=c("corine_east", "corine_east"),
                     crs=CORINE_CRS)
  new_points <- st_transform(points, crs=newCRS)
  out <- new_points
  return(out)
}



plot_curlew_obs_year_month <- function(curlew_long_lat_year_month, 
                                       obs_year, 
                                       obs_month, title=NULL){
  Month <- Month_strings[obs_month]
  if (is.null(title)) title <- paste("Curlew observations \n",Month, obs_year)
  out <- ggplot() + 
    geom_polygon(data=UKmap, aes(x=long, y=lat, group=group),
                 fill='gray90',
                 color='black') +
    geom_sf(data=curlew_long_lat_year_month, size=0.1, col="red") + 
    ggtitle(title) +
    theme_bw()
  return(out)
}



presabs <- function(r, pts){
  r2 <- r
  r2[] <- 0
  counts <- table(cellFromXY(r, pts))
  r2[as.numeric(names(counts))] <- ifelse(counts>=1,1,0)
  return(r2)
}

extract_coords_month_year <- function(observation_df, Year, Month){
  out <- filter(observation_df, obs_year==Year & obs_month==Month)[, c("corine_east", "corine_north")]
  return(out)
}

create_presence_absence_raster_month_year <- function(observation_df, Year, Month){
  coords <- extract_coords_month_year(observation_df, Year, Month)
  out <- presabs(corine_2006[[1]], coords)
  return(out)
}


corine_Year.Month_vector <- function(corine_to_use, Month_numeric){
  if (corine_to_use==2006){
    out <- paste0(2003:2008, ".", Month_numeric)
  }
  if (corine_to_use==2012){
    out <- paste0(2009:2014, ".", Month_numeric)
  }
  if (corine_to_use==2018){
    if (Month_numeric <= 8){
      out <- paste0(2015:2019, ".", Month_numeric)
    }
    if (Month_numeric > 8){
      out <- paste0(2015:2018, ".", Month_numeric)
    }
  }
  return(out)
}

corine_Year.Month_curlew_obs_index <- function(corine_to_use, Month_numeric,
                                               curlew_obs_list){
  L <- names(curlew_obs_list)
  if (corine_to_use == 2006){
    index <- seq(Month_numeric, Month_numeric + 60, by=12)
  }
  if (corine_to_use == 2012){
    index <- seq(Month_numeric + 72, Month_numeric+132, by=12)
  }
  if (corine_to_use == 2018){
    if (Month_numeric <= 8){
      index <- seq(Month_numeric + 144, Month_numeric + 192, by=12)
    }
    if (Month_numeric > 8){
      index <- seq(Month_numeric + 144, Month_numeric + 180, by=12)
    }
  }
  return(index)
}


build_exp_vars_rasterstack <- function(corine_to_use, Month_numeric,
                                       climate_vars_CORINE_proj_filestem){
  # CORINE landcover rasters for CORINE YEAR=corine_to_use
  corine <- get(paste0("corine_",corine_to_use))
  
  # Five climate variables month by month 
  # for CORINE window for Month=Month_numeric and CORINE YEAR=corine_to_use
  window <- corine_Year.Month_vector(corine_to_use, Month_numeric)
  drydaysSub_0_05mm_filenames <- paste0(climate_vars_CORINE_proj_filestem, "drydaysSub_0_05mm_",
                                        window, ".tif")
  drydaysSub_0_5mm_filenames <- paste0(climate_vars_CORINE_proj_filestem, "drydaysSub_0_5mm_",
                                       window, ".tif")
  frostdays_filenames <- paste0(climate_vars_CORINE_proj_filestem, "frostdays_",
                                window, ".tif")
  sub5Cdays_filenames <- paste0(climate_vars_CORINE_proj_filestem, "sub5Cdays_",
                                window, ".tif")
  totPrec_filenames <- paste0(climate_vars_CORINE_proj_filestem, "totPrec_",
                              window, ".tif")
  
  climate_vars_stack <- stack(c(drydaysSub_0_05mm_filenames,
                                drydaysSub_0_5mm_filenames,
                                frostdays_filenames,
                                sub5Cdays_filenames,
                                totPrec_filenames))
  exp_vars_stack <- stack(corine, climate_vars_stack, 
                          DTM_1km_OSGB_corine, Slope_percent_corine)  
  N_CORINE <- dim(corine)[3]
  
  N_exp_vars <- dim(exp_vars_stack)[3]
  names(exp_vars_stack)[1:N_CORINE] <- gsub(paste0("_",corine_to_use),"",names(corine))
  window <- corine_Year.Month_vector(corine_to_use, Month_numeric)
  L <- length(window)
  names(exp_vars_stack)[N_CORINE+1:L] <- paste0("X101.", 1:L)
  names(exp_vars_stack)[N_CORINE+L+1:L] <- paste0("X102.", 1:L)
  names(exp_vars_stack)[N_CORINE+2*L+1:L] <- paste0("X103.",1:L)
  names(exp_vars_stack)[N_CORINE+3*L+1:L] <- paste0("X104.",1:L)
  names(exp_vars_stack)[N_CORINE+4*L+1:L] <- paste0("X105.",1:L)
  names(exp_vars_stack)[(N_exp_vars-1):N_exp_vars] <- c("X201", "X202")
  
  
  return(exp_vars_stack)
}

month_year_string_to_year <- function(month_year_string){
  out <- as.numeric(strsplit(month_year_string, "_")[[1]][2])
  return(out)
}

corine_window_years <- function(corine_to_use, Month_numeric){
  if (corine_to_use == 2006) out <- 2003:2008
  if (corine_to_use == 2012) out <- 2009:2014
  if (corine_to_use == 2018){
    if (Month_numeric <= 8) out <- 2015:2019
    if (Month_numeric > 8) out <- 2015:2018
  }
  return(out)
}

build_response_vars_stack <- function(corine_to_use, Month_numeric,
                                      curlew_obs_list){
  L <- names(curlew_obs_list)
  Month_string <- Month_strings[Month_numeric]
  curlew_obs_years <- sapply(L, month_year_string_to_year)
  corine_window <- corine_window_years(corine_to_use, Month_numeric)
  index <- intersect(grep(Month_string, L), 
                     which(curlew_obs_years %in% corine_window))
  out <- stack(curlew_obs_list[index])
  N <- dim(out)[3]
  names(out) <- paste0("Y", 1:N)
  return(out)
}


save_exp_vars_and_response_vars_stacks <- function(output_filenamebase, 
                                                   corine_to_use, 
                                                   Month_numeric,
                                climate_vars_CORINE_proj_filestem){
  Month_string <- Month_strings[Month_numeric]
  exp_vars_stack <- build_exp_vars_rasterstack(corine_to_use, Month_numeric,
                                               climate_vars_CORINE_proj_filestem = 
                                                 climate_vars_CORINE_proj_filestem)
  response_vars_stack <- build_response_vars_stack(corine_to_use, 
                                                   Month_numeric,
                                                   curlew_obs_list)
  filename_exp_vars <- paste0(output_filenamebase, "exp_vars_CORINE", 
                              corine_to_use, "_", Month_string, "s",".tif")
  filename_response_vars <- paste0(output_filenamebase, "response_vars_CORINE", 
                                   corine_to_use, "_", Month_string,"s", ".tif")
  writeRaster(exp_vars_stack, filename_exp_vars, overwrite=TRUE)
  writeRaster(response_vars_stack, filename_response_vars, overwrite=TRUE)
}


rename_exp_vars_tif <- function(exp_vars_tif_file){
  exp_vars_brick <- brick(exp_vars_tif_file)
}

save_x_and_y_dfs_for_glmnet <- function(corine_to_use, Month_numeric,
                                        exp_and_response_vars_filestem,
                                        output_filestem){
  corine <- get(paste0("corine_", corine_to_use))
  Month <- Month_strings[Month_numeric]
  exp_vars_filename <- paste0(exp_and_response_vars_filestem,"exp_vars_CORINE", corine_to_use,"_",Month,"s.tif")
  response_vars_filename <- paste0(exp_and_response_vars_filestem,"response_vars_CORINE", corine_to_use,"_",Month,"s.tif")
  exp_vars_stack <- stack(exp_vars_filename)
  response_vars_stack <- stack(response_vars_filename)
  N_CORINE <- dim(corine)[3] # no. of land cover variables
  
  N_exp_vars <- dim(exp_vars_stack)[3]
  names(exp_vars_stack)[1:N_CORINE] <- gsub(paste0("_",corine_to_use),"",names(corine))
  window <- corine_Year.Month_vector(corine_to_use, Month_numeric)
  L <- length(window)
  names(exp_vars_stack)[N_CORINE+1:L] <- paste0("X101.", 1:L)
  names(exp_vars_stack)[N_CORINE+L+1:L] <- paste0("X102.", 1:L)
  names(exp_vars_stack)[N_CORINE+2*L+1:L] <- paste0("X103.",1:L)
  names(exp_vars_stack)[N_CORINE+3*L+1:L] <- paste0("X104.",1:L)
  names(exp_vars_stack)[N_CORINE+4*L+1:L] <- paste0("X105.",1:L)
  names(exp_vars_stack)[(N_exp_vars-1):N_exp_vars] <- c("X201", "X202")

# Averaging the five climate variables across years in the CORINE window  
  
  X101_index <- grep("X101", names(exp_vars_stack))
  X101 <- calc(exp_vars_stack[[X101_index]], mean)
  X102_index <- grep("X102", names(exp_vars_stack))
  X102 <- calc(exp_vars_stack[[X102_index]], mean)
  X103_index <- grep("X103", names(exp_vars_stack))
  X103 <- calc(exp_vars_stack[[X103_index]], mean)
  X104_index <- grep("X104", names(exp_vars_stack))
  X104 <- calc(exp_vars_stack[[X104_index]], mean)
  X105_index <- grep("X105", names(exp_vars_stack))
  X105 <- calc(exp_vars_stack[[X105_index]], mean)
  exp_vars_stack <- stack(exp_vars_stack, X101, X102, X103, X104, X105)
  N <- dim(exp_vars_stack)[3]
  names(exp_vars_stack)[(N-4):N] <- c("X101", "X102", "X103", "X104", "X105")
  
  N_response_vars <- dim(response_vars_stack)[3]
  names(response_vars_stack) <- paste0("Y", 1:N_response_vars)
  curlew_prop <- calc(response_vars_stack, mean)
  response_vars_stack <- stack(response_vars_stack, curlew_prop)
  N <- dim(response_vars_stack)[3]
  names(response_vars_stack)[N] <- "curlew_prop"
  
  exp_vars_and_response_vars_stack <- stack(exp_vars_stack, response_vars_stack)
  # writeRaster(exp_vars_and_response_vars_stack, paste0(filenamebase, "overall_CORINE",corine_to_use,"_",Month,"s.grd"), format="raster")
  
  exp_vars_and_response_vars.df <- data.frame(na.omit(values(exp_vars_and_response_vars_stack)))
  
  Y_index <- which(colnames(exp_vars_and_response_vars.df) %in% paste0("Y", 1:N_response_vars)) 
  weights <- rep(N_response_vars, dim(exp_vars_and_response_vars.df)[1])
  
  cols_to_use_x <- c(colnames(exp_vars_and_response_vars.df)[1:(N_CORINE-1)], 
                     "X101", "X102", "X103", "X104", "X105",
                     "X201", "X202")
  
  x <- exp_vars_and_response_vars.df[, cols_to_use_x]
  y <- cbind(1-exp_vars_and_response_vars.df[, "curlew_prop"],exp_vars_and_response_vars.df[,"curlew_prop"])
  y <- cbind(y, weights)
  colnames(y) <- c("no_sighting_prop", "sighting_prop", "weights")
  
  out <- list(x=as.matrix(x), y=as.matrix(y))
  saveRDS(out, file=paste0(output_filestem, "matrices_for_lasso_CORINE",corine_to_use,"_",Month,"s.rds"))
  
}


save_dfs_and_tifs_for_monthly_models <- function(index,
                                      exp_and_response_vars_filestem,
                                      output_filestem){
  corine_to_use <- corine_to_use_rep[index]
  Month_numeric <- Month_numeric_vector[index]
  corine <- get(paste0("corine_", corine_to_use))
  Month <- Month_strings[Month_numeric]
  exp_vars_filename <- paste0(exp_and_response_vars_filestem,
                              "exp_vars_CORINE", corine_to_use,"_",
                              Month,"s.tif")
  response_vars_filename <- paste0(exp_and_response_vars_filestem,
                                   "response_vars_CORINE", corine_to_use,"_",
                                   Month,"s.tif")
  exp_vars_stack <- stack(exp_vars_filename)
  response_vars_stack <- stack(response_vars_filename)
  N_CORINE <- dim(corine)[3]
  
  N_exp_vars <- dim(exp_vars_stack)[3]
  names(exp_vars_stack)[1:N_CORINE] <- gsub(paste0("_",corine_to_use),"",names(corine))
  window <- corine_Year.Month_vector(corine_to_use, Month_numeric)
  L <- length(window)
  names(exp_vars_stack)[N_CORINE+1:L] <- paste0("X101.", 1:L)
  names(exp_vars_stack)[N_CORINE+L+1:L] <- paste0("X102.", 1:L)
  names(exp_vars_stack)[N_CORINE+2*L+1:L] <- paste0("X103.",1:L)
  names(exp_vars_stack)[N_CORINE+3*L+1:L] <- paste0("X104.",1:L)
  names(exp_vars_stack)[N_CORINE+4*L+1:L] <- paste0("X105.",1:L)
  names(exp_vars_stack)[(N_exp_vars-1):N_exp_vars] <- c("X201", "X202")
  
  for (i in 1:L){
    print(paste0("i=",i))
    exp_vars_temp <- exp_vars_stack[[c(1:N_CORINE,
                                       N_CORINE + i,
                                       N_CORINE+L+i,
                                       N_CORINE+2*L+i,
                                       N_CORINE+3*L+i,
                                       N_CORINE+4*L+i,
                                       (N_exp_vars-1):N_exp_vars)]]
    names(exp_vars_temp) <- c(names(exp_vars_stack)[1:N_CORINE],
                              paste0("X10",1:5),
                              c("X201", "X202"))
    response_temp <- response_vars_stack[[i]]
    names(response_temp) <- "Y"
    exp_vars_and_response <- stack(exp_vars_temp, response_temp)
    #    tif_filename <- paste0(output_filestem,"month_",window[i], ".tif")
    #    writeRaster(exp_vars_and_response, tif_filename)
    exp_vars_and_response.df <- data.frame(na.omit(values(exp_vars_and_response)))
    filename <- paste0(output_filestem,"month_",window[i],".df.rds")
    saveRDS(exp_vars_and_response.df, filename)
  }
}



read_in_x_y_df_save_cv.glmnet <- function(index, random_seed,
                                          x_y_df_filestem,
                                          cv.glmnet_filestem,
                                          corine_to_use_rep,
                                          Month_numeric_vector,
                                          Month_strings){
  set.seed(random_seed)
  corine_to_use <- corine_to_use_rep[index]
  corine <- get(paste0("corine_", corine_to_use))
  N_CORINE <- dim(corine)[3]
  Month_numeric <- Month_numeric_vector[index]
  Month <- Month_strings[Month_numeric]
  x_y_df_filename <- paste0(x_y_df_filestem, "matrices_for_lasso_CORINE",
                            corine_to_use,"_", Month, "s.rds")
  x_y_df <- readRDS(x_y_df_filename)
  x <- x_y_df$x
  y <- x_y_df$y
  temp <- cv.glmnet(x, y[,1:2], family="binomial", weights=y[,3],
                    penalty.factor = c(rep(1, N_CORINE-1), rep(0, dim(x)[2]-(N_CORINE-1))))
  cv.glmnet_filename <- paste0(cv.glmnet_filestem,"CORINE",corine_to_use,
                               "_",Month,"s_cv.glmnet.rds")
  saveRDS(temp, cv.glmnet_filename)
}



correct_y_matrix <- function(corine_to_use, Month_numeric){
  Month <- Month_strings[Month_numeric]
  filename <- paste0(filenamebase, "matrices_for_lasso_CORINE",
                     corine_to_use,"_", Month,"s.rds")
  x_y_df <- readRDS(filename)
  y <- x_y_df$y
  temp <- matrix(as.numeric(y[,1:2]), ncol=2)
  window <- corine_Year.Month_vector(corine_to_use, Month_numeric)
  L <- length(window)
  weights <- rep(L, dim(y)[1])
  temp <- cbind(temp, weights)
  colnames(temp) <- c("F_prop", "S_prop", "weights")
  x_y_df$y <- temp
  saveRDS(x_y_df, filename)
}



find_match_index <- function(rownames_vector, overall_rownames,
                             inc_climate=FALSE){
  X101_index <- which(rownames_vector=="X101")
  L <- length(rownames_vector)
  if (inc_climate==FALSE) N <- X101_index-1
  if (inc_climate==TRUE) N <- length(rownames_vector)
  out <- numeric(N)
  for (i in 1:length(out)){
    out[i] <- which(overall_rownames == rownames_vector[i])
  }
  return(out)
}


x_all_zero <- function(x){
  temp <- sum(abs(x)==0)
  out <- ifelse(temp==length(x), 1, 0)
  return(out)
}


generate_month <- function(month_year_string){
  out <- strsplit(as.character(month_year_string), "20")[[1]][1]
  return(out)
}

generate_year <- function(month_year_string){
  out <- paste0("20",strsplit(as.character(month_year_string), "20")[[1]][2])
  return(out)
}

generate_year_long <- function(year_seed_string){
  temp <- strsplit(as.character(year_seed_string), "CORINE")[[1]][2]
  out <- strsplit(temp, "_")[[1]][1]
  return(out)
}

generate_glm_formula <- function(coef_matrix, col_index){
  vars_to_use <- rownames(coef_matrix)[which(coef_matrix[,col_index]!=0)]
  fmla <- as.formula(paste0("Y~", paste0(vars_to_use[1:length(vars_to_use)], 
                                         collapse="+")))
  out <- fmla
  return(out)
}



corine_to_use_fn <- function(Year){
  if (Year >= 2003 & Year <= 2008) corine_to_use <- 2006
  if (Year >= 2009 & Year <= 2014) corine_to_use <- 2012
  if (Year >=2015 & Year <= 2020) corine_to_use <- 2018
  out <- corine_to_use
  return(out)
}

fit_and_predict_GLM <- function(Year, Month_numeric, filenamebase, coef_matrix,
        monthly_dfs_filenamebase = "data/cleaned_curlew_obs_and_exp_vars_CORINE/monthly_dfs_and_tifs/"
){
  corine_to_use <- corine_to_use_fn(Year)
  Month <- Month_strings[Month_numeric]
  corine <- get(paste0("corine_", corine_to_use))
  exp_vars_filename <- paste0(filenamebase,"exp_vars_CORINE", corine_to_use,"_",Month,"s.tif")
  response_vars_filename <- paste0(filenamebase,"response_vars_CORINE", corine_to_use,"_",Month,"s.tif")
  exp_vars_stack <- stack(exp_vars_filename)
  response_vars_stack <- stack(response_vars_filename)
  N_CORINE <- dim(corine)[3]
  
  N_exp_vars <- dim(exp_vars_stack)[3]
  names(exp_vars_stack)[1:N_CORINE] <- gsub(paste0("_",corine_to_use),"",names(corine))
  window <- corine_Year.Month_vector(corine_to_use, Month_numeric)
  L <- length(window)
  names(exp_vars_stack)[N_CORINE+1:L] <- paste0("X101.", 1:L)
  names(exp_vars_stack)[N_CORINE+L+1:L] <- paste0("X102.", 1:L)
  names(exp_vars_stack)[N_CORINE+2*L+1:L] <- paste0("X103.",1:L)
  names(exp_vars_stack)[N_CORINE+3*L+1:L] <- paste0("X104.",1:L)
  names(exp_vars_stack)[N_CORINE+4*L+1:L] <- paste0("X105.",1:L)
  names(exp_vars_stack)[(N_exp_vars-1):N_exp_vars] <- c("X201", "X202")
  
  find_index <- function(window, Year){
    L <- length(window)
    temp <- numeric(L)
    for (i in 1:L){
      temp[i] <- as.numeric(strsplit(window[i],"[.]")[[1]][1])
    }
    index <- which(temp==Year)
    return(index)
  }
  
  i <- find_index(window, Year)

  exp_vars_temp <- exp_vars_stack[[c(1:N_CORINE,
                                     N_CORINE + i,
                                     N_CORINE+L+i,
                                     N_CORINE+2*L+i,
                                     N_CORINE+3*L+i,
                                     N_CORINE+4*L+i,
                                     (N_exp_vars-1):N_exp_vars)]]
  names(exp_vars_temp) <- c(names(exp_vars_stack)[1:N_CORINE],
                            paste0("X10",1:5),
                            c("X201", "X202"))
  response_temp <- response_vars_stack[[i]]
  names(response_temp) <- "Y"
  exp_vars_and_response <- stack(exp_vars_temp, response_temp)
  # This is a rasterstack
  
  filename <- paste0(monthly_dfs_filenamebase,"month_",window[i],".df.rds")
  exp_vars_and_response.df <- readRDS(filename)
  # This is the df (read in from memory)
  fmla <- generate_glm_formula(coef_matrix, col_index=Month_numeric)
  monthly.glm <- glm(formula = fmla, family=binomial,
                     data=exp_vars_and_response.df)
  #    predictions <- predict(exp_vars_and_response, monthly.glm, type="response")
  out <- list(monthly.glm = monthly.glm, rasterstack_for_preds = exp_vars_and_response)
  return(out)
}

generate_training_and_test_datasets <- function(exp_vars_and_response.df,
                                                random_seed){
  set.seed(random_seed)
  N <- dim(exp_vars_and_response.df)[1]
  presence_index <- which(exp_vars_and_response.df$Y == 1)
  absence_index <- which(exp_vars_and_response.df$Y == 0)
  N_presence <- length(presence_index)
  N_absence <- length(absence_index)
  N_presence_training <- ceiling(N_presence/2)
  N_absence_training <- ceiling (N_absence/2)
  training_index_presence <- sample(presence_index, N_presence_training)
  training_index_absence <- sample(absence_index, N_absence_training)
  training_index <- c(training_index_presence, training_index_absence)
  test_index <- (1:N)[-training_index]
  training.df <- exp_vars_and_response.df[training_index,]
  test.df <- exp_vars_and_response.df[test_index,]
  out <- list(training.df = training.df, test.df = test.df)
  return(out)
}


correct_prob_vector <- function(probability_vector, fitted_GLM){
  X <- t(coef(fitted_GLM))
  V_X <- vcov(fitted_GLM)
  XVXt <- X %*% V_X %*% t(X)
  C <- (0.5 - probability_vector) * probability_vector * (1-probability_vector) * c(XVXt)
  
  out <- probability_vector + C
  return(out)
}

correct_prob_raster <- function(prob_raster, fitted_GLM){
  X <- t(coef(fitted_GLM))
  V_X <- vcov(fitted_GLM)
  XVXt <- X %*% V_X %*% t(X)
  out <- (0.5-prob_raster)*prob_raster*(1-prob_raster)*c(XVXt)
  return(out)
}

fit_GLM_train_test <- function(Year, Month_numeric, 
                               exp_and_response_vars_filestem, coef_matrix,
                                monthly_dfs_filestem,
                               random_seed){
  corine_to_use <- corine_to_use_fn(Year)
  Month <- Month_strings[Month_numeric]
  corine <- get(paste0("corine_", corine_to_use))
  exp_vars_filename <- paste0(exp_and_response_vars_filestem,
                              "exp_vars_CORINE", corine_to_use,
                              "_",Month,"s.tif")
  response_vars_filename <- paste0(exp_and_response_vars_filestem,
                                   "response_vars_CORINE", corine_to_use,
                                   "_",Month,"s.tif")
  exp_vars_stack <- stack(exp_vars_filename)
  response_vars_stack <- stack(response_vars_filename)
  N_CORINE <- dim(corine)[3]
  
  N_exp_vars <- dim(exp_vars_stack)[3]
  names(exp_vars_stack)[1:N_CORINE] <- gsub(paste0("_",corine_to_use),"",names(corine))
  window <- corine_Year.Month_vector(corine_to_use, Month_numeric)
  L <- length(window)
  names(exp_vars_stack)[N_CORINE+1:L] <- paste0("X101.", 1:L)
  names(exp_vars_stack)[N_CORINE+L+1:L] <- paste0("X102.", 1:L)
  names(exp_vars_stack)[N_CORINE+2*L+1:L] <- paste0("X103.",1:L)
  names(exp_vars_stack)[N_CORINE+3*L+1:L] <- paste0("X104.",1:L)
  names(exp_vars_stack)[N_CORINE+4*L+1:L] <- paste0("X105.",1:L)
  names(exp_vars_stack)[(N_exp_vars-1):N_exp_vars] <- c("X201", "X202")
  
  find_index <- function(window, Year){
    L <- length(window)
    temp <- numeric(L)
    for (i in 1:L){
      temp[i] <- as.numeric(strsplit(window[i],"[.]")[[1]][1])
    }
    index <- which(temp==Year)
    return(index)
  }
  
  i <- find_index(window, Year)
  
  exp_vars_temp <- exp_vars_stack[[c(1:N_CORINE,
                                     N_CORINE + i,
                                     N_CORINE+L+i,
                                     N_CORINE+2*L+i,
                                     N_CORINE+3*L+i,
                                     N_CORINE+4*L+i,
                                     (N_exp_vars-1):N_exp_vars)]]
  names(exp_vars_temp) <- c(names(exp_vars_stack)[1:N_CORINE],
                            paste0("X10",1:5),
                            c("X201", "X202"))
  response_temp <- response_vars_stack[[i]]
  names(response_temp) <- "Y"
  exp_vars_and_response <- stack(exp_vars_temp, response_temp)
  # This is a rasterstack
  
  filename <- paste0(monthly_dfs_filestem,"month_",window[i],".df.rds")
  exp_vars_and_response.df <- readRDS(filename)
  # This is the df (read in from memory)
  
  train_test <- generate_training_and_test_datasets(exp_vars_and_response.df,
                                                    random_seed = random_seed)
  training.df <- train_test$training.df
  test.df <- train_test$test.df
  
  # find which column of coef_matrix corresponds to Month
  col_index <- grep(Month, colnames(coef_matrix))
  
  fmla <- generate_glm_formula(coef_matrix, col_index=col_index)
  training.glm <- glm(formula = fmla, family=binomial,
                     data=training.df)
  #    predictions <- predict(exp_vars_and_response, monthly.glm, type="response")
  out <- list(training.glm = training.glm, 
              training.df = training.df,
              test.df = test.df)
  return(out)
}




ggplot_of_predicted_prob_surface <- function(rasterstack_for_preds, glm_object, title=NULL){
  temp <-predict(rasterstack_for_preds, glm_object, type="response")
  temp2 <- projectRaster(temp, crs=geo.prj)
  temp2_masked <- mask(temp2, polygons_spdf)
  temp2_pts <- rasterToPoints(temp2_masked, spatial=TRUE)
  temp2_df <- data.frame(temp2_pts)
  probs <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1)
  temp2_df$col <- colour_codes_prob_vector_fn(temp2_df$layer, probs=probs)
  
  out <-  ggplot() +
    geom_polygon(data=UKmap, aes(x=long, y=lat, group=group),
                 fill='gray90',
                 color='black') +
    geom_tile(data = filter(temp2_df, col>1) , aes(x = x, y = y, fill = as.factor(col))) + 
    ggtitle(title) +
    labs(x = "long", y="lat", fill="Model probability class") +
    scale_fill_manual(values=cbPalette[1:(length(probs)-1)],
                      labels=c("below 20th percentile", 
                               "20th to 40th percentile", 
                               "40th to 60th percentile", 
                               "60th to 80th percentile",
                               "80th to 90th percentile",
                               "90th to 95th percentile",
                               "95th to 99th percentile",
                               "above 99th percentile")) +
    theme_void()
  return(out)
}



obs_preds_coefs_plot <- function(curlew_long_lat_sub, obs_year, 
                                 obs_month, coef_matrix,
                                 monthly_dfs_filenamebase = "data/cleaned_curlew_obs_and_exp_vars_CORINE/monthly_dfs_and_tifs/"
){
  Month <- Month_strings[obs_month]
  curlew_obs_plot <- plot_curlew_obs_year_month(curlew_long_lat_sub,
                                                obs_year=obs_year,
                                                obs_month=obs_month)
  GLM_year.month <- fit_and_predict_GLM(obs_year, obs_month,
                                        filenamebase, coef_matrix,
                          monthly_dfs_filenamebase=monthly_dfs_filenamebase)
  predicted_probs_plot <- ggplot_of_predicted_prob_surface(GLM_year.month$rasterstack_for_preds,
                                                           GLM_year.month$monthly.glm,
                                                           title=paste("Predicted probabilities for ",
                                                                       Month,obs_year))
  coefplot_month <- create_coefplot(GLM_year.month$monthly.glm)
  file <- paste0("figs/",Month,"_",obs_year,"_obs_preds_coefs.pdf")
  pdf(file, width=10)
  grid.arrange(curlew_obs_plot,
               predicted_probs_plot,
               coefplot_month, nrow=2)
  dev.off()
  out <- list(curlew_obs_plot=curlew_obs_plot,
              predicted_probs_plot = predicted_probs_plot,
              coefplot_month = coefplot_month)
  return(out)
}


check_var_coefs_plot <- function(variable_no){
  out <-  ggplot(trans_CORINE2006_df, aes(x=Month, y=get(paste0("X", variable_no)))) +
    geom_line() +
    facet_wrap(~Year_factor) +
    labs(y=clc_labels[variable_no]) +
    scale_x_continuous(breaks=1:12) +
    ggtitle(clc_labels[variable_no])
  return(out)
}


stripX <- function(Xstring){
  out <- as.numeric(strsplit(Xstring, "X")[[1]][2])
  return(out)
}


create_coefplot <- function(monthly.glm){
  oldlabels <- names(coefficients(monthly.glm))
  N <- length(oldlabels)
  oldlabels <- oldlabels[2:N]
  index_last_lcc <- which(oldlabels=="X101")-1
  newNames_vector <- numeric(length(oldlabels))
  for (i in 1:(index_last_lcc)){
    newNames_vector[i] <- clc_labels[stripX(oldlabels[i])]
  }
  newNames_vector[index_last_lcc+1] <- "Dry sub 0.05"
  newNames_vector[index_last_lcc+2] <- "Dry sub 0.5"
  newNames_vector[index_last_lcc+3] <- "Frost days"
  newNames_vector[index_last_lcc+4] <- "sub 5C days"
  newNames_vector[index_last_lcc+5] <- "Total precipitation"
  newNames_vector[index_last_lcc+6] <- "DTM"
  newNames_vector[index_last_lcc+7] <- "Slope percent"
  
  out <- coefplot(monthly.glm, coefficients = oldlabels) +
    scale_y_discrete(labels=newNames_vector)
  return(out)
}



breaks_prob_fn <- function(probability, threshold_vector){
  threshold_vector[1] <- -1
  N <- length(threshold_vector)
  for (i in 1:(N-1)){
    if (probability > threshold_vector[i] & probability <= threshold_vector[i+1]) out <- i
  }
  return(out)
}

colour_codes_prob_vector_fn <- function(prob_vector, probs=seq(0.25, 1, 0.25)){
  min_prob <- min(prob_vector)
  non_min_probs <- prob_vector[which(prob_vector>min_prob)]
  prob_quantiles <- quantile(non_min_probs, probs)
  threshold_vector <- c(min_prob, prob_quantiles)
  out <- sapply(prob_vector, breaks_prob_fn, threshold_vector=threshold_vector)
  return(out)
}

calculate_auc_plus_boyce_index <- function(GLM_train_test_object){
  GLM_test_pred <- predict(GLM_train_test_object$training.glm,
                                  newdata = GLM_train_test_object$test.df,
                                  type = "response")
  corrected_probs <- correct_prob_vector(GLM_test_pred, 
                                         GLM_train_test_object$training.glm)
  plotID <- 1:length(corrected_probs)
  df <- data.frame(plotID, obs = GLM_train_test_object$test.df$Y, 
                   pred = corrected_probs)
  area_under_curve <- auc(df, st.dev=FALSE)
  boyce_index <- boyce_calc(df$pred, df$obs)
  out <- list(AUC = area_under_curve, boyce_index = boyce_index)
  return(out)
}  
  
extract_coefs_sds_GLM <- function(GLM_test_train_object){
  GLM <- GLM_test_train_object$training.glm
  coefs <- summary(GLM)$coefficients[,1]
  sds <- summary(GLM)$coefficients[,2]
  out <- list(coefs = coefs, sds = sds)
  return(out)
}

index_GLM_coefs_rownames <- function(extracted_coefs_sds){
  index <- which(rownames(GLM_coefs) %in% names(extracted_coefs_sds$coefs))
  out <- index
  return(out)
}

plot_monthly_coefs_all_CORINE <- function(coef_df, 
                                          Month_string, 
                                          filename,
                                          ylim = c(-0.25, 0.25)){
  df <- coef_df
  
  test1 <- df %>% filter(Month==Month_string) %>%
    filter(coef!=0) %>%
    filter(Year>=2003 & Year<=2008) %>%
    filter(variable!="(Intercept)")
  
  test2 <- df%>% filter(Month==Month_string) %>%
    filter(coef!=0) %>%
    filter(Year >=2009 & Year<=2014) %>%
    filter(variable!="(Intercept)")
  
  test3 <- df%>% filter(Month==Month_string) %>%
    filter(coef!=0) %>%
    filter(Year >=2015 & Year <= 2019) %>%
    filter(variable!="(Intercept)")
  
  p1 <- ggplot(test1, aes(x=label, y=coef)) +
    geom_point(size=0.8) + 
    coord_flip() +
    facet_wrap(~Year, ncol=6) +
    geom_errorbar(aes(ymin = coef-sd, ymax = coef + sd),
                  width=0.5) +
    ylim(ylim) +
    ggtitle(paste0("CORINE2006 ", Month_string, " coefficients"))
  
if (!(Month_string %in% c("Sep", "Oct"))){
    
  p2 <- ggplot(test2, aes(x=label, y=coef)) +
    geom_point(size=0.8) + 
    coord_flip() +
    facet_wrap(~Year, ncol=6) +
    geom_errorbar(aes(ymin = coef-sd, ymax = coef + sd),
                  width=0.5) +
    ylim(ylim) +
    ggtitle(paste0("CORINE2012 ", Month_string, " coefficients"))
  
  p3 <- ggplot(test3, aes(x=label, y=coef)) +
    geom_point(size=0.8) + 
    coord_flip() +
    facet_wrap(~Year, ncol=6) +
    geom_errorbar(aes(ymin = coef-sd, ymax = coef + sd),
                  width=0.5) +
   ylim(ylim) +
    ggtitle(paste0("CORINE2018 ", Month_string, " coefficients"))
  
pdf(filename, width=10)
grid.arrange(grobs= list(p1, p2, p3), ncol=1)
dev.off()
}
  if (Month_string %in% c("Sep", "Oct")){
    pdf(filename, width=10)
    p1
    dev.off()
  }
}

plot_monthly_coefs_individual_CORINE_windows <- function(coef_df, 
                                          Month_string, 
                                          filename,
                                          ylim = c(-0.25, 0.25)){
  df <- coef_df
  
  test1 <- df %>% filter(Month==Month_string) %>%
    filter(coef!=0) %>%
    filter(Year>=2003 & Year<=2008) %>%
    filter(variable!="(Intercept)")
  
  test2 <- df%>% filter(Month==Month_string) %>%
    filter(coef!=0) %>%
    filter(Year >=2009 & Year<=2014) %>%
    filter(variable!="(Intercept)")
  
  test3 <- df%>% filter(Month==Month_string) %>%
    filter(coef!=0) %>%
    filter(Year >=2015 & Year <= 2019) %>%
    filter(variable!="(Intercept)")
  
  p1 <- ggplot(test1, aes(x=label, y=coef)) +
    geom_point(size=0.8) + 
    coord_flip(ylim=ylim) +
    facet_wrap(~Year, ncol=6) +
    geom_errorbar(aes(ymin = coef-sd, ymax = coef + sd),
                  width=0.5) +
    ggtitle(paste0("CORINE2006 ", Month_string, " coefficients"))
  

    p2 <- ggplot(test2, aes(x=label, y=coef)) +
      geom_point(size=0.8) + 
      coord_flip(ylim=ylim) +
      facet_wrap(~Year, ncol=6) +
      geom_errorbar(aes(ymin = coef-sd, ymax = coef + sd),
                    width=0.5) +
      ggtitle(paste0("CORINE2012 ", Month_string, " coefficients"))
    
    p3 <- ggplot(test3, aes(x=label, y=coef)) +
      geom_point(size=0.8) + 
      coord_flip(ylim=ylim) +
      facet_wrap(~Year, ncol=6) +
      geom_errorbar(aes(ymin = coef-sd, ymax = coef + sd),
                    width=0.5) +
      ggtitle(paste0("CORINE2018 ", Month_string, " coefficients"))
    
    pdf(filename, width=10)
    print(p1)
    print(p2)
    print(p3)
    dev.off()

}

plot_monthly_coefs_separate_pdfs_CORINE_windows <- function(coef_df, 
                                                         Month_string, 
                                                         filename2006,
                                                         filename2012,
                                                         filename2018,
                                                         ylim = c(-0.25, 0.25)){
  df <- coef_df
  
  test1 <- df %>% filter(Month==Month_string) %>%
    filter(coef!=0) %>%
    filter(Year>=2003 & Year<=2008) %>%
    filter(variable!="(Intercept)")
  
  test2 <- df%>% filter(Month==Month_string) %>%
    filter(coef!=0) %>%
    filter(Year >=2009 & Year<=2014) %>%
    filter(variable!="(Intercept)")
  
  test3 <- df%>% filter(Month==Month_string) %>%
    filter(coef!=0) %>%
    filter(Year >=2015 & Year <= 2019) %>%
    filter(variable!="(Intercept)")
  
  p1 <- ggplot(test1, aes(x=label, y=coef)) +
    geom_point(size=0.8) + 
    coord_flip(ylim=ylim) +
    facet_wrap(~Year, ncol=6) +
    geom_errorbar(aes(ymin = coef-sd, ymax = coef + sd),
                  width=0.5) +
    ggtitle(paste0("CORINE2006 ", Month_string, " coefficients"))
  
  
  p2 <- ggplot(test2, aes(x=label, y=coef)) +
    geom_point(size=0.8) + 
    coord_flip(ylim=ylim) +
    facet_wrap(~Year, ncol=6) +
    geom_errorbar(aes(ymin = coef-sd, ymax = coef + sd),
                  width=0.5) +
    ggtitle(paste0("CORINE2012 ", Month_string, " coefficients"))
  
  p3 <- ggplot(test3, aes(x=label, y=coef)) +
    geom_point(size=0.8) + 
    coord_flip(ylim=ylim) +
    facet_wrap(~Year, ncol=6) +
    geom_errorbar(aes(ymin = coef-sd, ymax = coef + sd),
                  width=0.5) +
    ggtitle(paste0("CORINE2018 ", Month_string, " coefficients"))
  
  pdf(filename2006, width=10)
  print(p1)
  dev.off()
  
  pdf(filename2012, width=10)
  print(p2)
  dev.off()
  
  pdf(filename2018, width=10)
  print(p3)
  dev.off()
  
}


exp_vars_raster_for_GLM_train_test_probs <- function(Year, 
                                                     Month_numeric, 
                                                     exp_and_response_vars_filestem){
  corine_to_use <- corine_to_use_fn(Year)
  Month <- Month_strings[Month_numeric]
  corine <- get(paste0("corine_", corine_to_use))
  exp_vars_filename <- paste0(exp_and_response_vars_filestem,
                              "exp_vars_CORINE", corine_to_use,
                              "_",Month,"s.tif")
  response_vars_filename <- paste0(exp_and_response_vars_filestem,
                                   "response_vars_CORINE", corine_to_use,
                                   "_",Month,"s.tif")
  exp_vars_stack <- stack(exp_vars_filename)
  response_vars_stack <- stack(response_vars_filename)
  N_CORINE <- dim(corine)[3]
  
  N_exp_vars <- dim(exp_vars_stack)[3]
  names(exp_vars_stack)[1:N_CORINE] <- gsub(paste0("_",corine_to_use),"",names(corine))
  window <- corine_Year.Month_vector(corine_to_use, Month_numeric)
  L <- length(window)
  names(exp_vars_stack)[N_CORINE+1:L] <- paste0("X101.", 1:L)
  names(exp_vars_stack)[N_CORINE+L+1:L] <- paste0("X102.", 1:L)
  names(exp_vars_stack)[N_CORINE+2*L+1:L] <- paste0("X103.",1:L)
  names(exp_vars_stack)[N_CORINE+3*L+1:L] <- paste0("X104.",1:L)
  names(exp_vars_stack)[N_CORINE+4*L+1:L] <- paste0("X105.",1:L)
  names(exp_vars_stack)[(N_exp_vars-1):N_exp_vars] <- c("X201", "X202")
  
  find_index <- function(window, Year){
    L <- length(window)
    temp <- numeric(L)
    for (i in 1:L){
      temp[i] <- as.numeric(strsplit(window[i],"[.]")[[1]][1])
    }
    index <- which(temp==Year)
    return(index)
  }
  
  i <- find_index(window, Year)
  
  exp_vars_temp <- exp_vars_stack[[c(1:N_CORINE,
                                     N_CORINE + i,
                                     N_CORINE+L+i,
                                     N_CORINE+2*L+i,
                                     N_CORINE+3*L+i,
                                     N_CORINE+4*L+i,
                                     (N_exp_vars-1):N_exp_vars)]]
  names(exp_vars_temp) <- c(names(exp_vars_stack)[1:N_CORINE],
                            paste0("X10",1:5),
                            c("X201", "X202"))
  #    predictions <- predict(exp_vars_and_response, monthly.glm, type="response")
  out <- exp_vars_temp
  return(out)
}

probs_corrected_raster <- function(Year, Month_numeric,
            filenamebase = "data/cleaned_curlew_obs_and_exp_vars_CORINE/"){
  exp_vars_raster <- exp_vars_raster_for_GLM_train_test_probs(Year, Month_numeric,
  filenamebase = filenamebase)
  GLMs_train_test_filenamebase <- paste0(filenamebase, "GLMs_by_month_train_test/")
  GLM_train_test_filename <- paste0(GLMs_train_test_filenamebase,
                         "GLM_", Year, ".", Month_numeric,
                         "_test_train.rds")
  GLM_train_test_obj <- readRDS(GLM_train_test_filename)
  GLM_train_test <- GLM_train_test_obj$training.glm
    probs_uncorrected <- predict(exp_vars_raster, 
                               GLM_train_test,
                               type="response")
    probs_corrected <- correct_prob_raster(probs_uncorrected,
                                           GLM_train_test)
out <- probs_corrected
return(out)
}


probs_uncorrected_raster <- function(Year, Month_numeric,
                                   exp_and_response_vars_filestem,
                                   GLM_train_test_filestem){
  exp_vars_raster <- exp_vars_raster_for_GLM_train_test_probs(Year, Month_numeric,
                    exp_and_response_vars_filestem = exp_and_response_vars_filestem)
  GLM_train_test_filename <- paste0(GLM_train_test_filestem,
                                    "GLM_", Year, ".", Month_numeric,
                                    "_train_test.rds")
  GLM_train_test_obj <- readRDS(GLM_train_test_filename)
  GLM_train_test <- GLM_train_test_obj$training.glm
  probs_uncorrected <- predict(exp_vars_raster, 
                               GLM_train_test,
                               type="response")
  out <- probs_uncorrected
  return(out)
}

cut_probs_raster <- function(probs_raster,
                            quantiles = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)){
  A <- values(probs_raster)
  M <- min(A, na.rm=TRUE)
  A <- A[A > M]
  breaks <- quantile(A, probs=quantiles, na.rm=TRUE)
  breaks[1] <- 0
  out <- cut(probs_raster, breaks=breaks)
  return(out)
}


model_prob_classes_raster <- function(probs_raster, 
            quantiles = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)){
prob_classes_raster <- cut_probs_raster(probs_raster, quantiles=quantiles)
prob_classes_raster <- as.factor(prob_classes_raster)
L <- levels(prob_classes_raster)[[1]]
L[["Model prob class"]] <- c("Less than 25th percentile","25th-50th percentile", "50th to 75th percentile","75th to 90th percentile", "90th to 95th percentile", "95th to 99th percentile", "Above 99th percentile")
levels(prob_classes_raster) <- L
out <- prob_classes_raster
return(out)
}

model_sd_prob_classes_raster <- function(sd_probs_raster, 
                                      quantiles = c(0, 0.2, 0.4, 0.6, 0.8, 1)){
  sd_probs_classes_raster <- cut_probs_raster(sd_probs_raster, 
                                              quantiles=quantiles)
  sd_probs_classes_raster <- as.factor(sd_probs_classes_raster)
  L <- levels(sd_probs_classes_raster)[[1]]
  L[["Model prob class"]] <- c("Less than 20th percentile",
                               "20th to 40th percentile",
                               "40th to 60th percentile",
                               "60th to 80th percentile",
                               "Above 80th percentile")
                               
  levels(sd_probs_classes_raster) <- L
  out <- sd_probs_classes_raster
  return(out)
}


levelplot_prob_classes_raster <- function(Year, Month_numeric,
                                          filenamebase = "data/cleaned_curlew_obs_and_exp_vars_CORINE/GLM_prob_rasters/probabilities_rasters/"){
  Month_string <- Month_strings[Month_numeric]
  filename <- paste0(filenamebase, Month_string, "_", Year, "_probs.rds")
  probs_raster <- readRDS(filename)
  probs_classes_raster <- model_prob_classes_raster(probs_raster)
  
  out <- levelplot(probs_classes_raster, col.regions=cbSeq, 
            scales=list(draw=FALSE), xlab="", ylab="",
            colorkey = list(height=0.6),
            main=paste0(Month_string, " ", Year,
            " fitted probs"))
return(out)
}

levelplot_mean_prob_raster <- function(mean_probs_raster, Month_string,
                                       col.regions){
  probs_classes_raster <- model_prob_classes_raster(mean_probs_raster)
  
  out <- levelplot(probs_classes_raster, col.regions=col.regions, 
                   scales=list(draw=FALSE), xlab="", ylab="",
                   colorkey = FALSE,
                   main=paste0(Month_string, " (means) "))
  return(out)
}

levelplot_sd_prob_raster <- function(sd_probs_raster, Month_string,
                                     col.regions){
  probs_classes_raster <- model_sd_prob_classes_raster(sd_probs_raster)
  
  out <- levelplot(probs_classes_raster, col.regions=col.regions, 
                   scales=list(draw=FALSE), xlab="", ylab="",
                   colorkey = FALSE,
                   main=paste0(Month_string, " (sds) "))
  return(out)
}

addColorTable <- function(inRstName, outRstName, rat.df){
  library(rgdal)
  r<- readGDAL(inRstName)
  rat.df$color<- as.character(rat.df$color)
  rat.df$attribute<- as.character(rat.df$attribute)
  outRst <- writeGDAL(r, outRstName, type="Byte", 
                      colorTable=list(rat.df$color), 
                      catNames=list(rat.df$attribute), mvFlag=11L)
  return(raster(outRst))
}


  
write_tiff <- function(Year, Month_numeric,
                       GLM_prob_rasters_filenamebase,
                       tiff_filenamebase,
                       quantiles = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1)){
  Month_code <- Month_codes[Month_numeric]
  Day_code <- "01"
  filename <- paste0(GLM_prob_rasters_filenamebase, 
                     Month_strings[Month_numeric],
                     "_", Year, "_probs.rds")
  probs_raster <- readRDS(filename)
  probs_classes_raster <- model_prob_classes_raster(probs_raster)
  r <- ratify(probs_classes_raster)
  levels(r)[[1]]$NAME <- c("Less than 20th percentile","20th-40th percentile", "40th-60th percentile","60th-80th percentile", "80th to 90th percentile", "90th to 95th percentile", "95th to 99th percentile", "Above 99th percentile")
  # check results with levels(r)[[1]]
  
    tiff_filename <- paste0(tiff_filenamebase,
                          "prob_surface_",
                          Year, Month_code, Day_code,
                          ".tif")
  writeRaster(r,
              tiff_filename, 
              options=c('TFW=YES'), overwrite=TRUE)
  # This defines the values, the color and the attribute
#  valT <- 0:7
#  colT <-  cbSeq
#  attT <- c("Less than 20th percentile","20th-40th percentile", "40th-60th percentile","60th-80th percentile", "80th to 90th percentile", "90th to 95th percentile", "95th to 99th percentile", "Above 99th percentile")
#  rat.df <- data.frame(value=valT,color=colT,attribute=attT)

#  tiff_symbology_filename <- paste0(tiff_filenamebase,
#                          "prob_surface_symbology_",
#                          Year, Month_code, Day_code,
#                          ".tif")
  
    
  # apply the magic function
#  rnew <- addColorTable(tiff_filename, tiff_symbology_filename, rat.df)
#  out <- rnew
#  return(out)
}

means_and_sds_prob_rasters <- function(probs_rasterstack,
                                       Months_in_study){
N_months <- length(Months_in_study)
raster_list_mean_probs_out <- vector("list", N_months)
raster_list_sd_probs_out <- vector("list", N_months)
index <- 1:N_months
    for (i in index){
      assign(paste0("index",i), grep(Months_in_study[i], 
                                     names(probs_rasterstack)))
    rasterstack_sub <- probs_rasterstack[[get(paste0("index",i))]]  
    A <- calc(rasterstack_sub, mean)
    B <- calc(rasterstack_sub, sd)
    raster_list_mean_probs_out[[i]] <- A
    raster_list_sd_probs_out[[i]] <- B
    }
  out <- list(mean_probs_rasters = raster_list_mean_probs_out,
              sd_probs_rasters = raster_list_sd_probs_out)
return(out)
}

save_means_and_sds_probs_rasters <- function(corine_to_use, Months_in_study,
                                            GLM_prob_rasters_filestem){
  if (corine_to_use == 2006) years <- 2003:2008
  if (corine_to_use == 2012) years <- 2009:2014
  if (corine_to_use == 2018) years <- 2015:2019
  Month_strings <- c("Jan", "Feb", "Mar",
                     "Apr", "May", "Jun",
                     "Jul", "Aug", "Sep",
                     "Oct", "Nov", "Dec")
  N_years <- length(years)
  N_months_in_study <- length(Months_in_study)
  Years_Months_df <- filter(Years_Months_df, year %in% years)
  probs_rasters_filenames_short <- paste0(Month_strings[Years_Months_df$month],
                                          "_", 
                                          Years_Months_df$year,
                                          "_probs.rds")
  probs_rasters_filenames_long <- paste0(GLM_prob_rasters_filestem, 
                                         probs_rasters_filenames_short)
  Month_years <- paste0(Month_strings[Years_Months_df$month],
                        "_",
                        Years_Months_df$year)
  for (i in 1:length(probs_rasters_filenames_long)){
    filename <- probs_rasters_filenames_long[i]
    assign(paste0("r", i), readRDS(filename))
  }
  raster_list <- mget(paste0("r", 1:length(probs_rasters_filenames_long)))
  probs_rasterstack <- stack(raster_list)
  names(probs_rasterstack) <- Month_years
  means_and_sds_probs_rasters <- means_and_sds_prob_rasters(probs_rasterstack,
                                                            Months_in_study)
  output_file <- paste0(GLM_prob_rasters_filestem, "CORINE", corine_to_use, 
                        "_means_and_sds_probs.rds")
  saveRDS(means_and_sds_probs_rasters,
          output_file)
}

plot_means_and_sds_prob_rasters <- function(corine_to_use,
                                            GLM_prob_rasters_filestem,
                                            Months_in_study){
  N_months <- length(Months_in_study)
  means_and_sds_probs_filename <- paste0(GLM_prob_rasters_filestem, "CORINE", 
                                         corine_to_use, "_means_and_sd_probs.rds")
  means_and_sds_prob_rasters <- readRDS(means_and_sds_probs_filename)
  means_sds <- means_and_sds_prob_rasters
  means <- means_sds$mean_probs_rasters
  sds <- means_sds$sd_probs_rasters
  for(i in 1:N_months){
    Month_string <- Month_strings[i]
    assign(paste0("p", i, "a"), levelplot_mean_prob_raster(means[[i]], Month_string))
    assign(paste0("p", i, "b"), levelplot_sd_prob_raster(sds[[i]], Month_string))
  }
  means_plots <- mget(paste0("p", 1:N_months, "a"))
  sds_plots <- mget(paste0("p", 1:N_months, "b"))
  overall_plots <- c(means_plots, sds_plots)
  means_plots_indices <- 1:N_months
  sds_plots_indices <- (N_months + 1):(2*N_months)
  means_sds_plots_indices_interlaced <- numeric(2*N_months)
  means_sds_plots_indices_interlaced[seq(1, (2*N_months)-1, by=2)] <- 
    months_plots_indices
  means_sds_plots_indices_interlaced[seq(2, (2*N_months), by=2)] <-
    sds_plots_indices
  output_file <- paste("figs/CORINE", corine_to_use, "means_and_sds_plots.pdf")
  pdf(output_file, width=10)
  grid.arrange(grobs=overall_plots[means_sds_plots_indices_interlaced], 
               ncol=4, 
               top=textGrob(paste0("CORINE", corine_to_use, " window - means and sds of fitted probabilities"), 
                            gp=gpar(fontsize=15, fontface="bold")))
  dev.off()
}  

x_axis_fn <- function(Month_numeric){
  if (Month_numeric == 1 | Month_numeric == 2) out <- Month_numeric
  if (Month_numeric ==11) out <- 3
  if (Month_numeric == 12) out <- 4
  return(out)
}


save_probs_raster <- function(index,
                              exp_and_response_vars_filestem,
                              GLM_train_test_filestem,
                              GLM_prob_rasters_filestem){
  Year <- Years[index]
  Month_numeric <- Months[index]
  probs_uncorrected <- probs_uncorrected_raster(Year, Month_numeric,
                                                exp_and_response_vars_filestem,
                                                GLM_train_test_filestem)
  filename <- paste0(GLM_prob_rasters_filestem, 
                     Month_strings[Month_numeric],
                     "_", Year, "_probs.rds")
  saveRDS(probs_uncorrected, filename)
}


plot_means_sds_prob_rasters <- function(corine_to_use,
            Months_in_study,
            GLM_prob_rasters_filestem,
            monthly_means_and_sds_filestem,
            col.regions_means,
            col.regions_sds,
            use_title = TRUE){
  if (corine_to_use == 2006) years_CORINE <- 2003:2008
  if (corine_to_use == 2012) years_CORINE <- 2009:2014
  if (corine_to_use == 2018) years_CORINE <- 2015:2019
  N_years <- length(years_CORINE)
  N_months <- length(Months_in_study)
  probs_rasters_files_short <- paste0(rep(Months_in_study, N_years),
                                          "_", rep(years_CORINE, 
                                                   rep(N_months, N_years)), 
                                          "_probs.rds")
  probs_rasters_files_long <- paste0(GLM_prob_rasters_filestem,
                                     probs_rasters_files_short)
  Month_years <- paste0(rep(Months_in_study, N_years),
                        "_", rep(years_CORINE, rep(N_months, N_years)))
  L <- length(probs_rasters_files_short)
  # The following is because we only have January and February
  # for 2019 (no November or December data)
  if(corine_to_use == 2018){
    probs_rasters_files_short <- probs_rasters_files_short[-c(L-1, L)]
    probs_rasters_files_long <- probs_rasters_files_long[-c(L-1, L)]
    Month_years <- Month_years[-c(L-1, L)]
  }

for (i in 1:length(probs_rasters_files_long)){
  filename <- probs_rasters_files_long[i]
  assign(paste0("r",i), readRDS(filename))
}

raster_list <- mget(paste0("r",1:length(probs_rasters_files_long)))

probs_rasterstack <- stack(raster_list)
names(probs_rasterstack) <- Month_years

means_sds_rasters <- means_and_sds_prob_rasters(probs_rasterstack,
                                                            Months_in_study)

output_file <- paste0(monthly_means_and_sds_filestem, "CORINE", corine_to_use,
                      "_means_and_sd_probs.rds")

mns_sds <- means_sds_rasters
mns <- mns_sds$mean_probs_rasters
sds <- mns_sds$sd_probs_rasters
for(i in 1:N_months){
  Month_string <- Months_in_study[i]
  assign(paste0("p", i, "a"), levelplot_mean_prob_raster(mns[[i]], 
                                            Month_string,
                                            col.regions = col.regions_means))
  assign(paste0("p", i, "b"), levelplot_sd_prob_raster(sds[[i]], 
                                            Month_string,
                                            col.regions=col.regions_sds))
}
means_plots <- mget(paste0("p", 1:N_months, "a"))
sds_plots <- mget(paste0("p", 1:N_months, "b"))
overall_plots <- c(means_plots,
                         sds_plots)

figs_output_file <- paste0("figs/CORINE", corine_to_use, 
                           "_means_and_sds_plots.pdf")
title <- paste0("CORINE", corine_to_use, 
                " window - means and sds of fitted probabilities")
# We plot November, December, January, February in that order

# With title bar
if (use_title==TRUE){
p <- arrangeGrob(grobs=overall_plots[c(3, 3+N_months, 4, 4+N_months,
                                        1, 1+N_months, 2, 2+N_months)], 
                  ncol=4,
                  top=textGrob(title, 
                               gp=gpar(fontsize=15, fontface="bold")))
}
# Without title bar
if (use_title==FALSE){
p <- arrangeGrob(grobs=overall_plots[c(3, 3+N_months, 4, 4+N_months,
                                       1, 1+N_months, 2, 2+N_months)], 
                 ncol=4)
}
ggsave(figs_output_file, p, width=10)
}

coef_matrix_5seeds <- function(corine_to_use, Months_in_study,
                               corine_to_use_rep,
                               Month_numeric_vector,
                               Month_strings,
                               cv.glmnet_filestem){
  if(corine_to_use == 2006){
    N_rows <- 45
    index <- 1:length(Months_in_study)
  }
  if(corine_to_use == 2012){
    N_rows <- 43
    index <- length(Months_in_study) + 1:(length(Months_in_study))
  }
  if(corine_to_use == 2018){
    N_rows <- 43
    index <- (2*length(Months_in_study) + (1:length(Months_in_study)))
  }
  coef_matrix_output_5seeds <- matrix(0, N_rows, length(Months_in_study)*5)
count <- 0
    for (i in index){
    count <- count + 1
    corine_to_use <- corine_to_use_rep[i]
    Month_numeric <- Month_numeric_vector[i]
    Month <- Month_strings[Month_numeric]
    cv.glmnet_filename <- paste0(cv.glmnet_filestem,"CORINE",corine_to_use,"_",
                                 Month,"s_cv.glmnet.rds")
    cv.glmnet_filename_seed111022 <- paste0(cv.glmnet_filestem, "seed_111022/",
                                            "CORINE", corine_to_use, "_",
                                            Month,"s_cv.glmnet.rds")
    cv.glmnet_filename_seed121022 <- paste0(cv.glmnet_filestem, "seed_121022/",
                                            "CORINE", corine_to_use, "_",
                                            Month,"s_cv.glmnet.rds")
    cv.glmnet_filename_seed131022 <- paste0(cv.glmnet_filestem, "seed_131022/",
                                            "CORINE", corine_to_use, "_",
                                            Month,"s_cv.glmnet.rds")
    cv.glmnet_filename_seed141022 <- paste0(cv.glmnet_filestem, "seed_141022/",
                                            "CORINE", corine_to_use, "_",
                                            Month,"s_cv.glmnet.rds")
    temp1 <- readRDS(cv.glmnet_filename)
    A <- as.matrix(coef(temp1, s="lambda.1se"))
    temp2 <- readRDS(cv.glmnet_filename_seed111022)
    B <- as.matrix(coef(temp2, s="lambda.1se"))
    temp3 <- readRDS(cv.glmnet_filename_seed121022)
    C <- as.matrix(coef(temp3, s="lambda.1se"))
    temp4 <- readRDS(cv.glmnet_filename_seed131022)
    D <- as.matrix(coef(temp2, s="lambda.1se"))
    temp5 <- readRDS(cv.glmnet_filename_seed141022)
    E <- as.matrix(coef(temp3, s="lambda.1se"))
    coef_matrix_output_5seeds[, count] <- A[,1]
    coef_matrix_output_5seeds[, count + length(Months_in_study)] <- B[,1]
    coef_matrix_output_5seeds[, count + 2*length(Months_in_study)] <- C[,1]
    coef_matrix_output_5seeds[, count + 3*length(Months_in_study)] <- D[,1]
    coef_matrix_output_5seeds[, count + 4*length(Months_in_study)] <- E[,1]
  }
rownames(coef_matrix_output_5seeds) <- rownames(A)
random_seeds <- c(101022, 111022, 121022, 131022, 141022)
colnames(coef_matrix_output_5seeds) <- paste0(rep(Months_in_study, 5), "_seed", 
                                      rep(random_seeds, rep(length(Months_in_study), 5)))
out <- coef_matrix_output_5seeds
return(out)
}

coef_matrix_multiple_seeds <- function(corine_to_use, Months_in_study,
                               corine_to_use_rep,
                               Month_numeric_vector,
                               Month_strings,
                               cv.glmnet_filestem,
                               random_seeds){
  N_random_seeds <- length(random_seeds)
  if(corine_to_use == 2006){
    N_rows <- 45
    index <- 1:length(Months_in_study)
  }
  if(corine_to_use == 2012){
    N_rows <- 43
    index <- length(Months_in_study) + 1:(length(Months_in_study))
  }
  if(corine_to_use == 2018){
    N_rows <- 43
    index <- (2*length(Months_in_study) + (1:length(Months_in_study)))
  }
  coef_matrix_output_multiple_seeds <- matrix(0, N_rows, 
                                              length(Months_in_study)*N_random_seeds)
  count <- 0
  for (i in index){
    count <- count + 1
    corine_to_use <- corine_to_use_rep[i]
    Month_numeric <- Month_numeric_vector[i]
    Month <- Month_strings[Month_numeric]
    for (j in 1:N_random_seeds){
      random_seed <- random_seeds[j]
      cv.glmnet_filestem_random_seed <- paste0(cv.glmnet_filestem, "seed_", random_seed, 
                                               "/CORINE", corine_to_use, "_",
                                               Month, "s_cv.glmnet.rds")
      temp <- readRDS(cv.glmnet_filestem_random_seed)
      A <- as.matrix(coef(temp, s="lambda.1se"))
      coef_matrix_output_multiple_seeds[, count + (j-1)*length(Months_in_study)] <- A[,1]
    }
  }
    rownames(coef_matrix_output_multiple_seeds) <- rownames(A)
    colnames(coef_matrix_output_multiple_seeds) <- paste0(rep(Months_in_study, 
                                                              N_random_seeds), 
                                                "_seed", 
                        rep(random_seeds, rep(length(Months_in_study), N_random_seeds)))
    out <- coef_matrix_output_multiple_seeds
    return(out)
  }    
    
    


generate_seed <- function(seed_string){
  out <- strsplit(as.character(seed_string), "seed")[[1]][2]
  return(out)
}


plot_coefs_multiple_seeds <- function(Month_of_interest,
            coef_matrix_2006_multiple_seeds,
            coef_matrix_2012_multiple_seeds,
            coef_matrix_2018_multiple_seeds,
            random_seeds,
            clc_labels,
            plot_output_file,
            colours = c("#8C510A", "#F5F5F5", "#01665E")){
            #colours = c("#FFF5F0", "#FB6A4A", "#67000D")){
            #colours = c("#075AFF", "#FFFFCC", "#FF0000")){
  N_random_seeds <- length(random_seeds)
  col_index_2006 <- grep(Month_of_interest, 
                         colnames(coef_matrix_2006_multiple_seeds))
  col_index_2012 <- grep(Month_of_interest,
                         colnames(coef_matrix_2012_multiple_seeds))
  col_index_2018 <- grep(Month_of_interest,
                         colnames(coef_matrix_2018_multiple_seeds))
  # Remove (Intercept) here
  coef_matrix_2006 <- coef_matrix_2006_multiple_seeds[2:dim(coef_matrix_2006_multiple_seeds)[1] ,col_index_2006]
  coef_matrix_2012 <- coef_matrix_2012_multiple_seeds[2:dim(coef_matrix_2012_multiple_seeds)[1], col_index_2012]
  coef_matrix_2018 <- coef_matrix_2018_multiple_seeds[2:dim(coef_matrix_2018_multiple_seeds)[1], col_index_2018]
  overall_set_of_rownames <- union(union(rownames(coef_matrix_2006),
                                         rownames(coef_matrix_2012)),
                                   rownames(coef_matrix_2018))
  vars_stripX <- sapply(overall_set_of_rownames, stripX)
  new_order <- sort(vars_stripX, index.return=TRUE)$ix
  overall_set_of_rownames_ordered <- overall_set_of_rownames[new_order]
  overall_coefs_matrix_2006_12_18 <- matrix(0, length(overall_set_of_rownames_ordered), 
                                            3*N_random_seeds)
  colnames(overall_coefs_matrix_2006_12_18) <- paste0(Month_of_interest, "s_CORINE", 
                                                     rep(c(2006, 2012, 2018), 
                                                         rep(N_random_seeds, 3)), 
                                                     "_seed", rep(random_seeds, 3))
  rownames(overall_coefs_matrix_2006_12_18) <- overall_set_of_rownames_ordered
  index_2006 <- find_match_index(rownames(coef_matrix_2006), 
                                 overall_set_of_rownames_ordered)
  index_2012 <- find_match_index(rownames(coef_matrix_2012), 
                                 overall_set_of_rownames_ordered)
  index_2018 <- find_match_index(rownames(coef_matrix_2018), 
                                 overall_set_of_rownames_ordered)
  X101_index_2006 <- which(rownames(coef_matrix_2006)=="X101")
  X101_index_2012 <- which(rownames(coef_matrix_2012)=="X101")
  X101_index_2018 <- which(rownames(coef_matrix_2018)=="X101")

  
  overall_coefs_matrix_2006_12_18[index_2006,1:N_random_seeds] <- coef_matrix_2006[1:(X101_index_2006-1),]
  overall_coefs_matrix_2006_12_18[index_2012,N_random_seeds+(1:N_random_seeds)] <- coef_matrix_2012[1:(X101_index_2012-1),]
  overall_coefs_matrix_2006_12_18[index_2018,2*N_random_seeds+(1:N_random_seeds)] <- coef_matrix_2018[1:(X101_index_2018-1),]
  overall_coefs_matrix_2006_12_18 <- cbind(overall_coefs_matrix_2006_12_18, 0)
  overall_coefs_matrix_2006_12_18[,3*N_random_seeds+1] <- apply(overall_coefs_matrix_2006_12_18[,1:(3*N_random_seeds)], 
                                                          1, x_all_zero)
  colnames(overall_coefs_matrix_2006_12_18)[3*N_random_seeds + 1] <- "zero_row"
  df <- data.frame(overall_coefs_matrix_2006_12_18)
  filtered_df <- filter(df, zero_row==0)
  X101_index <- which(rownames(filtered_df)=="X101")
  overall_coef_matrix_for_geom_tile <- sign(as.matrix(filtered_df[2:dim(filtered_df)[1],1:(3*N_random_seeds)]))
  A <- melt(overall_coef_matrix_for_geom_tile)
  A$seed <- sapply(A$X2, generate_seed)
  A$seed <- ordered(A$seed, levels=random_seeds)
  A$year <- sapply(A$X2, generate_year_long)
  vars_to_use <- sapply(rownames(overall_coef_matrix_for_geom_tile), stripX)
 
   p <- ggplot(A, aes(x = X1, y = seed, fill = factor(value))) +
   geom_tile(color="black") +
   scale_fill_manual(values=colours, labels=c("Negative", "Null", "Positive")) + 
   scale_x_discrete(labels=clc_labels[vars_to_use]) +
   coord_flip() +
   labs(x="CORINE land cover class", y="Random seed", fill= "Influence") +
   ggtitle(paste0("Influence of land cover classes on log-odds of curlew sightings for ", 
                  Month_of_interest))+
   facet_wrap(~year, nrow=2) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
   
   ggsave(filename = plot_output_file, plot = p, width=10)
   
   out <- overall_coef_matrix_for_geom_tile
   return(out)
}

generate_coef_sgns_inc_climate <- function(Month_of_interest,
                                      coef_matrix_2006_multiple_seeds,
                                      coef_matrix_2012_multiple_seeds,
                                      coef_matrix_2018_multiple_seeds,
                                      random_seeds){  
  N_random_seeds <- length(random_seeds)
  col_index_2006 <- grep(Month_of_interest, 
                         colnames(coef_matrix_2006_multiple_seeds))
  col_index_2012 <- grep(Month_of_interest,
                         colnames(coef_matrix_2012_multiple_seeds))
  col_index_2018 <- grep(Month_of_interest,
                         colnames(coef_matrix_2018_multiple_seeds))
  # remove (Intercept) row
  coef_matrix_2006 <- coef_matrix_2006_multiple_seeds[2:dim(coef_matrix_2006_multiple_seeds)[1], 
                                                      col_index_2006]
  coef_matrix_2012 <- coef_matrix_2012_multiple_seeds[2:dim(coef_matrix_2012_multiple_seeds)[1], 
                                                      col_index_2012]
  coef_matrix_2018 <- coef_matrix_2018_multiple_seeds[2:dim(coef_matrix_2018_multiple_seeds)[1], 
                                                      col_index_2018]
  overall_set_of_rownames <- union(union(rownames(coef_matrix_2006),
                                         rownames(coef_matrix_2012)),
                                   rownames(coef_matrix_2018))
  vars_stripX <- sapply(overall_set_of_rownames, stripX)
  new_order <- sort(vars_stripX, index.return=TRUE)$ix
  overall_set_of_rownames_ordered <- overall_set_of_rownames[new_order]
  overall_coefs_matrix_2006_12_18 <- matrix(0, length(overall_set_of_rownames_ordered), 
                                            3*N_random_seeds)
  colnames(overall_coefs_matrix_2006_12_18) <- paste0(Month_of_interest, "s_CORINE", 
                                                      rep(c(2006, 2012, 2018), 
                                                          rep(N_random_seeds, 3)), 
                                                      "_seed", rep(random_seeds, 3))
  rownames(overall_coefs_matrix_2006_12_18) <- overall_set_of_rownames_ordered
  index_2006 <- find_match_index(rownames(coef_matrix_2006), 
                                 overall_set_of_rownames_ordered,
                                 inc_climate=TRUE)
  index_2012 <- find_match_index(rownames(coef_matrix_2012), 
                                 overall_set_of_rownames_ordered,
                                 inc_climate=TRUE)
  index_2018 <- find_match_index(rownames(coef_matrix_2018), 
                                 overall_set_of_rownames_ordered,
                                 inc_climate=TRUE)

  
  overall_coefs_matrix_2006_12_18[index_2006,1:N_random_seeds] <- coef_matrix_2006  
  overall_coefs_matrix_2006_12_18[index_2012,N_random_seeds+(1:N_random_seeds)] <- coef_matrix_2012
  overall_coefs_matrix_2006_12_18[index_2018,2*N_random_seeds+(1:N_random_seeds)] <- coef_matrix_2018
  out <- sign(overall_coefs_matrix_2006_12_18)
  return(out)
}


plot_coef_sgns_matrix <- function(coef_sgns_matrix,
                                  Months_in_study,
                                  clc_labels,
                                  plot_output_file,
                                  colours = c("#8C510A", "#F5F5F5", "#01665E")){
                                  #colours = c("#FFF5F0", "#FB6A4A", "#67000D")){
                                  #colours = c("#075AFF", "#FFFFCC", "#FF0000")){
  # colnames in form Monthyear, e.g. Jan2006

A <- melt(coef_sgns_matrix)
A$month <- sapply(A$X2, generate_month)
A$month <- ordered(A$month, levels=Months_in_study[c(3, 4, 1, 2)])
A$year <- sapply(A$X2, generate_year)
vars_to_use <- sapply(rownames(coef_sgns_matrix), stripX)

p <- ggplot(A, aes(x = X1, y = month, fill = factor(value))) +
  geom_tile(color="black") +
  scale_fill_manual(values=colours, labels=c("Negative", "Null", "Positive")) + 
  scale_x_discrete(labels=clc_labels[vars_to_use]) +
  coord_flip() +
  labs(x="CORINE land cover class", y="Month", fill= "Influence") +
  ggtitle("Influence of land cover classes on log-odds of curlew sightings by month") +
  facet_wrap(~year, nrow=2) 

ggsave(plot_output_file, p, width=10)
}

generate_quantile_labels <- function(data_vector, N_quantiles = 4){
  splits <- (-N_quantiles):(N_quantiles)/N_quantiles
  labels <- character(length(splits))
  for (i in 1:N_quantiles){
    labels[i] <- paste0("[", splits[i], " to ", splits[i+1], ")")
  }
  for (i in (N_quantiles+2):(length(splits))){
    labels[i] <- paste0("(", splits[i-1], " to ", splits[i],"]")
  }
  labels[N_quantiles+1] <- "0"

  out <- numeric(length(data_vector))
    for (i in 1:length(data_vector)){
      if (data_vector[i] == 0) out[i] <- labels[N_quantiles+1]
      if (data_vector[i] == 1) out[i] <- labels[length(splits)]
      if (data_vector[i]!=0 & data_vector[i] < 1){
        index <- max(which(splits <= data_vector[i]))
        if (data_vector[i] < 0) out[i] <- labels[index]
        if (data_vector[i] > 0) out[i] <- labels[index+1]      
      }
    }
  out <- factor(out, levels=labels, ordered=TRUE)
  return(out)
}


plot_average_coef_sgns_matrix <- function(average_coef_sgns_matrix,
                                  Months_in_study,
                                  clc_labels,
                                  plot_output_file,
                                  colours = c("#8C510A", "#F5F5F5", "#01665E"),
                                  #colours = c("#FFF5F0", "#FB6A4A", "#67000D"),
                                  #colours = c("#0000FF", "#FFFFCC", "#FF0000"),
                                  N_quantiles=4,
                                  random_seeds){
  # colnames in form Monthyear, e.g. Jan2006
  
  N_random_seeds <- length(random_seeds)
  A <- melt(average_coef_sgns_matrix)
  A$month <- sapply(A$X2, generate_month)
  A$month <- ordered(A$month, levels=Months_in_study[c(3, 4, 1, 2)])
  A$year <- sapply(A$X2, generate_year)
  vars_to_use <- sapply(rownames(average_coef_sgns_matrix), stripX)
  A$quantile_lab <- generate_quantile_labels(A$value, N_quantiles=N_quantiles)
  colours_for_plot <- colorRampPalette(colours)(2*N_quantiles+1)
                        
  
  p <- ggplot(A, aes(x = X1, y = month, fill = quantile_lab)) +
    geom_tile(color="black") +
    scale_fill_manual(values=colours_for_plot) + 
    scale_x_discrete(labels=clc_labels[vars_to_use]) +
    coord_flip() +
    labs(x="CORINE land cover class", y="Month", fill= paste0("Average sign of coefs ", "(",N_random_seeds, " seeds)")) +
    ggtitle("Influence of land cover classes on log-odds of curlew sightings by month") +
    facet_wrap(~year, nrow=2) +
    guides(fill = guide_legend(reverse=TRUE))
  
  ggsave(plot_output_file, p, width=10, height=8)
}





choose_consensus_coefficient_sgn <- function(row_index, coef_matrix_Month_corine_5seeds){
  R <- coef_matrix_Month_corine_5seeds[row_index,]
  tab <- table(R)
  out <- as.numeric(names(tab)[which(tab == max(tab))])
  return(out)
}

sign_threshold_sum <- function(row_index, coef_sgn_matrix,
                                 threshold = 0.64){
  temp <- coef_sgn_matrix[row_index,]
  S_sgn <- sum(temp)
  prop_threshold <- abs(S_sgn)/length(temp)
  indicator <- (prop_threshold >= threshold)
  sign_out <- sign(S_sgn)*indicator
  out <- sign_out
  return(out)
}

consensus_coefficient_sgn_threshold <- function(coef_sgn_matrix_Month_corine_multipleseeds,
                                                threshold = 0.64){
  coef_sgns <- coef_sgn_matrix_Month_corine_multipleseeds
  N_row <- dim(coef_sgns)[1]
  sign_threshold <- sapply(1:N_row, sign_threshold_sum, coef_sgn_matrix=coef_sgns, 
                           threshold=threshold)
  out <- sign_threshold
  return(out)  
}

find_consensus_coefs <- function(coef_matrix_Month_2006_12_18_5seeds){
  A <- coef_matrix_Month_2006_12_18_5seeds
  N_row <- dim(coef_matrix_Month_2006_12_18_5seeds)[1]
  A_2006 <- A[,1:5]
  A_2012 <- A[,6:10]
  A_2018 <- A[,11:15]
  sgns_2006 <- sapply(1:N_row, choose_consensus_coefficient_sgn, 
                      coef_matrix_Month_corine_5seeds = A_2006)
  sgns_2012 <- sapply(1:N_row, choose_consensus_coefficient_sgn,
                      coef_matrix_Month_corine_5seeds = A_2012)
  sgns_2018 <- sapply(1:N_row, choose_consensus_coefficient_sgn,
                      coef_matrix_Month_corine_5seeds = A_2018)
  temp <- cbind(sgns_2006, sgns_2012, sgns_2018)
  rownames(temp) <- rownames(A)
  colnames(temp) <- c("2006", "2012", "2018")
  temp <- cbind(temp, 0)
  temp[,4] <- apply(temp[,1:3], 1, x_all_zero)
  colnames(temp)[4] <- "zero_row"
  df <- data.frame(temp)
  filtered_df <- filter(df, zero_row==0)
  out <- as.matrix(filtered_df[,1:3])
  return(out)
}

find_consensus_coef_sgns_threshold <- function(coef_sgn_matrix_Month_2006_12_18_multipleseeds,
                                           random_seeds,
                                           threshold=0.64){
  coef_sgns <- coef_sgn_matrix_Month_2006_12_18_multipleseeds
  N_random_seeds <-length(random_seeds)
  N_row <- dim(coef_sgns)[1]
  if(N_random_seeds*3 != ncol(coef_sgns)) stop("N_random_seeds*3 != ncol(coefs)")
  coef_sgns_2006 <- coef_sgns[,1:N_random_seeds]
  coef_sgns_2012 <- coef_sgns[, (N_random_seeds+1):(2*N_random_seeds)]
  coef_sgns_2018 <- coef_sgns[, ((2*N_random_seeds)+1):(3*N_random_seeds)]
  coef_sgns_threshold_2006 <- consensus_coefficient_sgn_threshold(coef_sgns_2006,
                                                                  threshold=threshold)
  coef_sgns_threshold_2012 <- consensus_coefficient_sgn_threshold(coef_sgns_2012,
                                                                  threshold=threshold)
  coef_sgns_threshold_2018 <- consensus_coefficient_sgn_threshold(coef_sgns_2018,
                                                                  threshold=threshold)
  temp <- cbind(coef_sgns_threshold_2006, coef_sgns_threshold_2012, coef_sgns_threshold_2018)
  rownames(temp) <- rownames(coef_sgns)
  colnames(temp) <- c("2006", "2012", "2018")
  temp <- cbind(temp, 0)
  temp[,4] <- apply(temp[,1:3], 1, x_all_zero)
  colnames(temp)[4] <- "zero_row"
  df <- data.frame(temp)
  filtered_df <- filter(df, zero_row==0)
  out <- as.matrix(filtered_df[,1:3])
  return(out)
}

  

month_average_sgns <- function(month_sgn_coefs_matrix_multiple_seeds,
                         corine_years = c(2006, 2012, 2018)){
N_CORINE <- length(corine_years)
coefs <- month_sgn_coefs_matrix_multiple_seeds
  out <- matrix(0, nrow(coefs), length(corine_years))
rownames(out) <- rownames(coefs)
colnames(out) <- paste0("CORINE", corine_years)
for (i in 1:N_CORINE){
  corine_to_use <- corine_years[i]
  index <- grep(paste0("CORINE", corine_to_use), colnames(coefs))
  temp <- coefs[, index]
  out[, i] <- apply(temp, 1, mean)
}
return(out)
}


generate_overall_average_sgns_matrix <- function(list_of_month_sgns_coef_matrices_multiple_seeds,
            Months_in_study,
            random_seeds,
            corine_years = c(2006, 2012, 2018)){
  coefs_list <- list_of_month_sgns_coef_matrices_multiple_seeds
  N_months <- length(Months_in_study)
  N_matrices <- length(coefs_list)
  N_random_seeds <- length(random_seeds)
  N_CORINE <- length(corine_years)
  if (N_months!=N_matrices) stop("different number of months to matrices")
  ncol_matrix <- function(data_matrix){
    out <- dim(data_matrix)[2]
    return(out)
  }
  matrix_ncols <- lapply(coefs_list, ncol_matrix)
  if (any(matrix_ncols!=N_random_seeds*N_CORINE)) stop("incorrect number of columns")
  rownames_list <- lapply(coefs_list, rownames)
  combined_rownames <- union(rownames_list[[1]], rownames_list[[2]])
  if (N_matrices>=3){
    for (i in 3:N_matrices){
      combined_rownames <- union(combined_rownames, rownames_list[[i]])
    }
  }
  vars_stripX <- sapply(combined_rownames, stripX)
  new_order <- sort(vars_stripX, index.return=TRUE)$ix
  combined_rownames_sorted <- combined_rownames[new_order]
  overall_Nrow <- length(combined_rownames_sorted)
  overall_average_sgns <- matrix(0, overall_Nrow, N_months*N_CORINE)
  rownames(overall_average_sgns) <- combined_rownames_sorted
  colnames(overall_average_sgns) <- paste0(rep(Months_in_study, N_CORINE),
                                           rep(corine_years, rep(N_months, N_CORINE)))
  for (i in 1:N_matrices){
    M <- coefs_list[[i]]
    index <- seq(i, i+(N_CORINE-1)*N_months, by=N_months)
    overall_average_sgns[rownames(M), index] <- month_average_sgns(M, corine_years=corine_years)
  }
out <- overall_average_sgns
return(out)
}



strip_year <- function(Month_year){
  out <- strsplit(Month_year, split="_")[[1]][1]
  return(out)
}

strip_month <- function(Month_year){
  out <- strsplit(Month_year, split="_")[[1]][2]
  return(out)
}

generate_label <- function(X_var, clc_climate_lookup){
  if (X_var=="(Intercept)") out <- NA
  else{
    index <- which(clc_climate_lookup[,1]==X_var)
    out <- clc_climate_lookup[index, 2]}
  return(out)
}



