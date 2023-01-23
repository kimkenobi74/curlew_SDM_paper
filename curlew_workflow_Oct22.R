library(ggplot2)
library(coefplot)
library(glmnet)
library(gridExtra)
library(plyr)
library(dplyr)
library(ggcorrplot)
library(raster)
library(sp)
library(sf)
library(sjmisc)
library(sjlabelled)
library(boot)
library(maps)
library(maptools)
library(ncdf4) # package for netcdf manipulation
library(rgdal) # package for geospatial analysis
library(RColorBrewer)
library(magrittr)
library(reshape)
library(scales)
library(PresenceAbsence)
library(ecospat)
library(rasterVis)
library(grid)

source("scripts/lasso_workflow_functions.R")

# read in aggregated presence data ------
# for GB and Ireland --------------------

GB_IRL_agg_presence <- readRDS("data/british_isles_agg_presence_month_corine_proj_sqkm.rds")
GB_IRL_agg_presence <- GB_IRL_agg_presence[,c("obs_year", "obs_month", "corine_east", "corine_north")]

# read in CORINE landcover classes
# for 2006, 2012 and 2018

corine_2006 <- readRDS("data/corine_2006_land_cover_proportions_from_centroids.rds")
corine_2012 <- readRDS("data/corine_2012_land_cover_proportions_from_centroids.rds")
corine_2018 <- readRDS("data/corine_2018_land_cover_proportions_from_centroids.rds")

# append CORINE year to names of variables

names(corine_2006) <- paste0(names(corine_2006),"_2006")
names(corine_2012) <- paste0(names(corine_2012),"_2012")
names(corine_2018) <- paste0(names(corine_2018),"_2018")

# read in legends for land cover classes

clc_legend <- read.csv("data/clc_legend.csv")
clc_labels <- clc_legend$LABEL3
# clc_labels[21] <- "Principally agriculture, with natural vegetation"
clc_labels <- revalue(clc_labels, 
                      c("Land principally occupied by agriculture, with significant areas of natural vegetation" = "Principally agriculture"))

# define months and years in study
# Jan, Feb, Nov, Dec for 2003 to 2018
# and Jan, Feb for 2019 
# (66 month-year combinations)

Years <- rep(2003:2019, c(rep(4, length(2003:2018)), 2))
Months <- c(rep(c(1,2,11,12), length(2003:2018)), 1, 2)
Month_strings <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

observation_df <- data.frame(GB_IRL_agg_presence)
Month_Year_strings <- paste(Month_strings[c(1,2,11,12)], Years, sep="_")

Years_Months_df <- data.frame(year = Years, month = Months)

# creating presence/absence rasters
# from curlew observation data

 t1 <- Sys.time()
 curlew_obs_list <- vector("list", length=length(Months))
 for (i in 1:length(Month_Year_strings)){
  print(i)
  Month <- Months[i]
  Year <- Years[i]
  curlew_obs_raster <- create_presence_absence_raster_month_year(observation_df, Year, Month)
  curlew_obs_list[[i]] <- curlew_obs_raster
 }
 names(curlew_obs_list) <- Month_Year_strings
 curlew_obs_complete <- brick(curlew_obs_list)
 t2 <- Sys.time()
 time_taken <- t2-t1
time_taken # about 6 seconds

saveRDS(curlew_obs_complete, file="output/curlew_obs_complete.rds")
curlew_obs_complete <- readRDS("output/curlew_obs_complete.rds")

# read in climate variables

drydaysSub_0_05mm_2003_2019 <- brick("data/climate_vars_IRL_UK/Merge_IRL_UK_drydaysSub_0_05mm_2003_2019.tif")
drydaysSub_0_5mm_2003_2019 <- brick("data/climate_vars_IRL_UK/Merge_IRL_UK_drydaysSub_0_5mm_2003_2019.tif")
frostdays_2003_2019 <- brick("data/climate_vars_IRL_UK/Merge_IRL_UK_frostdays_2003_2019.tif")
sub5Cdays_2003_2019 <- brick("data/climate_vars_IRL_UK/Merge_IRL_UK_sub5Cdays_2003_2019.tif")
totPrec_2003_2019 <- brick("data/climate_vars_IRL_UK/Merge_IRL_UK_totPrec_2003_2019.tif")

# read in DTM and slope percent variables

DTM_1km_OSGB <- brick("data/DTM_slope_vars/DTM_1km_OSGB.tif")
Slope_percent <- brick("data/DTM_slope_vars/Slope_Percent.tif")

# projecting DTM and slope percent
# to CORINE native grid

DTM_1km_OSGB_corine <- projectRaster(DTM_1km_OSGB, corine_2006)
Slope_percent_corine <- projectRaster(Slope_percent, corine_2006)

# Writing individual rasters for the climate variables
# Note, we have climate variables available
# for all 12 months each year 2003-2018
# and for January to August 2019.
# This is why we define 
# index_of_months_in_study.

index_of_months_in_study <- 
  c(c(1, 2, 11, 12) + 12*rep(0:15, rep(4, length(0:15))),
    12*16+c(1,2))

                              

filestem <- "output/climate_vars_CORINE_projection/"

# t1 <- Sys.time()
# count <- 0
# for (i in index_of_months_in_study){
#  count <- count+1
#  print(count)
#  Year <- Years[count]
#  Month <- Months[count]
##### dry days sub 0.05mm
#  temp <- projectRaster(drydaysSub_0_05mm_2003_2019[[i]],
#                        corine_2006[[1]])
#  filename <- paste0(filestem, "drydaysSub_0_05mm_", Year, ".", Month, ".tif")
#  writeRaster(temp, filename, overwrite=TRUE)
##### dry days sub 0.5mm
#  temp <- projectRaster(drydaysSub_0_5mm_2003_2019[[i]],
#                        corine_2006[[1]])
#  filename <- paste0(filestem, "drydaysSub_0_5mm_", Year, ".", Month, ".tif")
#  writeRaster(temp, filename, overwrite=TRUE)
##### frost days
#  temp <- projectRaster(frostdays_2003_2019[[i]],
#                        corine_2006[[1]])
#  filename <- paste0(filestem, "frostdays_", Year, ".", Month, ".tif")
#  writeRaster(temp, filename, overwrite=TRUE)
##### sub 5 deg C days
#    temp <- projectRaster(sub5Cdays_2003_2019[[i]],
#                        corine_2006[[1]])
#  filename <- paste0(filestem, "sub5Cdays_", Year, ".", Month, ".tif")
#  writeRaster(temp, filename, overwrite=TRUE)
##### Total precipitation
#    temp <- projectRaster(totPrec_2003_2019[[i]],
#                        corine_2006[[1]])
#  filename <- paste0(filestem, "totPrec_", Year, ".", Month, ".tif")
#  writeRaster(temp, filename, overwrite=TRUE)
#}
#t2 <- Sys.time()
#time_taken <- t2-t1
#time_taken
#### About 27 minutes

lookup_table_with_clc_plus_climate_plus_DTM <-
  data.frame(X_var = paste0("X", c(1:43, 101:105, 201:202)),
             label = c(clc_labels[1:43], "Dry sub 0.05mm", 
                       "Dry sub 0.5mm",
                       "Frost days",
                       "Days sub 5 degrees C",
                       "Total precipitation",
                       "DTM",
                       "Slope percent"))

A <- lookup_table_with_clc_plus_climate_plus_DTM
A[4,2] <- "Road and rail networks"
A[11,2] <- "Sport and leisure"
A[19, 2]<- "Annual-permanent crop mix"
clc_climate_lookup <- A

corine_to_use_rep <- rep(c(2006, 2012, 2018), rep(4, 3))
Month_numeric_vector <- rep(c(1, 2, 11, 12), 3)


# Creating the rasterstacks of explanatory variables and response variables
## NOTE: The function save_exp_vars_and_respnose_vars_stacks makes use of
# the function build_exp_vars_rasterstack, which is where the
# DTM and slope percent variables (from Crona) are included in the model.

output_filenamebase <- "output/cleaned_curlew_obs_and_exp_vars_CORINE/"
climate_vars_CORINE_proj_filestem <- "output/climate_vars_CORINE_projection/"

# t1 <- Sys.time()
# for (i in 1:length(corine_to_use_rep)){
#  print(i)
#  corine_to_use <- corine_to_use_rep[i]
#  Month_numeric <- Month_numeric_vector[i]

# save_exp_vars_and_response_vars_stacks(output_filenamebase, 
#                                       corine_to_use, 
#                                       Month_numeric,
#                                       climate_vars_CORINE_proj_filestem)
#}
#t2 <- Sys.time()
#time_taken <- t2-t1
#time_taken
# 1 minute 10 seconds 

# saving the data frames to use with glmnet
# The DTM and slope percent variables feature in the 
# save_x_and_y_dfs_for_glmnet function.

# t1 <- Sys.time()
# for (i in 1:length(corine_to_use_rep)){
#  print(i)
#  corine_to_use <- corine_to_use_rep[i]
#  Month_numeric <- Month_numeric_vector[i]
#  save_x_and_y_dfs_for_glmnet(corine_to_use, Month_numeric,
#                            exp_and_response_vars_filestem =
#                              "output/cleaned_curlew_obs_and_exp_vars_CORINE/",
#                            output_filestem = "output/x_and_y_dfs_for_glmnet/")
#}
# t2 <- Sys.time()
# time_taken <- t2-t1
# time_taken
# About 10 minutes

# The code here is getting the dataframes (on a monthly basis) ready 
# for the lasso step

exp_and_response_vars_filestem <- "output/cleaned_curlew_obs_and_exp_vars_CORINE/"
output_filestem <- "output/monthly_dfs/"
  
# t1 <- Sys.time()
# for (index in 1:length(corine_to_use_rep)){
#  print(index)
#  save_dfs_and_tifs_for_monthly_models(index,
#     exp_and_response_vars_filestem,
#     output_filestem)
# }
#t2 <- Sys.time()
#time_taken <- t2-t1
#time_taken
# 27 minutes

x_y_df_filestem <- "output/x_and_y_dfs_for_glmnet/"
cv.glmnet_filestem <- "output/cv.glmnet_fits/"

random_seeds <- as.numeric(paste0(1:31, "1022"))

# roughly 10 minutes per random seed = about 5-6 hours for 31 seeds
for(i in c(1:31)){
  print(i)
  random_seed <- random_seeds[i]
  cv.glmnet_filestem_seed <- paste0(cv.glmnet_filestem, "seed_",random_seed,"/")
  numCores <- parallel::detectCores()
  system.time(
    parallel::mclapply(1:length(corine_to_use_rep), read_in_x_y_df_save_cv.glmnet, mc.cores = numCores,
                       random_seed = random_seed, x_y_df_filestem = x_y_df_filestem,
                       cv.glmnet_filestem = cv.glmnet_filestem_seed,
                       corine_to_use_rep = corine_to_use_rep,
                       Month_numeric_vector = Month_numeric_vector,
                       Month_strings = Month_strings)
  )
}


 
 
  
 
 ####################################################################
 
 # Creating matrices for CORINE2006, CORINE2012 and CORINE2018
 # of the coefficients obtained by applying the cross-validated lasso
 # (Including the DTM and slope percent variables at this stage)
 
 
 cv.glmnet_filestem <- "output/cv.glmnet_fits/"
 
 Months_in_study <- c("Jan", "Feb", "Nov", "Dec")

#### coef_matrix_2006 
  
coef_matrix_2006_multiple_seeds <- coef_matrix_multiple_seeds(2006, Months_in_study,
                                              corine_to_use_rep,
                                              Month_numeric_vector,
                                              Month_strings,
                                              cv.glmnet_filestem,
                                              random_seeds)
coef_matrix_2012_multiple_seeds <- coef_matrix_multiple_seeds(2012, Months_in_study,
                                              corine_to_use_rep,
                                              Month_numeric_vector,
                                              Month_strings,
                                              cv.glmnet_filestem,
                                              random_seeds)
coef_matrix_2018_multiple_seeds <- coef_matrix_multiple_seeds(2018, Months_in_study,
                                              corine_to_use_rep,
                                              Month_numeric_vector,
                                              Month_strings,
                                              cv.glmnet_filestem,
                                              random_seeds)

## The following generates average sign matrix of coefficients
## (not including climate variables) and
## plots overall heatmap (all random seeds)

for (i in 1:length(Months_in_study)){
  Month_of_interest <- Months_in_study[i]
  plot_output_file <- paste0("figs/coefs_",Month_of_interest,"_multiple_seeds.pdf")
  assign(paste0("A", i), plot_coefs_multiple_seeds(Month_of_interest,
                    coef_matrix_2006_multiple_seeds,
                    coef_matrix_2012_multiple_seeds,
                    coef_matrix_2018_multiple_seeds,
                    random_seeds,
                    clc_labels,
                    plot_output_file))
}

L <- list(A1, A2, A3, A4)

average_sgns_matrix <- generate_overall_average_sgns_matrix(L, Months_in_study, 
                                  random_seeds)

plot_average_coef_sgns_matrix(average_sgns_matrix,
                              Months_in_study,
                              clc_labels,
                              plot_output_file = "figs/average_coef_sgns_matrix.pdf",
                              random_seeds=random_seeds)


## From here on we include the climate variables

for (i in 1:length(Months_in_study)){
  Month_of_interest <- Months_in_study[i]
  assign(paste0("B", i), generate_coef_sgns_inc_climate(Month_of_interest,
                                                   coef_matrix_2006_multiple_seeds,
                                                   coef_matrix_2012_multiple_seeds,
                                                   coef_matrix_2018_multiple_seeds,
                                                   random_seeds))
}

LB <- list(B1, B2, B3, B4)

                                                            
average_sgns_matrix_inc_climate <- generate_overall_average_sgns_matrix(LB,
                                                                  Months_in_study,
                                                                  random_seeds)

Jan_consensus <- find_consensus_coef_sgns_threshold(B1, random_seeds,
                                                              threshold=0.64)
Feb_consensus <- find_consensus_coef_sgns_threshold(B2, random_seeds,
                                                              threshold=0.64)
Nov_consensus <- find_consensus_coef_sgns_threshold(B3, random_seeds,
                                                              threshold=0.64)
Dec_consensus <- find_consensus_coef_sgns_threshold(B4, random_seeds,
                                                              threshold=0.64)

 
combined_row_names <- union(union(union(rownames(Jan_consensus), rownames(Feb_consensus)),
                            rownames(Nov_consensus)),
                            rownames(Dec_consensus))
vars_stripX <- sapply(combined_row_names, stripX)
new_order <- sort(vars_stripX, index.return=TRUE)$ix
combined_row_names_sorted <- combined_row_names[new_order]

overall_N_row <- length(combined_row_names_sorted)
overall_consensus_sgns <- matrix(0, overall_N_row, 12)
rownames(overall_consensus_sgns) <- combined_row_names_sorted
colnames(overall_consensus_sgns) <- paste0(Months_in_study, corine_to_use_rep)

overall_consensus_sgns[rownames(Jan_consensus), c(1, 5, 9)] <- Jan_consensus
overall_consensus_sgns[rownames(Feb_consensus), c(2, 6, 10)] <- Feb_consensus
overall_consensus_sgns[rownames(Nov_consensus), c(3, 7, 11)] <- Nov_consensus
overall_consensus_sgns[rownames(Dec_consensus), c(4, 8, 12)] <- Dec_consensus

saveRDS(overall_consensus_sgns, "output/coef_matrices/overall_consensus_sgns.rds")
overall_consensus_sgns <- readRDS("output/coef_matrices/overall_consensus_sgns.rds")

# The plot of consensus signs of landcover variables in models across all three CORINE windows

plot_coef_sgns_matrix(overall_consensus_sgns[-which(rownames(overall_consensus_sgns) %in% c("X101", "X102", "X103", "X104", "X105", "X201", "X202")),],
                      Months_in_study,
                      clc_labels,
                      plot_output_file = "figs/consensus_coefs.pdf")


overall_consensus_sgns_2006 <- overall_consensus_sgns[, 1:4]
overall_consensus_sgns_2012 <- overall_consensus_sgns[, 5:8]
overall_consensus_sgns_2018 <- overall_consensus_sgns[, 9:12]

Years_CORINE2006 <- rep(2003:2008, rep(4, length(2003:2008)))
Months_numeric_CORINE2006 <- rep(c(1, 2, 11, 12), length(2003:2008))
Years_CORINE2012 <- rep(2009:2014, rep(4, length(2009:2014)))
Months_numeric_CORINE2012 <- rep(c(1, 2, 11, 12), length(2009:2014))
Years_CORINE2018 <- c(rep(2015:2018, rep(4, length(2015:2018))), rep(2019,2))
Months_numeric_CORINE2018 <- c(rep(c(1, 2, 11, 12), 
                                   length(2015:2018)),
                               1:2)

# Fitting the GLMs (50% training data, 50% test data)

exp_and_response_vars_filestem <- "output/cleaned_curlew_obs_and_exp_vars_CORINE/"
monthly_dfs_filestem <- "output/monthly_dfs/"
GLMs_by_month_train_test_filestem <- "output/GLMs_by_month_train_test/"



t1 <- Sys.time()
for (i in 1:length(Years_CORINE2006)){
  print(i)
  Year <- Years_CORINE2006[i]
  Month_numeric <- Months_numeric_CORINE2006[i]
  GLM_train_test_object <- fit_GLM_train_test(Year, Month_numeric,
                                              exp_and_response_vars_filestem,
                                              overall_consensus_sgns_2006,
                                              monthly_dfs_filestem,
                                              random_seed = 2006)
  file_to_write_stem <- GLMs_by_month_train_test_filestem
  file_to_write <- paste0(file_to_write_stem, "GLM_", Year, ".", Month_numeric, "_train_test.rds")
  saveRDS(GLM_train_test_object, file_to_write)
  rm(GLM_train_test_object)
}
t2 <- Sys.time()
time_taken <- t2-t1
time_taken # 2.4 minutes


t1 <- Sys.time()
for (i in 1:length(Years_CORINE2012)){
  print(i)
  Year <- Years_CORINE2012[i]
  Month_numeric <- Months_numeric_CORINE2012[i]
  GLM_train_test_object <- fit_GLM_train_test(Year, Month_numeric,
                                              exp_and_response_vars_filestem, 
                                              overall_consensus_sgns_2012,
                                              monthly_dfs_filestem,
                                              random_seed = 2012)
  file_to_write_stem <- GLMs_by_month_train_test_filestem
  file_to_write <- paste0(file_to_write_stem, "GLM_", Year, ".", Month_numeric, "_train_test.rds")
  saveRDS(GLM_train_test_object, file_to_write)
  rm(GLM_train_test_object)
}
t2 <- Sys.time()
time_taken <- t2-t1
time_taken # 2.6 minutes


t1 <- Sys.time()
for (i in 1:length(Years_CORINE2018)){
  print(i)
  Year <- Years_CORINE2018[i]
  Month_numeric <- Months_numeric_CORINE2018[i]
  GLM_train_test_object <- fit_GLM_train_test(Year, Month_numeric,
                                              exp_and_response_vars_filestem, 
                                              overall_consensus_sgns_2018,
                                              monthly_dfs_filestem,
                                              random_seed = 2018)
  file_to_write_stem <- GLMs_by_month_train_test_filestem
  file_to_write <- paste0(file_to_write_stem, "GLM_", Year, ".", Month_numeric, "_train_test.rds")
  saveRDS(GLM_train_test_object, file_to_write)
  rm(GLM_train_test_object)
}
t2 <- Sys.time()
time_taken <- t2-t1
time_taken # 1.9 minutes

# Looking at the coefficients of the models
vars_stripX <- sapply(overall_set_of_rownames, stripX)
new_order <- sort(vars_stripX, index.return=TRUE)$ix
overall_set_of_rownames_ordered <- overall_set_of_rownames[new_order]


N_month_Year_combs <- length(Years)
GLM_coefs <- matrix(0, length(overall_set_of_rownames_ordered), N_month_Year_combs)
GLM_sd_coefs <- matrix(0, length(overall_set_of_rownames_ordered), N_month_Year_combs)
colnames(GLM_coefs) <- paste0(Month_strings[Months], "_", Years)
colnames(GLM_sd_coefs) <- colnames(GLM_coefs)
rownames(GLM_coefs) <- overall_set_of_rownames_ordered
rownames(GLM_sd_coefs) <- rownames(GLM_coefs)

t1 <- Sys.time()
for (i in 1:N_month_Year_combs){
  print(i)
  col_index <- i
  Year <- Years[i]
  Month_numeric <- Months[i]
  file_to_read_base <- "output/GLMs_by_month_train_test/"
  file_to_read <- paste0(file_to_read_base, 
                         "GLM_", Year, ".", 
                         Month_numeric, "_train_test.rds")
  GLM_train_test_object <- readRDS(file_to_read)
  coefs_sds <- extract_coefs_sds_GLM(GLM_train_test_object)
  index <- index_GLM_coefs_rownames(coefs_sds)
  # ignoring (Intercept) here
  GLM_coefs[index, col_index] <- coefs_sds$coefs[2:length(coefs_sds$coefs)]
  GLM_sd_coefs[index, col_index] <- coefs_sds$sds[2:length(coefs_sds$sds)]
}
t2 <- Sys.time()
time_taken <- t2-t1
time_taken # 34 seconds

saveRDS(GLM_coefs, "output/GLMs_by_month_train_test/GLM_coefs.rds")
saveRDS(GLM_sd_coefs, "output/GLMs_by_month_train_test/GLM_sd_coefs.rds")

GLM_coefs <- readRDS("output/GLMs_by_month_train_test/GLM_coefs.rds")
GLM_sd_coefs <- readRDS("output/GLMs_by_month_train_test/GLM_sd_coefs.rds")

var_IDs <- rownames(GLM_coefs)
df <- data.frame(variable = rep(var_IDs,
                                dim(GLM_coefs)[2]),
                 coef = as.numeric(GLM_coefs),
                 sd = as.numeric(GLM_sd_coefs),
                 Month = rep(sapply(colnames(GLM_coefs), 
                                    strip_year),
                             rep(dim(GLM_coefs)[1],
                                 dim(GLM_coefs)[2])),
                 Year = rep(sapply(colnames(GLM_coefs),
                                   strip_month),
                            rep(dim(GLM_coefs)[1],
                                dim(GLM_coefs)[2])))

df$variable <- factor(df$variable,
                      ordered=TRUE,
                      levels=gtools::mixedsort(rownames(GLM_coefs)))
labels <-sapply(unique(df$variable), generate_label,
                clc_climate_lookup = clc_climate_lookup)

df$label <- factor(rep(labels, dim(df)[1]/length(unique(df$variable))),
                   ordered=TRUE,
                   levels=labels)

# Plotting coefficients of models (including climate variables)
# Firstly a three-page pdf for each month (2006, 2012, 2018) ...

for (i in 1:length(Months_in_study)){
  Month_string <- Months_in_study[i]
  filename <- paste0("figs/", Month_string, "s_coefs_separate_CORINE.pdf")
  plot_monthly_coefs_individual_CORINE_windows(df, Month_string, filename, 
                                               ylim=c(-0.25, 0.25))
}

# then separate pdfs for each month-CORINE window combination

for (i in 1:length(Months_in_study)){
  Month_string <- Months_in_study[i]
  filename2006 <- paste0("figs/", Month_string, "s_coefs_separate_CORINE2006.pdf")
  filename2012 <- paste0("figs/", Month_string, "s_coefs_separate_CORINE2012.pdf")
  filename2018 <- paste0("figs/", Month_string, "s_coefs_separate_CORINE2018.pdf")
  plot_monthly_coefs_separate_pdfs_CORINE_windows(df, Month_string, 
                                                  filename2006,
                                                  filename2012,
                                                  filename2018,
                                                  ylim=c(-0.25, 0.25))
}

# Boyce and AUC values on test data (50%)

Years_CORINE_combined <- c(Years_CORINE2006, Years_CORINE2012, Years_CORINE2018)
Months_numeric_CORINE_combined <- c(Months_numeric_CORINE2006, 
                                    Months_numeric_CORINE2012,
                                    Months_numeric_CORINE2018)
Years_and_Months_numeric <- data.frame(Year=Years_CORINE_combined,
                                       Month_numeric = Months_numeric_CORINE_combined)
N_models <- dim(Years_and_Months_numeric)[1]

boyce_indices <- numeric(N_models)
auc_values <- numeric(N_models)

GLMs_filestem <- "output/GLMs_by_month_train_test/"

t1 <- Sys.time()
for (i in 1:N_models){
  print(i)
  Year <- Years_and_Months_numeric[i, "Year"]
  Month_numeric <- Years_and_Months_numeric[i, "Month_numeric"]
  file_to_read <- paste0(GLMs_filestem, "GLM_", Year, ".", Month_numeric, "_train_test.rds")
  GLM_train_test_object <- readRDS(file_to_read)
  auc_boyce <- calculate_auc_plus_boyce_index(GLM_train_test_object)
  auc_values[i] <- auc_boyce$AUC
  boyce_indices[i] <- auc_boyce$boyce_index
}
t2 <- Sys.time()
time_taken <- t2-t1
time_taken # 40 seconds

model_evaluations_df <- data.frame(Year=Years_and_Months_numeric$Year,
                                   Month_numeric = Years_and_Months_numeric$Month_numeric,
                                   AUC = auc_values,
                                   boyce_index = boyce_indices)

model_evaluations_output_file <- paste0(GLMs_filestem, 
                                        "model_evaluations_df.rds")
saveRDS(model_evaluations_df, model_evaluations_output_file)
model_evaluations_df <- readRDS(model_evaluations_output_file)



model_evaluations_df$x_axis <- sapply(model_evaluations_df$Month_numeric, 
                                      x_axis_fn)



AUC_plot <- ggplot(model_evaluations_df, aes(x=x_axis, y=AUC)) +
  geom_point() +
  ggtitle("AUC for all models") +
  facet_wrap(~Year, ncol=6) +
  scale_x_continuous(breaks=1:4, labels=Month_strings[c(1,2,11,12)]) +
  scale_y_continuous(breaks = seq(0, 1, by=0.25), labels=seq(0,1, by=0.25)) +
  labs(x="Month") +
  theme(axis.text.x = element_text(angle=90)) +
  coord_cartesian(ylim=c(0, 1.025)) +
  geom_hline(yintercept=0.85, lty=2) +
  theme_bw()

pdf("figs/AUC_all_models_plot.pdf")
AUC_plot
dev.off()

Boyce_plot <- ggplot(model_evaluations_df, aes(x=x_axis, y=boyce_index)) +
  geom_point() +
  ggtitle("Boyce index for all models") +
  facet_wrap(~Year, ncol=6) +
  scale_x_continuous(breaks=1:4, labels=Month_strings[c(1,2,11,12)]) +
  scale_y_continuous(breaks = seq(0, 1, by=0.25)) +
  labs(x="Month", y="Boyce index") +
  theme(axis.text.x = element_text(angle=90)) +
  coord_cartesian(ylim=c(0, 1.025)) +
  geom_hline(yintercept=0.75, lty=2) +
  theme_bw()

pdf("figs/Boyce_index_all_models_plot.pdf")
Boyce_plot
dev.off()





# The following is needed for the probability rasters generated below 

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbSeq <- brewer.pal(8, "YlOrRd")
cbSeqBlues <- brewer.pal(5, "Blues")


worldmap = map_data('world')
UKmap <- worldmap %>%
  filter(region %in% c("Ireland", "Isle of Man", "UK"))
coords <- as.matrix(UKmap[,c("long", "lat")])
groups <- unique(UKmap$group)
polygon_list <- vector("list", length(groups))
for (i in 1:length(groups)){
  print(i)
  group <- groups[i]
  index <- which(UKmap$group==group)
  polygon_list[[i]] <- Polygon(coords[index,])
}
polygons <- Polygons(polygon_list, "group")


# Convert UKmap dataframe to SpatialPolygonsDataFrame
geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
polygons_sp <- SpatialPolygons(list(polygons), proj4string = CRS(geo.prj))
df <- data.frame(temp = 1)
rownames(df) <- "group"
polygons_spdf <- SpatialPolygonsDataFrame(polygons_sp, df)

GLM_prob_rasters_filestem <- "output/GLM_prob_rasters/"
GLM_train_test_filestem <- "output/GLMs_by_month_train_test/"


# With parallelisation
numCores <- parallel::detectCores()
system.time(
  parallel::mclapply(1:N_models, save_probs_raster,
          mc.cores = numCores,
          exp_and_response_vars_filestem = exp_and_response_vars_filestem,
          GLM_train_test_filestem = GLM_train_test_filestem,
          GLM_prob_rasters_filestem = GLM_prob_rasters_filestem)
)
# About 7 minutes


Months_in_study <- c("Jan", "Feb", "Nov", "Dec") #(already defined earlier)

t1 <- Sys.time()
save_means_and_sds_probs_rasters(2006, Months_in_study,
                              GLM_prob_rasters_filestem)
save_means_and_sds_probs_rasters(2012, Months_in_study,
                              GLM_prob_rasters_filestem)
save_means_and_sds_probs_rasters(2018, Months_in_study,
                              GLM_prob_rasters_filestem)
t2 <- Sys.time()
time_taken <- t2-t1
time_taken

monthly_means_and_sds_filestem <- "output/GLM_prob_rasters/monthly_means_and_sds/"

t1 <- Sys.time()
CORINE2006_means_sds_plot_no_title <- plot_means_sds_prob_rasters(2006,
                            Months_in_study,
                            GLM_prob_rasters_filestem,
                            monthly_means_and_sds_filestem,
                            col.regions_means = cbSeq,
                            col.regions_sds = cbSeqBlues,
                            use_title=FALSE)
CORINE2012_means_sds_plot_no_title <- plot_means_sds_prob_rasters(2012,
                            Months_in_study,
                            GLM_prob_rasters_filestem,
                            monthly_means_and_sds_filestem,
                            col.regions_means = cbSeq,
                            col.regions_sds = cbSeqBlues,
                            use_title=FALSE)
CORINE2018_means_sds_plot_no_title <- plot_means_sds_prob_rasters(2018,
                            Months_in_study,
                            GLM_prob_rasters_filestem,
                            monthly_means_and_sds_filestem,
                            col.regions_means = cbSeq,
                            col.regions_sds = cbSeqBlues,
                            use_title=FALSE)
t2 <- Sys.time()
time_taken <- t2-t1 # 1.5 minutes


