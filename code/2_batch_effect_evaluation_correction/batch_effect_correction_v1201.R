
# Batch Effect Analysis and Correction Script
# Can be run directly in RStudio; intended for a 64-core server with large memory.

# Load required packages
# BiocManager::install("sva",force = TRUE )
# BiocManager::install("locfit",force = TRUE )
library(preprocessCore)  # for normalize.quantiles
library(sva)            # for ComBat
library(umap)           # for UMAP
library(ggplot2)        # for plotting
library(randomForest)   # for random forest
library(parallel)       # for parallel processing
library(BiocParallel)   # for BiocParallel
library(readxl)
library(magrittr)
library(dplyr)
setwd(r"(\\172.16.13.136\tphp\code\2_batch_effect_evaluation_correction\)")
source("../source/src_claude_v1.R",encoding = "UTF-8")
# Set parallel backend (use 64 cores)
register(MulticoreParam(workers = 64))

#===================== Main workflow =====================
  
data_matrix <- read.csv("../1_NA_filter_imputation/output/quantiled_log2_imputated_matrix_2957_12797_imputed.csv", row.names = 1,check.names = FALSE )
head(data_matrix[,1:10])
sample_info <- read.csv("../0_3_sample_info_mapping_data_matrix/output/batch_info_tab_batchserver.csv")
dim(data_matrix)  #[1] 12797  2957
dim(sample_info)   #[1] 2957   35
str(sample_info)
data_matrix<-data_matrix[,sample_info$FileName]
all(colnames(data_matrix) == as.character(sample_info$FileName)) #T

# 5. Run ComBat batch effect correction
message("Running ComBat batch effect correction...")
# Diagnostic 1: check batch and biological variable relationships
# Examine the relationship between batches and biological variables
# sample_info<-info
table(sample_info$batch, sample_info$tissue_type_major)
table(sample_info$instrument, sample_info$sample_type)
# If it is found that the batch is highly correlated with biological variables, this will lead to correction problems
library(limma)

# Method1: removeBatchEffect
corrected_data <- removeBatchEffect(
  x = data_matrix,
  batch = sample_info$new,  #new
  design = model.matrix(~ imputation_group, data = sample_info)
)
head(corrected_data[,1:10])
dim(corrected_data)
write.csv(corrected_data,"output/batch_corrected_data_limma3.csv",row.names = T)

# Method2: use ComBat (only for new)
# dim(data_matrix)
combat_result <- sva::ComBat(
  dat = data_matrix,
  batch = sample_info$new,  #date
  mod = model.matrix(~ imputation_group, data = sample_info),
  par.prior = TRUE
)

# Get corrected data
View(combat_result[,1:10])
dim(combat_result)
# corrected_data <- t(combat_result)

# Save corrected data
write.csv(combat_result, "output/batch_corrected_data_combat3.csv", quote = FALSE)


# Diagnostic2: different correction strategies
# Strategy A: stepwise correction (instrument first, then new)

# Step 1: correct instrument, protect other variables
step1_corrected <- removeBatchEffect(
  x = data_matrix,
  batch = sample_info$new,
  design = model.matrix(~ instrument + sample_type + tissue_type_major+imputation_group, 
                        data = sample_info)
)

# Step 2: correct new in already corrected data
step2_corrected <- removeBatchEffect(
  x = step1_corrected,
  batch = sample_info$instrument,
  design = model.matrix(~ sample_type + tissue_type_major+imputation_group, 
                        data = sample_info)
)
dim(step2_corrected)
View(step2_corrected[,1:10])
write.csv(step2_corrected, "output/batch_corrected_data_step2_limma3.csv", quote = FALSE)


# Strategy B: Use SVA (Surrogate Variable Analysis) to identify unknown batches
# library(sva)

# Define biological model
# sample_type + detailed/major tissue types failed
# mod_bio <- model.matrix(~ sample_type, data = sample_info)
# mod0 <- model.matrix(~ 1, data = sample_info)

# Identify surrogate variables
# sv_obj <- sva(data_matrix, mod_bio, mod0)
# dim(data_matrix)
# dim(mod_bio)
# dim(mod0)

# Use identified surrogate variables for correction
# corrected_sva <- removeBatchEffect(
#   x = data_matrix,
#   covariates = sv_obj$sv,
#   design = mod_bio
# )
# dim(corrected_sva)
# write.csv(corrected_sva, "batch_corrected_data_sva_obj.csv")

# Strategy C: Conservative correction - remove only part of the batch effect
# Use ComBat with mean.only parameter
combat_conservative <- sva::ComBat(
  dat = data_matrix,
  batch = sample_info$new,
  mod = model.matrix(~  imputation_group, data = sample_info),
  mean.only = TRUE  # correct mean only, not variance
)
dim(combat_conservative)
write.csv(combat_conservative, "output/batch_corrected_data_combat_mean_only3.csv", quote = FALSE)


# Strategy D: Use random forest to analyze pool matrix and remove batch-effect-related proteins
# message("Running random forest analysis...")
# dim(data_matrix)
# pool_matrix <-t(data_matrix[,sample_info$tissue_type_major == "pool" ])
# rf_data <- data.frame(pool_matrix, label = as.factor(pool_info$instrument))
# rf_results <- myRF(rf_data, ntree = 500, nodesize = 5)
# rf_data2 <- data.frame(pool_matrix, label = as.factor(pool_info$new))
# rf_results2 <- myRF(rf_data2, ntree = 500, nodesize = 5)
# sample_matrix<-t(data_matrix[,sample_info$tissue_type_major != "pool"])
# sample_info2 <- sample_info[sample_info$tissue_type_major != "pool", ]
# rf_data3 <- data.frame(sample_matrix, label = sample_info2$sample_type)
# rf_results3 <- myRF(rf_data3, ntree = 500, nodesize = 5)

# print(str(rf_data3))
# print(table(rf_data3$label))

# Check if there are non-numeric columns
# non_numeric <- !sapply(rf_data3[,-ncol(rf_data3)], is.numeric)
# if(any(non_numeric)) {
#   print(paste("Non-numeric columns:", names(rf_data)[non_numeric]))
# }

# rf_results3 <- myRF_v2(rf_data3, ntree = 500, nodesize = 5)


# rf_data4 <- data.frame(sample_matrix, label = as.factor(sample_info2$tissue_type_major))
# rf_results4 <- myRF(rf_data4, ntree = 500, nodesize = 5)

# rf_data5 <- data.frame(sample_matrix, label = as.factor(sample_info2$batch))
# rf_results5 <- myRF(rf_data5, ntree = 500, nodesize = 5)

# Get top N important features
# top_n <- 100
# top_features <- rf_results$features[1:top_n]
# message(paste("Top", top_n, "batch-related features identified"))
# top_features2<- rf_results2$features[1:top_n]
# intersect(top_features, top_features2 )
# top_features5<- rf_results5$features[1:top_n]
# # View(rf_results5)
# batch_proteins0<-union(top_features, top_features2)
# batch_proteins0<-union(batch_proteins0, top_features5)
# length(batch_proteins0)  #278

# top_n_b <-1000
# top_features3 <- rf_results3$features[1:top_n_b]
# intersect(top_features3, batch_proteins0 )
# top_features4 <- rf_results4$features[1:top_n_b]    
# intersect(top_features4, batch_proteins0)
# features_to_add<-union(intersect(top_features3, batch_proteins0), intersect(top_features4, batch_proteins0))      
# writeClipboard(features_to_add )
# features_to_add<-features_to_add[-which(features_to_add=="P69905")]
# features_to_add
# batch_proteins <-setdiff( batch_proteins0,features_to_add)
# length(batch_proteins)  

# Save feature importance lists
# # write.csv(data.frame(
# #   Feature = c(rf_results$features, rf_results2$features),
# #   Importance = c(rf_results$importance[rf_results$features,], rf_results2$importance[rf_results2$features,] ),
# #   RF_Model = c(rep("RF_Model1", length(rf_results$features), rep("RF_Model2", length(rf_results2$features))))
# # ), "rf_feature_importance_20250730.csv", quote = FALSE, row.names = FALSE)
# writexl::write_xlsx(list(
#   RF_Model1 = data.frame(Feature = rf_results$features, Importance = rf_results$importance[rf_results$features,]),
#   RF_Model2 = data.frame(Feature = rf_results2$features, Importance = rf_results2$importance[rf_results2$features,]),
#   RF_Model3 = data.frame(Feature = rf_results3$features, Importance = rf_results3$importance[rf_results3$features,]),
#   RF_Model4 = data.frame(Feature = rf_results4$features, Importance = rf_results4$importance[rf_results4$features,]),
#   RF_Model5 = data.frame(Feature = rf_results5$features, Importance = rf_results5$importance[rf_results5$features,])
# ), "rf_feature_importance_20250812.xlsx")

# # Create dataset after removing top features
# features_to_keep <- setdiff(colnames(data_matrix), batch_proteins)
# filtered_data <- data_matrix[, features_to_keep]
# dim(filtered_data)
# write.csv(t(filtered_data), "batch_corrected_data_RF_pool_removal.csv", quote = FALSE)

# #Strategy E: Use pool
