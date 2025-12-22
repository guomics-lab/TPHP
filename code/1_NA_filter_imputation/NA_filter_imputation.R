# Batch Effect Analysis and Correction Script
# 

# add libarary
# BiocManager::install("sva",force = TRUE )
# BiocManager::install("locfit",force = TRUE )
library(preprocessCore)  # for normalize.quantiles
library(sva)            # for ComBat
library(umap)           # for UMAP
library(ggplot2)        # for plotting
library(randomForest)   # for random forest
library(parallel)       # for parallel processing
library(BiocParallel)   # for BiocParallel
setwd(r"(\\172.16.13.136\tphp\code\1_NA_filter_imputation)")
source("../source/src_claude_v1.R")
# Set parallel processing parameters (to utilize 64 cores)
register(MulticoreParam(workers = 64))

#===================== Main analysis workflow  =====================#
  
  # 1. imput data
  message("Data input...")
# Change to your file path
data_file <- "input/all_mapped_pg_matrix_13609_2957_batchserver.csv"  # 或 .xlsx
sample_info_file <- "input/batch_info_tab_batchserver.csv"  # 或 .xlsx

# read data
if(grepl("\\.xlsx$", data_file)) {
  myd <- openxlsx::read.xlsx(data_file)
} else {
  myd <- read.csv(data_file, row.names = 1, check.names = FALSE)
}

# read sample information
if(grepl("\\.xlsx$", sample_info_file)) {
  sample_info <- openxlsx::read.xlsx(sample_info_file)
} else {
  sample_info <- read.csv(sample_info_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
}

# 2. Data preprocessing
message("Data preprocessing...")
data_matrix <- t(myd)
dim(data_matrix)  #2957 13609
dim(sample_info) #2957   34
# Missing value handling
na_filter_reult <- removeHighNAProteins(data_matrix, sample_info,
                                        group_columns = c("sample_type", "tissue_type_detailed"),
                                        data_start_col = 1,  #
                                        na_threshold = 0.5,
                                        group_threshold = 1.0)
#  sample_type, tissue_type_detailed 
# total 13609 prots，del 812 prots
# NA rate >= 50% 
# del prots(top 10)：A0A1B0GVR7, Q69YU5, Q96M34, Q99801, P01100, P21731, Q6UX46, A6NGC4, K7EJ46, O43184...  

data_matrix <- na_filter_reult$filtered_data
dim(data_matrix)   #[1]  2957 12797
write.csv(data_matrix, "output/data_matrix_filtered_05NA_2957_12797.csv", row.names = TRUE )


# Optional: Quantile normalization
data_matrix_q <- normalize.quantiles(as.matrix(t(data_matrix)))
View(head(data_matrix_q))
dim(data_matrix_q)
# Calculate the correlation per protein (column) before and after quantile normalization
cor_quantile <- sapply(1:nrow(data_matrix), function(i) {
  cor(data_matrix[,i], data_matrix_q[i,], use = "pairwise.complete.obs", method = "spearman")
})

cor_df <- data.frame(
  Protein = rownames(data_matrix),
  Correlation = cor_quantile
)
summary(cor_df$Correlation)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -0.4049  0.5659  0.7172  0.6700  0.8191  0.9712
write.csv(cor_df, "output/quantile_correlation.csv", row.names = FALSE)


# Optional:log2 transformation
data_q_log<- log2(data_matrix_q + 1)  
sample_na<-apply(data_q_log, 2, function(x) sum(is.na(x))) 
which(sample_na==nrow(data_q_log))  #0
pro_na<-apply(data_q_log, 1, function(x) sum(is.na(x))) 
which(pro_na==ncol(data_q_log)) #0
write.csv(data_q_log, "output/data_matrix_filtered_05NA_2957_12797.csv", row.names = TRUE )

#NA impute
#jitter_sd：default 1e-3, If your data are log2-transformed values and you want the noise to be smaller/larger, adjust this parameter (e.g., to 1e-2).
ncol(data_matrix)
min=min(data_q_log, na.rm = TRUE  )
min   #[1] 7.769049
# setseed(2025)
imputation<- log2(0.5) + min
imputation   #[1] 6.769049

# quantiled matrix， imputate with 0.5*min, log2
expr_imputed <- data_q_log
expr_imputed[is.na(expr_imputed)] <- imputation    
dim(expr_imputed)   #[1] 12797  2957
colnames(expr_imputed) <-  rownames(data_matrix)
rownames(expr_imputed) <- colnames(data_matrix)
write.csv(expr_imputed, "output/quantiled_log2_imputated_matrix_2957_12797_imputed.csv", row.names = TRUE )


#####NO quantile normalization, imputate with 0.5*min, log2
min <- min(data_matrix, na.rm = TRUE  )
min   #[1] 12.77158
max(data_matrix, na.rm = TRUE)
# setseed(2025)
imputation<- 0.5*min
imputation   #[1] 6.38579
expr_imputed <- data_matrix
expr_imputed[is.na(expr_imputed)] <- imputation   
dim(expr_imputed)
colnames(expr_imputed) <-  colnames(data_matrix)
rownames(expr_imputed) <- rownames(data_matrix)
write.csv(t(log2(expr_imputed)), "output/imputed_raw_12754proteins_2957samples.csv", row.names = TRUE )
