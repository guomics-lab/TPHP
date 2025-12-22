#20240921
rm(list = ls())
library(readxl)
setwd("//172.16.13.136/tphp/code/3_normal_unsupervised_clustering_heterogeneity/")
stringAsFactor =F

#0.read input (normal sample matrix)
pm<-read.csv("../2_batch_effect_evaluation_correction/output/batch_corrected_data_step2_limma3.csv",row.names = 1,check.names = F, stringsAsFactors = F)
pm<-as.data.frame(pm)
# pm_inf<-pm[,c(1:15)]
info<-read_xlsx("../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V6.xlsx")
# info<-read_xlsx("input/20240620_PUH_sample_information_1781files_info_edited_v8.xlsx")
info$FileName<-gsub(pattern = ".d$","",info$FileName)
# sample_info<-info[(info$sample_type=="N" &  info$Low_protein_IDs == FALSE),] #549
sample_info<-info[(info$sample_type=="N"),]#560
# sample_info<-info[(info$sample_type=="N" & info$Tissue_heterogeneity_abnormal == FALSE & info$Low_protein_IDs == FALSE),]
data<-apply(t(pm[,sample_info$FileName]),c(1,2),as.numeric)
dim(data)
sample_info[sample_info$anatomical_classification =="pancreas",]
#1.pre-process of data
# Data normalization
# After proteomics data have already undergone quantile normalization and log2 transformation,
# Z-score normalization is typically not required
# data_scaled <- scale(data)

# # Fill missing values (minimum-value imputation plus random noise)
# min_value <- min(data, na.rm = TRUE)
# 
# # Impute missing values as the minimum value plus random noise
# set.seed(42)  # Ensure reproducibility
# data_filled <- data
# data_filled[is.na(data)] <- log2(2^min_value * runif(sum(is.na(data)), min = 0.9, max = 1.1))
# data_scaled<-data_filled
# Sample information processing
rownames(sample_info)<-sample_info$FileName
group_labels<-data.frame(major_category = sample_info$anatomical_classification,
                         sub_category = sample_info$tissue_name)
rownames(group_labels)<-sample_info$FileName
# # 2. Calculate similarity and heterogeneity
# # Inspect the normalized data (mean should be 0 and standard deviation should be 1)
# summary(as.numeric(unlist(data_filled)))
# # summary(as.numeric(unlist(data_scaled)))
# 
# # Compute Pearson correlation coefficients between samples
# cor_matrix <- cor(t(data_filled), method = "pearson")
# rownames(cor_matrix)=colnames(cor_matrix)<-sample_info1$FileName
# # Show the first few rows of the correlation matrix
# head(cor_matrix)
# summary(as.numeric(unlist(cor_matrix)))
# #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.2336  0.6150  0.6780  0.6747  0.7419  1.0000 
# 
# # Compute coefficient of variation (CV = standard deviation / mean)
# cv <- apply(data_filled,1, function(x) sd(2^x) / mean(2^x))
# 
# # Inspect the distribution of CV
# summary(cv)
# 
# # Load required package
# library(stats)
# 
# # Compute Euclidean distance matrix
# dist_matrix <- dist(data_filled, method = "euclidean")
# 
# # Convert to matrix format
# dist_matrix <- as.matrix(dist_matrix)
# 
# # Show the first few rows of the distance matrix
# head(dist_matrix)

#3. Assess whether subcategories deviate from major categories
library(dplyr)
library(parallel)  # For parallel computing


# Count the number of samples in each sub_category
subcat_counts <- group_labels %>%
  group_by(sub_category) %>%
  tally()

# Filter out sub_categories with fewer than 2 samples
valid_subcats <- subcat_counts %>%
  filter(n >= 2) %>%
  pull(sub_category)

# Filtered data and labels
data_filled<-data
group_labels_filtered <- group_labels[group_labels$sub_category %in% valid_subcats, ]
data_filtered <- data_filled[rownames(group_labels_filtered), ]
dim(data_filtered)
group_labels_filtered_unique<-unique(group_labels_filtered )
rownames(group_labels_filtered_unique)<-group_labels_filtered_unique$sub_category

# Get the number of CPU cores
num_cores <- detectCores() - 1
# Create a parallel cluster
cl <- makeCluster(num_cores)

# Export the filtered data objects to the parallel cluster
clusterExport(cl, varlist = c("data_filtered", "group_labels_filtered"))

# Compute Euclidean distances and correlations between a specified tissue and other tissues
# Inputs:
#   mat: expression matrix (genes x samples)
#   sample_info: sample information data frame; must include the 'tissue' column
#   target_tissue: name of the target tissue
# Output:
#   distance summary matrix
# Compute pairwise Euclidean distances between all tissue types
calculate_pairwise_tissue_distances <- function(mat, index, sample_info) {
  # Input validation
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("mat must be a matrix or a data frame")
  }
  
  if (!is.character(index) || length(index) != 1) {
    stop("index must be a single string")
  }
  
  if (!(index %in% colnames(sample_info))) {
    stop(paste0("Sample information must include column '", index, "'"))
  }
  
  # Ensure matrix column names match sample information row names
  if (is.null(colnames(mat)) || is.null(rownames(sample_info))) {
    stop("The matrix must have column names and sample_info must have row names")
  }
  
  common_samples <- intersect(colnames(mat), rownames(sample_info))
  
  if (length(common_samples) == 0) {
    stop("No matches between matrix column names and sample_info row names")
  }
  
  if (length(common_samples) < ncol(mat) * 0.5) {
    warning(paste("Only", length(common_samples), "samples matched out of", 
                  ncol(mat), "samples"))
  }
  
  # Keep only matched samples
  mat <- mat[, common_samples, drop = FALSE]
  sample_info <- sample_info[common_samples, , drop = FALSE]
  
  # Remove genes containing NA values
  na_genes <- apply(mat, 1, function(x) any(is.na(x)))
  if (sum(na_genes) > 0) {
    warning(paste("Removing", sum(na_genes), "genes containing NA values"))
    mat <- mat[!na_genes, , drop = FALSE]
  }
  
  if (nrow(mat) == 0) {
    stop("No genes remain after removing NA values")
  }
  
  # Transpose matrix (samples as rows)
  mat_t <- t(mat)
  
  # Get all tissue types and their samples
  tissue_vector <- sample_info[[index]]
  tissues <- sort(unique(tissue_vector))  # Sort for reproducibility
  n_tissues <- length(tissues)
  
  if (n_tissues < 2) {
    stop("At least two tissue types are required for comparison")
  }
  
  # Compute distance matrix across all samples
  all_dist <- as.matrix(dist(mat_t, method = "euclidean"))
  
  # Compute correlation matrix across all samples
  all_cor <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  
  # Prepare sample indices for each tissue
  tissue_samples <- list()
  for (tissue in tissues) {
    samples <- which(tissue_vector == tissue)
    if (length(samples) > 0) {
      tissue_samples[[tissue]] <- samples
    } else {
      warning(paste0("Tissue '", tissue, "' has no samples"))
    }
  }
  
  # Initialize results list
  results <- list()
  
  # Iterate only over the upper triangle to avoid redundant calculations
  for (i in seq_along(tissues)) {
    tissue1 <- tissues[i]
    
    if (!(tissue1 %in% names(tissue_samples))) next
    
    idx1 <- tissue_samples[[tissue1]]
    samples1 <- rownames(sample_info)[idx1]
    n1 <- length(idx1)
    
    # Start from i (not i+1) to include within-tissue distances (diagonal blocks)
    for (j in i:n_tissues) {
      tissue2 <- tissues[j]
      
      if (!(tissue2 %in% names(tissue_samples))) next
      
      idx2 <- tissue_samples[[tissue2]]
      samples2 <- rownames(sample_info)[idx2]
      n2 <- length(idx2)
      
      # Extract distance/correlation submatrices between the two tissue sample sets
      if (i == j) {
        # Within-tissue distances/correlations (upper triangle only)
        if (n1 > 1) {
          dist_values <- all_dist[samples1, samples2][upper.tri(all_dist[samples1, samples2])]
          cor_values <- all_cor[samples1, samples2][upper.tri(all_cor[samples1, samples2])]
          n_comp <- n1 * (n1 - 1) / 2
        } else {
          next  # Cannot compute within-group distance with only one sample
        }
      } else {
        # Between-tissue distances/correlations (full block)
        dist_values <- as.vector(all_dist[samples1, samples2])
        cor_values <- as.vector(all_cor[samples1, samples2])
        n_comp <- n1 * n2
      }
      
      if (length(dist_values) == 0) next
      
      # Distance summary statistics
      mean_dist <- mean(dist_values, na.rm = TRUE)
      median_dist <- median(dist_values, na.rm = TRUE)
      sd_dist <- sd(dist_values, na.rm = TRUE)
      min_dist <- min(dist_values, na.rm = TRUE)
      max_dist <- max(dist_values, na.rm = TRUE)
      q25_dist <- quantile(dist_values, 0.25, na.rm = TRUE)
      q75_dist <- quantile(dist_values, 0.75, na.rm = TRUE)
      
      # Correlation summary statistics
      mean_cor <- mean(cor_values, na.rm = TRUE)
      median_cor <- median(cor_values, na.rm = TRUE)
      sd_cor <- sd(cor_values, na.rm = TRUE)
      min_cor <- min(cor_values, na.rm = TRUE)
      max_cor <- max(cor_values, na.rm = TRUE)
      q25_cor <- quantile(cor_values, 0.25, na.rm = TRUE)
      q75_cor <- quantile(cor_values, 0.75, na.rm = TRUE)
      
      # Store results
      results[[length(results) + 1]] <- data.frame(
        tissue1 = tissue1,
        tissue2 = tissue2,
        is_within_tissue = (i == j),
        # Distance-related results
        mean_distance = mean_dist,
        median_distance = median_dist,
        sd_distance = sd_dist,
        min_distance = min_dist,
        max_distance = max_dist,
        q25_distance = q25_dist,
        q75_distance = q75_dist,
        # Correlation-related results
        mean_correlation = mean_cor,
        median_correlation = median_cor,
        sd_correlation = sd_cor,
        min_correlation = min_cor,
        max_correlation = max_cor,
        q25_correlation = q25_cor,
        q75_correlation = q75_cor,
        # Sample metadata
        n_samples_tissue1 = n1,
        n_samples_tissue2 = n2,
        n_comparisons = n_comp,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Merge results
  if (length(results) == 0) {
    stop("Unable to compute distances/correlations for any tissue pairs")
  }
  
  summary_mat <- do.call(rbind, results)
  
  # Sort between-tissue and within-tissue results separately
  between_tissue <- summary_mat[!summary_mat$is_within_tissue, ]
  within_tissue <- summary_mat[summary_mat$is_within_tissue, ]
  
  # Sort by mean distance (alternatively, could sort by correlation)
  between_tissue <- between_tissue[order(between_tissue$mean_distance), ]
  within_tissue <- within_tissue[order(within_tissue$mean_distance), ]
  
  # Combine sorted results
  summary_mat <- rbind(between_tissue, within_tissue)
  rownames(summary_mat) <- NULL
  
  # Build full distance matrix (using mean distance)
  dist_matrix <- matrix(NA, nrow = n_tissues, ncol = n_tissues)
  rownames(dist_matrix) <- colnames(dist_matrix) <- tissues
  
  # Build full correlation matrix (using mean correlation)
  cor_matrix <- matrix(NA, nrow = n_tissues, ncol = n_tissues)
  rownames(cor_matrix) <- colnames(cor_matrix) <- tissues
  
  for (i in 1:nrow(summary_mat)) {
    t1 <- summary_mat$tissue1[i]
    t2 <- summary_mat$tissue2[i]
    d <- summary_mat$mean_distance[i]
    c <- summary_mat$mean_correlation[i]
    
    dist_matrix[t1, t2] <- d
    cor_matrix[t1, t2] <- c
    
    if (t1 != t2) {
      dist_matrix[t2, t1] <- d  # Symmetric fill
      cor_matrix[t2, t1] <- c   # Symmetric fill
    }
  }
  
  # Attach attributes
  attr(summary_mat, "distance_matrix") <- dist_matrix
  attr(summary_mat, "correlation_matrix") <- cor_matrix
  attr(summary_mat, "n_genes") <- nrow(mat)
  attr(summary_mat, "n_samples") <- ncol(mat)
  attr(summary_mat, "n_tissues") <- n_tissues
  attr(summary_mat, "distance_method") <- "euclidean"
  attr(summary_mat, "correlation_method") <- "pearson"
  
  return(summary_mat)
}

# Compute pairwise distances among all sub-tissue types
result <- calculate_pairwise_tissue_distances(t(data_filtered), "sub_category", group_labels_filtered)
result$tissue1_major_tissue<-group_labels_filtered_unique[result$tissue1,]$major_category
result$tissue2_major_tissue<-group_labels_filtered_unique[result$tissue2,]$major_category  
result$is_within_major_tissue<-ifelse(result$tissue1_major_tissue==result$tissue2_major_tissue,"TRUE","FALSE")
write.csv(result, "output/Distance_20251208_NOremoving_samples_with_cors_v6_step2limma.csv",row.names = F)

result2<-calculate_pairwise_tissue_distances(t(data_filtered), "major_category", group_labels_filtered)
write.csv(result2, "output/Distance_20251208_NOremoving_samples_with_cor_major_tissue_v6_step2limma.csv",row.names = F)

# Helper: get only between-tissue distances (exclude within-tissue)
get_between_tissue_distances <- function(result) {
  result[!result$is_within_tissue, ]
}

# Helper: get only within-tissue distances
get_within_tissue_distances <- function(result) {
  result[result$is_within_tissue, ]
}

# Helper: print distance matrix
print_distance_matrix <- function(result, digits = 2) {
  dist_mat <- attr(result, "distance_matrix")
  if (!is.null(dist_mat)) {
    cat("\nDistance matrix:\n")
    print(round(dist_mat, digits))
  }
}

# Helper: get the most similar tissue pairs
get_closest_tissues <- function(result, n = 5, exclude_within = TRUE) {
  if (exclude_within) {
    result <- result[!result$is_within_tissue, ]
  }
  head(result[, c("tissue1", "tissue2", "mean_distance", "median_distance")], n)
}

# Helper: plot boxplot-style comparisons

result$pair <- ifelse(result$is_within_tissue,
                          paste0(result$tissue1, " (within)"),
                          paste0(result$tissue1, "-", result$tissue2))
    
    # Build boxplot-style summary data
plot_data <- data.frame(
      pair = result$pair,
      type = ifelse(as.logical(result$is_within_tissue), "Within", "Between"),
      mean = result$mean_distance,
      median = result$median_distance,
      q25 = result$q25_distance,
      q75 = result$q75_distance,
      min = result$min_distance,
      max = result$max_distance
 )


# Plot
p <- ggplot(plot_data, aes(x = reorder(pair, median), color = type)) +
      geom_pointrange(aes(y = median, ymin = q25, ymax = q75), 
                      position = position_dodge(0.5)) +
      geom_point(aes(y = median), size = 3) +
      coord_flip() +
      scale_color_manual(values = c("Between" = "#CCCCCC", "Within" = "black")) +  # Transparent gray and black
      labs(title = "Tissue Distance Comparison",
           x = "Tissue Pair",
           y = "",  # Remove y-axis label
           color = "Type") +
      theme_classic() +  # White background
      theme(legend.position = "top",
            axis.text.y =  element_blank(),  # Enlarge x-axis labels (due to coord_flip, this is the y-axis)
            axis.text.x = element_blank(),  # Remove axis tick text
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA))


pdf("output/all_distance_subtype_20251208_v6_step2limma_NOremoving.pdf",height=4,width=4)
print(p)
dev.off()


library(ggstatsplot)
library(palmerpenguins)
library(tidyverse)
library(readxl)

#read data
#example
# data("penguins", package = "palmerpenguins")
# penguins <- drop_na(penguins)
# plt <- ggbetweenstats(
#   data = penguins,
#   x = species,
#   y = bill_length_mm
# )
# plt

##read tissue heterogeneity data
# d1<-read_xlsx("E:/YUEL/code/Figure_S5_tissue_heterogeneity/TableS3_normal_autopsy_heterogeneity_specificity.xlsx",sheet=2)
# d2<-read_xlsx("E:/YUEL/code/Figure_S5_tissue_heterogeneity/TableS3_normal_autopsy_heterogeneity_specificity.xlsx",sheet=3)

# Detailed tissue labels
result$detailed_tissue<-ifelse(result$is_within_tissue, "within_detailed_tissue","between_detailed_tissue")
# Major tissue labels   
result$major_tissue<-ifelse(result$is_within_major_tissue, "within_major_tissue","between_major_tissue")
# Major tissue label
result2$major_tissue<-ifelse(result2$is_within_tissue, "within_major_tissue","between_major_tissue")

plt <- ggbetweenstats(
  data = result,
  x = detailed_tissue,
  y = median_distance
)
# plt

plt2<-ggbetweenstats(
  data = result,
  x = major_tissue,
  y = median_distance,
  # ggsignif.args    = list(textsize = 4, tip_length = 0.01),
  # p.adjust.method  = "bonferroni"
)
# plt2

plt3<-ggbetweenstats(
  data = result2,
  x = major_tissue,
  y = median_distance,
  outlier.tagging = FALSE,
  outlier.label = NULL
)
# plt3

plt4 <- ggbetweenstats(
  data = result,
  x = detailed_tissue,
  y = median_correlation
)
# plt4

plt5<-ggbetweenstats(
  data = result,
  x = major_tissue,
  y =  median_correlation
  # ggsignif.args    = list(textsize = 4, tip_length = 0.01),
  # p.adjust.method  = "bonferroni"
)

plt6<-ggbetweenstats(
  data = result2,
  x = major_tissue,
  y = median_correlation,
  outlier.tagging = FALSE,
  outlier.label = NULL
)
# plt6

# Combine plots
library(cowplot)
plt7 <- plot_grid(plt, plt2, plt3, plt4, plt5, plt6, ncol = 3)
# plt7

pdf("output/distance_correlation_ggstatplot_removing_samples_20251208_v6_step2_limma_NOremoving.pdf", width = 8 , height = 5)
plt7
dev.off()



# 
# ##### Plot hclust for major tissues
# # A more complete version, ensuring distances from each tissue to itself (0) are included
# library(tidyr)
# library(dplyr)
# 
# # Assume the original data
# df <- data.frame(tissue1=result2$tissue1,
#                  tissue2=result2$tissue2,
#                  distance=result2$median_distance)
# 
# # Get all unique tissue names
# all_tissues <- unique(c(df$tissue1, df$tissue2))
# 
# # Step 1: complete symmetry
# reverse_df <- df %>%
#   rename(tissue1 = tissue2, tissue2 = tissue1)
# 
# complete_df <- bind_rows(df, reverse_df) %>%
#   distinct(tissue1, tissue2, .keep_all = TRUE)
# 
# # Step 2: add the diagonal (distance from a tissue to itself is 0)
# diagonal_df <- data.frame(
#   tissue1 = all_tissues,
#   tissue2 = all_tissues,
#   distance = 0
# )
# 
# final_df <- bind_rows(complete_df, diagonal_df) %>%
#   distinct(tissue1, tissue2, .keep_all = TRUE)
# 
# # Step 3: convert to wide format
# distance_matrix <- final_df %>%
#   pivot_wider(
#     names_from = tissue2,
#     values_from = distance
#   ) %>%
#   rename(tissue = tissue1) %>%
#   # Ensure column order matches row order
#   select(tissue, all_of(all_tissues))
# 
# 
# # Optional: convert to matrix format
# matrix_format <- as.matrix(distance_matrix[, -1])
# rownames(matrix_format) <- distance_matrix$tissue
# 
# # Plot hclust
# 
# distance_matrix_based_hclust <- function(
#   Xmed, # distance matrix (tissue_type x tissue_type)
#   main = NULL,
#   label_converter = NULL,  # can be a *named character vector* or a *function*
#   dist_cut = NULL,         # distance cutoff for clustering (replaces rho_min)
#   linkage = "complete",    # hclust linkage method
#   plot = TRUE
# ){
#   if (!requireNamespace("dendextend", quietly = TRUE))
#     stop("Package 'dendextend' is required.")
#   if (!requireNamespace("colorspace", quietly = TRUE))
#     stop("Package 'colorspace' is required.")
#   require(dendextend)
#   
#   # Ensure the input is a distance-matrix format
#   # If input is a data frame, convert it to a matrix
#   if (is.data.frame(Xmed)) {
#     # Assume the first column is tissue names and the rest are distance values
#     rownames(Xmed) <- Xmed[, 1]
#     Xmed <- as.matrix(Xmed[, -1])
#   }
#   
#   # Ensure it is a symmetric distance matrix
#   if (!isSymmetric(Xmed)) {
#     warning("Distance matrix is not symmetric. Making it symmetric by taking the upper triangle.")
#     Xmed[lower.tri(Xmed)] <- t(Xmed)[lower.tri(Xmed)]
#   }
#   
#   # Convert to a dist object
#   res.dist <- as.dist(Xmed)
#   
#   # Perform hierarchical clustering
#   res.hcls <- stats::hclust(res.dist, method = linkage)
#   res.dend <- as.dendrogram(res.hcls)
#   
#   # Determine cut height and number of clusters
#   if (is.null(dist_cut)) {
#     # If no cut height is provided, use the mean distance
#     dist_cut <- mean(res.dist, na.rm = TRUE)
#   }
#   k <- length(unique(stats::cutree(res.hcls, h = dist_cut)))
#   
#   # Label conversion
#   old_lab <- labels(res.dend)
#   if (is.null(label_converter)) {
#     new_lab <- old_lab
#   } else if (is.function(label_converter)) {
#     new_lab <- label_converter(old_lab)
#   } else {
#     # Assume a named character vector; for any NA, keep original labels
#     new_lab <- unname(label_converter[old_lab])
#     new_lab[is.na(new_lab)] <- old_lab
#   }
#   
#   # Generate colors
#   set.seed(1)
#   branch_colors <- sample(colorspace::rainbow_hcl(max(10, k), c = 70, l = 50))
#   
#   # Plotting
#   if (plot) {
#     # Increase the right margin to fit labels
#     op <- par(no.readonly = TRUE)
#     on.exit(par(op), add = TRUE)
#     cur <- par("mar") # c(bottom, left, top, right)
#     par(mar = c(cur[1], cur[2], cur[3], 8))
#     par(xpd = NA)
#     
#     # Create a colored dendrogram
#     p_dend <- res.dend %>%
#       set("labels", new_lab) %>%
#       set("labels_cex", 0.7) %>%
#       set("branches_k_color", value = branch_colors, k = k) %>%
#       color_branches(h = dist_cut, col = branch_colors[1:k])
#     
#     # Draw dendrogram
#     plot(p_dend, horiz = TRUE, axes = TRUE,
#          xlab = "Euclidean Distance", main = main)
#     
#     # Add cutoff line
#     abline(v = dist_cut, lty = 2, col = "red")
#     text(x = dist_cut, y = -0.5, labels = paste("Cutoff =", round(dist_cut, 3)), 
#          pos = 4, cex = 0.8, col = "red")
#     
#     # Optional: add colored bars indicating cluster assignments
#     cluster_assignments <- cutree(res.hcls, h = dist_cut)
#     colored_bars(colors = branch_colors[cluster_assignments], 
#                  dend = p_dend, 
#                  rowLabels = "Cluster",
#                  horiz = TRUE,
#                  cex.rowLabels = 0.7)
#   }
#   
#   # Return results
#   return(list(
#     res = list(
#       distance_matrix = Xmed,
#       dist_object = res.dist,
#       hclust_result = res.hcls,
#       dendrogram = res.dend,
#       leaf_labels = new_lab,
#       cluster_assignments = cutree(res.hcls, h = dist_cut)
#     ),
#     params = list(
#       dist_cut = dist_cut,
#       k = k,
#       linkage = linkage
#     )
#   ))
# }
# 
# 
# 
# # pdf('output/TPHP_normal_tissue_cluster_20251026.pdf', width = 3, height = 7.5)
# # hclust_list <- lapply(names(mat_med_list) %>% setNames(names(mat_med_list)), function(nm){
# #   message('Analysis of ', nm, '...\n')
# #   Xmed <- mat_med_list[[nm]]
# #   spearman_based_hclust(Xmed, main = nm, label_converter = label_converter)
# # })
# # graphics.off()
# 
# # pdf('output/TPHP_normal_tissue_cluster_20251026.pdf', width = 3, height = 7.5)
# pdf("output/tissue_clster_V6_step2_limma.pdf")
# res<-distance_matrix_based_hclust(Xmed = matrix_format,
#                                main = "Tissue Clustering based on Euclidean Distance",
#                                dist_cut = 300,
#                                linkage = "complete",plot = T)
# dev.off()
# 
# 
# ########### Plot normal tSNE
# ####### PCA-based t-SNE --------
# library(RSpectra)
# library(Rtsne)
# normal_tsne_list <- lapply(names(mat_list), function(nm){
#   message('PCA-based t-SNE of matrix ', nm, '...\n')
#   X <- mat_list[[nm]]
#   dim(X) # 12754 560
#   DF <- X %>% t()
#   M <- apply(X, 1, function(v) {
#     (v - mean(v, na.rm = T)) / sd(v, na.rm = T)
#   })
#   dim(M) # 506 12754
#   m1 <- prcomp(M)
#   # summary(m1)
#   # m1$x; # PCA scores
#   # m1$rotation # eigenvectors
#   # m1$sdev ^ 2 # eigenvalues
#   
#   # plot(cumsum(m1$sdev ^ 2)/sum(m1$sdev ^ 2)) # 500+
#   # (cumsum(m1$sdev ^ 2)/sum(m1$sdev ^ 2))[94] # 0.7687
#   # which(cumsum(m1$sdev ^ 2)/sum(m1$sdev ^ 2) <= 0.763) %>% max() # 90
#   # which(cumsum(m1$sdev ^ 2)/sum(m1$sdev ^ 2) <= 0.9) %>% max() # 246
#   # (cumsum(m1$sdev ^ 2)/sum(m1$sdev ^ 2))[420] # 0.97
#   
#   k <- which(cumsum(m1$sdev ^ 2)/sum(m1$sdev ^ 2) <= 0.9) %>% max()
#   
#   y <- as.factor(colnames(DF))
#   pca <- prcomp_svds(DF, retx = T, center = F, scale. = T, tol = NULL, k = k)
#   
#   
#   # calculate tSNE
#   set.seed(2020)
#   res <- Rtsne(pca$x,pca=F,max_iter=2500,theta=0,verbose=T, perplexity = 20
#                # perplexity must be <= floor((nrow(pca$x) - 1) / 3)
#   )
#   
#   tsne <- res$Y %>% as.data.frame
#   colnames(tsne) <- c('tsne1', 'tsne2')
#   hc <- hclust(dist(tsne))
#   
#   res.hc <- tsne %>% 
#     mutate(hclust = factor(cutree(hc, 20)),
#            hclust = str_c('Cluster ', hclust),
#            FileName = colnames(X)) %>% 
#     inner_join(info1) %>% 
#     set_rownames(.$FileName)
#   
#   hc.cent <- res.hc %>%
#     group_by(hclust) %>%
#     dplyr::select(tsne1, tsne2) %>%
#     summarize_all(mean)
#   
#   return(list(res.pca = pca,
#               pca.k = k,
#               res.tsne = res,
#               hc = hc,
#               res.hc = res.hc,
#               hc.cent = hc.cent))
# })
# 
# names(normal_tsne_list) <- names(mat_list)
# saveRDS(normal_tsne_list, 'normal_tsne_list.rds')
# 
# 
# for(nm in names(normal_tsne_list)){
#   ntsne <- normal_tsne_list[[nm]]
#   res.hc <- ntsne$res.hc
#   my_colors <- c('#368FC6', '#55BC6D', '#D8894E', '#7C0823', '#9D7DDB', '#BCBD22', '#CF4E9C', '#000000', '#888888')
#   my_shapes <- c(
#     1, 2, 6, 3, 4, 15, 17, 18, 19
#   )
#   
#   shape_color_pairs <- expand.grid(my_shapes, my_colors)
#   
#   # anatomical classification
#   p1 <- ggplot(res.hc) + 
#     geom_point(aes(x = tsne1, y = tsne2, shape = anatomical_classification, color = anatomical_classification), size = 2)+
#     scale_shape_manual(values = shape_color_pairs[, 1])+
#     scale_color_manual(values = as.character(shape_color_pairs[, 2]))+
#     theme_bw() +
#     theme(text = element_text(size = 10))
#   p2 <- p1 + theme(legend.position = 'none')
#   p <- ggarrange(p1, p2, nrow = 1, ncol = 2)
#   ggsave(str_c('output/20251025_57normalAll_tissue_pca_based_tsne_class_', nm, '.pdf'), p, width = 12, height = 6)
#   
#   # patient
#   p1 <- ggplot(res.hc) + 
#     geom_point(aes(x = tsne1, y = tsne2, shape = patient_ID, color = patient_ID), size = 2)+
#     scale_shape_manual(values = shape_color_pairs[, 1])+
#     scale_color_manual(values = as.character(shape_color_pairs[, 2]))+
#     theme_bw() +
#     theme(text = element_text(size = 10))
#   p2 <- p1 + theme(legend.position = 'none')
#   p <- ggarrange(p1, p2, nrow = 1, ncol = 2)
#   ggsave(str_c('output/20251025_57normalAll_tissue_pca_based_tsne_patient_', nm, '.pdf'), p, width = 12, height = 6)
#   
#   # Batchg
#   res.hc %<>%
#     arrange(date, instrument) %>% 
#     mutate(Batch = str_c(instrument, str_sub(date, 1, 6)),
#            Batch = factor(Batch, levels = unique(Batch)))
#   p1 <- ggplot(res.hc) + 
#     geom_point(aes(x = tsne1, y = tsne2, color = Batch), size = 2)+
#     scale_color_manual(values = RColorBrewer::brewer.pal(n = length(unique(res.hc$Batch)), name = 'Spectral'))+
#     theme_bw() +
#     theme(text = element_text(size = 10))
#   p2 <- p1 + theme(legend.position = 'none')
#   p <- ggarrange(p1, p2, nrow = 1, ncol = 2)
#   ggsave(str_c('output/20251027_57normalAll_tissue_pca_based_tsne_Batch_', nm, '.pdf'), p, width = 12, height = 6)
# }
# 
# rio::export(lapply(normal_tsne_list, function(X) X$res.hc),
#             'output/20251027_57normalAll_tissue_pca_based_tsne.xlsx')
# 
# 
# 
