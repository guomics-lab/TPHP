# Load required libraries
library(tidyverse)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(limma)
library(VennDiagram)
library(factoextra)
library(ggpubr)
library(magrittr)
library(dplyr)
library(purrr)
library(statmod)

# Set working directory and create output folders
setwd(r"(\\172.16.13.136\tphp\code\5_2_compare_with_previous)")
# dir.create("output", showWarnings = FALSE)
dir.create("output/plots", showWarnings = FALSE)
dir.create("output/tables", showWarnings = FALSE)

# Load sample information files

# =============================================================================
# 1. DATA LOADING AND PREPROCESSING
# =============================================================================

# Load sample information files
sample_info_proprietary <- read.csv("output/20251210_tphp_mapping_previous_study_info.xlsx", stringsAsFactors = FALSE,,row.names = 1)
sample_info_proprietary$tissue_type<-sample_info_proprietary$anatomical_classification
sample_info_proprietary<-sample_info_proprietary[-grep("TPHP_DIA_002|TPHP_DIA_2_rep",sample_info_proprietary$FileName),]
dim(sample_info_proprietary)
sample_info_proprietary$sample_id<-sample_info_proprietary$FileName
sample_info_published1 <- read.csv("output/2020cell_individual_age_sex.csv", stringsAsFactors = FALSE,row.names = 1)
sample_info_published1$sample_id<-sample_info_published1$GTEx.Sample_ID
# Note: sample_info_published2 not needed as expr_published2 is already tissue means

# standardize sample info   
# Function to standardize sample info columns
standardize_sample_info <- function(sample_info, dataset_name) {
  cat(paste("Standardizing", dataset_name, "sample info...\n"))
  cat(paste("Original columns:", paste(colnames(sample_info), collapse = ", "), "\n"))
  
  # Ensure required columns exist
  required_cols <- c("sample_id", "tissue_type")
  
  # Check for age column (Age or age)
  age_col <- NULL
  if("Age" %in% colnames(sample_info)) {
    age_col <- "Age"
  } else if("age" %in% colnames(sample_info)) {
    age_col <- "age"
  }
  
  if(is.null(age_col)) {
    stop(paste("No age column found in", dataset_name, "sample info"))
  }
  
  # Create standardized data frame with consistent column names
  standardized_info <- data.frame(
    sample_id = sample_info$sample_id,
    tissue_type = sample_info$tissue_type,
    age_original = sample_info[[age_col]],
    stringsAsFactors = FALSE
  )
  
  # Add gender if available
  if("gender" %in% colnames(sample_info)) {
    standardized_info$gender <- sample_info$gender
  } else if("Gender" %in% colnames(sample_info)) {
    standardized_info$gender <- sample_info$Gender
  } else if("sex" %in% colnames(sample_info)) {
    standardized_info$gender <- sample_info$sex
  } else {
    standardized_info$gender <- "Unknown"
  }
  
  # Convert age to numeric
  standardized_info$age <- convert_age_to_numeric(standardized_info$age_original, dataset_name)
  
  cat(paste("Standardized columns:", paste(colnames(standardized_info), collapse = ", "), "\n"))
  cat(paste("Age range:", min(standardized_info$age, na.rm = TRUE), "to", 
            max(standardized_info$age, na.rm = TRUE), "\n\n"))
  
  return(standardized_info)
}

# Function to convert age strings to numeric values
convert_age_to_numeric <- function(age_values, dataset_name) {
  # If already numeric, return as is
  if(is.numeric(age_values)) {
    return(age_values)
  }
  
  # Convert to character for string processing
  age_values <- as.character(age_values)
  
  # Handle age ranges like "30-39", "40-49", etc.
  numeric_ages <- sapply(age_values, function(x) {
    if(is.na(x) || x == "" || x == "Unknown") {
      return(NA)
    }
    
    # Check if it's a range (contains hyphen)
    if(grepl("-", x)) {
      # Extract numbers from range and take the midpoint
      range_parts <- strsplit(x, "-")[[1]]
      if(length(range_parts) == 2) {
        start_age <- as.numeric(range_parts[1])
        end_age <- as.numeric(range_parts[2])
        if(!is.na(start_age) && !is.na(end_age)) {
          return((start_age + end_age) / 2)  # Use midpoint
        }
      }
    }
    
    # Try to extract number directly
    numeric_value <- as.numeric(x)
    if(!is.na(numeric_value)) {
      return(numeric_value)
    }
    
    # If all else fails, return NA
    return(NA)
  })
  
  cat(paste("Age conversion for", dataset_name, ":\n"))
  cat(paste("Original format examples:", paste(head(unique(age_values), 5), collapse = ", "), "\n"))
  cat(paste("Converted range:", min(numeric_ages, na.rm = TRUE), "to", 
            max(numeric_ages, na.rm = TRUE), "\n"))
  
  # Display age group distribution based on Published2 groupings
  age_groups_pub2 <- cut(numeric_ages, 
                        breaks = c(0, 30, 45, 68, 100),
                        labels = c("Young_14-30", "Adult_30-45", "Middle_45-68", "Older_68-93"),
                        include.lowest = TRUE, right = FALSE)
  cat("Age group distribution (Published2 compatible):\n")
  print(table(age_groups_pub2, useNA = "ifany"))
  
  return(numeric_ages)
}
# Standardize all sample info files
# Standardize all sample info files
sample_info_proprietary <- standardize_sample_info(sample_info_proprietary, "Proprietary")
sample_info_published1 <- standardize_sample_info(sample_info_published1, "Published1")
head(sample_info_proprietary)
head(sample_info_published1)
sample_info_proprietary$gender<-str_to_lower(sample_info_proprietary$gender)

# Load proteomics expression data
expr_proprietary0 <- read.csv("output/tphp_protein_normalized_pm_10169_ensemblID.csv", check.names = FALSE)
expr_published10 <- read.csv("output/cell2020_protein_normalized_pm_ensemblID.csv",check.names = FALSE)  # ratio matrix [0.002, 442]
expr_published20 <- read.csv("output/liugh2025_protein_normalized_pm_ensemblID.csv",check.names = FALSE)  # tissue mean matrix

expr_proprietary <- apply(expr_proprietary0[,-c(1:2)],c(1,2),as.numeric)
rownames(expr_proprietary) <- expr_proprietary0$`Gene stable ID` 
expr_published1 <- apply(expr_published10[,-c(1:2)],c(1,2),as.numeric)
rownames(expr_published1) <- expr_published10$gene.id   
expr_published2 <-  apply(expr_published20[,-c(1:4)],c(1,2),as.numeric)
rownames(expr_published2) <- expr_published20$`Gene stable ID`

# 1.1 GENE DEDUPLICATION - Remove duplicate rownames, keep row with highest sum
# =============================================================================

remove_duplicate_genes <- function(expr_matrix) {
  cat(paste("Original matrix dimensions:", nrow(expr_matrix), "x", ncol(expr_matrix), "\n"))
  
  # Check for duplicated row names
  duplicated_genes <- rownames(expr_matrix)[duplicated(rownames(expr_matrix))]
  
  if(length(duplicated_genes) > 0) {
    cat(paste("Found", length(duplicated_genes), "duplicated gene names\n"))
    
    # Calculate row sums for all genes
    row_sums <- rowSums(expr_matrix, na.rm = TRUE)
    
    # Create a data frame with gene names and row sums
    gene_data <- data.frame(
      gene = rownames(expr_matrix),
      row_sum = row_sums,
      row_index = 1:nrow(expr_matrix),
      stringsAsFactors = FALSE
    )
    
    # For each gene, keep only the row with maximum sum
    gene_data <- gene_data %>%
      group_by(gene) %>%
      slice_max(row_sum, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    # Extract the filtered matrix
    filtered_matrix <- expr_matrix[gene_data$row_index, ]
    rownames(filtered_matrix) <- gene_data$gene
    
    cat(paste("After deduplication:", nrow(filtered_matrix), "x", ncol(filtered_matrix), "\n"))
    cat(paste("Removed", nrow(expr_matrix) - nrow(filtered_matrix), "duplicate rows\n\n"))
    
    return(filtered_matrix)
  } else {
    cat("No duplicated gene names found\n\n")
    return(expr_matrix)
  }
}

# Apply gene deduplication to all datasets
cat("Deduplicating proprietary dataset...\n")
expr_proprietary <- remove_duplicate_genes(expr_proprietary)

cat("Deduplicating published1 dataset...\n")
expr_published1 <- remove_duplicate_genes(expr_published1)

cat("Deduplicating published2 dataset...\n")
expr_published2 <- remove_duplicate_genes(expr_published2)
colnames(expr_published2)<-tolower(colnames(expr_published2))
colnames(expr_published2)<-gsub(pattern = "Lymph","lymph node",colnames(expr_published2))
colnames(expr_published2)<-gsub(pattern = "Adrenal","adrenal gland",colnames(expr_published2))
colnames(expr_published2)<-gsub(pattern = "Muscle","skeletal muscle",colnames(expr_published2))

# Define overlapping proteins and tissues
overlapping_proteins <- intersect(intersect(expr_proprietary0$`Gene stable ID`,expr_published10$gene.id),
                                   expr_published20$`Gene stable ID`)# 8398 proteins
overlapping_tissues <- intersect(intersect(sample_info_proprietary$tissue_type, sample_info_published1$tissue_type),
                                 colnames(expr_published2))# 6 tissue types
overlapping_tissues
# Alternative: If you have the lists in files
# overlapping_proteins <- read.table("overlapping_proteins.txt", stringsAsFactors = FALSE)$V1
# overlapping_tissues <- read.table("overlapping_tissues.txt", stringsAsFactors = FALSE)$V1

# Filter proprietary data
filter_intensity_data <- function(expr_data, sample_info, proteins, tissues) {
  # Filter samples by tissue type
#   tissues<-overlapping_tissues
#   sample_info<-sample_info_proprietary 
  valid_samples <- sample_info[sample_info$tissue_type %in% tissues, ]
#   dim(valid_samples)
#   expr_data<-expr_proprietary
#   View(expr_data[,1:10])
#   proteins<-overlapping_proteins                
  # Filter expression data
  expr_filtered <- expr_data[rownames(expr_data) %in% proteins, 
                            colnames(expr_data) %in% valid_samples$sample_id]
  dim(expr_filtered)
  # Ensure sample info matches expression data
  sample_info_filtered <- valid_samples[valid_samples$sample_id %in% colnames(expr_filtered), ]
  sample_info_filtered <- sample_info_filtered[match(colnames(expr_filtered), sample_info_filtered$sample_id), ]
  
  return(list(expression = expr_filtered, sample_info = sample_info_filtered))
}


# Filter published2 data (tissue mean matrix)
filter_published2_data <- function(expr_data, proteins, tissues) {
  # Filter proteins
  expr_filtered <- expr_data[rownames(expr_data) %in% proteins, ]
  
  # Filter tissues (assuming column names are tissue types)
  tissue_cols <- intersect(colnames(expr_filtered), tissues)
  expr_filtered <- expr_filtered[, tissue_cols]
  
  return(expr_filtered)
}

# Apply filtering
data_proprietary <- filter_intensity_data(expr_proprietary, sample_info_proprietary, 
                                           overlapping_proteins, overlapping_tissues)
data_published1 <- filter_intensity_data(expr_published1, sample_info_published1, 
                                         overlapping_proteins, overlapping_tissues)
data_published2 <- filter_published2_data(expr_published2, overlapping_proteins, overlapping_tissues)
dim(data_proprietary$expression)
dim(data_published1$expression)
dim(data_published2)
# > dim(data_proprietary$expression)
# [1] 8398   62
# > dim(data_published1$expression)
# [1] 8398   64
# > dim(data_published2)
# [1] 8398    7
# =============================================================================
# 2. DATA NORMALIZATION AND PREPROCESSING
# =============================================================================

cat("Data ranges before normalization:\n")
cat(paste("Proprietary (intensity):", round(min(data_proprietary$expression, na.rm = TRUE), 3), 
          "to", round(max(data_proprietary$expression, na.rm = TRUE), 3), "\n"))
cat(paste("Published1 (ratio):", round(min(data_published1$expression, na.rm = TRUE), 3), 
          "to", round(max(data_published1$expression, na.rm = TRUE), 3), "\n"))
cat(paste("Published2 (intensity means):", round(min(data_published2, na.rm = TRUE), 3), 
          "to", round(max(data_published2, na.rm = TRUE), 3), "\n\n"))

# Normalization strategy:
# 1. For intensity data (proprietary & published2): Log2 transform + Z-score normalization
# 2. For ratio data (published1): Log2 transform + median centering + Z-score normalization

# Function for Z-score normalization
zscore_normalize <- function(expr_matrix) {
  expr_matrix <- as.matrix(expr_matrix)
  expr_matrix[is.na(expr_matrix)]<-min(expr_matrix,na.rm = T) *0.5
  # Z-score normalization: (x - mean) / sd for each protein across samples
  normalized <- t(scale(t(expr_matrix), center = TRUE, scale = TRUE))
  
  # Handle proteins with zero variance
  # normalized[is.na(normalized)] <- 0
  
  rownames(normalized) <- rownames(expr_matrix)
  colnames(normalized) <- colnames(expr_matrix)
  
  return(normalized)
}

zscore_normalize_noImpu <- function(expr_matrix) {
  expr_matrix <- as.matrix(expr_matrix)
  
  # Z-score normalization: (x - mean) / sd for each protein across samples
  normalized <- t(scale(t(expr_matrix), center = TRUE, scale = TRUE))
  
  # Handle proteins with zero variance
#   normalized[is.na(normalized)] <- 0
  
  rownames(normalized) <- rownames(expr_matrix)
  colnames(normalized) <- colnames(expr_matrix)
  
  return(normalized)
}

# Function for median centering (for ratio data)
median_center <- function(expr_matrix) {
  expr_matrix <- as.matrix(expr_matrix)
  
  # Subtract median from each sample (column-wise centering)
  medians <- apply(expr_matrix, 2, median, na.rm = TRUE)
  centered <- sweep(expr_matrix, 2, medians, "-")
  
  return(centered)
}

# Process proprietary data (intensity matrix)
cat("Processing proprietary data (intensity matrix)...\n")
# Log2 transform intensity data
data_proprietary$expression <- log2(data_proprietary$expression + 1)
# Z-score normalization across proteins
data_proprietary$expression_normalized <- zscore_normalize(data_proprietary$expression)
head(data_proprietary$expression_normalized[,1:10])
# Process published1 data (ratio matrix)
cat("Processing published1 data (ratio matrix)...\n")
# Log2 transform ratio data
data_published1$expression <- 2^data_published1$expression   ##inf created if log2, use 2^
# Median center each sample (remove sample-specific bias)
data_published1$expression <- median_center(data_published1$expression)
# Z-score normalization across proteins
data_published1$expression_normalized <- zscore_normalize(data_published1$expression)

# Process published2 data (intensity means matrix)
cat("Processing published2 data (intensity means matrix)...\n")
# Log2 transform intensity means
# data_published2_log <- 2^data_published2
# Z-score normalization across proteins (row-wise)
data_published2_normalized <- zscore_normalize(data_published2)

cat("Data ranges after normalization (Z-scores):\n")
cat(paste("Proprietary:", round(min(data_proprietary$expression_normalized, na.rm = TRUE), 3), 
          "to", round(max(data_proprietary$expression_normalized, na.rm = TRUE), 3), "\n"))
cat(paste("Published1:", round(min(data_published1$expression_normalized, na.rm = TRUE), 3), 
          "to", round(max(data_published1$expression_normalized, na.rm = TRUE), 3), "\n"))
cat(paste("Published2:", round(min(data_published2_normalized, na.rm = TRUE), 3), 
          "to", round(max(data_published2_normalized, na.rm = TRUE), 3), "\n\n"))

# =============================================================================
# 3. CALCULATE TISSUE-SPECIFIC MEANS FOR COMPARISON
# =============================================================================

# Calculate tissue means for proprietary data
calc_tissue_means <- function(expr_data, sample_info) {
  tissue_means <- data.frame(row.names = rownames(expr_data))
  
  for(tissue in overlapping_tissues) {
    tissue_samples <- sample_info[sample_info$tissue_type == tissue, "sample_id"]
    tissue_samples <- intersect(tissue_samples, colnames(expr_data))
    
    if(length(tissue_samples) > 0) {
      if(length(tissue_samples) == 1) {
        tissue_means[[tissue]] <- expr_data[, tissue_samples]
      } else {
        tissue_means[[tissue]] <- rowMeans(expr_data[, tissue_samples], na.rm = TRUE)
      }
    }
  }
  return(tissue_means)
}

# Calculate tissue means for proprietary and published1 data using normalized data
proprietary_tissue_means <- calc_tissue_means(data_proprietary$expression_normalized, data_proprietary$sample_info)
published1_tissue_means <- calc_tissue_means(data_published1$expression_normalized, data_published1$sample_info)


# =============================================================================
# 4. DATASET SIMILARITY ANALYSIS BY TISSUE TYPE
# =============================================================================
analyze_tissue_similarity <- function(tissue_name) {
    # tissue_name <-"muscle"
  cat(paste("Analyzing tissue:", tissue_name, "\n"))
  
  # Get tissue-specific means (using normalized data)
  prop_mean <- proprietary_tissue_means[,tissue_name]
  pub1_mean <- published1_tissue_means[,tissue_name]
  names(prop_mean) <-rownames(proprietary_tissue_means)
  names(pub1_mean) <-rownames(published1_tissue_means)

  # Check if tissue exists in published2
  if(!tissue_name %in% colnames(data_published2_normalized)) {
    cat(paste("WARNING: Tissue", tissue_name, "not found in Published2 dataset\n"))
    pub2_mean <- NULL
  } else {
    pub2_mean <- data_published2_normalized[, tissue_name]
    names(pub2_mean) <-rownames(data_published2_normalized)
  }
  
  # Debug: Check data availability
  cat(paste("  Proprietary means available:", !is.null(prop_mean), "- length:", ifelse(is.null(prop_mean), 0, length(prop_mean)), "\n"))
  cat(paste("  Published1 means available:", !is.null(pub1_mean), "- length:", ifelse(is.null(pub1_mean), 0, length(pub1_mean)), "\n"))
  cat(paste("  Published2 means available:", !is.null(pub2_mean), "- length:", ifelse(is.null(pub2_mean), 0, length(pub2_mean)), "\n"))
  
  # Find common proteins for each comparison
  common_proteins_prop_pub1 <- c()
  common_proteins_prop_pub2 <- c()
  common_proteins_pub1_pub2 <- c()
  
  if(!is.null(prop_mean) && !is.null(pub1_mean)) {
    common_proteins_prop_pub1 <- intersect(names(prop_mean), names(pub1_mean))
    # Remove NA values
    valid_prop <- !is.na(prop_mean[common_proteins_prop_pub1])
    valid_pub1 <- !is.na(pub1_mean[common_proteins_prop_pub1])
    common_proteins_prop_pub1 <- common_proteins_prop_pub1[valid_prop & valid_pub1]
  }
  
  if(!is.null(prop_mean) && !is.null(pub2_mean)) {
    common_proteins_prop_pub2 <- intersect(names(prop_mean), names(pub2_mean))
    valid_prop <- !is.na(prop_mean[common_proteins_prop_pub2])
    valid_pub2 <- !is.na(pub2_mean[common_proteins_prop_pub2])
    common_proteins_prop_pub2 <- common_proteins_prop_pub2[valid_prop & valid_pub2]
  }
  
  if(!is.null(pub1_mean) && !is.null(pub2_mean)) {
    common_proteins_pub1_pub2 <- intersect(names(pub1_mean), names(pub2_mean))
    valid_pub1 <- !is.na(pub1_mean[common_proteins_pub1_pub2])
    valid_pub2 <- !is.na(pub2_mean[common_proteins_pub1_pub2])
    common_proteins_pub1_pub2 <- common_proteins_pub1_pub2[valid_pub1 & valid_pub2]
  }
  
  cat(paste("  Common proteins Prop-Pub1:", length(common_proteins_prop_pub1), "\n"))
  cat(paste("  Common proteins Prop-Pub2:", length(common_proteins_prop_pub2), "\n"))
  cat(paste("  Common proteins Pub1-Pub2:", length(common_proteins_pub1_pub2), "\n"))
  
  # Calculate correlations
  corr_prop_pub1 <- NA
  corr_prop_pub2 <- NA
  corr_pub1_pub2 <- NA
  
  if(length(common_proteins_prop_pub1) >= 10) {
    corr_prop_pub1 <- cor(prop_mean[common_proteins_prop_pub1], 
                         pub1_mean[common_proteins_prop_pub1], use = "complete.obs", method = "pearson")
  }
  
  if(length(common_proteins_prop_pub2) >= 10) {
    corr_prop_pub2 <- cor(prop_mean[common_proteins_prop_pub2], 
                         pub2_mean[common_proteins_prop_pub2], use = "complete.obs", method = "pearson")
  }
  
  if(length(common_proteins_pub1_pub2) >= 10) {
    corr_pub1_pub2 <- cor(pub1_mean[common_proteins_pub1_pub2], 
                         pub2_mean[common_proteins_pub1_pub2], use = "complete.obs", method = "pearson")
  }
  
  cat(paste("  Correlations: Prop-Pub1 =", round(corr_prop_pub1, 3), 
            "Prop-Pub2 =", round(corr_prop_pub2, 3), 
            "Pub1-Pub2 =", round(corr_pub1_pub2, 3), "\n\n"))
  
  # Create comparison plots only if sufficient data
  plots_list <- list()
  
  # Proprietary vs Published1
  if(!is.na(corr_prop_pub1) && length(common_proteins_prop_pub1) >= 10) {
    df1 <- data.frame(
      proprietary = prop_mean[common_proteins_prop_pub1],
      published1 = pub1_mean[common_proteins_prop_pub1]
    )
    
    plots_list$p1 <- ggplot(df1, aes(x = proprietary, y = published1)) +
      geom_point(alpha = 0.6, size = 1) +
      geom_smooth(method = "lm", color = "red", se = TRUE) +
      labs(title = paste(tissue_name, ": Proprietary vs Published1"),
           x = "Proprietary (Z-score normalized)",
           y = "Published1 (Z-score normalized)") +
      annotate("text", x = Inf, y = (-Inf), 
               label = paste("r =", round(corr_prop_pub1, 3), "\nn =", length(common_proteins_prop_pub1)),
               hjust = 1.1, vjust = -0.5, size = 3) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))
  }
  
  # Proprietary vs Published2
  if(!is.na(corr_prop_pub2) && length(common_proteins_prop_pub2) >= 10) {
    df2 <- data.frame(
      proprietary = prop_mean[common_proteins_prop_pub2],
      published2 = pub2_mean[common_proteins_prop_pub2]
    )
    
    plots_list$p2 <- ggplot(df2, aes(x = proprietary, y = published2)) +
      geom_point(alpha = 0.6, size = 1) +
      geom_smooth(method = "lm", color = "red", se = TRUE) +
      labs(title = paste(tissue_name, ": Proprietary vs Published2"),
           x = "Proprietary (Z-score normalized)",
           y = "Published2 (Z-score normalized)") +
      annotate("text", x = Inf, y = -Inf, 
               label = paste("r =", round(corr_prop_pub2, 3), "\nn =", length(common_proteins_prop_pub2)),
               hjust = 1.1, vjust = -0.5, size = 3) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))
  }
  
  # Published1 vs Published2
  if(!is.na(corr_pub1_pub2) && length(common_proteins_pub1_pub2) >= 10) {
    df3 <- data.frame(
      published1 = pub1_mean[common_proteins_pub1_pub2],
      published2 = pub2_mean[common_proteins_pub1_pub2]
    )
    
    plots_list$p3 <- ggplot(df3, aes(x = published1, y = published2)) +
      geom_point(alpha = 0.6, size = 1) +
      geom_smooth(method = "lm", color = "red", se = TRUE) +
      labs(title = paste(tissue_name, ": Published1 vs Published2"),
           x = "Published1 (Z-score normalized)",
           y = "Published2 (Z-score normalized)") +
      annotate("text", x = Inf, y = -Inf, 
               label = paste("r =", round(corr_pub1_pub2, 3), "\nn =", length(common_proteins_pub1_pub2)),
               hjust = 1.1, vjust = -0.5, size = 3) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))
  }
  
  # Arrange and save plots
  if(length(plots_list) > 0) {
    if(length(plots_list) == 3) {
      combined_plot <- grid.arrange(plots_list$p1, plots_list$p2, plots_list$p3, ncol = 3)
    } else if(length(plots_list) == 2) {
      combined_plot <- grid.arrange(plots_list[[1]], plots_list[[2]], ncol = 2)
    } else {
      combined_plot <- plots_list[[1]]
    }
    
    ggsave(paste0("results/plots/", tissue_name, "_dataset_similarity.png"), 
           combined_plot, width = 15, height = 5, dpi = 300)
  }
  
  return(list(
    tissue = tissue_name,
    n_proteins_prop_pub1 = length(common_proteins_prop_pub1),
    n_proteins_prop_pub2 = length(common_proteins_prop_pub2),
    n_proteins_pub1_pub2 = length(common_proteins_pub1_pub2),
    corr_prop_pub1 = corr_prop_pub1,
    corr_prop_pub2 = corr_prop_pub2,
    corr_pub1_pub2 = corr_pub1_pub2
  ))
}

# Run similarity analysis for all tissues
similarity_results <- map(overlapping_tissues, analyze_tissue_similarity)
names(similarity_results) <- overlapping_tissues



# Create summary table
correlation_summary <- data.frame(
  Tissue = overlapping_tissues,
  N_Proteins_Prop_Pub1 = map_dbl(similarity_results, ~.x$n_proteins_prop_pub1),
  N_Proteins_Prop_Pub2 = map_dbl(similarity_results, ~.x$n_proteins_prop_pub2),
  N_Proteins_Pub1_Pub2 = map_dbl(similarity_results, ~.x$n_proteins_pub1_pub2),
  Corr_Prop_Pub1 = map_dbl(similarity_results, ~ifelse(is.na(.x$corr_prop_pub1), NA, .x$corr_prop_pub1)),
  Corr_Prop_Pub2 = map_dbl(similarity_results, ~ifelse(is.na(.x$corr_prop_pub2), NA, .x$corr_prop_pub2)),
  Corr_Pub1_Pub2 = map_dbl(similarity_results, ~ifelse(is.na(.x$corr_pub1_pub2), NA, .x$corr_pub1_pub2))
)

write.csv(correlation_summary, "output/tables/dataset_correlations_summary_pearson_20250822.csv", row.names = FALSE)

###########################without imputation############################

# data_proprietary$expression_normalized_noImp <- zscore_normalize_noImpu(data_proprietary$expression)
# head(data_proprietary$expression_normalized_noImp[,1:10])

# # Process published1 data (ratio matrix)
# cat("Processing published1 data (ratio matrix)...\n")
# # Median center each sample (remove sample-specific bias)
# # data_published1$expression <- median_center(data_published1$expression)
# # Z-score normalization across proteins
# data_published1$expression_normalized_noImp <- zscore_normalize_noImpu(data_published1$expression)


# # Process published2 data (intensity means matrix)
# cat("Processing published2 data (intensity means matrix)...\n")
# # Log2 transform intensity means
# data_published2_log <- log2(data_published2 + 1)
# # Z-score normalization across proteins (row-wise)
# data_published2_normalized <- zscore_normalize_noImpu(data_published2_log)

# cat("Data ranges after normalization (Z-scores):\n")
# cat(paste("Proprietary:", round(min(data_proprietary$expression_normalized_noImp, na.rm = TRUE), 3), 
#           "to", round(max(data_proprietary$expression_normalized_noImp, na.rm = TRUE), 3), "\n"))
# cat(paste("Published1:", round(min(data_published1$expression_normalized_noImp, na.rm = TRUE), 3), 
#           "to", round(max(data_published1$expression_normalized_noImp, na.rm = TRUE), 3), "\n"))
# cat(paste("Published2:", round(min(data_published2_normalized, na.rm = TRUE), 3), 
#           "to", round(max(data_published2_normalized, na.rm = TRUE), 3), "\n\n"))
# #Calculate tissue means for proprietary and published1 data using normalized data
# proprietary_tissue_means <- calc_tissue_means(data_proprietary$expression_normalized_noImp, data_proprietary$sample_info)%>%as.matrix()
# published1_tissue_means <- calc_tissue_means(data_published1$expression_normalized_noImp, data_published1$sample_info)%>%as.matrix()
# head(proprietary_tissue_means)
# head(data_published2_normalized)
# sum(is.na(data_published2_normalized))
# proprietary_tissue_means[is.nan(proprietary_tissue_means)] <- NA
# published1_tissue_means[is.nan(published1_tissue_means)] <- NA
# similarity_results <- map(overlapping_tissues, analyze_tissue_similarity)
# names(similarity_results) <- overlapping_tissues



# Create summary table
# correlation_summary <- data.frame(
#   Tissue = overlapping_tissues,
#   N_Proteins_Prop_Pub1 = map_dbl(similarity_results, ~.x$n_proteins_prop_pub1),
#   N_Proteins_Prop_Pub2 = map_dbl(similarity_results, ~.x$n_proteins_prop_pub2),
#   N_Proteins_Pub1_Pub2 = map_dbl(similarity_results, ~.x$n_proteins_pub1_pub2),
#   Corr_Prop_Pub1 = map_dbl(similarity_results, ~ifelse(is.na(.x$corr_prop_pub1), NA, .x$corr_prop_pub1)),
#   Corr_Prop_Pub2 = map_dbl(similarity_results, ~ifelse(is.na(.x$corr_prop_pub2), NA, .x$corr_prop_pub2)),
#   Corr_Pub1_Pub2 = map_dbl(similarity_results, ~ifelse(is.na(.x$corr_pub1_pub2), NA, .x$corr_pub1_pub2))
# )

# write.csv(correlation_summary, "output/tables/dataset_correlations_summary_spearman_no_imp_complete_obs.csv", row.names = FALSE)


###########################imputated, age-matched############################
data_proprietary$age_matched_sample_info <- data_proprietary$sample_info %>%
  filter(!is.na(age) & age >= 50 & age <= 69) 
dim(data_proprietary$age_matched_sample_info )
data_proprietary$age_matched_sample_matrix <- data_proprietary$expression_normalized[,data_proprietary$age_matched_sample_info$sample_id]
dim(data_proprietary$age_matched_sample_matrix )
head(data_proprietary$age_matched_sample_matrix[,1:10])

data_published1$age_matched_sample_info <- data_published1$sample_info %>%
  filter(!is.na(age) & age >= 50 & age <= 69)
dim(data_published1$age_matched_sample_info )
data_published1$age_matched_sample_matrix <- data_published1$expression_normalized[,data_published1$age_matched_sample_info$sample_id]  
head(data_published1$age_matched_sample_matrix[,1:10])

proprietary_tissue_means <- calc_tissue_means(data_proprietary$age_matched_sample_matrix, data_proprietary$age_matched_sample_info)%>%as.matrix()
published1_tissue_means <- calc_tissue_means(data_published1$age_matched_sample_matrix, data_published1$age_matched_sample_info)%>%as.matrix()
head(proprietary_tissue_means)
head(data_published2_normalized)
# sum(is.na(data_published2_normalized))
# proprietary_tissue_means[is.nan(proprietary_tissue_means)] <- NA
# published1_tissue_means[is.nan(published1_tissue_means)] <- NA
similarity_results <- map(overlapping_tissues, analyze_tissue_similarity)
names(similarity_results) <- overlapping_tissues

# Create summary table
correlation_summary <- data.frame(
  Tissue = overlapping_tissues,
  N_Proteins_Prop_Pub1 = map_dbl(similarity_results, ~.x$n_proteins_prop_pub1),
  N_Proteins_Prop_Pub2 = map_dbl(similarity_results, ~.x$n_proteins_prop_pub2),
  N_Proteins_Pub1_Pub2 = map_dbl(similarity_results, ~.x$n_proteins_pub1_pub2),
  Corr_Prop_Pub1 = map_dbl(similarity_results, ~ifelse(is.na(.x$corr_prop_pub1), NA, .x$corr_prop_pub1)),
  Corr_Prop_Pub2 = map_dbl(similarity_results, ~ifelse(is.na(.x$corr_prop_pub2), NA, .x$corr_prop_pub2)),
  Corr_Pub1_Pub2 = map_dbl(similarity_results, ~ifelse(is.na(.x$corr_pub1_pub2), NA, .x$corr_pub1_pub2))
)

write.csv(correlation_summary, "output/tables/dataset_correlations_summary_pearson_imp_age_matched_20250822.csv", row.names = FALSE)

# Data visualization    
library(ggplot2)
library(tidyr)

# Create data
# data <- data.frame(
#   Tissue = c("muscle", "heart", "liver", "lung", "skin", "spleen", "adrenal gland", "pancreas"),
#   Corr_Prop_Pub1 = c(0.484754, 0.628893, 0.624479, 0.633895, 0.493675, 0.392046, 0.48181, 0.320341),
#   Corr_Prop_Pub1_age_matched = c(0.467568, 0.611374, 0.574954, 0.617598, 0.40118, 0.368458, 0.449298, 0.304392),
#   Corr_Prop_Pub2 = c(0.374592, 0.499011, 0.488376, 0.43302, 0.207761, 0.218186, 0.334601, 0.265759),
#   Corr_Pub1_Pub2 = c(0.574162, 0.644819, 0.615748, 0.541057, 0.440689, 0.595942, 0.51635, 0.502963)
# )
cor1<-read.csv("output/tables/dataset_correlations_summary_pearson_20250822.csv")
cor2<-read.csv("output/tables/dataset_correlations_summary_pearson_imp_age_matched_20250822.csv")
data<-cbind(cor1[,c("Tissue", "Corr_Prop_Pub1","Corr_Prop_Pub2","Corr_Pub1_Pub2")],
            Corr_Prop_Pub1_age_matched = cor2[,c("Corr_Prop_Pub1")])
data
# Reshape data to long format
data_long <- pivot_longer(data, cols = -Tissue, names_to = "Correlation_Type", values_to = "Value")
cols02<-c("#005496","#e8c559","#ea9c9d","#606f8a")
pdf("output/plots/tissue_correlation_comparison.pdf", width = 10, height = 5)
# 创建分组柱状图
ggplot(data_long, aes(x = Tissue, y = Value, fill = Correlation_Type)) +
  geom_bar( width= .5 , stat='identity', position=position_dodge( .7 )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = cols02)+ 
  labs(y = "Correlation Value", x = "Tissue", fill = "Correlation Type")+
  theme(  panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	     plot.title = element_text(size=15),
	    axis.line.x = element_line(color="black", size = 0.5),
	    axis.line.y = element_line(color="black", size = 0.5),
	    panel.background = element_blank())
dev.off()

##########################tissue heterogeneity############################
published1_tissue_means <- calc_tissue_means(data_published1$expression_normalized, data_published1$sample_info)
data_cor<-do.call(cbind,list(data_proprietary$expression_normalized, published1=published1_tissue_means,published2=data_published2_normalized))
dim(data_cor)                    
tail(colnames(data_cor),n=16)                    
corr_matrix <- cor(data_cor,  use = "complete.obs", method = "pearson")
write.csv(corr_matrix, "output/tables/dataset_correlation_matrix_pearson.csv", row.names = TRUE)
write.csv(data_proprietary$expression_normalized, "output/tables/dataset_proprietary_expression_normalized.csv", row.names = TRUE   )
write.csv(data_proprietary$sample_info, "output/tables/dataset_proprietary_sample_info.csv", row.names = FALSE)
# =============================================================================
# 5. AGE-RELATED PROTEIN EXPRESSION ANALYSIS (PROPRIETARY + PUBLISHED1 ONLY)
# =============================================================================

analyze_age_effects <- function(tissue_name) {
  cat(paste("Analyzing age effects in tissue:", tissue_name, "\n"))
  
  # Only use proprietary and published1 data (both have individual samples)
  prop_samples <- data_proprietary$sample_info[data_proprietary$sample_info$tissue_type == tissue_name, ]
  pub1_samples <- data_published1$sample_info[data_published1$sample_info$tissue_type == tissue_name, ]
  
  # Check if both datasets have samples for this tissue
  if(nrow(prop_samples) == 0 && nrow(pub1_samples) == 0) {
    return(list(tissue = tissue_name, status = "No samples available"))
  }
  
  # Ensure both datasets have the same columns before combining
  # Get common columns (should include: sample_id, tissue_type, age, gender, age_original)
  common_cols <- intersect(colnames(prop_samples), colnames(pub1_samples))
  
  if(!"age" %in% common_cols) {
    return(list(tissue = tissue_name, status = "Age column missing"))
  }
  
  # Select common columns for both datasets
  if(nrow(prop_samples) > 0) {
    prop_samples <- prop_samples[, common_cols]
    prop_samples$dataset <- "Proprietary"
  }
  
  if(nrow(pub1_samples) > 0) {
    pub1_samples <- pub1_samples[, common_cols]
    pub1_samples$dataset <- "Published1"
  }
  
  # Combine sample info (only if both have data)
  if(nrow(prop_samples) > 0 && nrow(pub1_samples) > 0) {
    combined_samples <- rbind(prop_samples, pub1_samples)
  } else if(nrow(prop_samples) > 0) {
    combined_samples <- prop_samples
  } else {
    combined_samples <- pub1_samples
  }
  
  # Get expression data (use normalized data for age analysis)
  prop_expr <- NULL
  pub1_expr <- NULL
  
  if(nrow(prop_samples) > 0) {
    available_prop_samples <- intersect(prop_samples$sample_id, colnames(data_proprietary$expression_normalized))
    if(length(available_prop_samples) > 0) {
      prop_expr <- data_proprietary$expression_normalized[, available_prop_samples, drop = FALSE]
    }
  }
  
  if(nrow(pub1_samples) > 0) {
    available_pub1_samples <- intersect(pub1_samples$sample_id, colnames(data_published1$expression_normalized))
    if(length(available_pub1_samples) > 0) {
      pub1_expr <- data_published1$expression_normalized[, available_pub1_samples, drop = FALSE]
    }
  }
  
  # Combine expression data
  if(!is.null(prop_expr) && !is.null(pub1_expr)) {
    combined_expr <- cbind(prop_expr, pub1_expr)
  } else if(!is.null(prop_expr)) {
    combined_expr <- prop_expr
  } else if(!is.null(pub1_expr)) {
    combined_expr <- pub1_expr
  } else {
    return(list(tissue = tissue_name, status = "No expression data available"))
  }
  
  # Ensure sample order matches between expression data and sample info
  available_samples <- intersect(combined_samples$sample_id, colnames(combined_expr))
  combined_samples <- combined_samples[combined_samples$sample_id %in% available_samples, ]
  combined_expr <- combined_expr[, available_samples]
  
  # Reorder to match
  combined_samples <- combined_samples[match(colnames(combined_expr), combined_samples$sample_id), ]
  
  # Remove samples with missing age information
  valid_age_idx <- !is.na(combined_samples$age)
  combined_samples <- combined_samples[valid_age_idx, ]
  combined_expr <- combined_expr[, valid_age_idx]
  
  if(nrow(combined_samples) < 5) {
    return(list(tissue = tissue_name, status = "Insufficient samples with age data"))
  }
  
  # Age categorization - consistent with Published2 dataset groupings
  combined_samples$age_group <- cut(combined_samples$age, 
                                   breaks = c(0, 30, 45, 68, 100),
                                   labels = c("Young_14-30", "Adult_30-45", "Middle_45-68", "Older_68-93"),
                                   include.lowest = TRUE, right = FALSE)
  
  # Check data quality before analysis
  cat(paste("  Age range:", min(combined_samples$age, na.rm = TRUE), "to", 
            max(combined_samples$age, na.rm = TRUE), "\n"))
  cat("  Age group distribution (consistent with Published2):\n")
  age_table <- table(combined_samples$age_group, useNA = "ifany")
  for(i in 1:length(age_table)) {
    cat(paste("   ", names(age_table)[i], ":", age_table[i], "samples\n"))
  }
  cat(paste("  Gender distribution:", paste(table(combined_samples$gender), collapse = ", "), "\n"))
  cat(paste("  Dataset distribution:", paste(table(combined_samples$dataset), collapse = ", "), "\n"))
  
  # Primary focus on age analysis with gender as covariate when appropriate
  # Create design matrix with age as main factor of interest
  design_factors <- "age"
  
  # Add dataset as covariate if multiple datasets present
  dataset_table <- table(combined_samples$dataset)
  dataset_table <- dataset_table[dataset_table > 0]
  
  if(length(dataset_table) >= 2) {
    design_factors <- c(design_factors, "dataset")
    cat("  Including dataset as covariate\n")
  }
  
  # Add gender as covariate if there's sufficient variation (but not as main factor)
  if("gender" %in% colnames(combined_samples)) {
    gender_table <- table(combined_samples$gender)
    gender_table <- gender_table[gender_table > 0]
    
    # More lenient criteria for gender as covariate (not main effect)
    if(length(gender_table) >= 2 && min(gender_table) >= 1 && nrow(combined_samples) >= 6) {
      design_factors <- c(design_factors, "gender")
      cat("  Including gender as covariate\n")
    } else {
      cat("  Excluding gender (insufficient samples for reliable covariate adjustment)\n")
    }
  }
  
  cat(paste("  Design formula: ~", paste(design_factors, collapse = " + "), "\n"))
  cat("  Note: Primary interest is age effects; other factors serve as covariates\n")
  
  # Create design formula
  design_formula <- as.formula(paste("~", paste(design_factors, collapse = " + ")))
  
  # Clean the combined_samples data
  combined_samples$gender <- factor(combined_samples$gender)
  combined_samples$dataset <- factor(combined_samples$dataset)
  
  # Remove unused factor levels
  combined_samples$gender <- droplevels(combined_samples$gender)
  combined_samples$dataset <- droplevels(combined_samples$dataset)
  
  # Create design matrix
  design_matrix <- model.matrix(design_formula, data = combined_samples)
  
  # Check for singularities in design matrix
  qr_design <- qr(design_matrix)
  if(qr_design$rank < ncol(design_matrix)) {
    cat("  WARNING: Design matrix is singular, using only linearly independent columns\n")
    # Keep only linearly independent columns
    keep_cols <- qr_design$pivot[1:qr_design$rank]
    design_matrix <- design_matrix[, keep_cols]
    
    # Ensure age coefficient is still present
    if(!"age" %in% colnames(design_matrix)) {
      cat("  ERROR: Age coefficient lost in singularity resolution\n")
      return(list(tissue = tissue_name, status = "Age coefficient unavailable due to design singularity"))
    }
  }
  
  cat(paste("  Final design matrix dimensions:", nrow(design_matrix), "x", ncol(design_matrix), "\n"))
  cat(paste("  Design matrix columns:", paste(colnames(design_matrix), collapse = ", "), "\n"))
  
  # Remove proteins with zero or near-zero variance
  protein_vars <- apply(combined_expr, 1, var, na.rm = TRUE)
  valid_proteins <- !is.na(protein_vars) & protein_vars > 1e-10
  
  if(sum(valid_proteins) < 100) {
    cat(paste("  WARNING: Only", sum(valid_proteins), "proteins with sufficient variance\n"))
    if(sum(valid_proteins) < 10) {
      return(list(tissue = tissue_name, status = "Insufficient protein variance"))
    }
  }
  
  combined_expr_filtered <- combined_expr[valid_proteins, ]
  cat(paste("  Using", nrow(combined_expr_filtered), "proteins for age analysis\n"))
  
  # Fit linear model with error handling
  tryCatch({
    fit <- lmFit(combined_expr_filtered, design_matrix)
    fit <- eBayes(fit, robust = TRUE, trend = TRUE)  # Use robust estimation with trend
    
    # Get age-associated proteins (age is always the main coefficient of interest)
    age_results <- topTable(fit, coef = "age", number = Inf, sort.by = "P")
    age_results$protein <- rownames(age_results)
    age_results <- age_results[, c("protein", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
    
    cat(paste("  Found", sum(age_results$adj.P.Val < 0.05, na.rm = TRUE), "age-associated proteins (FDR < 0.05)\n"))
    
  }, error = function(e) {
    cat(paste("  ERROR in linear modeling:", e$message, "\n"))
    return(list(tissue = tissue_name, status = paste("Model fitting failed:", e$message)))
  })
  
  # If we reach here, the model fitting was successful
  
  # Significant age-associated proteins (FDR < 0.05)
  sig_age_proteins <- age_results[age_results$adj.P.Val < 0.05, ]
  
  # Create age correlation plot for top proteins
  top_proteins <- head(sig_age_proteins$protein, min(20, nrow(sig_age_proteins)))
  
  if(length(top_proteins) > 0) {
    age_expr_data <- combined_expr[top_proteins, , drop = FALSE]
    age_corr_data <- data.frame(
      age = rep(combined_samples$age, each = length(top_proteins)),
      expression = as.vector(age_expr_data),
      protein = rep(top_proteins, ncol(age_expr_data)),
      dataset = rep(combined_samples$dataset, each = length(top_proteins))
    )
    
    # Age correlation plot
    age_plot <- ggplot(age_corr_data, aes(x = age, y = expression, color = dataset)) +
      geom_point(alpha = 0.6, size = 0.8) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
      facet_wrap(~protein, scales = "free_y", ncol = 4) +
            labs(title = paste(tissue_name, ": Age-associated Protein Expression (Top", length(top_proteins), "proteins)"),
           x = "Age (years)", y = "Normalized Expression (Z-score)")+
      scale_color_manual(values = c("Proprietary" = "#E31A1C", "Published1" = "#1F78B4")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(size = 8))
    
    ggsave(paste0("results/plots/", tissue_name, "_age_effects_top_proteins.png"), 
           age_plot, width = 16, height = 12, dpi = 300)
  }
  
  # Dataset-specific age effects analysis
  dataset_age_effects <- list()
  
  for(dataset in c("Proprietary", "Published1")) {
    dataset_samples <- combined_samples[combined_samples$dataset == dataset, ]
    
      if(nrow(dataset_samples) > 3) {
        dataset_expr <- combined_expr[, dataset_samples$sample_id]
        
        # Remove proteins with zero variance in this dataset
        protein_vars <- apply(dataset_expr, 1, var, na.rm = TRUE)
        dataset_expr <- dataset_expr[!is.na(protein_vars) & protein_vars > 0, ]
        
        if(nrow(dataset_expr) > 10) {  # Ensure enough proteins for analysis
          design_factors_ds <- "age"
          if("gender" %in% colnames(dataset_samples) && length(unique(dataset_samples$gender)) > 1) {
            design_factors_ds <- c(design_factors_ds, "gender")
          }
          
          design_formula_ds <- as.formula(paste("~", paste(design_factors_ds, collapse = " + ")))
          design_ds <- model.matrix(design_formula_ds, data = dataset_samples)
          
          fit_ds <- lmFit(dataset_expr, design_ds)
          fit_ds <- eBayes(fit_ds)
          
          age_results_ds <- topTable(fit_ds, coef = "age", number = Inf, sort.by = "P")
          age_results_ds$dataset <- dataset
          dataset_age_effects[[dataset]] <- age_results_ds
        }
      }
  }
  
  return(list(
    tissue = tissue_name,
    n_samples_total = nrow(combined_samples),
    n_samples_proprietary = sum(combined_samples$dataset == "Proprietary"),
    n_samples_published1 = sum(combined_samples$dataset == "Published1"),
    age_range = paste(min(combined_samples$age, na.rm = TRUE), "to", max(combined_samples$age, na.rm = TRUE)),
    age_group_distribution = table(combined_samples$age_group, useNA = "ifany"),
    combined_age_results = age_results,
    significant_proteins = sig_age_proteins,
    dataset_specific_results = dataset_age_effects,
    sample_info = combined_samples,
    status = "Complete"
  ))
}

# Run age analysis for all tissues
age_results <- map(overlapping_tissues, analyze_age_effects)
names(age_results) <- overlapping_tissues

# Create summary of age effects across tissues
age_summary <- data.frame(
  Tissue = overlapping_tissues,
  Status = map_chr(age_results, ~ifelse("status" %in% names(.x), .x$status, "Complete")),
  Total_Samples = map_dbl(age_results, ~ifelse("n_samples_total" %in% names(.x), .x$n_samples_total, 0)),
  Samples_Proprietary = map_dbl(age_results, ~ifelse("n_samples_proprietary" %in% names(.x), .x$n_samples_proprietary, 0)),
  Samples_Published1 = map_dbl(age_results, ~ifelse("n_samples_published1" %in% names(.x), .x$n_samples_published1, 0)),
  Significant_Age_Proteins = map_dbl(age_results, ~ifelse("significant_proteins" %in% names(.x), nrow(.x$significant_proteins), 0)),
  Total_Proteins_Tested = map_dbl(age_results, ~ifelse("combined_age_results" %in% names(.x), nrow(.x$combined_age_results), 0))
)

age_summary$Percent_Significant <- round(age_summary$Significant_Age_Proteins / age_summary$Total_Proteins_Tested * 100, 2)
age_summary$Percent_Significant[is.na(age_summary$Percent_Significant)] <- 0

write.csv(age_summary, "output/tables/age_effects_summary.csv", row.names = FALSE)

# Save detailed age results for each tissue
for(tissue in overlapping_tissues) {
  if("combined_age_results" %in% names(age_results[[tissue]])) {
    write.csv(age_results[[tissue]]$combined_age_results,
              paste0("output/tables/", tissue, "_age_effects_all_proteins.csv"),
              row.names = FALSE)
    
    write.csv(age_results[[tissue]]$significant_proteins,
              paste0("output/tables/", tissue, "_age_effects_significant.csv"),
              row.names = FALSE)
  }
}

# =============================================================================
# 6. COMPREHENSIVE VISUALIZATION AND REPORTING
# =============================================================================

# Overall correlation heatmap across tissues (only available correlations)
corr_cols <- c("Corr_Prop_Pub1", "Corr_Prop_Pub2", "Corr_Pub1_Pub2")
corr_matrix <- correlation_summary[, corr_cols]
rownames(corr_matrix) <- correlation_summary$Tissue

# Remove rows/columns with all NAs
corr_matrix_clean <- corr_matrix[rowSums(!is.na(corr_matrix)) > 0, 
                                colSums(!is.na(corr_matrix)) > 0]
dim(corr_matrix_clean)
if(nrow(corr_matrix_clean) > 0 && ncol(corr_matrix_clean) > 0) {
  png("results/plots/tissue_correlation_heatmap.png", width = 12, height = 10, units = "in", res = 300)
  pheatmap(as.matrix(corr_matrix_clean), 
           cluster_rows = TRUE, cluster_cols = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Dataset Correlations Across Tissues",
           display_numbers = TRUE, number_format = "%.3f",
           na_col = "grey90")
  dev.off()
}

# Age effects summary plot (only tissues with data)
age_summary_complete <- age_summary[age_summary$Status == "Complete", ]

if(nrow(age_summary_complete) > 0) {
  age_summary_plot <- ggplot(age_summary_complete, 
                            aes(x = reorder(Tissue, Percent_Significant), y = Percent_Significant)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Percentage of Age-Associated Proteins by Tissue",
         subtitle = "Based on Proprietary and Published1 datasets only",
         x = "Tissue Type", y = "% Proteins with Significant Age Effects") +
    theme_minimal()
  
  ggsave("results/plots/age_effects_summary.png", age_summary_plot, width = 10, height = 8, dpi = 300)
}

# Sample size summary plot
sample_summary_plot <- ggplot(age_summary_complete, aes(x = Tissue)) +
  geom_bar(aes(y = Samples_Proprietary, fill = "Proprietary"), stat = "identity", alpha = 0.8) +
  geom_bar(aes(y = Samples_Published1, fill = "Published1"), stat = "identity", alpha = 0.8, 
           position = position_stack()) +
  scale_fill_manual(values = c("Proprietary" = "#E31A1C", "Published1" = "#1F78B4")) +
  labs(title = "Sample Sizes by Tissue and Dataset",
       x = "Tissue Type", y = "Number of Samples", fill = "Dataset") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/plots/sample_sizes_summary.png", sample_summary_plot, width = 12, height = 6, dpi = 300)

# Create final summary report
sink("results/proteomics_analysis_summary.txt")
cat("=============================================================================\n")
cat("MULTI-TISSUE PROTEOMICS COMPARATIVE ANALYSIS SUMMARY\n")
cat("=============================================================================\n\n")

cat("Analysis Parameters:\n")
cat(paste("- Number of overlapping proteins:", length(overlapping_proteins), "\n"))
cat(paste("- Number of overlapping tissues:", length(overlapping_tissues), "\n"))
cat("- Dataset characteristics:\n")
cat("  * Proprietary: Individual samples\n")
cat("  * Published1: Ratio matrix [0.002, 442] - Individual samples\n")
cat("  * Published2: Tissue mean matrix - No individual samples\n\n")

cat("Dataset Similarity Results (based on tissue means):\n")
print(correlation_summary)

cat("\nAge Effects Analysis (Proprietary and Published1 only):\n")
print(age_summary)

cat("\nKey Findings:\n")

# Correlation findings
valid_corrs <- correlation_summary[, c("Corr_Prop_Pub1", "Corr_Prop_Pub2", "Corr_Pub1_Pub2")]
valid_corrs <- valid_corrs[!is.na(valid_corrs)]

if(length(valid_corrs) > 0) {
  cat(paste("- Highest dataset correlation:", round(max(valid_corrs), 3), "\n"))
  cat(paste("- Lowest dataset correlation:", round(min(valid_corrs), 3), "\n"))
  cat(paste("- Average correlation across all comparisons:", round(mean(valid_corrs), 3), "\n"))
}

# Age effect findings
if(nrow(age_summary_complete) > 0) {
  cat(paste("- Tissues with age analysis completed:", nrow(age_summary_complete), "out of", nrow(age_summary), "\n"))
  cat(paste("- Tissue with most age-associated proteins:", 
            age_summary_complete$Tissue[which.max(age_summary_complete$Significant_Age_Proteins)], "\n"))
  cat(paste("- Average % of age-associated proteins:", 
            round(mean(age_summary_complete$Percent_Significant), 2), "%\n"))
  cat(paste("- Total samples used for age analysis:", sum(age_summary_complete$Total_Samples), "\n"))
}

cat("\nData Quality Notes:\n")
cat("- Published1 data was log2 transformed from ratio values [0.002, 442]\n")
cat("- Published2 provides tissue-level means only, cannot be used for age analysis\n")
cat("- Age analysis limited to datasets with individual sample information\n")

sink()

cat("\n=============================================================================\n")
cat("Analysis complete! Check the 'results' folder for detailed outputs:\n")
cat("=============================================================================\n")
cat("TABLES:\n")
cat("- dataset_correlations_summary.csv: Correlation between datasets by tissue\n")
cat("- age_effects_summary.csv: Summary of age-related effects by tissue\n") 
cat("- [tissue]_age_effects_all_proteins.csv: Detailed age analysis per tissue\n")
cat("- [tissue]_age_effects_significant.csv: Significant age-associated proteins\n\n")

