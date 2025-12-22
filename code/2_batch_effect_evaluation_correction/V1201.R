#!/usr/bin/env Rscript
# ========================================================================
# Batch effect correction methods - batch evaluation script (connection-fix version)
# Version: 3.1
# ========================================================================

# Clean environment and close all connections
rm(list = ls())
gc()
setwd("E:/YUEL/code/qc_20251201/")
source("../src/src_claude_v1.R")
# Close all open connections
closeAllConnections()

# Check current connection status
message("Current connection status:")
message(paste("Connections used:", length(showConnections())))
message(paste("Maximum connections:", getOption("connections")))

# ========================================================================
# 1. Environment setup and package loading
# ========================================================================

message("=== Batch effect correction batch evaluation started ===")
message(paste("Start time:", Sys.time()))

# Reduce parallel cores to avoid connection issues
library(parallel)
max_cores <- detectCores() - 1
# Limit max cores to 16 to avoid connection issues
n_cores <- min(16, max_cores)  
message(paste("Cores used:", n_cores, "(capped at 16 to avoid connection issues)"))

# Load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  library(umap)
  library(randomForest)
  library(cluster)
  library(viridis)
  library(foreach)
  library(doParallel)
  library(pheatmap)
})

# Do not create a parallel cluster yet; create it only when needed
cl <- NULL

# Source evaluation functions
tryCatch({
  source("../20250803_claude_comprehensive_batch_effect_evaluation_v2.R")
  source("../src/src_claude_v1.R")
}, error = function(e) {
  message(paste0("Warning: Failed to load all source files; some functions may be unavailable. Details: ", e$message))
})

# ========================================================================
# 2. Configuration
# ========================================================================

config <- list(
  # Data files
  raw_data_file = "imputed_raw_12754proteins_2957samples.csv",
  sample_info_file = "batch_info_tab_batchserver.csv",
  
  # Corrected method file pattern
  corrected_data_pattern = "batch_corrected_data_(.*?)3.csv",
  
  # Biological variables
  biological_vars = c("tissue_type_major", "sample_type", "imputation_group"),
  technical_vars = c("new", "instrument", "month", "trans", "batch"),
  
  # PVCA includes all variables
  pvca_vars = c("imputation_group", "sample_type", "tissue_type_major", "instrument", "batch", "new", "month", "trans"),
  
  # Analysis parameters
  pool_label = "pool",
  pvca_threshold = 0.6,
  umap_n_neighbors = 15,
  umap_n_epochs = 200,
  rf_ntree = 500,
  rf_top_features = 50,
  
  # Output directory
  output_dir = paste0("batch_evaluation_", format(Sys.Date(), "%Y%m%d")),
  
  # Visualization parameters
  plot_width = 12,
  plot_height = 8,
  color_palette = "Set2",
  
  # Performance parameters
  use_parallel = FALSE,  # Do not use parallel by default to avoid connection issues
  max_cores = 16
)

# Create output directory
if (!dir.exists(config$output_dir)) {
  dir.create(config$output_dir, recursive = TRUE)
}

# ========================================================================
# 3. Safe data loading function
# ========================================================================

load_data <- function() {
  message("\n[Data loading] Reading raw data...")
  
  # Use tryCatch to ensure file connections are not leaked
  data_matrix <- NULL
  sample_info <- NULL
  
  # Read raw data
  tryCatch({
    con <- file(config$raw_data_file, "r")
    data_matrix <- read.csv(con, row.names = 1)
    close(con)
  }, error = function(e) {
    if (exists("con")) try(close(con), silent = TRUE)
    stop(paste("Failed to read data matrix:", e$message))
  })
  
  tryCatch({
    con <- file(config$sample_info_file, "r")
    sample_info <- read.csv(con, stringsAsFactors = FALSE)
    close(con)
  }, error = function(e) {
    if (exists("con")) try(close(con), silent = TRUE)
    stop(paste("Failed to read sample information:", e$message))
  })
  
  # Ensure correct data type
  data_matrix <- as.matrix(data_matrix)
  
  # Validate dimensions
  if (ncol(data_matrix) != nrow(sample_info)) {
    stop("Sample count mismatch!")
  }
  
  # Read all corrected datasets
  corrected_files <- list.files(pattern = config$corrected_data_pattern, full.names = TRUE)
  corrected_data_list <- list()
  
  if (length(corrected_files) > 0) {
    message(paste("[Data loading] Found", length(corrected_files), "corrected data files"))
    
    for (file in corrected_files) {
      method_name <- gsub("batch_corrected_data_|2.csv", "", basename(file))
      message(paste("  - Reading:", method_name))
      
      tryCatch({
        con <- file(file, "r")
        temp_data <- read.csv(con, row.names = 1)
        close(con)
        corrected_data_list[[method_name]] <- as.matrix(temp_data)
      }, error = function(e) {
        if (exists("con")) try(close(con), silent = TRUE)
        warning(paste("Failed to read", method_name, ":", e$message))
      })
    }
  } else {
    warning("No corrected data files found!")
  }
  
  # Split pool samples
  pool_indices <- which(sample_info$tissue_type_major == config$pool_label)
  
  result <- list(
    data_matrix = data_matrix,
    sample_info = sample_info,
    corrected_data_list = corrected_data_list,
    pool_indices = pool_indices,
    pool_matrix = if(length(pool_indices) > 0) data_matrix[, pool_indices, drop = FALSE] else NULL,
    pool_info = if(length(pool_indices) > 0) sample_info[pool_indices, , drop = FALSE] else NULL
  )
  
  message(paste("[Data loading] Total samples:", ncol(data_matrix)))
  message(paste("[Data loading] Pool samples:", length(pool_indices)))
  message(paste("[Data loading] Features:", nrow(data_matrix)))
  
  # Clean memory
  gc()
  
  return(result)
}

# ========================================================================
# 4. UMAP analysis function (reduced thread count)
# ========================================================================

perform_umap_analysis <- function(data_matrix, sample_info, color_vars, title_prefix = "", n_threads = 1) {
  message(paste("[UMAP]", title_prefix))
  
  # Ensure data is a numeric matrix
  data_matrix <- as.matrix(data_matrix)
  
  # Limit threads to 4 to avoid connection issues
  n_threads <- min(n_threads, 4)
  
  # Run UMAP
  umap_result <- umap(
    t(data_matrix),
    n_neighbors = min(config$umap_n_neighbors, ncol(data_matrix) - 1),
    metric = "euclidean",
    n_epochs = config$umap_n_epochs,
    min_dist = 0.1,
    n_threads = n_threads
  )
  
  # Create data frame
  umap_df <- data.frame(
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2]
  )
  
  # Add sample information
  for (var in color_vars) {
    if (var %in% names(sample_info)) {
      umap_df[[var]] <- as.factor(sample_info[[var]])
    }
  }
  
  # Generate UMAP plots
  plots <- list()
  for (var in color_vars) {
    if (var %in% names(umap_df)) {
      p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = .data[[var]])) +
        geom_point(size = 3, alpha = 0.7) +
        theme_minimal() +
        theme(legend.position = "right") +
        labs(title = paste(title_prefix, "- Colored by", var)) +
        scale_color_brewer(palette = config$color_palette)
      
      plots[[var]] <- p
    }
  }
  
  return(list(
    umap_result = umap_result,
    umap_df = umap_df,
    plots = plots
  ))
}

# ========================================================================
# 11. Main execution function (optimized version)
# ========================================================================

main <- function() {
  # Initialize
  cl <- NULL
  
  # Periodic connection cleanup
  cleanup_connections <- function() {
    cons <- showConnections()
    if (nrow(cons) > 100) {
      message("Cleaning connections...")
      closeAllConnections()
      gc()
    }
  }
  
  tryCatch({
    # Start timing
    start_time <- Sys.time()
    
    # 1. Load data
    data <- load_data()
    cleanup_connections()
    
    # 2. If parallel is needed, create cluster here
    if (config$use_parallel && length(data$corrected_data_list) > 2) {
      message("Creating parallel cluster...")
      cl <- makeCluster(min(4, config$max_cores))  # Limit to 4 cores
      registerDoParallel(cl)
    }
    
    # 3. Batch process all methods
    evaluation_results <- batch_process_methods(data, config$output_dir)
    cleanup_connections()
    
    # Check for valid results
    if (length(evaluation_results) == 0) {
      stop("No evaluations completed successfully")
    }
    
    # 4. Generate comparison results
    comparison <- compare_methods(evaluation_results)
    cleanup_connections()
    
    # 5. Create comparison plots
    if (length(comparison) > 0) {
      create_comparison_plots(comparison, config$output_dir)
    }
    cleanup_connections()
    
    # 6. Create combined visualization
    if (length(evaluation_results) > 1) {
      create_combined_visualization(data, evaluation_results, config$output_dir)
    }
    cleanup_connections()
    
    # 7. Generate comprehensive report
    generate_comprehensive_report(data, evaluation_results, comparison, config$output_dir)
    
    # 8. Save R data objects
    save(
      data,
      evaluation_results,
      comparison,
      file = file.path(config$output_dir, "evaluation_results.RData")
    )
    
    # End timing
    end_time <- Sys.time()
    time_elapsed <- difftime(end_time, start_time, units = "mins")
    
    # message("\n" * 2)
    # message("=" * 80)
    message("Batch effect correction evaluation completed!")
    message(paste("Total time:", round(time_elapsed, 2), "minutes"))
    message(paste("Results saved to:", config$output_dir))
    # message("=" * 80)
    
  }, error = function(e) {
    message("\nError occurred:")
    message(e$message)
    message("\nCall stack:")
    traceback()
  }, finally = {
    # Clean up resources
    if (!is.null(cl) && inherits(cl, "cluster")) {
      tryCatch({
        stopCluster(cl)
        message("Parallel computing resources released")
      }, error = function(e) {
        # Ignore cleanup errors
      })
    }
    
    # Final cleanup of all connections
    closeAllConnections()
    gc()
    
    message(paste("\nFinal connection status: used", length(showConnections()), "connections"))
  })
}

# ========================================================================
# Simplified serial version (if parallel version has issues)
# ========================================================================

main_serial <- function() {
  message("Running serial version (more stable)...")
  
  # Close all connections
  closeAllConnections()
  
  tryCatch({
    start_time <- Sys.time()
    
    # Load data
    data <- load_data()
    
    # Evaluate original data
    message("\nEvaluating original data...")
    orig_result <- comprehensive_evaluation(
      data$data_matrix,
      data$sample_info,
      "Original",
      FALSE
    )
    
    # Save results
    results <- list(orig_result)
    
    # If pool samples exist
    if (!is.null(data$pool_matrix)) {
      pool_result <- comprehensive_evaluation(
        data$pool_matrix,
        data$pool_info,
        "Original",
        TRUE
      )
      results[[length(results) + 1]] <- pool_result
    }
    
    # Simple result output
    for (res in results) {
      if (!is.null(res$silhouette)) {
        message("\nSilhouette results:")
        for (var in names(res$silhouette)) {
          message(paste("  ", var, ":", round(res$silhouette[[var]]$average, 4)))
        }
      }
      
      if (!is.null(res$pool_cv)) {
        message(paste("\nPool CV:", round(res$pool_cv$mean_cv, 4)))
      }
    }
    
    # Save basic results
    saveRDS(results, file.path(config$output_dir, "basic_results.rds"))
    
    end_time <- Sys.time()
    message(paste("\nDone! Time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
    
  }, error = function(e) {
    message(paste("\nError:", e$message))
  }, finally = {
    closeAllConnections()
  })
}

# ========================================================================
# Load other required functions (from previous scripts)
# ========================================================================

# This section needs to include all functions defined previously
# Such as perform_pvca_analysis, perform_enhanced_pca, calculate_silhouette_scores, etc.
# Omitted here due to length; please copy from the previous script
# 4.1 PVCA analysis function
perform_pvca_analysis <- function(data_matrix, sample_info, effect_names, threshold = 0.6, label = "") {
  message(paste("[PVCA]", label))
  
  # Check whether pvca function exists
  if (!exists("pvcaBF")) {
    message("  - PVCA function does not exist, skipping")
    return(list(success = FALSE, error = "pvcaBF function not found"))
  }
  
  tryCatch({
    pvca_results <- pvcaBF(data_matrix, sample_info, effect_names, threshold = threshold)
    
    pvca_df <- data.frame(
      Effect = names(pvca_results),
      Variance = as.numeric(unlist(pvca_results)),
      stringsAsFactors = FALSE
    )
    
    return(list(
      results = pvca_results,
      data = pvca_df,
      success = TRUE
    ))
  }, error = function(e) {
    warning(paste("[PVCA] Error:", e$message))
    return(list(success = FALSE, error = e$message))
  })
}

# 4.2 PCA analysis function
perform_enhanced_pca <- function(data_matrix, sample_info, color_vars, title_prefix = "") {
  message(paste("[PCA]", title_prefix))
  
  # Ensure data is a numeric matrix
  data_matrix <- as.matrix(data_matrix)
  
  # Run PCA
  pca_result <- prcomp(t(data_matrix), scale. = TRUE, center = TRUE)
  
  # Compute variance explained
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  
  # Create data frame
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    PC3 = pca_result$x[, 3],
    PC4 = pca_result$x[, 4]
  )
  
  # Add sample information
  for (var in color_vars) {
    if (var %in% names(sample_info)) {
      pca_df[[var]] <- as.factor(sample_info[[var]])
    }
  }
  
  # Generate multiple PCA plots
  plots <- list()
  for (var in color_vars) {
    if (var %in% names(pca_df)) {
      p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[var]])) +
        geom_point(size = 3, alpha = 0.7) +
        theme_minimal() +
        theme(legend.position = "right") +
        labs(
          title = paste(title_prefix, "- Colored by", var),
          x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
          y = paste0("PC2 (", round(var_explained[2], 1), "%)")
        ) +
        scale_color_brewer(palette = config$color_palette)
      
      plots[[var]] <- p
    }
  }
  
  return(list(
    pca_result = pca_result,
    pca_df = pca_df,
    var_explained = var_explained,
    plots = plots
  ))
}

# 4.3 UMAP analysis function
perform_umap_analysis <- function(data_matrix, sample_info, color_vars, title_prefix = "", n_threads = 1) {
  message(paste("[UMAP]", title_prefix))
  
  # Ensure data is a numeric matrix
  data_matrix <- as.matrix(data_matrix)
  
  # Run UMAP
  umap_result <- umap(
    t(data_matrix),
    n_neighbors = min(config$umap_n_neighbors, ncol(data_matrix) - 1),
    metric = "euclidean",
    n_epochs = config$umap_n_epochs,
    min_dist = 0.1,
    n_threads = n_threads
  )
  
  # Create data frame
  umap_df <- data.frame(
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2]
  )
  
  # Add sample information
  for (var in color_vars) {
    if (var %in% names(sample_info)) {
      umap_df[[var]] <- as.factor(sample_info[[var]])
    }
  }
  
  # Generate UMAP plots
  plots <- list()
  for (var in color_vars) {
    if (var %in% names(umap_df)) {
      p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = .data[[var]])) +
        geom_point(size = 3, alpha = 0.7) +
        theme_minimal() +
        theme(legend.position = "right") +
        labs(title = paste(title_prefix, "- Colored by", var)) +
        scale_color_brewer(palette = config$color_palette)
      
      plots[[var]] <- p
    }
  }
  
  return(list(
    umap_result = umap_result,
    umap_df = umap_df,
    plots = plots
  ))
}

# 4.4 Silhouette analysis function
calculate_silhouette_scores <- function(data_matrix, sample_info, grouping_vars) {
  message("[Silhouette] Calculating silhouette coefficients...")
  
  # Ensure data is a numeric matrix
  data_matrix <- as.matrix(data_matrix)
  
  # Check data validity
  if (ncol(data_matrix) < 2) {
    message("  - Not enough samples, skipping silhouette analysis")
    return(list())
  }
  
  scores <- list()
  
  # Compute distance matrix
  tryCatch({
    dist_matrix <- dist(t(data_matrix))
    
    for (var in grouping_vars) {
      if (var %in% names(sample_info)) {
        # Get grouping information
        groups <- as.factor(sample_info[[var]])
        
        # Check whether there are enough groups
        if (length(unique(groups)) < 2) {
          message(paste("  -", var, "has only one level, skipping"))
          next
        }
        
        # Compute silhouette
        tryCatch({
          clusters <- as.numeric(groups)
          sil <- silhouette(clusters, dist_matrix)
          avg_sil <- mean(sil[, 3], na.rm = TRUE)
          
          scores[[var]] <- list(
            silhouette = sil,
            average = avg_sil,
            by_cluster = tapply(sil[, 3], clusters, mean, na.rm = TRUE)
          )
          
          message(paste("  -", var, "average silhouette:", round(avg_sil, 4)))
          
        }, error = function(e) {
          message(paste("  -", var, "failed:", e$message))
        })
      }
    }
    
  }, error = function(e) {
    message(paste("  - Distance matrix computation failed:", e$message))
  })
  
  return(scores)
}

batch_process_methods <- function(data, output_dir) {
  message("\n[Batch processing] Starting evaluation of all methods...")
  
  all_results <- list()
  
  # 1. Evaluate original data
  message("\n--- Evaluating original data ---")
  
  # All samples
  orig_all <- comprehensive_evaluation(
    data$data_matrix,
    data$sample_info,
    "Original",
    FALSE
  )
  all_results[[length(all_results) + 1]] <- orig_all
  
  # Pool samples
  if (!is.null(data$pool_matrix) && ncol(data$pool_matrix) > 0) {
    orig_pool <- comprehensive_evaluation(
      data$pool_matrix,
      data$pool_info,
      "Original",
      TRUE
    )
    all_results[[length(all_results) + 1]] <- orig_pool
  }
  
  # 2. Evaluate all correction methods
  for (method_name in names(data$corrected_data_list)) {
    message(paste("\n--- Evaluating", method_name, "---"))
    
    corrected_data <- data$corrected_data_list[[method_name]]
    
    # Validate dimensions
    if (!all(dim(corrected_data) == dim(data$data_matrix))) {
      warning(paste("Data dimensions do not match, skipping:", method_name))
      next
    }
    
    # All samples
    method_all <- comprehensive_evaluation(
      corrected_data,
      data$sample_info,
      method_name,
      FALSE
    )
    all_results[[length(all_results) + 1]] <- method_all
    
    # Pool samples
    if (!is.null(data$pool_indices) && length(data$pool_indices) > 0) {
      pool_corrected <- corrected_data[, data$pool_indices]
      method_pool <- comprehensive_evaluation(
        pool_corrected,
        data$pool_info,
        method_name,
        TRUE
      )
      all_results[[length(all_results) + 1]] <- method_pool
    }
    
    # Save detailed results for this method
    method_dir <- file.path(output_dir, method_name)
    if (!dir.exists(method_dir)) {
      dir.create(method_dir)
    }
    
    # Save PCA plots
    if (!is.null(method_all$pca$plots)) {
      for (var_name in names(method_all$pca$plots)) {
        ggsave(
          file.path(method_dir, paste0("pca_", var_name, ".pdf")),
          method_all$pca$plots[[var_name]],
          width = 8, height = 6
        )
      }
    }
    
    # Save UMAP plots
    if (!is.null(method_all$umap$plots)) {
      for (var_name in names(method_all$umap$plots)) {
        ggsave(
          file.path(method_dir, paste0("umap_", var_name, ".pdf")),
          method_all$umap$plots[[var_name]],
          width = 8, height = 6
        )
      }
    }
  }
  
  return(all_results)
}

create_combined_visualization <- function(data, evaluation_results, output_dir) {
  message("\n[Combined visualization] Generating comparison plots...")
  
  # Extract results for all samples
  all_samples_results <- evaluation_results[!sapply(evaluation_results, function(x) x$is_pool)]
  
  # 1. Create PCA comparison (Original vs best method)
  if (length(all_samples_results) >= 2) {
    original_result <- all_samples_results[[1]]
    best_result <- all_samples_results[[length(all_samples_results)]]  # Assume the last is the best
    
    if (!is.null(original_result$pca) && !is.null(best_result$pca)) {
      # Combine PCA plots
      combined_pca <- plot_grid(
        original_result$pca$plots[[1]] + ggtitle("Original Data"),
        best_result$pca$plots[[1]] + ggtitle(best_result$method),
        ncol = 2,
        labels = c("A", "B")
      )
      
      ggsave(
        file.path(output_dir, "pca_comparison_combined.pdf"),
        combined_pca,
        width = 14, height = 6
      )
    }
    
    # 2. Create UMAP comparison
    if (!is.null(original_result$umap) && !is.null(best_result$umap)) {
      combined_umap <- plot_grid(
        original_result$umap$plots[[1]] + ggtitle("Original Data"),
        best_result$umap$plots[[1]] + ggtitle(best_result$method),
        ncol = 2,
        labels = c("A", "B")
      )
      
      ggsave(
        file.path(output_dir, "umap_comparison_combined.pdf"),
        combined_umap,
        width = 14, height = 6
      )
    }
  }
  
  # 3. Create heatmap comparison
  create_heatmap_comparison(evaluation_results, output_dir)
}

# ========================================================================
# 6. Method comparison function
# ========================================================================

compare_methods <- function(evaluation_results) {
  message("\n[Method comparison] Generating comparison report...")
  
  comparison <- list()
  
  # Extract all method names
  all_methods <- unique(sapply(evaluation_results, function(x) x$method))
  all_samples_results <- evaluation_results[!sapply(evaluation_results, function(x) x$is_pool)]
  pool_results <- evaluation_results[sapply(evaluation_results, function(x) x$is_pool)]
  
  # 1. PVCA comparison
  pvca_comparison <- data.frame()
  for (result in all_samples_results) {
    if (result$pvca$success) {
      df <- result$pvca$data
      df$Method <- result$method
      pvca_comparison <- rbind(pvca_comparison, df)
    }
  }
  comparison$pvca <- pvca_comparison
  
  # 2. Silhouette comparison
  sil_comparison <- data.frame()
  for (result in all_samples_results) {
    for (var in names(result$silhouette)) {
      sil_comparison <- rbind(sil_comparison, data.frame(
        Method = result$method,
        Variable = var,
        Silhouette = result$silhouette[[var]]$average,
        stringsAsFactors = FALSE
      ))
    }
  }
  comparison$silhouette <- sil_comparison
  
  # 3. k-BET comparison
  kbet_comparison <- data.frame()
  for (result in all_samples_results) {
    for (var in names(result$kbet)) {
      if (result$kbet[[var]]$success) {
        kbet_comparison <- rbind(kbet_comparison, data.frame(
          Method = result$method,
          Variable = var,
          RejectionRate = result$kbet[[var]]$rejection_rate,
          PValue = result$kbet[[var]]$p_value,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  comparison$kbet <- kbet_comparison
  
  # 4. Batch strength comparison
  batch_comparison <- data.frame()
  for (result in all_samples_results) {
    for (var in names(result$batch_strength)) {
      batch_comparison <- rbind(batch_comparison, data.frame(
        Method = result$method,
        Variable = var,
        SignificantProportion = result$batch_strength[[var]]$significant_proportion,
        MedianP = result$batch_strength[[var]]$median_p,
        stringsAsFactors = FALSE
      ))
    }
  }
  comparison$batch_strength <- batch_comparison
  
  # 5. Pool CV comparison
  if (length(pool_results) > 0) {
    pool_cv_comparison <- data.frame()
    for (result in pool_results) {
      if (!is.null(result$pool_cv)) {
        pool_cv_comparison <- rbind(pool_cv_comparison, data.frame(
          Method = result$method,
          MeanCV = result$pool_cv$mean_cv,
          MedianCV = result$pool_cv$median_cv,
          stringsAsFactors = FALSE
        ))
      }
    }
    comparison$pool_cv <- pool_cv_comparison
  }
  
  return(comparison)
}

# ========================================================================
# 7. Visualization functions
# ========================================================================

create_comparison_plots <- function(comparison, output_dir) {
  message("\n[Visualization] Generating comparison plots...")
  
  # 1. PVCA comparison plot
  if (!is.null(comparison$pvca) && nrow(comparison$pvca) > 0) {
    p_pvca <- ggplot(comparison$pvca, aes(x = Effect, y = Variance, fill = Method)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "PVCA Comparison Across Methods",
           y = "Proportion of Variance",
           x = "Effect") +
      scale_fill_brewer(palette = config$color_palette)
    
    ggsave(file.path(output_dir, "pvca_comparison.pdf"), p_pvca, 
           width = config$plot_width, height = config$plot_height)
  }
  
  # 2. Silhouette comparison plot
  if (!is.null(comparison$silhouette) && nrow(comparison$silhouette) > 0) {
    # Separate biological and technical variables
    bio_vars <- comparison$silhouette[comparison$silhouette$Variable %in% config$biological_vars, ]
    tech_vars <- comparison$silhouette[comparison$silhouette$Variable %in% config$technical_vars, ]
    
    if (nrow(bio_vars) > 0) {
      p_sil_bio <- ggplot(bio_vars, aes(x = Method, y = Silhouette, fill = Variable)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Silhouette Scores - Biological Variables",
             y = "Average Silhouette Score") +
        scale_fill_brewer(palette = "Set1") +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)
      
      ggsave(file.path(output_dir, "silhouette_biological.pdf"), p_sil_bio,
             width = config$plot_width, height = config$plot_height * 0.75)
    }
    
    if (nrow(tech_vars) > 0) {
      p_sil_tech <- ggplot(tech_vars, aes(x = Method, y = Silhouette, fill = Variable)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Silhouette Scores - Technical Variables",
             y = "Average Silhouette Score") +
        scale_fill_brewer(palette = "Set2") +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)
      
      ggsave(file.path(output_dir, "silhouette_technical.pdf"), p_sil_tech,
             width = config$plot_width, height = config$plot_height * 0.75)
    }
  }
  
  # 3. k-BET comparison plot
  if (!is.null(comparison$kbet) && nrow(comparison$kbet) > 0) {
    p_kbet <- ggplot(comparison$kbet, aes(x = Method, y = 1 - RejectionRate, fill = Variable)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "k-BET Acceptance Rate Comparison",
           y = "Acceptance Rate (1 - Rejection Rate)") +
      scale_fill_brewer(palette = config$color_palette) +
      ylim(0, 1)
    
    ggsave(file.path(output_dir, "kbet_comparison.pdf"), p_kbet,
           width = config$plot_width, height = config$plot_height * 0.75)
  }
  
  # 4. Batch strength comparison plot
  if (!is.null(comparison$batch_strength) && nrow(comparison$batch_strength) > 0) {
    p_batch <- ggplot(comparison$batch_strength, 
                      aes(x = Method, y = 1 - SignificantProportion, fill = Variable)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Batch Effect Reduction Comparison",
           y = "Proportion of Non-significant Features") +
      scale_fill_brewer(palette = config$color_palette) +
      ylim(0, 1)
    
    ggsave(file.path(output_dir, "batch_strength_comparison.pdf"), p_batch,
           width = config$plot_width, height = config$plot_height * 0.75)
  }
  
  # 5. Pool CV comparison plot
  if (!is.null(comparison$pool_cv) && nrow(comparison$pool_cv) > 0) {
    p_cv <- ggplot(comparison$pool_cv, aes(x = Method, y = MeanCV)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Pool Sample CV Comparison",
           y = "Mean CV") +
      geom_text(aes(label = round(MeanCV, 3)), vjust = -0.5)
    
    ggsave(file.path(output_dir, "pool_cv_comparison.pdf"), p_cv,
           width = config$plot_width * 0.75, height = config$plot_height * 0.75)
  }
  
  message("[Visualization] Comparison plots saved")
}

# ========================================================================
# 8. Report generation function
# ========================================================================

generate_comprehensive_report <- function(data, evaluation_results, comparison, output_dir) {
  message("\n[Report] Generating comprehensive evaluation report...")
  
  report_file <- file.path(config$output_dir, "comprehensive_evaluation_report.txt")
  
  sink(report_file)
  
  # cat("=" * 80, "\n")
  cat("Comprehensive Evaluation Report for Batch Effect Correction Methods\n")
  # cat("=" * 80, "\n\n")
  
  cat("Generation time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("R version:", R.version.string, "\n\n")
  
  # Data overview
  cat("Data overview\n")
  # cat("-" * 40, "\n")
  cat("Total samples:", ncol(data$data_matrix), "\n")
  cat("Pool samples:", length(data$pool_indices), "\n")
  cat("Features:", nrow(data$data_matrix), "\n")
  cat("Number of evaluated methods:", length(unique(sapply(evaluation_results, function(x) x$method))), "\n")
  cat("Biological variables:", paste(config$biological_vars, collapse = ", "), "\n")
  cat("Technical variables:", paste(config$technical_vars, collapse = ", "), "\n\n")
  
  # Method list
  cat("Evaluated correction methods\n")
  # cat("-" * 40, "\n")
  methods <- unique(sapply(evaluation_results, function(x) x$method))
  for (i in seq_along(methods)) {
    cat(sprintf("%d. %s\n", i, methods[i]))
  }
  cat("\n")
  
  # PVCA summary
  if (!is.null(comparison$pvca)) {
    cat("PVCA analysis results\n")
    # cat("-" * 40, "\n")
    
    # Summarize by method
    for (method in methods) {
      method_pvca <- comparison$pvca[comparison$pvca$Method == method, ]
      if (nrow(method_pvca) > 0) {
        cat(sprintf("\n%s:\n", method))
        tech_var <- sum(method_pvca$Variance[method_pvca$Effect %in% config$technical_vars])
        bio_var <- sum(method_pvca$Variance[method_pvca$Effect %in% config$biological_vars])
        cat(sprintf("  Total variance of technical variables: %.2f%%\n", tech_var ))
        cat(sprintf("  Total variance of biological variables: %.2f%%\n", bio_var))
        # cat(sprintf("  Bio/tech ratio: %.2f\n", bio_var / max(tech_var, 0.001)))
      }
    }
    cat("\n")
  }
  
  # Silhouette summary
  if (!is.null(comparison$silhouette)) {
    cat("Silhouette coefficient analysis\n")
    # cat("-" * 40, "\n")
    
    # Compute improvement
    original_method <- "Original"
    if (original_method %in% comparison$silhouette$Method) {
      original_sil <- comparison$silhouette[comparison$silhouette$Method == original_method, ]
      
      for (method in setdiff(methods, original_method)) {
        method_sil <- comparison$silhouette[comparison$silhouette$Method == method, ]
        cat(sprintf("\n%s vs Original:\n", method))
        
        for (var in config$biological_vars) {
          orig_val <- original_sil$Silhouette[original_sil$Variable == var]
          method_val <- method_sil$Silhouette[method_sil$Variable == var]
          if (length(orig_val) > 0 && length(method_val) > 0) {
            improvement <- ((method_val - orig_val) / abs(orig_val)) * 100
            cat(sprintf("  %s: %.4f -> %.4f (improvement: %+.1f%%)\n", 
                        var, orig_val, method_val, improvement))
          }
        }
      }
    }
    cat("\n")
  }
  
  # k-BET summary
  if (!is.null(comparison$kbet)) {
    cat("k-BET batch effect testing\n")
    # cat("-" * 40, "\n")
    
    for (method in methods) {
      method_kbet <- comparison$kbet[comparison$kbet$Method == method, ]
      if (nrow(method_kbet) > 0) {
        cat(sprintf("\n%s:\n", method))
        for (i in 1:nrow(method_kbet)) {
          cat(sprintf("  %s: acceptance rate=%.3f, p value=%.4f\n",
                      method_kbet$Variable[i],
                      1 - method_kbet$RejectionRate[i],
                      method_kbet$PValue[i]))
        }
      }
    }
    cat("\n")
  }
  
  # Batch strength summary
  if (!is.null(comparison$batch_strength)) {
    cat("Batch effect strength analysis\n")
    # cat("-" * 40, "\n")
    
    for (method in methods) {
      method_batch <- comparison$batch_strength[comparison$batch_strength$Method == method, ]
      if (nrow(method_batch) > 0) {
        cat(sprintf("\n%s:\n", method))
        for (i in 1:nrow(method_batch)) {
          cat(sprintf("  %s: proportion of non-significant features=%.3f, median p value=%.4f\n",
                      method_batch$Variable[i],
                      1 - method_batch$SignificantProportion[i],
                      method_batch$MedianP[i]))
        }
      }
    }
    cat("\n")
  }
  
  # Pool CV summary
  if (!is.null(comparison$pool_cv)) {
    cat("Pool sample coefficient of variation analysis\n")
    # cat("-" * 40, "\n")
    
    # Sort by CV ascending
    comparison$pool_cv <- comparison$pool_cv[order(comparison$pool_cv$MeanCV), ]
    
    cat("\nMethod ranking (smaller CV is better):\n")
    for (i in 1:nrow(comparison$pool_cv)) {
      cat(sprintf("%d. %s: mean CV=%.4f, median CV=%.4f\n",
                  i,
                  comparison$pool_cv$Method[i],
                  comparison$pool_cv$MeanCV[i],
                  comparison$pool_cv$MedianCV[i]))
    }
    
    # Relative improvement
    if ("Original" %in% comparison$pool_cv$Method) {
      original_cv <- comparison$pool_cv$MeanCV[comparison$pool_cv$Method == "Original"]
      cat("\nImprovement relative to original data:\n")
      for (i in 1:nrow(comparison$pool_cv)) {
        if (comparison$pool_cv$Method[i] != "Original") {
          improvement <- ((original_cv - comparison$pool_cv$MeanCV[i]) / original_cv) * 100
          cat(sprintf("  %s: %.1f%%\n", comparison$pool_cv$Method[i], improvement))
        }
      }
    }
    cat("\n")
  }
  
  # Overall scoring and recommendation
  # cat("=" * 80, "\n")
  cat("Overall evaluation and recommendation\n")
  # cat("=" * 80, "\n\n")
  
  # Compute overall score
  method_scores <- data.frame(Method = methods, stringsAsFactors = FALSE)
  method_scores$Score <- 0
  
  # Score based on each metric
  if (!is.null(comparison$pvca)) {
    for (method in methods) {
      method_pvca <- comparison$pvca[comparison$pvca$Method == method, ]
      if (nrow(method_pvca) > 0) {
        bio_var <- sum(method_pvca$Variance[method_pvca$Effect %in% config$biological_vars])
        tech_var <- sum(method_pvca$Variance[method_pvca$Effect %in% config$technical_vars])
        score <- bio_var / max(tech_var, 0.001)
        method_scores$Score[method_scores$Method == method] <- 
          method_scores$Score[method_scores$Method == method] + score
      }
    }
  }
  
  if (!is.null(comparison$silhouette)) {
    for (method in methods) {
      method_sil <- comparison$silhouette[comparison$silhouette$Method == method, ]
      bio_sil <- mean(method_sil$Silhouette[method_sil$Variable %in% config$biological_vars], na.rm = TRUE)
      method_scores$Score[method_scores$Method == method] <- 
        method_scores$Score[method_scores$Method == method] + bio_sil * 2
    }
  }
  
  if (!is.null(comparison$pool_cv) && "Original" %in% comparison$pool_cv$Method) {
    original_cv <- comparison$pool_cv$MeanCV[comparison$pool_cv$Method == "Original"]
    for (method in methods) {
      if (method %in% comparison$pool_cv$Method) {
        method_cv <- comparison$pool_cv$MeanCV[comparison$pool_cv$Method == method]
        cv_improvement <- (original_cv - method_cv) / original_cv
        method_scores$Score[method_scores$Method == method] <- 
          method_scores$Score[method_scores$Method == method] + cv_improvement
      }
    }
  }
  
  # Sort and output
  method_scores <- method_scores[order(method_scores$Score, decreasing = TRUE), ]
  
  cat("Overall method score ranking:\n")
  # cat("-" * 40, "\n")
  for (i in 1:nrow(method_scores)) {
    cat(sprintf("%d. %s: %.3f\n", i, method_scores$Method[i], method_scores$Score[i]))
  }
  
  cat("\nRecommendation:\n")
  # cat("-" * 40, "\n")
  
  if (nrow(method_scores) > 1) {
    best_method <- method_scores$Method[1]
    second_method <- method_scores$Method[2]
    
    score_diff <- (method_scores$Score[1] - method_scores$Score[2]) / method_scores$Score[2] * 100
    
    if (score_diff > 20) {
      cat(sprintf("Strongly recommend %s, as its overall performance is clearly better than others (%.1f%% ahead of the runner-up)\n", 
                  best_method, score_diff))
    } else if (score_diff > 10) {
      cat(sprintf("Recommend %s, with good overall performance\n", best_method))
      cat(sprintf("Alternative: %s\n", second_method))
    } else {
      cat(sprintf("%s and %s perform similarly; choose based on specific needs\n", best_method, second_method))
    }
    
    # Special note
    if (!is.null(comparison$pool_cv)) {
      best_cv_method <- comparison$pool_cv$Method[1]
      if (best_cv_method != best_method) {
        cat(sprintf("\nNote: If technical reproducibility is the priority, %s performs best on pool samples\n", best_cv_method))
      }
    }
  }
  
  cat("\n")
  # cat("=" * 80, "\n")
  cat("End of report\n")
  # cat("=" * 80, "\n")
  
  sink()
  
  message(paste("[Report] Comprehensive evaluation report saved to:", report_file))
}

# ========================================================================

# Heatmap function
create_heatmap_comparison <- function(evaluation_results, output_dir) {
  message("[Heatmap] Generating evaluation metrics heatmap...")
  
  # Prepare data matrix
  all_samples_results <- evaluation_results[!sapply(evaluation_results, function(x) x$is_pool)]
  
  if (length(all_samples_results) < 2) return()
  
  # Create scoring matrix
  methods <- sapply(all_samples_results, function(x) x$method)
  metrics <- c()
  scores <- matrix(NA, nrow = 0, ncol = length(methods))
  colnames(scores) <- methods
  
  # Add Silhouette scores
  for (var in config$biological_vars) {
    metric_name <- paste("Silhouette", var)
    metrics <- c(metrics, metric_name)
    row_scores <- sapply(all_samples_results, function(x) {
      if (!is.null(x$silhouette[[var]])) {
        return(x$silhouette[[var]]$average)
      }
      return(NA)
    })
    scores <- rbind(scores, row_scores)
  }
  
  # Add batch strength
  for (var in config$technical_vars[1:2]) {
    metric_name <- paste("BatchStrength", var)
    metrics <- c(metrics, metric_name)
    row_scores <- sapply(all_samples_results, function(x) {
      if (!is.null(x$batch_strength[[var]])) {
        return(1 - x$batch_strength[[var]]$significant_proportion)
      }
      return(NA)
    })
    scores <- rbind(scores, row_scores)
  }
  
  rownames(scores) <- metrics
  
  # Normalize scores (0-1)
  scores_normalized <- t(apply(scores, 1, function(x) {
    if (all(is.na(x))) return(x)
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }))
  
  # Create heatmap
  library(pheatmap)
  
  pdf(file.path(output_dir, "evaluation_heatmap.pdf"), width = 10, height = 8)
  pheatmap(
    scores_normalized,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    number_format = "%.2f",
    main = "Batch Correction Methods Evaluation Heatmap",
    color = colorRampPalette(c("red", "white", "blue"))(100),
    na_col = "grey",
    angle_col = 45
  )
  dev.off()
  
  message("[Heatmap] Evaluation metrics heatmap saved")
}

# ========================================================================
# 10. Create combined comparison plots
# ========================================================================

create_combined_visualization <- function(data, evaluation_results, output_dir) {
  message("\n[Combined visualization] Generating comparison plots...")
  
  # Extract results for all samples
  all_samples_results <- evaluation_results[!sapply(evaluation_results, function(x) x$is_pool)]
  
  # 1. Create PCA comparison (Original vs best method)
  if (length(all_samples_results) >= 2) {
    original_result <- all_samples_results[[1]]
    best_result <- all_samples_results[[length(all_samples_results)]]  # Assume the last is the best
    
    if (!is.null(original_result$pca) && !is.null(best_result$pca)) {
      # Combine PCA plots
      combined_pca <- plot_grid(
        original_result$pca$plots[[1]] + ggtitle("Original Data"),
        best_result$pca$plots[[1]] + ggtitle(best_result$method),
        ncol = 2,
        labels = c("A", "B")
      )
      
      ggsave(
        file.path(output_dir, "pca_comparison_combined.pdf"),
        combined_pca,
        width = 14, height = 6
      )
    }
    
    # 2. Create UMAP comparison
    if (!is.null(original_result$umap) && !is.null(best_result$umap)) {
      combined_umap <- plot_grid(
        original_result$umap$plots[[1]] + ggtitle("Original Data"),
        best_result$umap$plots[[1]] + ggtitle(best_result$method),
        ncol = 2,
        labels = c("A", "B")
      )
      
      ggsave(
        file.path(output_dir, "umap_comparison_combined.pdf"),
        combined_umap,
        width = 14, height = 6
      )
    }
  }
  
  # 3. Create heatmap comparison
  create_heatmap_comparison(evaluation_results, output_dir)
}

# Heatmap function
create_heatmap_comparison <- function(evaluation_results, output_dir) {
  message("[Heatmap] Generating evaluation metrics heatmap...")
  
  # Prepare data matrix
  all_samples_results <- evaluation_results[!sapply(evaluation_results, function(x) x$is_pool)]
  
  if (length(all_samples_results) < 2) return()
  
  # Create scoring matrix
  methods <- sapply(all_samples_results, function(x) x$method)
  metrics <- c()
  scores <- matrix(NA, nrow = 0, ncol = length(methods))
  colnames(scores) <- methods
  
  # Add Silhouette scores
  for (var in config$biological_vars) {
    metric_name <- paste("Silhouette", var)
    metrics <- c(metrics, metric_name)
    row_scores <- sapply(all_samples_results, function(x) {
      if (!is.null(x$silhouette[[var]])) {
        return(x$silhouette[[var]]$average)
      }
      return(NA)
    })
    scores <- rbind(scores, row_scores)
  }
  
  # Add batch strength
  for (var in config$technical_vars[1:2]) {
    metric_name <- paste("BatchStrength", var)
    metrics <- c(metrics, metric_name)
    row_scores <- sapply(all_samples_results, function(x) {
      if (!is.null(x$batch_strength[[var]])) {
        return(1 - x$batch_strength[[var]]$significant_proportion)
      }
      return(NA)
    })
    scores <- rbind(scores, row_scores)
  }
  
  rownames(scores) <- metrics
  
  # Normalize scores (0-1)
  scores_normalized <- t(apply(scores, 1, function(x) {
    if (all(is.na(x))) return(x)
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }))
  
  # Create heatmap
  library(pheatmap)
  
  pdf(file.path(output_dir, "evaluation_heatmap.pdf"), width = 10, height = 8)
  pheatmap(
    scores_normalized,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    number_format = "%.2f",
    main = "Batch Correction Methods Evaluation Heatmap",
    color = colorRampPalette(c("red", "white", "blue"))(100),
    na_col = "grey",
    angle_col = 45
  )
  dev.off()
  
  message("[Heatmap] Evaluation metrics heatmap saved")
}

# message("=" * 60)
message("Connection-issue fixed version loaded")
message("Main changes:")
message("1. Limit parallel cores to 16")
message("2. Limit UMAP threads to 4")
message("3. Add connection cleanup mechanism")
message("4. Provide serial version as an alternative")
# message("=" * 60)
message("")
message("Usage:")
message("1. Parallel version: main()")
message("2. Serial version: main_serial() (more stable)")
# message("=" * 60)

#4.5 k-BET analysis function (fixed version)
perform_kbet_analysis <- function(data_matrix, sample_info, batch_var, k0 = NULL) {
  message(paste("[k-BET] Analyzing batch effect -", batch_var))
  
  # Check if kBET package is available
  if (!requireNamespace("kBET", quietly = TRUE)) {
    message("  - kBET package not installed, skipping")
    return(list(success = FALSE, error = "kBET not installed"))
  }
  
  # k-BET should only be used for technical variables
  if (!batch_var %in% config$technical_vars) {
    message(paste("  -", batch_var, "is not a technical variable, skipping k-BET analysis"))
    return(list(success = FALSE, error = "Not a technical variable"))
  }
  
  tryCatch({
    batch <- as.factor(sample_info[[batch_var]])
    
    # Check number of batches
    if (length(unique(batch)) < 2) {
      message(paste("  -", batch_var, "has only one batch, skipping k-BET"))
      return(list(success = FALSE, error = "Only one batch"))
    }
    
    # Set k0 automatically
    if (is.null(k0)) {
      min_batch_size <- min(table(batch))
      k0 <- max(5, min(10, floor(min_batch_size * 0.5)))
    }
    
    # Check sample count
    if (ncol(data_matrix) < k0 + 1) {
      message(paste("  - Insufficient samples, skipping k-BET"))
      return(list(success = FALSE, error = "Insufficient samples"))
    }
    
    # Ensure data format is correct
    data_for_kbet <- t(as.matrix(data_matrix))
    
    # Run k-BET test
    kbet_result <- kBET::kBET(
      data_for_kbet,
      batch,
      k0 = k0,
      plot = FALSE,
      do.pca = TRUE,
      verbose = FALSE
    )
    
    # Safely extract results
    rejection_rate <- NA
    p_value <- NA
    
    if (is.list(kbet_result)) {
      if ("summary" %in% names(kbet_result)) {
        if (is.data.frame(kbet_result$summary) || is.matrix(kbet_result$summary)) {
          if ("kBET.observed" %in% colnames(kbet_result$summary)) {
            rejection_rate <- as.numeric(kbet_result$summary[1, "kBET.observed"])
            if (nrow(kbet_result$summary) > 1) {
              p_value <- as.numeric(kbet_result$summary[2, "kBET.observed"])
            }
          }
        }
      }
    }
    
    if (!is.na(rejection_rate)) {
      message(paste("  - Rejection rate:", round(rejection_rate, 3)))
      return(list(
        success = TRUE,
        rejection_rate = rejection_rate,
        p_value = p_value,
        result = kbet_result
      ))
    } else {
      return(list(
        success = FALSE,
        error = "Could not extract kBET results"
      ))
    }
    
  }, error = function(e) {
    message(paste("  - k-BET error:", e$message))
    return(list(success = FALSE, error = e$message))
  })
}

# 4.6 Batch strength calculation
calculate_batch_strength <- function(data_matrix, sample_info, batch_var) {
  message(paste("[Batch strength]", batch_var))
  
  # Compute batch strength only for technical variables
  if (!batch_var %in% config$technical_vars) {
    return(NULL)
  }
  
  # Ensure data is a numeric matrix
  data_matrix <- as.matrix(data_matrix)
  
  # Use ANOVA to compute batch effect strength
  batch <- as.factor(sample_info[[batch_var]])
  
  if (length(unique(batch)) < 2) {
    return(NULL)
  }
  
  batch_strength <- apply(data_matrix, 1, function(x) {
    tryCatch({
      model <- aov(x ~ batch)
      summary(model)[[1]][1, 5]  # p-value
    }, error = function(e) {
      return(1)
    })
  })
  
  # Proportion of significant features
  sig_proportion <- mean(batch_strength < 0.05, na.rm = TRUE)
  
  return(list(
    p_values = batch_strength,
    significant_proportion = sig_proportion,
    median_p = median(batch_strength, na.rm = TRUE)
  ))
}
# 
# # ========================================================================
# # 5. Comprehensive evaluation function (fixed version)
# # ========================================================================

comprehensive_evaluation <- function(data_matrix, sample_info, method_name = "Unknown", is_pool = FALSE) {
  message(paste("\n[Comprehensive evaluation]", method_name, ifelse(is_pool, "- Pool samples", "- All samples")))
  
  results <- list(
    method = method_name,
    is_pool = is_pool,
    timestamp = Sys.time()
  )
  
  # Ensure data is a numeric matrix
  data_matrix <- as.matrix(data_matrix)
  
  # Check data validity
  if (ncol(data_matrix) < 2) {
    message("  - Insufficient samples, skipping most analyses")
    results$error <- "Insufficient samples"
    return(results)
  }
  
  # 1. PVCA analysis
  if (ncol(data_matrix) >= 10) {
    results$pvca <- perform_pvca_analysis(
      data_matrix,
      sample_info,
      config$pvca_vars,
      config$pvca_threshold,
      paste(method_name, ifelse(is_pool, "Pool", "All"))
    )
  } else {
    message("  - Too few samples, skipping PVCA")
    results$pvca <- list(success = FALSE, error = "Insufficient samples for PVCA")
  }
  
  # 2. PCA analysis
  if (ncol(data_matrix) >= 3) {
    # Select key visualization variables
    pca_vars <- c("batch", "instrument","sample_type")
    pca_vars <- pca_vars[pca_vars %in% names(sample_info)]
    
    if (length(pca_vars) > 0) {
      results$pca <- perform_enhanced_pca(
        data_matrix,
        sample_info,
        pca_vars,
        paste(method_name, ifelse(is_pool, "Pool", "All"))
      )
    }
  } else {
    message("  - Too few samples, skipping PCA")
    results$pca <- list(success = FALSE)
  }
  
  # 3. UMAP analysis
  if (ncol(data_matrix) >= config$umap_n_neighbors + 1) {
    # Select key variables
    umap_vars <- c("batch", "tissue_type_major", "sample_type")
    umap_vars <- umap_vars[umap_vars %in% names(sample_info)]
    
    if (length(umap_vars) > 0) {
      results$umap <- perform_umap_analysis(
        data_matrix,
        sample_info,
        umap_vars,
        paste(method_name, ifelse(is_pool, "Pool", "All")),
        n_threads = min(8, n_cores)
      )
    }
  } else {
    message("  - Too few samples, skipping UMAP")
    results$umap <- list(success = FALSE)
  }
  
  # 4. Silhouette analysis
  if (is_pool) {
    # Pool samples: analyze only valid variables
    valid_vars <- c()
    # Only analyze technical variables (biological variables are usually identical for pool samples)
    for (var in config$technical_vars) {
      if (var %in% names(sample_info) && length(unique(sample_info[[var]])) > 1) {
        valid_vars <- c(valid_vars, var)
      }
    }
    
    if (length(valid_vars) > 0) {
      results$silhouette <- calculate_silhouette_scores(
        data_matrix,
        sample_info,
        valid_vars
      )
    } else {
      message("  - No valid grouping variables for pool samples, skipping silhouette")
      results$silhouette <- list()
    }
  } else {
    # Non-pool samples: analyze all variables (except imputation_group if groups are too small)
    all_vars <- c(config$biological_vars, config$technical_vars)
    results$silhouette <- calculate_silhouette_scores(
      data_matrix,
      sample_info,
      all_vars
    )
  }
  
  # 5. k-BET analysis (technical variables only)
  results$kbet <- list()
  for (batch_var in config$technical_vars) {
    if (batch_var %in% names(sample_info)) {
      if (length(unique(sample_info[[batch_var]])) > 1) {
        results$kbet[[batch_var]] <- perform_kbet_analysis(
          data_matrix,
          sample_info,
          batch_var
        )
      } else {
        message(paste("  -", batch_var, "has only one level, skipping k-BET"))
      }
    }
  }
  
  # 6. Batch strength analysis (technical variables only)
  results$batch_strength <- list()
  for (batch_var in config$technical_vars) {
    if (batch_var %in% names(sample_info)) {
      if (length(unique(sample_info[[batch_var]])) > 1) {
        strength_result <- calculate_batch_strength(
          data_matrix,
          sample_info,
          batch_var
        )
        if (!is.null(strength_result)) {
          results$batch_strength[[batch_var]] <- strength_result
        }
      }
    }
  }
  
  # 7. Pool sample CV
  if (is_pool && nrow(data_matrix) > 0 && ncol(data_matrix) > 1) {
    cv_values <- apply(data_matrix, 1, function(x) {
      x_mean <- mean(x, na.rm = TRUE)
      x_sd <- sd(x, na.rm = TRUE)
      if (is.na(x_sd) || is.na(x_mean) || x_mean == 0) return(NA)
      abs(x_sd / x_mean)
    })
    
    cv_values <- cv_values[!is.na(cv_values) & !is.infinite(cv_values)]
    
    if (length(cv_values) > 0) {
      results$pool_cv <- list(
        mean_cv = mean(cv_values, na.rm = TRUE),
        median_cv = median(cv_values, na.rm = TRUE),
        cv_distribution = cv_values
      )
      message(paste("  - Pool mean CV:", round(results$pool_cv$mean_cv, 4)))
    } else {
      message("  - No valid CV values")
      results$pool_cv <- NULL
    }
  }
  
  return(results)
}

# # ========================================================================
# # Subsequent functions unchanged or use the previous fixed version...
# # ========================================================================
# 
# # Other functions (compare_methods, create_comparison_plots, etc.) are omitted here
# # Please use the previous version
# 
# message("Full fixed version loaded")
# message("Main fixes:")
# message("1. Properly distinguish biological vs technical variables")
# message("2. Fix k-BET result extraction")
# message("3. Ensure all data are numeric matrices")
# message("4. Special handling for pool samples")
# message("Run main() to start analysis")
