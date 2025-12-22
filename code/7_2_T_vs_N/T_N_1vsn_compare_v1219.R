# 0. Packages ----
# rm(list = ls())
# pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Reproducibility & parallel
set.seed(1L)
options(stringsAsFactors = FALSE)
n_cores <- min(64L, parallel::detectCores())
BiocParallel::register(BiocParallel::SnowParam(workers = max(1L, n_cores - 2L)))

# Working dir & helpers
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('../source/source_code.R')
library(data.table)


# 1.Read data --------
info.all <- rio::import('../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V6.xlsx') %>%
  filter(!Tissue_heterogeneity_abnormal, !Low_protein_IDs, !FileName %in% pancancer_low_pearson)
dim(info.all) # 2921 34

datall <- rio::import('../2_batch_effect_evaluation_correction/output/batch_corrected_data_step2_limma3.csv') %>% column_to_rownames('V1')
df.all <- datall %>% t() %>% as.data.frame() %>% rownames_to_column('FileName') %>% 
  inner_join(info.all, .)

# T and N data
meta <- info.all %>% filter(sample_type %in% c('T', 'N'))
expr <- datall[, meta$FileName]

dim(meta) # 1584   34
dim(expr) # 12797  1584

# meta <- rio::import('../6_3paired_tumor_nontumor/T_NT_info_check.xlsx')
# df_dep <- rio::import('../6_3paired_tumor_nontumor/output/tumor_nontumor_analysis_reports_v1202.xlsx', sheet = 3)
# prot_meta <- rio::import('../6_3paired_tumor_nontumor/T_NT_shared_DEP_1_22.xlsx')


# 2.DEP ----
meta %<>%
  mutate(cancer_subtype = cancer_abbr, .after = class_abbr) %>% 
  mutate(cancer_abbr = ifelse(str_detect(cancer_abbr, '^BRCA'), 'BRCA', cancer_abbr),
         cancer_abbr = ifelse(str_detect(cancer_abbr, 'CCOC|HGSOC'), 'OC', cancer_abbr)) %>% 
  mutate(cancer = ifelse(str_detect(cancer, '^Breast carcinoma'), 'Breast carcinoma', cancer),
         cancer = ifelse(str_detect(cancer, 'ovarian cancer$'), 'Ovarian cancer', cancer)) %>% 
  mutate(Dataset = c('data1' = 'DS1', 'data2' = 'DS2', 'data3' = 'DS3')[new],
         .before = new) %>%
  mutate(sample_type = factor(sample_type, levels = c('N', 'T')),  # N as reference
         Dataset = factor(Dataset),
         cancer_abbr = ifelse(sample_type == 'N', 'NORMAL', cancer_abbr),
         class_abbr = str_replace_all(class_abbr, '\\-', '.'))

# df <- meta %>% select(FileName, patient_ID, sample_type, cancer_abbr, cancer_subtype, Group) %>% 
#   mutate(sample_type = factor(as.character(sample_type), levels = c("N", "T"))) %>% 
#   inner_join(expr %>% t() %>% as.data.frame() %>% rownames_to_column('FileName'))

vec_ca <- sort(unique(meta$cancer_abbr))

## method1 T vs N (1 vs n) -----
### helper -----
# Differential expression function with adjusted Hedges'g added
perform_de_analysis <- function(cancer_type, metadata, expression_matrix) {
  
  # Filter samples for comparison (unchanged)
  if (cancer_type == 'ALL_TUMORS') {
    samples_meta <- metadata
  } else {
    samples_meta <- metadata %>%
      dplyr::filter(cancer_abbr %in% c(cancer_type, 'NORMAL') | 
                      (cancer_abbr == cancer_type & sample_type == 'T') |
                      sample_type == 'N')
  }
  
  # Get expression data (unchanged)
  samples_expr <- expression_matrix[, samples_meta$FileName]
  
  # Check sample requirements (unchanged)
  n_tumor <- sum(samples_meta$sample_type == 'T')
  n_normal <- sum(samples_meta$sample_type == 'N')
  
  if (n_tumor < 3 || n_normal < 3) {
    return(NULL)
  }
  
  # Helper: Hedges' small-sample bias correction factor J for adjusted effects (as in Code1)
  J_df <- function(df) {
    ifelse(
      is.finite(df) && df > 1,
      exp(lgamma(df/2) - 0.5*log(df/2) - lgamma((df-1)/2)),
      NA_real_
    )
  }
  
  # DE analysis for each protein (unchanged loop body structure)
  de_results <- parallel::mclapply(rownames(samples_expr), function(protein) {
    cat(sprintf("Analyzing %s...\r", protein))
    
    # Prepare model data (unchanged)
    model_data <- data.frame(
      expression = as.numeric(samples_expr[protein, ]),
      sample_type = samples_meta$sample_type,
      Dataset = samples_meta$Dataset
    )
    
    # Remove NA (unchanged)
    model_data <- model_data[!is.na(model_data$expression), ]
    
    # Check groups after NA removal (unchanged)
    if (length(unique(model_data$sample_type)) < 2 || nrow(model_data) < 6) {
      return(NULL)
    }
    
    # Fit model with Dataset covariate if needed (unchanged)
    if (length(unique(model_data$Dataset)) > 1) {
      fit <- stats::lm(expression ~ sample_type + Dataset, data = model_data)
    } else {
      fit <- stats::lm(expression ~ sample_type, data = model_data)
    }
    
    # Extract summaries and coefficients (extended to get sigma and df as in Code1)
    sm <- summary(fit)
    coef_summary <- sm$coefficients
    
    if (!'sample_typeT' %in% rownames(coef_summary)) {
      return(NULL)
    }
    
    # Calculate metrics (existing calculations retained)
    tumor_expr  <- model_data$expression[model_data$sample_type == 'T']
    normal_expr <- model_data$expression[model_data$sample_type == 'N']
    
    # Log2FC (unchanged; from linear scale given log2 data)
    log2FC <- log2(mean(2^tumor_expr) / mean(2^normal_expr))
    
    # Cohen's d (unchanged)
    pooled_sd <- sqrt(((length(tumor_expr)-1)*stats::var(tumor_expr) + 
                         (length(normal_expr)-1)*stats::var(normal_expr)) / 
                        (length(tumor_expr) + length(normal_expr) - 2))
    cohens_d <- (mean(tumor_expr) - mean(normal_expr)) / pooled_sd
    
    # Statistical tests (unchanged)
    t_test <- stats::t.test(tumor_expr, normal_expr)
    wilcox_test <- stats::wilcox.test(tumor_expr, normal_expr)
    
    # === Added from Code1: raw Hedges' g for independent samples ===
    n1 <- length(tumor_expr)
    n2 <- length(normal_expr)
    d_smd <- cohens_d
    J_unadj <- if (is.finite(n1 + n2) && (4*(n1 + n2) - 9) != 0) J_df(n1 + n2) else NA_real_
    hedges_g <- if (is.finite(d_smd) && is.finite(J_unadj)) J_unadj * d_smd else NA_real_
    
    # === Added from Code1: covariate-adjusted standardized effect and unbiased g ===
    beta_T <- coef_summary['sample_typeT', 'Estimate']
    sigma_hat <- sm$sigma
    df_resid  <- stats::df.residual(fit)
    es_adj <- if (is.finite(sigma_hat) && sigma_hat > 0) beta_T / sigma_hat else NA_real_
    g_adj  <- J_df(df_resid) * es_adj  # unbiased standardized adjusted effect size
    
    # Return results (original fields preserved, new fields appended)
    data.frame(
      protein = protein,
      cancer = cancer_type,
      n_tumor = length(tumor_expr),
      n_normal = length(normal_expr),
      mean_tumor = mean(tumor_expr),
      mean_normal = mean(normal_expr),
      log2FC = log2FC,
      cohens_d = cohens_d,
      hedges_g = hedges_g,            # (raw, unadjusted)
      es_adj = es_adj,                # (covariate-adjusted standardized effect)
      g_adj = g_adj,                  # (bias-corrected adjusted g)
      sigma = sigma_hat,              # (residual SD)
      df_resid = df_resid,            # (residual degrees of freedom)
      beta = coef_summary['sample_typeT', 'Estimate'],
      se = coef_summary['sample_typeT', 'Std. Error'],
      t_statistic = coef_summary['sample_typeT', 't value'],
      p_value_lm = coef_summary['sample_typeT', 'Pr(>|t|)'],
      p_value_ttest = t_test$p.value,
      p_value_wilcox = wilcox_test$p.value,
      stringsAsFactors = FALSE
    )
    
  }, mc.cores = 1)
  
  # Combine results (unchanged)
  de_results <- dplyr::bind_rows(de_results)
  
  # Multiple testing correction (unchanged)
  if (nrow(de_results) > 0) {
    de_results <- de_results %>%
      dplyr::mutate(
        p_adj_lm = stats::p.adjust(p_value_lm, method = 'BH'),
        p_adj_t = stats::p.adjust(p_value_ttest, method = 'BH'),
        p_adj_wilcox = stats::p.adjust(p_value_wilcox, method = 'BH')
      )
  }
  
  return(de_results)
}

#### main ----
analysis_df <- meta

# Run analysis for all cancer types
cancer_types <- unique(analysis_df$cancer_abbr[analysis_df$sample_type == 'T'])

de_results_list <- lapply(cancer_types, function(cancer) {
  cat(sprintf("\nAnalyzing %s...\n", cancer))
  perform_de_analysis(cancer, analysis_df, expr)
})
# saveRDS(de_results_list, 'de_results_list1.rds')
# de_results_list <- readRDS('de_results_list1.rds')

# Combine results
de_results_all <- dplyr::bind_rows(de_results_list[!sapply(de_results_list, is.null)]) %>% 
  inner_join(dfprot %>% rename(protein = Protein.Group))


de_results_all %<>%
  rename(cancer_abbr = cancer) %>% 
  mutate(
    significant.lm = p_adj_lm < 0.05 & abs(g_adj) > 0.5,
    significant.wilcox = p_adj_wilcox < 0.05 & abs(log2FC) > log2(2)
    # direction = sum(p_adj_lm < 0.05 & g_adj > 0.5, na.rm = TRUE),
    # n_downregulated.lm = sum(p_adj_lm < 0.05 & g_adj < -0.5, na.rm = TRUE),
    # n_upregulated = sum(p_adj_lm < 0.05 & log2FC > log2(2), na.rm = TRUE),
    # n_downregulated = sum(p_adj_lm < 0.05 & log2FC < -log2(2), na.rm = TRUE),
  )

pdf('T_N_volcano.pdf', width = 10, height =10)
ggplot(de_results_all, aes(x = g_adj, y = -log10(p_adj_lm))) +
  facet_wrap(~ cancer_abbr, scale = 'free') +
  geom_point(aes(color = (p_adj_lm < 0.05) & (abs(g_adj) >= 0.5))) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(size = 10), legend.position = 'top')
ggplot(de_results_all, aes(x = log2FC, y = -log10(p_adj_wilcox))) +
  facet_wrap(~ cancer_abbr, scale = 'free') +
  geom_point(aes(color = (p_adj_wilcox < 0.05) & (abs(log2FC) >= log2(2)))) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  theme(text = element_text(size = 10), legend.position = 'top')
graphics.off()


dep1 <- de_results_all %>% 
  filter(significant.lm) %>% 
  mutate(direction.lm = ifelse(g_adj > 0, 'Up', 'Down')) %>%
  count(cancer_abbr, direction.lm)

dep2 <- de_results_all %>% 
  filter(significant.wilcox) %>% 
  mutate(direction.wilcox = ifelse(log2FC > 0, 'Up', 'Down')) %>%
  count(cancer_abbr, direction.wilcox)

dep.stat <- rbind(
  dep1 %>% mutate(direction.lm = str_c(direction.lm, '.lm')) %>% rename(direction = direction.lm),
  dep2 %>% mutate(direction.wilcox = str_c(direction.wilcox, '.wilcox')) %>% rename(direction = direction.wilcox)
) %>% 
  pivot_wider(id_cols = cancer_abbr, names_from = direction, values_from = n)
rio::export(dep.stat, 'T_N_1vsn_compare_stat.xlsx')

# 
# # Effect size distribution
# effect_size_plot <- de_results_all %>%
#   dplyr::filter(p_adj_lm < 0.05) %>%
#   ggplot2::ggplot(ggplot2::aes(x = cancer, y = g_adj, fill = cancer)) +
#   ggplot2::geom_violin(alpha = 0.7) +
#   ggplot2::geom_boxplot(width = 0.2, alpha = 0.8) +
#   ggplot2::coord_flip() +
#   ggplot2::labs(
#     title = "Effect Sizes (Hedges'g) for Significant proteins",
#     x = 'Cancer Type',
#     y = "Hedges'g"
#   ) +
#   ggplot2::theme_minimal() +
#   ggplot2::theme(legend.position = 'none')
# 
# # P-value distribution
# pval_dist_plot <- de_results_all %>%
#   ggplot2::ggplot(ggplot2::aes(x = p_value_lm)) +
#   ggplot2::geom_histogram(bins = 50, fill = 'steelblue', alpha = 0.7) +
#   ggplot2::facet_wrap(~cancer, scales = 'free_y') +
#   ggplot2::labs(
#     title = 'P-value Distribution Across Cancer Types',
#     x = 'P-value',
#     y = 'Frequency'
#   ) +
#   ggplot2::theme_minimal()



## XX method2 T vs N1/N2/... (1 vs 1/1/...) XX -----
#
# This method treats Tumor ("T") as one level and each Normal subtype comes from `class_abbr`.

# Outputs:
#   - ANOVA table per protein for the grp factor (with BH-adjusted p-values across proteins)
#   - Post-hoc pairwise comparisons (Dunnett): T vs each Normal subtype (multiplicity via Dunnett; BH also available within protein)


# Returns a list with two tibbles: list(anova = <per-protein>, posthoc = <per-contrast>)
#
# Notes:
#   * Requires `emmeans` for post-hoc comparisons; `car` is used if available for Type II/III ANOVA.
#   * Normals are defined as samples with sample_type == 'N' and non-NA `class_abbr`.
#   * Tumors are samples with sample_type == 'T'; for a specific cancer_type, only those tumors are included.
#   * BH adjustment for ANOVA is done within-cancer across proteins handled by this function call.
#   * BH adjustment for post-hoc is done within each protein across its T-vs-Normal contrasts (optional; see code).
#   * two-sided Dunnett p-values: `p_value_dunnett`.
### helper ------
perform_de_analysis_method2 <- function(
    cancer_type,
    metadata,
    expression_matrix,
    group_col = "class_abbr",
    # covariates = c("Dataset"),   # kept for interface compatibility; not used
    proteins = rownames(expression_matrix),
    min_per_group = 1,
    # anova_type = c("II", "III"), # kept for interface compatibility; not used
    adjust_method = "BH",
    n_cores = 1,
    max_proteins = Inf   # set to a number to limit proteins (refines prior 1:10 stub)
) {
  has_mvtnorm <- requireNamespace("mvtnorm", quietly = TRUE)
  if (!has_mvtnorm) stop("Package 'mvtnorm' is required for exact Dunnett adjustment.")
  
  # --- Select samples ---
  if (cancer_type == "ALL_TUMORS") {
    keep <- (metadata$sample_type == 'T') | (metadata$sample_type == 'N' & !is.na(metadata[[group_col]]))
  } else {
    keep <- (metadata$cancer_abbr == cancer_type & metadata$sample_type == 'T') |
      (metadata$sample_type == 'N' & !is.na(metadata[[group_col]]))
  }
  samples_meta <- metadata[keep, , drop = FALSE]
  
  # Intersect expression columns by FileName
  common_files <- intersect(colnames(expression_matrix), samples_meta$FileName)
  if (length(common_files) < 6) return(list(anova = dplyr::tibble(), posthoc = dplyr::tibble()))
  samples_meta <- samples_meta[match(common_files, samples_meta$FileName), , drop = FALSE]
  samples_expr <- expression_matrix[, common_files, drop = FALSE]
  
  # Sanity counts
  n_tumor <- sum(samples_meta$sample_type == 'T')
  n_norm  <- sum(samples_meta$sample_type == 'N')
  if (n_tumor < min_per_group || n_norm < min_per_group) {
    return(list(anova = dplyr::tibble(), posthoc = dplyr::tibble()))
  }
  
  # Per-protein analysis
  protein_indices <- which(rownames(samples_expr) %in% proteins)
  if (length(protein_indices) == 0) return(list(anova = dplyr::tibble(), posthoc = dplyr::tibble()))
  
  res_list <- parallel::mclapply(protein_indices, function(i) {
    protein <- rownames(samples_expr)[i]
    cat(sprintf("Analyzing %s...
", protein))
    
    # Build analysis data (log2 scale as provided)
    grp <- ifelse(samples_meta$sample_type == 'T', 'T', as.character(samples_meta[[group_col]]))
    grp <- factor(grp)
    model_data <- data.frame(
      expression = as.numeric(samples_expr[i, ]),
      grp = grp,
      stringsAsFactors = FALSE
    )
    
    # Drop NAs in response
    keep_row <- !is.na(model_data$expression)
    model_data <- model_data[keep_row, , drop = FALSE]
    if (nrow(model_data) < 2*min_per_group) return(NULL)
    model_data$grp <- droplevels(model_data$grp)
    
    # Enforce minimum per-group size and presence of Tumor level
    tab_grp <- table(model_data$grp)
    small_lvls <- names(tab_grp)[tab_grp < min_per_group]
    if (length(small_lvls) > 0) {
      model_data <- model_data[!model_data$grp %in% small_lvls, , drop = FALSE]
      model_data$grp <- droplevels(model_data$grp)
    }
    if (!('T' %in% levels(model_data$grp)) || length(levels(model_data$grp)) < 2) return(NULL)
    
    # --- Sufficient statistics by group ---
    split_by <- split(model_data$expression, model_data$grp)
    n_g <- vapply(split_by, function(x) sum(!is.na(x)), numeric(1))
    mean_g <- vapply(split_by, function(x) mean(x, na.rm = TRUE), numeric(1))
    var_g <- vapply(split_by, function(x) stats::var(x, na.rm = TRUE), numeric(1))
    # If a group has n <= 1, var() returns NA; but (n-1)*var should be 0 -> coerce safely
    var_g[n_g <= 1] <- 0
    groups <- names(n_g)
    G <- length(n_g)
    N <- sum(n_g)
    if (!('T' %in% groups)) return(NULL)
    
    # One-way ANOVA (cell-means) from sums of squares
    grand_mean <- sum(n_g * mean_g) / N
    SSB <- sum(n_g * (mean_g - grand_mean)^2)
    SSE <- sum((n_g - 1) * var_g)  # finite by construction
    df1 <- G - 1
    df2 <- N - G
    
    F_val <- (SSB / df1) / (SSE / max(df2, 1))  # safe when df2>0; unused otherwise
    p_val <- if (df2 > 0 && SSE > 0) stats::pf(F_val, df1, df2, lower.tail = FALSE) else NA_real_
    anova_row <- dplyr::tibble(
      protein = protein,
      cancer = cancer_type,
      df1 = as.numeric(df1),
      df2 = as.numeric(df2),
      F_value = as.numeric(F_val),
      p_value_anova = as.numeric(p_val),
      n_total = N,
      n_groups = G
    )
    
    # --- Exact Dunnett (single-step, multivariate t) ---
    posthoc_tbl <- dplyr::tibble()
    normals <- setdiff(groups, 'T')
    finite_guard <- is.finite(df2) && df2 > 0 && is.finite(SSE) && SSE > 0
    if (length(normals) >= 1 && finite_guard) {
      n_T <- n_g['T']
      m_T <- mean_g['T']
      x_T <- model_data$expression[model_data$grp == 'T']
      # linear-scale means / medians for fold-change summaries
      m_lin_T <- mean(2^x_T, na.rm = TRUE)
      med_lin_T <- stats::median(2^x_T, na.rm = TRUE)
      median_T <- stats::median(x_T, na.rm = TRUE)
      
      # pooled variance and SEs
      s2_p <- SSE / df2
      
      # Build correlation matrix among all T - N_k contrasts
      v_vec <- numeric(length(normals)); names(v_vec) <- normals
      for (k in normals) v_vec[k] <- (1 / n_T) + (1 / n_g[k])
      m <- length(normals)
      R <- diag(m)
      if (m > 1) {
        for (a in seq_len(m)) for (b in seq_len(m)) if (a != b) {
          k1 <- normals[a]; k2 <- normals[b]
          R[a, b] <- (1 / n_T) / sqrt(v_vec[k1] * v_vec[k2])
        }
      }
      
      # Compute estimates, SEs, t-stats
      t_vals <- se_vals <- est_vals <- rep(NA_real_, m)
      names(t_vals) <- names(se_vals) <- names(est_vals) <- normals
      # Also store medians and GM-based log2FC
      mean_norm <- median_norm <- rep(NA_real_, m)
      names(mean_norm) <- names(median_norm) <- normals
      log2FC_lin <- log2_median_FC <- rep(NA_real_, m)
      names(log2FC_lin) <- names(log2_median_FC) <- normals
      
      for (k in normals) {
        x_N  <- model_data$expression[model_data$grp == k]
        m_N  <- mean_g[k]
        mean_norm[k] <- m_N
        median_norm[k] <- stats::median(x_N, na.rm = TRUE)
        med_lin_N <- stats::median(2^x_N, na.rm = TRUE)
        m_lin_N <- mean(2^x_N, na.rm = TRUE)
        log2FC_lin[k] <- tryCatch(log2(m_lin_T / m_lin_N), error = function(e) NA_real_)
        log2_median_FC[k] <- tryCatch(log2(med_lin_T / med_lin_N), error = function(e) NA_real_)
        
        est <- (m_T - m_N)
        se  <- sqrt(s2_p * v_vec[k])
        t   <- est / se
        est_vals[k] <- est
        se_vals[k]  <- se
        t_vals[k]   <- t
      }
      
      # Non-adjusted (per-contrast, two-sided) p-values
      p_raw <- 2 * stats::pt(-abs(t_vals), df = df2)
      names(p_raw) <- normals
      
      # # Single-step Dunnett adjusted p-values via multivariate t (Genz & Bretz)
      # p_adj <- rep(NA_real_, m)
      # for (j in seq_len(m)) {
      #   cval <- abs(t_vals[j])
      #   if (!is.finite(cval)) { p_adj[j] <- NA_real_; next }
      #   if (m == 1) {
      #     p_adj[j] <- 2 * stats::pt(-cval, df = df2)
      #   } else {
      #     prob <- as.numeric(mvtnorm::pmvt(
      #       lower = rep(-cval, m), upper = rep(cval, m),
      #       df = df2, corr = R,
      #       algorithm = mvtnorm::GenzBretz(maxpts = 25000, abseps = 1e-06, releps = 0),
      #       keepAttr = FALSE, seed = 123
      #     ))
      #     p_adj[j] <- 1 - prob
      #   }
      # }
      
      posthoc_tbl <- dplyr::tibble(
        protein = protein,
        cancer = cancer_type,
        normal_level = normals,
        n_tumor = as.numeric(n_T),
        n_normal = as.numeric(n_g[normals]),
        mean_tumor = as.numeric(m_T),
        mean_normal = as.numeric(mean_norm[normals]),
        median_tumor = as.numeric(median_T),
        median_normal = as.numeric(median_norm[normals]),
        log2FC = as.numeric(log2FC_lin[normals]),
        log2_median_FC = as.numeric(log2_median_FC[normals]),
        estimate_T_minus_N = as.numeric(est_vals[normals]),
        se_dunnett_pooled = as.numeric(se_vals[normals]),
        t_statistic_dunnett = as.numeric(t_vals[normals]),
        df = as.numeric(df2),
        p_value_dunnett = as.numeric(p_raw[normals])#,              # unadjusted
        # p_value_dunnett_exact = as.numeric(p_adj)                   # adjusted (single-step)
      )
    }
    
    list(anova = anova_row, posthoc = posthoc_tbl)
  }, mc.cores = n_cores)
  
  # Combine
  anova_res  <- dplyr::bind_rows(lapply(res_list, `[[`, "anova"))
  posthoc_res <- dplyr::bind_rows(lapply(res_list, `[[`, "posthoc"))
  
  # Adjustments across proteins / within protein
  if (nrow(anova_res) > 0) {
    anova_res <- anova_res %>%
      dplyr::mutate(p_adj_anova = stats::p.adjust(p_value_anova, method = adjust_method))
  }
  if (nrow(posthoc_res) > 0) {
    posthoc_res <- posthoc_res %>%
      dplyr::group_by(cancer, protein) %>%
      # Keep BH across Dunnett-adjusted p-values (optional)
      # dplyr::mutate(p_adj_within_protein = stats::p.adjust(p_value_dunnett_exact, method = adjust_method)) %>%
      dplyr::ungroup()
  }
  
  return(list(anova = anova_res, posthoc = posthoc_res))
}

### main ------
# analysis_df <- meta
# 
# # Tumor cancer types (same logic as Method 1)
# cancer_types <- unique(analysis_df$cancer_abbr[analysis_df$sample_type == "T"])
# 
# # Run analysis for all cancer types
# res_m2_list <- lapply(cancer_types, function(cancer) {
#   cat(sprintf("\nAnalyzing %s...\n", cancer))
#   perform_de_analysis_method2(
#     cancer_type       = cancer,
#     metadata          = analysis_df,
#     expression_matrix = expr,
#     group_col         = "class_abbr",     # Normal subtypes live here
#     covariates        = c("Dataset"),     # include if it varies; else ignored
#     proteins          = rownames(expr),   # or a subset
#     min_per_group     = 3,
#     anova_type        = "II",             # or "III"
#     adjust_method     = "BH",
#     n_cores           = 1
#   )
# })
# 
# # Drop NULL/empty results
# res_m2_list <- Filter(function(x) {
#   !is.null(x) && (nrow(x$anova) + nrow(x$posthoc) > 0)
# }, res_m2_list)
# 
# # Combine across cancers
# anova_all   <- dplyr::bind_rows(lapply(res_m2_list, `[[`, "anova"))
# posthoc_all <- dplyr::bind_rows(lapply(res_m2_list, `[[`, "posthoc"))
# 
# # Optional: global BH across all cancers (additive columns; do not overwrite)
# if (nrow(anova_all) > 0) {
#   anova_all <- anova_all %>%
#     dplyr::mutate(p_adj_anova_global = stats::p.adjust(p_value_anova, method = "BH"))
# }
# if (nrow(posthoc_all) > 0) {
#   posthoc_all <- posthoc_all %>%
#     dplyr::mutate(p_adj_global = stats::p.adjust(p_value, method = "BH"))
# }
# 
# # Examples: significance filters (choose one or both as needed)
# anova_sig        <- dplyr::filter(anova_all, p_adj_anova < 0.05)                 # within-cancer (per-function) BH
# posthoc_sig_protein <- dplyr::filter(posthoc_all, p_adj_within_protein < 0.05)         # within-protein BH
# posthoc_sig_glob <- dplyr::filter(posthoc_all, p_adj_global < 0.05)              # global BH across all contrasts


# 3.Intensity 2-fold -----
# expr %>% t() %>% as.data.frame() %>% 
#   rownames_to_column('FileName') %>% 
#   inner_join()
# meta %>% select(FileName, sample_type, cancer_abbr)

# vec_ca <- sort(unique(de_results_all$cancer_abbr))
# for(ca in vec_ca){
#   submeta <- meta %>% filter(sample_type == 'N' | cancer_abbr == ca) %>% 
#     mutate(cancer_abbr = ifelse(cancer_abbr == 'NORMAL', str_c(cancer_abbr, '.', class_abbr), cancer_abbr))
#   subexpr <- expr[, submeta$FileName]
#   dfsub <- subexpr %>% t() %>% as.data.frame() %>% 
#     rownames_to_column('FileName') %>% 
#     inner_join(submeta %>% select(FileName, cancer_abbr), .)
#   subexpr.med <- dfsub %>% select(-FileName) %>% group_by(cancer_abbr) %>% 
#     summarise_all(function(x) median(2^x))
#   
# }

meta.med <- meta %>%
  mutate(cancer_abbr = ifelse(cancer_abbr == 'NORMAL', str_c(cancer_abbr, '.', class_abbr), cancer_abbr))
df.med <- expr %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>%
  inner_join(meta.med %>% select(FileName, cancer_abbr), .)
expr.med <- df.med %>% select(-FileName) %>% group_by(cancer_abbr) %>%
  summarise_all(function(x) log2(median(2^x)))

# n.long.med <- expr.med %>% filter(str_detect(cancer_abbr, '^NORMAL\\.')) %>%
#   pivot_longer(cols = -cancer_abbr, names_to = 'protein', values_to = 'log2.n') %>% 
#   rename(N.class = cancer_abbr)
# t.long.med <- expr.med %>% filter(!str_detect(cancer_abbr, '^NORMAL\\.')) %>%
#   pivot_longer(cols = -cancer_abbr, names_to = 'protein', values_to = 'log2.t') %>% 
#   rename(T.class = cancer_abbr)
# 
# t.enrich <- split(t.long.med, t.long.med$T.class) %>%
#   lapply(function(X) {
#     X %>% inner_join(n.long.med) %>% mutate(log2.t.minus.n = log2.t - log2.n)
#   })
# df.t.enrich <- bind_rows(t.enrich)
# # saveRDS(df.t.enrich, 'df.t.enrich.rds')
# 
# t.enriched <- df.t.enrich %>%
#   filter(log2.t.minus.n > log2(2)) %>% 
#   group_by(T.class, protein) %>%
#   count(name = 'enriched') %>% 
#   filter(enriched == length(unique(df.t.enrich$N.class)))
long <- expr.med %>%
  pivot_longer(-cancer_abbr, names_to = "protein", values_to = "log2x")

# Count distinct normal classes overall
N_total <- long %>% filter(startsWith(cancer_abbr, "NORMAL.")) %>%
  distinct(cancer_abbr) %>% nrow()

# Per-protein n_max and its class
n_max_tbl <- long %>%
  filter(startsWith(cancer_abbr, "NORMAL.")) %>%
  group_by(protein) %>%
  summarise(n_max = max(log2x, na.rm = TRUE), .groups = "drop")

n_max_class_tbl <- long %>%
  filter(startsWith(cancer_abbr, "NORMAL.")) %>%
  group_by(protein) %>%
  slice_max(log2x, n = 1, with_ties = FALSE) %>%
  transmute(protein, N.class_max = cancer_abbr)

n_thr <- n_max_tbl %>%
  inner_join(n_max_class_tbl, by = "protein")

t.enrich <- long %>%
  filter(!startsWith(cancer_abbr, "NORMAL.")) %>%
  transmute(T.class = cancer_abbr, protein, log2.t = log2x) %>%
  inner_join(n_thr, by = "protein") %>%
  mutate(delta = log2.t - n_max) %>%
  arrange(T.class, desc(delta))
# saveRDS(t.enrich, 'output/t.enrich.rds')
# t.enrich <- readRDS('output/t.enrich.rds')

t.enriched <- t.enrich %>%
  filter(delta > log2(2)) %>%
  mutate(cancer_abbr = T.class, .before = 1)

t.de.en <- de_results_all %>% inner_join(t.enriched)



stat1 <- t.de.en %>% 
  filter(significant.lm) %>% 
  mutate(direction.lm = ifelse(g_adj > 0, 'Up', 'Down')) %>%
  filter(direction.lm == 'Up') %>% 
  count(cancer_abbr, direction.lm)

stat2 <- t.de.en %>% 
  filter(significant.wilcox) %>% 
  mutate(direction.wilcox = ifelse(log2FC > 0, 'Up', 'Down')) %>%
  filter(direction.wilcox == 'Up') %>% 
  count(cancer_abbr, direction.wilcox)

tn.en.stat <- rbind(
  stat1 %>% mutate(direction.lm = str_c(direction.lm, '.lm')) %>% rename(direction = direction.lm),
  stat2 %>% mutate(direction.wilcox = str_c(direction.wilcox, '.wilcox')) %>% rename(direction = direction.wilcox)
) %>% 
  pivot_wider(id_cols = cancer_abbr, names_from = direction, values_from = n)
rio::export(tn.en.stat, 'T_N_1vsn_compare_enrich_stat.xlsx')


# 4.combine with T-NT results -----
df_dep <- rio::import('../6_3paired_tumor_nontumor/output/tumor_nontumor_analysis_reports_v1202.xlsx', sheet = 3)


tntn.de.en <- t.de.en %>%
  filter(significant.lm, significant.wilcox, g_adj > 0.5, log2FC > log2(2)) %>% 
  inner_join(df_dep %>% filter(direction == 'Up') %>% select(cancer_abbr, protein))
top1 <- tntn.de.en %>%
  filter(p_adj_lm < 0.001, p_adj_t < 0.001, p_adj_wilcox < 0.001) %>% 
  group_by(cancer_abbr) %>%
  arrange(desc(delta)) %>% slice(1) %>% 
  ungroup() %>% 
  mutate(label = str_c(protein, '_', Genes))


dfbox <- datall[unique(top1$protein), ] %>% t() %>% as.data.frame() %>% 
  rownames_to_column('FileName') %>% 
  inner_join(
    info.all %>%
      filter(sample_type %in% c('T', 'NT', 'N')) %>%
      mutate(cancer_subtype = cancer_abbr, .after = class_abbr) %>% 
      mutate(cancer_abbr = ifelse(str_detect(cancer_abbr, '^BRCA'), 'BRCA', cancer_abbr),
             cancer_abbr = ifelse(str_detect(cancer_abbr, 'CCOC|HGSOC'), 'OC', cancer_abbr)) %>% 
      mutate(cancer = ifelse(str_detect(cancer, '^Breast carcinoma'), 'Breast carcinoma', cancer),
             cancer = ifelse(str_detect(cancer, 'ovarian cancer$'), 'Ovarian cancer', cancer)) %>% 
      mutate(Dataset = c('data1' = 'DS1', 'data2' = 'DS2', 'data3' = 'DS3')[new],
             .before = new) %>% 
      select(FileName, sample_type, class_abbr, cancer_abbr)
  , .) %>% 
  pivot_longer(cols = -c(FileName:cancer_abbr), names_to = 'protein', values_to = 'Log2Intensity') %>% 
  as.data.frame() %>% 
  left_join(top1 %>% distinct(protein, label))

box.plots <- list()
for(ca in setdiff(vec_ca, 'NORMAL')){
  sub.top1 <- top1 %>% filter(cancer_abbr %in% ca)
  subox <- dfbox %>% filter(sample_type == 'N' | cancer_abbr == ca) %>%
    filter(protein %in% sub.top1$protein) %>%
    mutate(cancer_abbr = ifelse(cancer_abbr == 'NORMAL', str_c(cancer_abbr, '.', class_abbr), cancer_abbr))
  sample_color <- c(T = '#CF3F54', NT = '#79AFD1', N ='#03A2B3')
  
  box.plots[[ca]] <- ggplot(subox)+
    aes(x = sample_type, y = Log2Intensity, color = sample_type) +
    # facet_wrap(~ cancer_abbr,
    #            ncol = 1, scales = "free_y", strip.position = "top") +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75)) +
    labs(x = "Sample type", y = "Log2Intensity", title = str_c(ca, ' - ', sub.top1$label)) +
    # stat_pvalue_manual(
    #   data = pvals, label = "p.label", size = 3,
    #   y.position = "y.position", xmin = "xmin", xmax = "xmax",
    #   hide.ns = FALSE, tip.length = 0,
    #   inherit.aes = FALSE
    # ) +
    scale_color_manual(values = sample_color) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      strip.text.y.right = element_text(size = 10, angle = 0, hjust = 0, margin = margin(l = 2)),
      # strip.background = element_blank(),
    )
}
length(box.plots)
# box.plots
# box.one <- wrap_plots(box.plots, ncol = 5)
box.one <- ggpubr::ggarrange(plotlist = box.plots, nrow = 5, ncol = 5, common.legend = T)

ggsave('T_NT_dep_T_N_dep_T_N_enrich_box_top1.pdf', box.one, width = 3*5, height = 3*5)

list(all = tntn.de.en, top1 = top1) %>%
  rio::export('T_NT_dep_T_N_dep_T_N_enrich_box.xlsx')


# 5.Combine results v2 ------
# T-N enriched results
t.enrich <- readRDS('output/t.enrich.rds')
t.enriched <- t.enrich %>%
  filter(delta > log2(2)) %>%
  mutate(cancer_abbr = T.class, .before = 1)

# T-N compare results
de_results_list <- readRDS('de_results_list1.rds')
de_results_all <- dplyr::bind_rows(de_results_list[!sapply(de_results_list, is.null)]) %>% 
  inner_join(dfprot %>% rename(protein = Protein.Group))

de_results_all %<>%
  rename(cancer_abbr = cancer) %>% 
  mutate(
    significant.lm = p_adj_lm < 0.05 & abs(g_adj) > 0.5,
    significant.wilcox = p_adj_wilcox < 0.05 & abs(log2FC) > log2(2)
    # direction = sum(p_adj_lm < 0.05 & g_adj > 0.5, na.rm = TRUE),
    # n_downregulated.lm = sum(p_adj_lm < 0.05 & g_adj < -0.5, na.rm = TRUE),
    # n_upregulated = sum(p_adj_lm < 0.05 & log2FC > log2(2), na.rm = TRUE),
    # n_downregulated = sum(p_adj_lm < 0.05 & log2FC < -log2(2), na.rm = TRUE),
  )

# T-NT compare results
df_dep <- rio::import('../6_3paired_tumor_nontumor/output/tumor_nontumor_analysis_reports_v1202.xlsx', sheet = 3)

# combine results
tntn.de.en <- de_results_all %>%
  inner_join(t.enriched) %>% 
  filter(p_adj_wilcox < 0.05, log2FC > log2(2)) %>% 
  inner_join(df_dep %>% filter(direction == 'Up') %>% select(cancer_abbr, protein))

# filter: normal.max < tumor.q1

### q1 ---------
# Tumor Q1 (type 7) per protein, per tumor group
tum <- df.med %>%
  select(-FileName) %>%
  filter(!startsWith(cancer_abbr, "NORMAL")) %>%
  group_by(cancer_abbr) %>%
  reframe({
    m <- as.matrix(pick(where(is.numeric)))
    tibble::as_tibble_row(setNames(
      as.numeric(matrixStats::colQuantiles(m, probs = 0.25, na.rm = TRUE, type = 7)),
      colnames(m)
    ))
  })

# Normal max per protein, per normal group
norm <- df.med %>%
  select(-FileName) %>%
  filter(startsWith(cancer_abbr, "NORMAL")) %>%
  group_by(cancer_abbr) %>%
  reframe({
    m <- as.matrix(pick(where(is.numeric)))
    tibble::as_tibble_row(setNames(
      matrixStats::colMaxs(m, na.rm = TRUE),
      colnames(m)
    ))
  })

# Long tables
tum_long  <- tum  %>% rename(tumor_abbr  = cancer_abbr)  %>%
  pivot_longer(-tumor_abbr,  names_to = "protein", values_to = "tumor_q1")
norm_long <- norm %>% rename(normal_abbr = cancer_abbr) %>%
  pivot_longer(-normal_abbr, names_to = "protein", values_to = "normal_max")

# Pairwise comparison: counts/proportion where tumor > normal
joint <- tum_long %>% inner_join(norm_long, by = "protein")

pair_counts <- joint %>%
  group_by(tumor_abbr, normal_abbr) %>%
  summarise(n = dplyr::n(),
            n_gt = sum(tumor_q1 > normal_max, na.rm = TRUE),
            p_gt = n_gt / n,
            .groups = "drop")
pair_counts

# proteins_per_pair <- joint %>%
#   filter(tumor_q1 > normal_max) %>%
#   group_by(tumor_abbr, normal_abbr) %>%
#   summarise(n_gt = dplyr::n(), proteins = list(protein), .groups = "drop")

# Collapse normals to one row per protein: the largest normal_max
normal_by_protein <- joint %>%
  distinct(protein, normal_abbr, normal_max) %>%                 # drop cross-join dupes
  group_by(protein) %>%
  summarise(any_normal_max = max(normal_max, na.rm = TRUE),
            n_normals      = sum(!is.na(normal_max)),
            .groups = "drop")

# Distinct tumor entries (one row per tumor_abbr-protein)
tumor_by_tp <- joint %>%
  distinct(tumor_abbr, protein, tumor_q1)

# Keep tumor-protein if tumor_q1 > max normal_max for that protein
keepers <- tumor_by_tp %>%
  inner_join(normal_by_protein, by = "protein") %>%
  filter(n_normals > 0, tumor_q1 > any_normal_max) %>%
  select(tumor_abbr, protein)

tq1.gt.nmax <- joint %>%
  semi_join(keepers)
tntn.final <- tntn.de.en %>%
  inner_join(tq1.gt.nmax, by = c(cancer_abbr = 'tumor_abbr', protein = 'protein'))
# saveRDS(tntn.final, 'tntn.final.q1.rds')
# saveRDS(tntn.final, 'tntn.final.q1.20251219.rds')
tntn.final <- readRDS('tntn.final.q1.20251219.rds')

top1 <- tntn.final %>%
  # filter(p_adj_lm < 0.001, p_adj_t < 0.001, p_adj_wilcox < 0.001) %>%
  # inner_join(df_dep %>% filter(t_test_paired.p < 0.05) %>% select(cancer_abbr, protein)) %>%
  group_by(cancer_abbr) %>%
  arrange(desc(delta)) %>% slice(1) %>%
  ungroup() %>%
  mutate(label = str_c(protein, '_', Genes)) %>% 
  inner_join(df_dep %>% select(cancer_abbr, protein, p_adj_BH))


#### boxplot -----
dfbox <- datall[unique(top1$protein), ] %>% t() %>% as.data.frame() %>% 
  rownames_to_column('FileName') %>% 
  inner_join(
    info.all %>%
      filter(sample_type %in% c('T', 'NT', 'N')) %>%
      mutate(cancer_subtype = cancer_abbr, .after = class_abbr) %>% 
      mutate(cancer_abbr = ifelse(str_detect(cancer_abbr, '^BRCA'), 'BRCA', cancer_abbr),
             cancer_abbr = ifelse(str_detect(cancer_abbr, 'CCOC|HGSOC'), 'OC', cancer_abbr)) %>% 
      mutate(cancer = ifelse(str_detect(cancer, '^Breast carcinoma'), 'Breast carcinoma', cancer),
             cancer = ifelse(str_detect(cancer, 'ovarian cancer$'), 'Ovarian cancer', cancer)) %>% 
      mutate(Dataset = c('data1' = 'DS1', 'data2' = 'DS2', 'data3' = 'DS3')[new],
             .before = new) %>% 
      select(FileName, sample_type, class_abbr, cancer_abbr)
    , .) %>% 
  pivot_longer(cols = -c(FileName:cancer_abbr), names_to = 'protein', values_to = 'Log2Intensity') %>% 
  as.data.frame() %>% 
  left_join(top1 %>% distinct(protein, label))

box.plots <- list()
for(ca in setdiff(vec_ca, 'NORMAL')){
  sub.top1 <- top1 %>% filter(cancer_abbr %in% ca)
  subox <- dfbox %>% filter(sample_type == 'N' | cancer_abbr == ca) %>%
    filter(protein %in% sub.top1$protein) %>%
    mutate(cancer_abbr = ifelse(cancer_abbr == 'NORMAL', str_c(cancer_abbr, '.', class_abbr), cancer_abbr)) %>% 
    mutate(sample_type = factor(sample_type, levels = c("N","NT","T")))
  sample_color <- c(T = '#CF3F54', NT = '#79AFD1', N ='#03A2B3')
  
  box.plots[[ca]] <- ggplot(subox)+
    aes(x = sample_type, y = Log2Intensity, color = sample_type) +
    # facet_wrap(~ cancer_abbr,
    #            ncol = 1, scales = "free_y", strip.position = "top") +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75)) +
    labs(x = "Sample type", y = "Log2Intensity", title = str_c(ca, ' - ', sub.top1$label)) +
    # stat_pvalue_manual(
    #   data = pvals, label = "p.label", size = 3,
    #   y.position = "y.position", xmin = "xmin", xmax = "xmax",
    #   hide.ns = FALSE, tip.length = 0,
    #   inherit.aes = FALSE
    # ) +
    ggpubr::stat_compare_means(
      comparisons = list(c("T","NT"), c("T","N"), c("NT","N")),
      method = "wilcox.test",           # or "t.test"
      p.adjust.method = "BH",
      label = "p.format",               # numeric p-values
      hide.ns = FALSE,
      tip.length = 0,
      step.increase = 0.15              # vertical spacing between brackets
    ) +
    scale_color_manual(values = sample_color) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      strip.text.y.right = element_text(size = 10, angle = 0, hjust = 0, margin = margin(l = 2)),
      # strip.background = element_blank(),
    )
}
# box.plots[[1]]
# box.one <- wrap_plots(box.plots, ncol = 5)
box.one <- ggpubr::ggarrange(plotlist = box.plots, nrow = 5, ncol = 5, common.legend = T)

ggsave('T_NT_dep_T_N_dep_T_N_enrich_box_q1_top1_20251202.pdf', box.one, width = 3*5, height = 3*5)


# 6.PRM overlap -----
# "Q15910" "Q02548"
df_prm <- rio::import('//192.168.99.100/share/members/jiangwenhao/TPHP/20220908/PRM/PUH_PRM_prot_info_20230313.xlsx')
prm <- readxl::read_excel('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/PRM/20230318_TPHP_PRM_dysregulated_filter50NAByOrgan_DEP.xlsx', sheet = 'all')
prm %<>% rename(class_abbr = organ) %>% inner_join(meta %>% distinct(class_abbr, cancer_abbr))

prm_prot <- intersect(tntn.final$protein, prm$prot)


## for DIA data -----
top_prm <- tntn.final %>%
  filter(protein %in% prm_prot) %>% 
  # filter(p_adj_lm < 0.001, p_adj_t < 0.001, p_adj_wilcox < 0.001) %>%
  # inner_join(df_dep %>% filter(t_test_paired.p < 0.05) %>% select(cancer_abbr, protein)) %>%
  mutate(label = str_c(protein, '_', Genes)) %>% 
  distinct(cancer_abbr, protein, Genes, label, p_adj_wilcox) %>% 
  inner_join(df_dep %>% select(cancer_abbr, protein, p_adj_BH) %>% mutate(cancer = cancer_abbr))

dfbox <- datall[prm_prot, ] %>% t() %>% as.data.frame() %>% 
  rownames_to_column('FileName') %>% 
  inner_join(
    info.all %>%
      filter(sample_type %in% c('T', 'NT', 'N')) %>%
      mutate(cancer_subtype = cancer_abbr, .after = class_abbr) %>% 
      mutate(cancer_abbr = ifelse(str_detect(cancer_abbr, '^BRCA'), 'BRCA', cancer_abbr),
             cancer_abbr = ifelse(str_detect(cancer_abbr, 'CCOC|HGSOC'), 'OC', cancer_abbr)) %>% 
      mutate(cancer = ifelse(str_detect(cancer, '^Breast carcinoma'), 'Breast carcinoma', cancer),
             cancer = ifelse(str_detect(cancer, 'ovarian cancer$'), 'Ovarian cancer', cancer)) %>% 
      mutate(Dataset = c('data1' = 'DS1', 'data2' = 'DS2', 'data3' = 'DS3')[new],
             .before = new) %>% 
      select(FileName, sample_type, class_abbr, cancer_abbr)
    , .) %>% 
  pivot_longer(cols = -c(FileName:cancer_abbr), names_to = 'protein', values_to = 'Log2Intensity') %>% 
  as.data.frame() %>% 
  left_join(top_prm %>% distinct(protein, label))

box.plots <- list()
for(i in 1:nrow(top_prm)){
  ca <- top_prm$cancer_abbr[i]
  sub.top1 <- top_prm %>% slice(i)
  subox <- dfbox %>% filter(sample_type == 'N' | cancer_abbr == ca) %>%
    filter(protein %in% sub.top1$protein) %>%
    mutate(cancer_abbr = ifelse(cancer_abbr == 'NORMAL', str_c(cancer_abbr, '.', class_abbr), cancer_abbr)) %>% 
    mutate(sample_type = factor(sample_type, levels = c("N","NT","T")))
  sample_color <- c(T = '#CF3F54', NT = '#79AFD1', N ='#03A2B3')
  
  ## --- add: manual p-value labels for T–NT and T–N ---
  yr <- range(subox$Log2Intensity, na.rm = TRUE)
  ymax <- yr[2]; yspan <- diff(yr); if(!is.finite(yspan) || yspan == 0) yspan <- 1
  pvals <- tibble::tibble(
    group1    = c("T","T"),
    group2    = c("NT","N"),
    p.label   = c(
      paste0(signif(sub.top1$p_adj_BH, 3)),
      paste0(signif(sub.top1$p_adj_wilcox, 3))
    ),
    y.position = ymax + c(0.15, 0.30) * yspan
  )
  
  box.plots[[i]] <- ggplot(subox)+
    aes(x = sample_type, y = Log2Intensity, color = sample_type) +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75)) +
    labs(x = "Sample type", y = "Log2Intensity", title = str_c(ca, ' - ', sub.top1$label)) +
    ## --- use the manual labels for the two requested pairs ---
    ggpubr::stat_pvalue_manual(
      data = pvals, label = "p.label", size = 3,
      y.position = "y.position", xmin = "group1", xmax = "group2",
      hide.ns = FALSE, tip.length = 0, inherit.aes = FALSE
    ) +
    ## --- keep an automatic test only for NT–N
    # ggpubr::stat_compare_means(
    #   comparisons = list(c("NT","N")),
    #   method = "t.test",
    #   p.adjust.method = "BH",
    #   label = "p.format",
    #   hide.ns = FALSE,
    #   tip.length = 0,
    #   step.increase = 0.15
    # ) +
    scale_color_manual(values = sample_color) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      strip.text.y.right = element_text(size = 10, angle = 0, hjust = 0, margin = margin(l = 2))
    )
}
# box.plots[[1]]
# box.one <- wrap_plots(box.plots, ncol = 5)
length(box.plots)
box.one <- ggpubr::ggarrange(plotlist = box.plots, nrow = 1, ncol = 3, common.legend = T)

ggsave('T_NT_dep_T_N_dep_T_N_enrich_box_q1_PRM_validated_DIA_20251202.pdf', box.one, width = 3*3, height = 3*1)


## for PRM data -----
# top_prm <- tntn.final %>%
#   filter(protein %in% prm_prot) %>%
#   # filter(p_adj_lm < 0.001, p_adj_t < 0.001, p_adj_wilcox < 0.001) %>%
#   # inner_join(df_dep %>% filter(t_test_paired.p < 0.05) %>% select(cancer_abbr, protein)) %>%
#   mutate(label = str_c(protein, '_', Genes)) %>%
#   distinct(cancer_abbr, protein, Genes, label, p_adj_wilcox) %>%
#   inner_join(prm %>% select(cancer_abbr, prot, p_value) %>% rename(protein = prot))

dfbox <- df_prm %>% filter(!is.na(sample_type), anatomical_classification == 'lymph node') %>% 
  mutate(sample_type = c(carcinoma = 'T', adjacent = 'NT', Normal = 'N')[sample_type]) %>% 
  select(sample_type, all_of(prm_prot)) %>% 
  mutate(cancer_abbr = 'DLBCL', .before = 1) %>% 
  pivot_longer(cols = all_of(prm_prot), names_to = 'protein', values_to = 'Log2Intensity') %>% 
  mutate(Log2Intensity = ifelse(!is.na(Log2Intensity), Log2Intensity, min(Log2Intensity, na.rm = T) + log2(0.5)))

box.plots <- list()
for(i in 1:nrow(top_prm)){
  ca <- top_prm$cancer_abbr[i]
  sub.top1 <- top_prm %>% slice(i)
  subox <- dfbox %>% filter(sample_type == 'N' | cancer_abbr == ca) %>%
    filter(protein %in% sub.top1$protein) %>%
    mutate(cancer_abbr = ifelse(cancer_abbr == 'NORMAL', str_c(cancer_abbr, '.', class_abbr), cancer_abbr)) %>% 
    mutate(sample_type = factor(sample_type, levels = c("N","NT","T")))
  sample_color <- c(T = '#CF3F54', NT = '#79AFD1', N ='#03A2B3')
  
  # ## --- add: manual p-value labels for T–NT and T–N ---
  # yr <- range(subox$Log2Intensity, na.rm = TRUE)
  # ymax <- yr[2]; yspan <- diff(yr); if(!is.finite(yspan) || yspan == 0) yspan <- 1
  # pvals <- tibble::tibble(
  #   group1    = c("T","T"),
  #   group2    = c("NT","N"),
  #   p.label   = c(
  #     paste0(signif(sub.top1$p_adj_BH, 3)),
  #     paste0(signif(sub.top1$p_adj_wilcox, 3))
  #   ),
  #   y.position = ymax + c(0.15, 0.30) * yspan
  # )
  
  box.plots[[i]] <- ggplot(subox)+
    aes(x = sample_type, y = Log2Intensity, color = sample_type) +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75)) +
    labs(x = "Sample type", y = "Log2Intensity", title = str_c(ca, ' - ', sub.top1$label)) +
    # ## --- use the manual labels for the two requested pairs ---
    # ggpubr::stat_pvalue_manual(
    #   data = pvals, label = "p.label", size = 3,
    #   y.position = "y.position", xmin = "group1", xmax = "group2",
    #   hide.ns = FALSE, tip.length = 0, inherit.aes = FALSE
    # ) +
    # # --- keep an automatic test only for NT–N
    # ggpubr::stat_compare_means(
    #   comparisons = list(c("NT","N")),
    #   method = "t.test",
    #   p.adjust.method = "BH",
    #   label = "p.format",
    #   hide.ns = FALSE,
    #   tip.length = 0,
    #   step.increase = 0.15
    # ) +
    ggpubr::stat_compare_means(
      comparisons = list(c('NT', 'T'), c('N', 'T')),
      method = "wilcox",
      p.adjust.method = "BH",
      label = "p.format",
      hide.ns = FALSE,
      tip.length = 0,
      step.increase = 0.15
    ) +
    scale_color_manual(values = sample_color) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      strip.text.y.right = element_text(size = 10, angle = 0, hjust = 0, margin = margin(l = 2))
    )
}
length(box.plots) # 3

box.one <- ggpubr::ggarrange(plotlist = box.plots, nrow = 1, ncol = 3, common.legend = T)

ggsave('T_NT_dep_T_N_dep_T_N_enrich_box_q1_PRM_validated_PRM_20251202.pdf', box.one, width = 3*3, height = 3*1)



# 7.Tumor-enriched mapped to others -------
## 7.1 Location -----
tntn.final <- readRDS('tntn.final.q1.20251219.rds')
hpa_loc <- rio::import('input/subcellular_location.tsv.HPA_v24.1_Ensembl_v109/subcellular_location.tsv')
# hpa_loc %>% filter(`Gene name` %in% tntn.final$Genes)
# writeClipboard(unique(unlist(str_split(unique(hpa_loc$`Main location`), ';'))))
# unique(unlist(str_split(str_subset(unique(hpa_loc$`Main location`), '[Mm]embrane', negate = T), ';')))
# str_subset(unique(hpa_loc$`Main location`), '[Mm]embrane')
# Plasma membrane; Endoplasmic reticulum; Golgi apparatus; Endosomes; Lysosomes; Peroxisomes; Vesicles; Mitochondria; Nuclear membrane
hpa_pm_loc <- hpa_loc %>% filter(str_detect(`Main location`, 'Plasma membrane'),
                                 `Gene name` %in% tntn.final$Genes)
tntn.final.pm <- tntn.final %>%
  select(-(normal_abbr:normal_max)) %>% distinct() %>% 
  filter(Genes %in% hpa_pm_loc$`Gene name`) %>% 
  mutate(label = str_c(protein, '_', Genes)) %>% 
  inner_join(df_dep %>% select(cancer_abbr, protein, p_adj_BH) %>% mutate(cancer = cancer_abbr)) %>% 
  arrange(protein, cancer_abbr)

intersect(prm_prot, tntn.final.pm$protein) # character(0)

dfbox <- datall[unique(tntn.final.pm$protein), ] %>% t() %>% as.data.frame() %>% 
  rownames_to_column('FileName') %>% 
  inner_join(
    info.all %>%
      filter(sample_type %in% c('T', 'NT', 'N')) %>%
      mutate(cancer_subtype = cancer_abbr, .after = class_abbr) %>% 
      mutate(cancer_abbr = ifelse(str_detect(cancer_abbr, '^BRCA'), 'BRCA', cancer_abbr),
             cancer_abbr = ifelse(str_detect(cancer_abbr, 'CCOC|HGSOC'), 'OC', cancer_abbr)) %>% 
      mutate(cancer = ifelse(str_detect(cancer, '^Breast carcinoma'), 'Breast carcinoma', cancer),
             cancer = ifelse(str_detect(cancer, 'ovarian cancer$'), 'Ovarian cancer', cancer)) %>% 
      mutate(Dataset = c('data1' = 'DS1', 'data2' = 'DS2', 'data3' = 'DS3')[new],
             .before = new) %>% 
      select(FileName, sample_type, class_abbr, cancer_abbr)
    , .) %>% 
  pivot_longer(cols = -c(FileName:cancer_abbr), names_to = 'protein', values_to = 'Log2Intensity') %>% 
  as.data.frame() %>% 
  left_join(tntn.final.pm %>% distinct(protein, label))

box.plots <- list()
for(i in 1:nrow(tntn.final.pm)){
  sub.pm <- tntn.final.pm %>% slice(i)
  prot <- sub.pm$protein
  ca <- sub.pm$cancer_abbr
  
  subox <- dfbox %>% filter(sample_type == 'N' | cancer_abbr == ca) %>%
    filter(protein %in% sub.pm$protein) %>%
    mutate(cancer_abbr = ifelse(cancer_abbr == 'NORMAL', str_c(cancer_abbr, '.', class_abbr), cancer_abbr)) %>% 
    mutate(sample_type = factor(sample_type, levels = c("N","NT","T")))
  sample_color <- c(T = '#CF3F54', NT = '#79AFD1', N ='#03A2B3')
  
  ## --- add: manual p-value labels for T–NT and T–N ---
  yr <- range(subox$Log2Intensity, na.rm = TRUE)
  ymax <- yr[2]; yspan <- diff(yr); if(!is.finite(yspan) || yspan == 0) yspan <- 1
  pvals <- tibble::tibble(
    group1    = c("T","T"),
    group2    = c("NT","N"),
    p.label   = c(
      paste0(signif(sub.pm$p_adj_BH, 3)),
      paste0(signif(sub.pm$p_adj_wilcox, 3))
    ),
    y.position = ymax + c(0.15, 0.30) * yspan
  )
  
  box.plots[[i]] <- ggplot(subox)+
    aes(x = sample_type, y = Log2Intensity, color = sample_type) +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75)) +
    labs(x = "Sample type", y = "Log2Intensity", title = str_c(ca, ' - ', sub.pm$label)) +
    ## --- use the manual labels for the two requested pairs ---
    ggpubr::stat_pvalue_manual(
      data = pvals, label = "p.label", size = 3,
      y.position = "y.position", xmin = "group1", xmax = "group2",
      hide.ns = FALSE, tip.length = 0, inherit.aes = FALSE
    ) +
    ## --- keep an automatic test only for NT–N
    # ggpubr::stat_compare_means(
    #   comparisons = list(c("NT","N")),
    #   method = "t.test",
    #   p.adjust.method = "BH",
    #   label = "p.format",
    #   hide.ns = FALSE,
    #   tip.length = 0,
    #   step.increase = 0.15
    # ) +
    scale_color_manual(values = sample_color) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      strip.text.y.right = element_text(size = 10, angle = 0, hjust = 0, margin = margin(l = 2))
    )
}
length(box.plots) # 56
# saveRDS(box.plots, 'box.plots.plasma.membrane.rds')
box.one <- ggpubr::ggarrange(plotlist = box.plots, nrow = 7, ncol = 8, common.legend = T)
ggsave('T_NT_dep_T_N_dep_T_N_enrich_box_q1_Plasma_membrane_20251202.pdf', box.one, width = 3*8, height = 3*7)
rio::export(tntn.final.pm, 'T_NT_dep_T_N_dep_T_N_enrich_box_q1_Plasma_membrane_20251202.xlsx')


## 7.2 Potiential drug targets -----
pdt <- rio::import('input/HPA_protein_class_Potential_drug_targets.tsv')

tntn.final.pdt <- tntn.final %>%
  select(-(normal_abbr:normal_max)) %>% distinct() %>% 
  filter(Genes %in% pdt$Gene) %>% 
  mutate(label = str_c(protein, '_', Genes)) %>% 
  inner_join(df_dep %>% select(cancer_abbr, protein, p_adj_BH) %>% mutate(cancer = cancer_abbr)) %>% 
  arrange(protein, cancer_abbr)
intersect(prm_prot, tntn.final.pdt$protein) # character(0)

dfbox <- datall[unique(tntn.final.pdt$protein), ] %>% t() %>% as.data.frame() %>% 
  rownames_to_column('FileName') %>% 
  inner_join(
    info.all %>%
      filter(sample_type %in% c('T', 'NT', 'N')) %>%
      mutate(cancer_subtype = cancer_abbr, .after = class_abbr) %>% 
      mutate(cancer_abbr = ifelse(str_detect(cancer_abbr, '^BRCA'), 'BRCA', cancer_abbr),
             cancer_abbr = ifelse(str_detect(cancer_abbr, 'CCOC|HGSOC'), 'OC', cancer_abbr)) %>% 
      mutate(cancer = ifelse(str_detect(cancer, '^Breast carcinoma'), 'Breast carcinoma', cancer),
             cancer = ifelse(str_detect(cancer, 'ovarian cancer$'), 'Ovarian cancer', cancer)) %>% 
      mutate(Dataset = c('data1' = 'DS1', 'data2' = 'DS2', 'data3' = 'DS3')[new],
             .before = new) %>% 
      select(FileName, sample_type, class_abbr, cancer_abbr)
    , .) %>% 
  pivot_longer(cols = -c(FileName:cancer_abbr), names_to = 'protein', values_to = 'Log2Intensity') %>% 
  as.data.frame() %>% 
  left_join(tntn.final.pdt %>% distinct(protein, label))

box.plots <- list()
for(i in 1:nrow(tntn.final.pdt)){
  sub.pdt <- tntn.final.pdt %>% slice(i)
  prot <- sub.pdt$protein
  ca <- sub.pdt$cancer_abbr

  subox <- dfbox %>% filter(sample_type == 'N' | cancer_abbr == ca) %>%
    filter(protein %in% sub.pdt$protein) %>%
    mutate(cancer_abbr = ifelse(cancer_abbr == 'NORMAL', str_c(cancer_abbr, '.', class_abbr), cancer_abbr)) %>%
    mutate(sample_type = factor(sample_type, levels = c("N","NT","T")))
  sample_color <- c(T = '#CF3F54', NT = '#79AFD1', N ='#03A2B3')

  ## --- add: manual p-value labels for T–NT and T–N ---
  yr <- range(subox$Log2Intensity, na.rm = TRUE)
  ymax <- yr[2]; yspan <- diff(yr); if(!is.finite(yspan) || yspan == 0) yspan <- 1
  pvals <- tibble::tibble(
    group1    = c("T","T"),
    group2    = c("NT","N"),
    p.label   = c(
      paste0(signif(sub.pdt$p_adj_BH, 3)),
      paste0(signif(sub.pdt$p_adj_wilcox, 3))
    ),
    y.position = ymax + c(0.15, 0.30) * yspan
  )

  box.plots[[i]] <- ggplot(subox)+
    aes(x = sample_type, y = Log2Intensity, color = sample_type) +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75)) +
    labs(x = "Sample type", y = "Log2Intensity", title = str_c(ca, ' - ', sub.pdt$label)) +
    ## --- use the manual labels for the two requested pairs ---
    ggpubr::stat_pvalue_manual(
      data = pvals, label = "p.label", size = 3,
      y.position = "y.position", xmin = "group1", xmax = "group2",
      hide.ns = FALSE, tip.length = 0, inherit.aes = FALSE
    ) +
    ## --- keep an automatic test only for NT–N
    # ggpubr::stat_compare_means(
    #   comparisons = list(c("NT","N")),
    #   method = "t.test",
    #   p.adjust.method = "BH",
    #   label = "p.format",
    #   hide.ns = FALSE,
    #   tip.length = 0,
    #   step.increase = 0.15
    # ) +
    scale_color_manual(values = sample_color) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      strip.text.y.right = element_text(size = 10, angle = 0, hjust = 0, margin = margin(l = 2))
    )
}
length(box.plots) # 180
# box.plots[[1]]
# box.one <- wrap_plots(box.plots, ncol = 5)
box.one <- ggpubr::ggarrange(plotlist = box.plots, nrow = 12, ncol = 15, common.legend = T)

ggsave('T_NT_dep_T_N_dep_T_N_enrich_box_q1_potential_drug_targets_20251202.pdf', box.one, width = 3*15, height = 3*12)
rio::export(tntn.final.pdt, 'T_NT_dep_T_N_dep_T_N_enrich_box_q1_potential_drug_targets_20251202.xlsx')


# Output ----
res.q1 <- tntn.final %>%
  # filter(p_adj_lm < 0.001, p_adj_t < 0.001, p_adj_wilcox < 0.001) %>%
  # inner_join(df_dep %>% filter(t_test_paired.p < 0.05) %>% select(cancer_abbr, protein)) %>%
  group_by(cancer_abbr, protein) %>%
  arrange(desc(normal_max)) %>% slice(1) %>%
  ungroup() %>%
  mutate(label = str_c(protein, '_', Genes)) %>% 
  inner_join(df_dep %>% select(cancer_abbr, protein, p_adj_BH))

list(Tumor.enriched = res.q1,
     PRM.validation = top_prm,
     Plasma.membrane.loc = res.q1 %>% semi_join(tntn.final.pm),
     Potiential.druggble = res.q1 %>% semi_join(tntn.final.pdt)) %>% 
  rio::export('output/tumor_enriched_q1Filter_20251219.xlsx')




# check ----
# rio::export(de_results_all, 'T_N_1vsn_compare.xlsx')
# 
# # data.box <- expr %>% t() %>% as.data.frame() %>% 
# #   rownames_to_column('FileName') %>% 
# #   inner_join(meta, .)
# 
# 
# ca <- 'ESCA'
# vec_dep <- de_results_all %>%
#   filter(cancer_abbr == ca, significant.lm) %>%
#   pull(protein)
# 
# vmax <- 5
# # hm_cols <- colorRampPalette(brewer.pal(11, "RdYlBu")[c(10:6, 5:1)])(101)  # red→blue diverging
# hm_cols   <- viridisLite::viridis(101)          # set direction = -1 to invert
# hm_breaks <- seq(-vmax, vmax, length.out = 101)
# pheatmap::pheatmap(
#   expr[vec_dep, meta %>% filter(sample_type == 'N' | cancer_abbr == ca) %>% pull(FileName)],
#   scale = 'row',
#   color = hm_cols, breaks = hm_breaks,
#   cluster_rows = T, cluster_cols = T,
#   clustering_method = 'complete',
#   show_rownames = F, show_colnames = F,
#   # fontsize_row = 6, fontsize_col = 10,
#   border_color = NA, na_col = "#CCCCCC",
#   # annotation_row = ann_row,
#   # annotation_row = prm_speprot %>% column_to_rownames('protein') %>% select(cancer_abbr),
#   annotation_col = meta %>% column_to_rownames('FileName') %>% select(cancer_abbr, sample_type),
#   annotation_colors = list(
#     cancer_color[-which(names(cancer_color) %in% c('CCOC', 'HGSOC'))] %>%
#       append(c(cancer_color['HGSOC'], '#000000') %>% setNames(c('OC', 'NORMAL'))),
#     sample_type = sample_color
#   ),
#   fontsize = 7,
#   main = str_c("T-N matrix ", ca, ' - lm', collapse = ' '),
#   filename = "ESCA_compare_to_N_heatmap_lm.pdf", width = 8, height = 4
# )
# 
# vec_dep <- de_results_all %>%
#   filter(cancer_abbr == ca, significant.wilcox) %>%
#   pull(protein)
# 
# vmax <- 5
# # hm_cols <- colorRampPalette(brewer.pal(11, "RdYlBu")[c(10:6, 5:1)])(101)  # red→blue diverging
# hm_cols   <- viridisLite::viridis(101)          # set direction = -1 to invert
# hm_breaks <- seq(-vmax, vmax, length.out = 101)
# pheatmap::pheatmap(
#   expr[vec_dep, meta %>% filter(sample_type == 'N' | cancer_abbr == ca) %>% pull(FileName)],
#   scale = 'row',
#   color = hm_cols, breaks = hm_breaks,
#   cluster_rows = T, cluster_cols = T,
#   clustering_method = 'complete',
#   show_rownames = F, show_colnames = F,
#   # fontsize_row = 6, fontsize_col = 10,
#   border_color = NA, na_col = "#CCCCCC",
#   # annotation_row = ann_row,
#   # annotation_row = prm_speprot %>% column_to_rownames('protein') %>% select(cancer_abbr),
#   annotation_col = meta %>% column_to_rownames('FileName') %>% select(cancer_abbr, sample_type),
#   annotation_colors = list(
#     cancer_color[-which(names(cancer_color) %in% c('CCOC', 'HGSOC'))] %>%
#       append(c(cancer_color['HGSOC'], '#000000') %>% setNames(c('OC', 'NORMAL'))),
#     sample_type = sample_color
#   ),
#   fontsize = 7,
#   main = str_c("T-N matrix ", ca, ' - wilcox', collapse = ' '),
#   filename = "ESCA_compare_to_N_heatmap_wilcox.pdf", width = 8, height = 4
# )
# 
# # LARCA − Q9GZR2
# dfbox %>% count(sample_type)
# dfbox %>% filter(!FileName %in% meta$FileName) %>% count(sample_type)
# View(dfbox %>% filter(sample_type == 'N', protein == 'Q9GZR2'))
# View(df.med %>% filter(cancer_abbr == 'NORMAL.EYE.i') %>% select(cancer_abbr, Q9GZR2))
# View(norm %>% filter(cancer_abbr == 'NORMAL.EYE.i') %>% select(cancer_abbr, Q9GZR2))
# 
# 
# df_dep %>% filter(cancer_abbr == 'THCA', protein == 'Q15434')
# df_dep %>% filter(cancer_abbr == 'THYM', protein == 'Q9H7Z6')

