# 0. Packages ----
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Reproducibility & parallel
set.seed(1L)
options(stringsAsFactors = FALSE)
n_cores <- min(64L, parallel::detectCores())
BiocParallel::register(BiocParallel::SnowParam(workers = max(1L, n_cores - 2L)))

# Working dir & helpers
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('../source/source_code.R')

## helper functions -----
get_scicol <- function(n){
  if (requireNamespace("Polychrome", quietly=TRUE)) {
    Polychrome::createPalette(n, seedcolors = c(
      "#000000","#E69F00","#56B4E9","#009E73",
      "#F0E442","#0072B2","#D55E00","#CC79A7"
    ))[1:n]
  } else if (requireNamespace("colorspace", quietly=TRUE)) {
    colorspace::qualitative_hcl(n, palette = "Dark 3")
  } else {
    hcl(h = seq(15, 375, length.out = n+1)[-1], c = 60, l = 65)
  }
}

# 1. Config ------
# User paths (ASCII-only)
info_pan <- rio::import('../6_2pancancer_clustering_survival/input/PUH_pancancer_sample_information_2datesets_1146files.xlsx')
info.all <- rio::import('../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V6.xlsx') %>% 
  filter(!Tissue_heterogeneity_abnormal, !Low_protein_IDs, !FileName %in% pancancer_low_pearson)
datall <- rio::import('../2_batch_effect_evaluation_correction/output/batch_corrected_data_step2_limma3.csv') %>% column_to_rownames('V1')

out_dir   <- "subtype_cluster_V1009_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# 2. Load data ----
## 2.1 Pan-cancer data -----
meta <- info.all %>% filter(sample_type == 'T')
expr <- datall[, meta$FileName]


dim(meta) # 1070   34
dim(expr) # 12797  1070


# 3.Pan-cancer with BECA ----
# mat.rm.pool <- log2(as.matrix(expr))
mat.rm.pool <- as.matrix(expr)

# 4. Molecular subtypes Clustering --------
## 4.1 input -------
# CV top 50%
prot_cv <- apply(mat.rm.pool, 1, function(x){
  sd(x, na.rm = T) / mean(x, na.rm = T)
})
mat_consensus <- mat.rm.pool[prot_cv >= quantile(prot_cv, 0.50), ]
dim(mat_consensus) # 6399 1070
quantile(mat_consensus)
# 0%       25%       50%       75%      100% 
# 3.161037  6.434478  6.808211 11.031736 24.635597 

ann_df <- meta %>%
  select(FileName, tissue_name, class_abbr, cancer_abbr, instrument, new) %>%
  mutate(across(-FileName, \(x) factor(x))) %>%
  column_to_rownames("FileName")

ann_colors <- lapply(ann_df, function(f){
  lv <- levels(f)
  setNames(get_scicol(length(lv)), lv)
})

# 5.Features vio -----
## 5.1 Pan-cancer features ------
df_sil1 <- readRDS('df_sil1.rds')
# df_sil1_5a <- readRDS('df_sil1_5a.rds')

df_fea1 <- rio::import('ssgsea_pancancer_LM22_heatmap_v2_addImmTypeCls_addESTIMATE.xlsx',
                       sheet = 1)
df_fea2 <- rio::import('ssgsea_pancancer_LM22_heatmap_v2_addImmTypeCls_addESTIMATE.xlsx',
                       sheet = 2)
df_psu <- rio::import('../4_F-T-NT-N_analysis/output_chh/trajectory/step2_limma3_allsample_2728samples_info.xlsx')
setdiff(df_fea2$FileName, df_psu$FileName) # character(0)

df_fea <- df_fea2 %>% inner_join(df_psu %>% select(FileName, pseudotime))# %>% 
  # inner_join(
  #   df_fea1 %>% column_to_rownames('LM22') %>% t() %>% as.data.frame() %>% rownames_to_column('FileName')
  # )

### a.cancers plot -----
dfvio <- df_fea %>%
  select(cancer_abbr, cluster, `TGF-Beta`:last_col()) %>%
  pivot_longer(cols = `TGF-Beta`:last_col(), names_to = 'Feature', values_to = 'value') %>% 
  arrange(cluster) %>% 
  mutate(cancer_abbr = factor(cancer_abbr, levels = unique(cancer_abbr)),
         cluster = factor(cluster),
         Feature = factor(Feature, levels = colnames(df_fea)[-(1:4)]))
vlins_x <- df_fea %>% distinct(cancer_abbr, cluster) %>% count(cluster) %>% pull(n) %>% cumsum() %>% 
  `+`(0.5)

plot.vio <- ggplot(dfvio) +
  # facet_wrap(~interaction(cluster, Feature), scales = 'free', ncol = nlevels(dfvio$cluster)) +
  facet_wrap(~ Feature, scales = 'free_y', ncol = 1) +
  aes(x = cancer_abbr, y = value, fill = cancer_abbr) +
  # geom_jitter(color = '#CCCCCC', size=0.1, alpha=0.3,
  #             position = position_jitterdodge(jitter.width = 1, seed = 0)) +
  geom_violin(aes(color = cancer_abbr), width = 0.9) +
  # geom_vline(xintercept = vlins_x, linetype = "dashed",    # or "solid", "dotdash", etc.
  #            color = "#AAAAAA", linewidth = 0.8) +
  scale_color_manual(values = cancer_color) +
  scale_fill_manual(values = cancer_color) +
  geom_boxplot(width=0.1, color = '#000000', outlier.size = 0.5, show.legend = FALSE) +
  labs(x = '', y = '') +
  theme_classic() +
  theme(text = element_text(size = 15, color = "black"),
        axis.text.x.bottom = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
ggsave('pancancer_subtype_features_v1202.pdf', plot.vio, width = 25, height = 15, limitsize = FALSE)


### b.clusters plot -----
# plot.vio.cls <- ggplot(dfvio) +
#   # facet_wrap(~interaction(cluster, Feature), scales = 'free', ncol = nlevels(dfvio$cluster)) +
#   facet_wrap(~ Feature, scales = 'free_y', ncol = 1) +
#   aes(x = cluster, y = value, fill = cancer_abbr) +
#   # geom_jitter(color = '#CCCCCC', size=0.1, alpha=0.3,
#   #             position = position_jitterdodge(jitter.width = 1, seed = 0)) +
#   geom_violin(aes(color = cancer_abbr), width = 0.9) +
#   # geom_vline(xintercept = vlins_x, linetype = "dashed",    # or "solid", "dotdash", etc.
#   #            color = "#AAAAAA", linewidth = 0.8) +
#   scale_color_manual(values = cluster_color) +
#   scale_fill_manual(values = cluster_color) +
#   geom_boxplot(width=0.1, color = '#000000', outlier.size = 0.5, show.legend = FALSE) +
#   labs(x = '', y = '') +
#   theme_classic() +
#   theme(text = element_text(size = 15, color = "black"),
#         axis.text.x.bottom = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
# ggsave('pancancer_subtype_features_cluster_stat.pdf', plot.vio.cls, width = 25, height = 15, limitsize = FALSE)

# Check normality assumption for ANOVA
normality_test <- dfvio %>%
  dplyr::group_by(Feature, cluster) %>%
  rstatix::shapiro_test(value) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()

# Check homogeneity of variance
levene_test_results <- dfvio %>%
  dplyr::group_by(Feature) %>%
  rstatix::levene_test(value ~ cluster)

# One-way ANOVA for each Feature
anova_results <- dfvio %>%
  dplyr::group_by(Feature) %>%
  rstatix::anova_test(value ~ cluster) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance() %>%
  as.data.frame() %>% 
  dplyr::mutate(eta_squared = ges)

# Post-hoc analysis using Tukey HSD with PROPER position calculation
tukey_results <- dfvio %>%
  dplyr::group_by(Feature) %>%
  rstatix::tukey_hsd(value ~ cluster) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()

# Summary statistics
summary_stats <- dfvio %>%
  dplyr::group_by(Feature, cluster) %>%
  rstatix::get_summary_stats(value, type = "mean_sd") %>%
  dplyr::select(Feature, cluster, n, mean, sd)

generate_cld <- function(df, feature_name) {
  tryCatch({
    aov_model <- stats::aov(value ~ cluster, data = df)
    tukey_cld <- multcomp::glht(aov_model, linfct = multcomp::mcp(cluster = "Tukey"))
    cld_result <- multcomp::cld(tukey_cld, level = 0.05)
    
    letters_df <- data.frame(
      Feature = feature_name,
      cluster = names(cld_result$mcletters$Letters),
      cld = cld_result$mcletters$Letters,
      stringsAsFactors = FALSE
    )
    return(letters_df)
  }, error = function(e) {
    print(paste("Error in CLD for", feature_name, ":", e$message))
    return(data.frame(Feature = feature_name, cluster = NA, cld = NA))
  })
}

# Apply CLD to all features
cld_results <- dfvio %>%
  dplyr::group_split(Feature) %>%
  purrr::map_dfr(~ generate_cld(.x, unique(.x$Feature)))

# Get y-range statistics
y_range <- dfvio %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(
    max_val = max(value, na.rm = TRUE),
    min_val = min(value, na.rm = TRUE),
    range_val = max_val - min_val,
    .groups = 'drop'
  )

# Add positions to CLD results
cld_results <- cld_results %>%
  dplyr::left_join(y_range, by = "Feature") %>%
  dplyr::mutate(y_pos = max_val + range_val * 0.1)

# Prepare ANOVA label data
anova_label_data <- anova_results %>%
  dplyr::left_join(y_range, by = "Feature") %>%
  dplyr::mutate(label = paste0("ANOVA: p", 
                               ifelse(p < 0.001, "<0.001", 
                                      paste0("=", round(p, 3)))))

# Create main faceted plot
plot.vio.cls <- dfvio %>%
  dplyr::left_join(cld_results %>% dplyr::select(Feature, cluster, cld), 
                   by = c("Feature", "cluster")) %>%
  ggplot2::ggplot() +
  ggplot2::facet_wrap(~ Feature, scales = 'free_y', ncol = 1) +
  ggplot2::aes(x = cluster, y = value, fill = cluster) +
  ggplot2::geom_violin(ggplot2::aes(color = cluster), width = 0.9, alpha = 0.7) +
  ggplot2::scale_color_manual(values = cluster_color) +
  ggplot2::scale_fill_manual(values = cluster_color) +
  ggplot2::geom_boxplot(width = 0.1, color = '#000000', 
                        outlier.size = 0.5, show.legend = FALSE) +
  ggplot2::geom_text(data = cld_results %>% dplyr::filter(!is.na(cld)),
                     ggplot2::aes(x = cluster, y = y_pos, label = cld),
                     vjust = 0, size = 3, fontface = "bold", inherit.aes = FALSE) +
  ggplot2::geom_text(data = anova_label_data,
                     ggplot2::aes(x = 1, y = max_val * 0.95, label = label),
                     hjust = 0, vjust = 1, size = 3, fontface = "italic", 
                     inherit.aes = FALSE) +
  ggplot2::labs(x = '', y = '') +
  ggplot2::theme_classic() +
  ggplot2::theme(text = ggplot2::element_text(size = 15, color = "black"),
                 axis.text.x.bottom = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 0.5),
                 strip.background = ggplot2::element_blank(),
                 legend.position = "right")
ggplot2::ggsave('pancancer_subtype_features_cluster_stat_v1202_v2.pdf', plot.vio.cls, width = 5, height = 20, limitsize = FALSE)

# Export statistical results
list(anova = anova_results,
     tukey = tukey_results,
     cld = cld_results,
     summary_stats = summary_stats) %>% 
  rio::export('pancancer_subtype_features_cluster_stat_v1202.xlsx')

## 5.2 T-NT paired data -----------

# Final check ------
dfvio %>% filter(Feature == 'pseudotime', cancer_abbr == 'PACA') %>% pull(value)
dfvio %>% filter(Feature == 'pseudotime', cancer_abbr == 'BRCA-HER2+') %>% pull(value)
dfvio %>% filter(Feature == 'TGF-Beta', cancer_abbr == 'BRCA-HER2+') %>% pull(value) %>% sort()
dfvio %>% 
  filter(Feature == 'TGF-Beta') %>% 
  ggplot() +
  # facet_wrap(~interaction(cluster, Feature), scales = 'free', ncol = nlevels(dfvio$cluster)) +
  facet_wrap(~ Feature, scales = 'free', ncol = 1) +
  aes(x = cancer_abbr, y = value, fill = cancer_abbr) +
  # geom_jitter(color = '#CCCCCC', size=0.1, alpha=0.3,
  #             position = position_jitterdodge(jitter.width = 1, seed = 0)) +
  geom_violin(aes(color = cancer_abbr), width = 0.9) +
  # geom_vline(xintercept = vlins_x, linetype = "dashed",    # or "solid", "dotdash", etc.
  #            color = "#AAAAAA", linewidth = 0.8) +
  scale_color_manual(values = cancer_color) +
  scale_fill_manual(values = cancer_color) +
  coord_cartesian(ylim = c(0.0003, 0.0104)) +
  geom_boxplot(width=0.1, color = '#000000', outlier.size = 0.5, show.legend = FALSE) +
  labs(x = '', y = '') +
  theme_classic() +
  theme(text = element_text(size = 15, color = "black"),
        axis.text.x.bottom = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

meta$FileName[meta$cancer_abbr == 'BRCA-HER2+']
df_fea2 %>% filter(FileName %in% meta$FileName[meta$cancer_abbr == 'BRCA-HER2+'])
length(meta$FileName) # 1070
length(df_fea2$FileName) # 969

