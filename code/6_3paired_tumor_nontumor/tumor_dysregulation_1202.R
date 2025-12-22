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


# 1.Config ------
## 1.1 Read data -----
info.all <- rio::import('../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V5.xlsx') %>% 
  mutate(Low_protein_IDs = ifelse(is.na(Low_protein_IDs), FALSE, Low_protein_IDs))

datall <- rio::import('../0_process_DIA-NN_data/input/mapped_pg_matrix_2856_13609.csv') %>% 
  column_to_rownames('V1') %>% t() %>% log2()

datapool <- rio::import('../0_process_DIA-NN_data/input/pool_pg_matrix_149_13609.csv') %>% 
  column_to_rownames('V1') %>% t() %>% log2()

info.qc <- rio::import('../0_QC/output/QC2856_identity_20251201_source_check.xlsx')

## 1.2 Filter data -----
# QC filter
# info.qc.filter <- info.qc %>% filter(!(Is.Lower.Ingroup.class & Is.Lower.Ingroup.tissue))
info.qc.filter <- info.qc %>% filter(!(Is.Lower.Ingroup.major & Is.Lower.Ingroup.detailed))

# T and NT
meta.all <- info.all %>%
  # filter(!HE.abnormal, !Low_protein_IDs) %>%
  filter(FileName %in% info.qc.filter$FileName) %>% 
  filter(sample_type %in% c('T', 'NT', 'p'), !FileName %in% pool_N) %>% 
  mutate(cancer_subtype = cancer_abbr, .after = class_abbr) %>% 
  mutate(cancer_abbr = ifelse(str_detect(cancer_abbr, '^BRCA'), 'BRCA', cancer_abbr),
         cancer_abbr = ifelse(str_detect(cancer_abbr, 'CCOC|HGSOC'), 'OC', cancer_abbr)) %>% 
  mutate(cancer = ifelse(str_detect(cancer, '^Breast carcinoma'), 'Breast carcinoma', cancer),
         cancer = ifelse(str_detect(cancer, 'ovarian cancer$'), 'Ovarian cancer', cancer)) %>% 
  mutate(Dataset = c('data1' = 'DS1', 'data2' = 'DS2', 'data3' = 'DS3')[new],
         .before = new)

# remove pooling
meta0 <- meta.all %>% filter(sample_type != 'p')
cls_ca_map <- meta0 %>% select(class_abbr, cancer_abbr, cancer) %>% distinct() %>% drop_na()
meta0.a <- meta0 %>% filter(!is.na(cancer))
meta0.b <- meta0 %>% filter(is.na(cancer)) %>%
  select(-cancer, -cancer_abbr) %>% inner_join(cls_ca_map) %>%
  mutate(cancer_subtype = ifelse(is.na(cancer_subtype), cancer_abbr, cancer_subtype)) %>% 
  select(all_of(colnames(meta0)))
meta0 <- rbind(meta0.a, meta0.b)

expr <- datall[rowSums(!is.na(datall)) != 0, meta0$FileName] %>% 
  removeRowsAllNa()

dim(meta0) # 2146   36
dim(expr) # 13427  2146


# 
# rio::export(meta0, 'output/T_NT_info_check0.xlsx')
# meta0 %>% select(cancer_abbr, sample_type, patient_ID, Gender, Age) %>% 
#   str()


## 1.3 T-NT data QC ------
# identity cutoff = max(q1 - 1.5 * iqr, 0.6 * q2);
#   but not used for filter in 2.1
# pearson'r cutoff = 0.6

### correlation calculate -----
cor.methods <- c("pearson", "spearman")[1]
cor.list <- list()
for(i in seq_along(cor.methods)){
  cat('Calculate correlation - method', cor.methods[i], '...\n')
  cor.list[[i]] <- cor(expr, method = cor.methods[i], use = "pairwise.complete.obs")
}
names(cor.list) <- cor.methods
# saveRDS(cor.list, 'output/T_NT_qc_cor_list.rds')
# cor.list <- readRDS('output/T_NT_qc_cor_list.rds')

### Choose Pearson correlation matrix for display -----
cor_mat_pearson <- cor.list$pearson
qc.meta <- meta0 %>% select(FileName, sample_type, cancer_abbr, `# proteins`)

#### a) Mark low intra-group `# proteins` via robust IQR fence ----------
# Tukey-style lower fence within each group: Q1 - 1.5 * IQR
prot_cutoffs <- qc.meta %>%
  group_by(sample_type, cancer_abbr) %>%
  summarise(
    q1 = quantile(`# proteins`, 0.25, na.rm = TRUE),
    q2 = quantile(`# proteins`, 0.25, na.rm = TRUE),
    q3 = quantile(`# proteins`, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    cutoff = max(q1 - 1.5 * iqr, 0.6 * q2),
    .groups = "drop"
  )

qc.meta <- qc.meta %>%
  left_join(prot_cutoffs, by = c("sample_type", "cancer_abbr")) %>%
  mutate(low_proteins_intragroup = `# proteins` < cutoff) %>%
  # Keep original columns + the two QC flags (2nd flag is added below)
  select(FileName, sample_type, cancer_abbr, `# proteins`, low_proteins_intragroup)

#### b) Mark low intra-group correlation with iterative pairwise filtering ----------
# Build a tidy pairwise table (Pearson only) within each (sample_type, cancer_abbr),
# compute an IQR-based *pairwise* lower fence (kept here but final filtering uses a
# fixed cutoff column as in your code), then iteratively remove the sample that
# participates in the *lowest* correlation pair until the minimum pairwise r >= cutoff.
# Removed samples are *flagged* (low_corr_intragroup=TRUE); nothing is deleted.

.get_pairwise <- function(cor_mat, method_name){
  stopifnot(is.matrix(cor_mat))
  qc.meta %>%
    group_by(sample_type, cancer_abbr) %>%
    group_modify(function(df, key){
      samples <- df$FileName
      if (length(samples) <= 1L) {
        tibble(File1 = character(0), File2 = character(0), corr = numeric(0), method = character(0))
      } else {
        S <- cor_mat[samples, samples, drop = FALSE]
        idx <- which(upper.tri(S), arr.ind = TRUE)
        tibble(
          File1  = rownames(S)[idx[,"row"]],
          File2  = colnames(S)[idx[,"col"]],
          corr   = as.numeric(S[idx]),
          method = method_name
        )
      }
    }) %>%
    ungroup()
}

# Pearson only for filtering
pairwise_corr_long <- .get_pairwise(cor_mat_pearson, "pearson")
pairwise_corr_long %>% filter(corr < 0.6) %>% count(cancer_abbr) %>% arrange(desc(n))

# Pairwise IQR fences per group & method (kept for reference); use `cutoff` below
pair_cutoffs <- pairwise_corr_long %>%
  group_by(method, sample_type, cancer_abbr) %>%
  summarise(
    q1 = quantile(corr, 0.25, na.rm = TRUE),
    q3 = quantile(corr, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    cutoff0 = q1 - 1.5 * iqr,
    cutoff = 0.6, # 0.6 will be a reasonable cutoff for correlation (as in your code)
    .groups = "drop"
  )

# --- Iterative filtering helper: remove one sample per iteration based on the
#     current *lowest* pairwise correlation, until min r >= cutoff.
.iterative_flag_low_corr <- function(S, cutoff){
  flagged <- character(0)
  current <- colnames(S)
  while(length(current) > 1L){
    M <- S[current, current, drop = FALSE]
    diag(M) <- NA
    min_val <- suppressWarnings(min(M, na.rm = TRUE))
    if (!is.finite(min_val) || is.na(min_val) || min_val >= cutoff) break
    idx <- which(M == min_val, arr.ind = TRUE)[1, , drop = FALSE]
    s1 <- rownames(M)[idx[1, "row"]]
    s2 <- colnames(M)[idx[1, "col"]]
    others1 <- setdiff(current, s1)
    others2 <- setdiff(current, s2)
    m1 <- mean(M[s1, others1], na.rm = TRUE)
    m2 <- mean(M[s2, others2], na.rm = TRUE)
    # remove the one with lower mean correlation to the rest; deterministic tie-breaker
    rm <- if (!is.na(m1) && !is.na(m2)) {
      if (m1 < m2) s1 else if (m2 < m1) s2 else sort(c(s1, s2))[1]
    } else if (is.na(m1) && !is.na(m2)) {
      s1
    } else if (!is.na(m1) && is.na(m2)) {
      s2
    } else {
      sort(c(s1, s2))[1]
    }
    flagged <- c(flagged, rm)
    current <- setdiff(current, rm)
  }
  unique(flagged)
}

# Compute flags per group using the iterative rule and the per-group cutoff (default 0.6)
.groups_df <- qc.meta %>% distinct(sample_type, cancer_abbr)

iter_flags <- purrr::pmap_dfr(list(.groups_df$sample_type, .groups_df$cancer_abbr), function(st, ca){
  g_samples <- qc.meta %>% filter(sample_type == st, cancer_abbr == ca) %>% pull(FileName)
  if (length(g_samples) <= 1L) return(tibble(sample_type = st, cancer_abbr = ca, FileName = character(0), low_corr_intragroup = logical(0)))
  submat <- cor_mat_pearson[g_samples, g_samples, drop = FALSE]
  r_cutoff <- pair_cutoffs %>% filter(method == "pearson", sample_type == st, cancer_abbr == ca) %>% pull(cutoff)
  if (length(r_cutoff) == 0L || is.na(r_cutoff[1])) r_cutoff <- 0.6
  flagged <- .iterative_flag_low_corr(submat, r_cutoff[1])
  tibble(sample_type = st, cancer_abbr = ca, FileName = flagged, low_corr_intragroup = TRUE)
})

# Set NA for groups with <=1 sample; otherwise FALSE if not flagged by the iterative rule
.group_sizes <- qc.meta %>% count(sample_type, cancer_abbr, name = "cs_n")
qc.meta <- qc.meta %>%
  left_join(.group_sizes, by = c("sample_type", "cancer_abbr")) %>%
  left_join(iter_flags, by = c("sample_type", "cancer_abbr", "FileName")) %>%
  mutate(low_corr_intragroup = dplyr::case_when(
    cs_n <= 1 ~ NA,
    is.na(low_corr_intragroup) ~ FALSE,
    TRUE ~ TRUE
  )) %>%
  select(-cs_n)
# saveRDS(qc.meta, 'output/T_NT_QC_meta.rds')
# qc.meta <- readRDS('output/T_NT_QC_meta.rds')
# rio::export(qc.meta, 'output/T_NT_QC_infov5_based.xlsx')


#### c) Heatmaps per (sample_type, cancer_abbr) group ----------
# Draw clustered heatmaps of the *sample-by-sample Pearson correlation* within
# each group. Annotate rows/cols with the two QC flags (kept as factors to show legends).

# Helper: safe filename
.safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# Output dir
out_dir <- "T_NT_qc"
dir.create(out_dir, showWarnings = FALSE)

# List of groups
groups_df <- qc.meta %>% distinct(sample_type, cancer_abbr) %>% arrange(sample_type, cancer_abbr)

ann_clr <- list(low_proteins_intragroup = c('TRUE' = 'red4', 'FALSE' = 'green4'),
                low_corr_intragroup = c('TRUE' = 'red4', 'FALSE' = 'green4'))

# Multi-page PDF with all groups
pdf(file = file.path(out_dir, "TNT_intragroup_correlation_heatmaps.pdf"), width = 10, height = 10)
for (i in seq_len(nrow(groups_df))) {
  st <- groups_df$sample_type[i]
  ca <- groups_df$cancer_abbr[i]
  g_samples <- qc.meta %>% filter(sample_type == st, cancer_abbr == ca) %>% pull(FileName)
  
  if (length(g_samples) <= 1L) {
    plot.new()
    title(main = paste0(st, " | ", ca, " (n=", length(g_samples), ")"))
    next
  }
  
  submat <- cor_mat_pearson[g_samples, g_samples, drop = FALSE]
  
  ann <- qc.meta %>%
    filter(sample_type == st, cancer_abbr == ca) %>%
    select(FileName, low_proteins_intragroup, low_corr_intragroup) %>%
    mutate(
      low_proteins_intragroup = factor(ifelse(is.na(low_proteins_intragroup), "NA", as.character(low_proteins_intragroup)),
                                       levels = c("FALSE", "TRUE", "NA")),
      low_corr_intragroup     = factor(ifelse(is.na(low_corr_intragroup), "NA", as.character(low_corr_intragroup)),
                                       levels = c("FALSE", "TRUE", "NA"))
    ) %>%
    column_to_rownames("FileName")
  
  pheatmap::pheatmap(
    submat,
    color = viridisLite::cividis(79),
    breaks = seq(0.22, 1, length.out = 79),
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "ward.D2",
    border_color = NA,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = paste0(st, " | ", ca, " (n=", length(g_samples), ")"),
    annotation_col = ann,
    annotation_row = ann,
    annotation_colors = ann_clr,
    legend = TRUE
  )
}
graphics.off()


#### d) Intragroup heatmaps after excluding low-identity OR low-correlation samples -----
# Treat NA flags as pass (FALSE). Exclude if flagged low on either metric.
qc_pass <- qc.meta %>%
  dplyr::mutate(
    low_proteins_intragroup = dplyr::coalesce(low_proteins_intragroup, FALSE),
    low_corr_intragroup     = dplyr::coalesce(low_corr_intragroup, FALSE)
  )

# Multi-page PDF with pass-only groups
pdf(file = file.path(out_dir, "TNT_intragroup_correlation_heatmaps_PASS_ONLY.pdf"), width = 10, height = 10)
for (i in seq_len(nrow(groups_df))) {
  st <- groups_df$sample_type[i]
  ca <- groups_df$cancer_abbr[i]
  
  g_samples <- qc_pass %>%
    dplyr::filter(sample_type == st, cancer_abbr == ca,
                  !low_proteins_intragroup, !low_corr_intragroup) %>%
    dplyr::pull(FileName)
  
  if (length(g_samples) <= 1L) {
    plot.new()
    title(main = paste0(st, " | ", ca, " (pass-only; n=", length(g_samples), ")"))
    next
  }
  
  submat <- cor_mat_pearson[g_samples, g_samples, drop = FALSE]
  
  ann <- qc.meta %>%
    dplyr::filter(FileName %in% g_samples) %>%
    dplyr::select(FileName, low_proteins_intragroup, low_corr_intragroup) %>%
    dplyr::mutate(
      low_proteins_intragroup = factor(ifelse(is.na(low_proteins_intragroup), "NA", as.character(low_proteins_intragroup)),
                                       levels = c("FALSE", "TRUE", "NA")),
      low_corr_intragroup     = factor(ifelse(is.na(low_corr_intragroup), "NA", as.character(low_corr_intragroup)),
                                       levels = c("FALSE", "TRUE", "NA"))
    ) %>%
    tibble::column_to_rownames("FileName")
  
  pheatmap::pheatmap(
    submat,
    color = viridisLite::cividis(79),
    breaks = seq(0.22, 1, length.out = 79),
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "ward.D2",
    border_color = NA,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = paste0(st, " | ", ca, " (pass-only; n=", length(g_samples), ")"),
    annotation_col = ann,
    annotation_row = ann,
    annotation_colors = ann_clr,
    legend = TRUE
  )
}
graphics.off()

# 2.Data preparation ----
## 2.1 QC filter -----
qc.meta <- readRDS('output/T_NT_QC_meta.rds')

meta <- meta0 %>% inner_join(qc.meta) %>% 
  filter(#!low_proteins_intragroup,
         !FileName %in% pancancer_low_pearson,
         !low_corr_intragroup)
expr <- datall[rowSums(!is.na(datall)) != 0, meta$FileName] %>% 
  removeRowsAllNa()

dim(meta) # 2100   38
dim(expr) # 13427  2100

# rio::export(meta, 'output/T_NT_info_check.xlsx')
# saveRDS(meta, 'output/T_NT_info_check.rds')
# saveRDS(expr, 'output/expr.rds')

## 2.2 NA impute after sample aggregation -----
### helper functions ------
# Keeps a consistent set of proteins across groups (intersection of kept columns).
preproc_fill_na_by_group <- function(df,
                                     group_col,          # e.g., "subtype" / "organ" / "sample_type"
                                     na_cutoff = 0.50,   # drop columns with NA ratio > 50% within ANY group
                                     noise_sd  = 0.001,  # tiny jitter around the floor
                                     seed      = 2022,   # reproducibility
                                     min_log2 = NULL) {
  stopifnot(is.data.frame(df), group_col %in% names(df))
  set.seed(seed)
  
  # Separate metadata (non-double) and expression (double, assumed log2)
  num_idx   <- vapply(df, is.double, logical(1))
  expr_cols <- names(df)[num_idx]
  meta_cols <- setdiff(names(df), expr_cols)
  if (length(expr_cols) == 0) return(df)
  
  grp      <- df[[group_col]]
  grp_lvls <- unique(grp)
  
  # 1) Per-group NA filtering -> intersect columns that pass in ALL groups
  keep_list <- lapply(grp_lvls, function(g) {
    rows_g <- which(grp == g & !is.na(grp))
    if (!length(rows_g)) return(rep(FALSE, length(expr_cols)))
    colMeans(is.na(df[rows_g, expr_cols, drop = FALSE])) <= na_cutoff
  })
  keep_mask <- Reduce("&", keep_list)                # intersection across groups
  keep_cols <- expr_cols[keep_mask]
  if (length(keep_cols) == 0) return(df[, c(meta_cols), drop = FALSE])
  
  # 2) Within each group, impute remaining NAs using group's min_log2 - 1 (+ tiny jitter)
  out <- df
  for (g in grp_lvls) {
    rows_g <- which(grp == g & !is.na(grp))
    if (!length(rows_g)) next
    sub <- as.matrix(out[rows_g, keep_cols, drop = FALSE])
    n_na <- sum(is.na(sub))
    if (n_na > 0) {
      min_log2 <- ifelse(!is.na(min_log2), min_log2, suppressWarnings(min(sub, na.rm = TRUE)))
      if (is.finite(min_log2)) {
        if (noise_sd != 0){
          jit <- rnorm(ceiling(n_na * 1.1), mean = 1, sd = noise_sd)  # >0 to keep log2 valid
        } else {
          jit <- rep(1, ceiling(n_na * 1.1))
        }
        jit <- jit[jit > 0]
        imp <- min_log2 + log2(0.5) + log2(jit[seq_len(n_na)])
        sub[is.na(sub)] <- imp
        out[rows_g, keep_cols] <- sub
      }
    }
  }
  
  # Return metadata + retained expression columns
  out[, c(meta_cols, keep_cols), drop = FALSE]
}

df <- meta %>% select(FileName, patient_ID, sample_type, cancer_abbr, cancer_subtype, Gender, Age, Dataset) %>% 
  mutate(sample_type = factor(as.character(sample_type), levels = c("NT", "T")),
         Age = as.integer(Age),
         Gender = factor(Gender), Dataset = factor(Dataset)) %>% 
  inner_join(expr %>% t() %>% as.data.frame() %>% rownames_to_column('FileName'))
df %>% select(FileName:Dataset) %>% vapply(is.double, logical(1)) %>% any() # FALSE

df1 <- df %>%
  select(-FileName) %>% 
  pivot_longer(cols = -(patient_ID:Dataset)) %>% 
  pivot_wider(id_cols = patient_ID:Dataset, values_fn = function(x) log2(median(2^x, na.rm = T)))
dim(df1) # 1916 13434
# saveRDS(df1, 'output/df1_patientMatrix_1916_13434.rds')
# df1 <- readRDS('output/df1_patientMatrix_1916_13434.rds')

# View(df1 %>% count(patient_ID, sample_type) %>% filter(n == 2) %>% semi_join(df1, .))
# only "PUH_CA 1"

## Strictly enforce paired design in df1b: T & NT must both be present per Dataset × cancer_abbr × patient_ID
paired_ids <- df1 %>%
  dplyr::count(Dataset, cancer_abbr, patient_ID, sample_type) %>%
  tidyr::pivot_wider(names_from = sample_type, values_from = n, values_fill = 0) %>%
  dplyr::filter(`T` > 0, NT > 0) %>%
  dplyr::select(Dataset, cancer_abbr, patient_ID)
dim(paired_ids) # 890 3

df1b <- df1 %>%
  dplyr::inner_join(paired_ids)
dim(df1b) # 1780 13434

# remove contaminants
con.lib <- rio::import('../../code.20251201.archieved/1_1_check_counterintuitive_protein_coverage/output/blk_DDA/599blk_protein_0.1FreqPeptide.xlsx')
df1b <- df1b[, setdiff(colnames(df1b), con.lib$Protein.ID)]
dim(df1b) # 1780 13214

### b. paired test (LM + limma) -------------
quantile(as.matrix(df1b %>% select(-(patient_ID:Dataset))), na.rm = T)
# 0%       25%       50%       75%      100% 
# 3.674865 10.192197 11.452714 12.996761 24.948503 

# # Hedges' small-sample bias correction factor (J)
# J <- function(df) ifelse(is.finite(df) & df > 1, 1 - 3/(4*df - 1), NA_real_)

# Hedges J (exact, stable)
J <- function(df) {
  ifelse(
    is.finite(df) && df > 1,
    exp(lgamma(df/2) - 0.5*log(df/2) - lgamma((df-1)/2)),
    NA_real_
  )
}
# safe_sd  <- function(x) {
#   s <- stats::sd(x, na.rm = TRUE); ifelse(is.finite(s) && s > 0, s, NA_real_)
# }
safe_div <- function(a, b) {
  ifelse(is.finite(a) & is.finite(b) & b != 0, a / b, NA_real_)
}

min_log2_tnt <- min(df1b %>% select(-(1:Dataset)), na.rm = T)

vec_ca <- sort(unique(df1b$cancer_abbr))

res.all.list <- lapply(vec_ca, function(nm) {
  # nm <- 'BRCA'
  dfsub <- df1b %>% dplyr::filter(cancer_abbr == nm)
  
  # NA imputation per cancer_abbr (left-censoring style)
  dfsub_imp <- preproc_fill_na_by_group(
    df = dfsub, group_col = "cancer_abbr", na_cutoff = 0.50, noise_sd = 0.001, seed = 2022, min_log2 = min_log2_tnt
  )
  
  # Long then paired-wide with T/NT
  long_dfsub <- dfsub_imp %>%
    tidyr::pivot_longer(
      cols = -(patient_ID:Dataset),
      names_to = "protein",
      values_to = "value"
    )
  
  wide_dfsub <- long_dfsub %>%
    tidyr::pivot_wider(
      id_cols = c(cancer_abbr, cancer_subtype, patient_ID, Gender, Age, Dataset, protein),
      names_from = sample_type,
      values_from = value
    ) %>%
    tidyr::drop_na(`T`, NT)
  
  cat("\nAnalysis of cancer ", nm, "...\n")
  vec_feat <- sort(unique(wide_dfsub$protein))
  lmer_df <- lapply(vec_feat, function(feat) {
    cat("Protein: ", feat, "...\r")
    d <- dplyr::filter(long_dfsub, protein == feat) %>%
      dplyr::mutate(
        sample_type = factor(as.character(sample_type), levels = c("NT","T")),
        Gender      = factor(Gender),
        Dataset     = factor(Dataset),
        cancer_subtype = factor(cancer_subtype),
        Age_c       = Age - mean(Age, na.rm = TRUE)
      )
    
    # need repeated measures to estimate (1|patient_ID)
    if (dplyr::n_distinct(d$patient_ID) < 3) return(data.frame())
    
    # --- make fixed-effect terms subset-suitable (mirrors code2 pattern) ---
    n <- nrow(d)
    st_ok <- length(unique(stats::na.omit(d$sample_type))) >= 2
    if (!st_ok) return(data.frame())
    
    g2         <- all(c("F", "M") %in% unique(stats::na.omit(d$Gender)))
    subtype_ok <- length(unique(stats::na.omit(d$cancer_subtype))) >= 2
    ds_ok      <- length(unique(stats::na.omit(d$Dataset))) >= 2
    age_ok <- {
      x <- d$Age_c
      x <- x[is.finite(x)]
      length(unique(stats::na.omit(x))) >= 2
    }
    
    opt_terms <- c(
      if (subtype_ok) "cancer_subtype" else NA_character_,
      if (g2)         "Gender"        else NA_character_,
      if (age_ok)     "Age_c"          else NA_character_,
      if (ds_ok)      "Dataset"        else NA_character_
    )
    opt_terms <- opt_terms[!is.na(opt_terms)]
    
    # Limit complexity by sample size
    k_opt <- if (n >= 12) 3 else if (n >= 8) 2 else 1
    opt_terms <- head(opt_terms, k_opt)
    
    fixed_terms <- c("sample_type", opt_terms)
    fml_chr <- paste0(
      "value ~ 1 + ", paste(fixed_terms, collapse = " + "),
      " + (1 | patient_ID)"
    )
    # --- end subset-suitable term logic ---
    
    fml <- as.formula(fml_chr)
    
    fit <- tryCatch(
      suppressWarnings(
        suppressMessages(
          lmerTest::lmer(fml, data = d, REML = TRUE)
        )
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(data.frame())
    
    sing <- lme4::isSingular(fit, tol = 1e-4)
    re_var <- tryCatch({
      vc <- as.data.frame(lme4::VarCorr(fit))
      vc$vcov[vc$grp == "patient_ID"][1]
    }, error = function(e) NA_real_)
    
    sm <- summary(fit)
    coef_name <- "sample_typeT"
    if (!coef_name %in% rownames(sm$coefficients)) return(data.frame())
    
    est <- sm$coefficients[coef_name, "Estimate"]
    se  <- sm$coefficients[coef_name, "Std. Error"]
    t   <- sm$coefficients[coef_name, "t value"]
    df  <- sm$coefficients[coef_name, "df"]
    p   <- sm$coefficients[coef_name, "Pr(>|t|)"]
    
    sig <- sigma(fit)
    es  <- safe_div(est, sig)
    g   <- J(df) * es
    
    data.frame(
      cancer_abbr   = nm,
      protein       = feat,
      effect   = as.numeric(est),
      se       = as.numeric(se),
      t        = as.numeric(t),
      df       = as.numeric(df),
      p        = as.numeric(p),
      sigma    = as.numeric(sig),
      es_adj   = as.numeric(es),
      g_adj    = as.numeric(g),
      formula  = fml_chr,
      is_singular   = as.logical(sing),
      re_var_patient= as.numeric(re_var),
      stringsAsFactors = FALSE
    )
  }) %>% plyr::ldply() %>%
    dplyr::mutate(p_adj_BH = p.adjust(p, method = "BH"), .after = p)
  
  ret <- list(res = lmer_df,
              impute.data = long_dfsub)
  return(ret)
})
names(res.all.list) <- vec_ca
# saveRDS(res.all.list, 'output/res.list_paired_v1202.rds')
# res.all.list <- readRDS('output/res.list_paired_v1202.rds')

res.df <- dplyr::bind_rows(lapply(res.all.list, function(D) D$res))
impute.data <- plyr::ldply(res.all.list, function(D) D$impute.data) %>%
  pivot_wider(id_cols = c(patient_ID, Dataset, sample_type,
                          cancer_abbr, cancer_subtype, Gender, Age),
              names_from = 'protein') %>%
  mutate(ID = str_c(patient_ID, Dataset, '_', sample_type), .before = 1)



### c. post-process ------
length(setdiff(res.df$protein, dfprot$Protein.Group)) == 0 # TRUE
res_main <- res.df %>% 
  mutate(direction = ifelse(g_adj > 0, "Up", "Down")) %>% 
  inner_join(dfprot %>% rename(protein = Protein.Group))

# filter
res_main_filter <- res_main %>%
  filter(#abs(Log2FC) >= log2(2), abs(effect) >= 0.5,
    abs(g_adj) >= 0.5, p_adj_BH < 0.05)

# write an explanation for the columns in a data frame
instru_cn <- data.frame(
  field = c(
    "cancer_abbr","protein","effect","se","t","df","p","p_adj_BH",
    "sigma","es_adj","g_adj","formula","direction","Genes"
  ),
  meaning = c(
    "Cancer abbreviation / cohort identifier",
    "UniProt accession ID",
    "Estimated T–NT effect (difference) in log2 units from the fitted model",
    "Standard error of the estimated effect",
    "t statistic testing effect = 0",
    "Degrees of freedom used for the t test (Satterthwaite/KR df for mixed models)",
    "Raw (unadjusted) p-value for the effect",
    "Benjamini–Hochberg adjusted p-value (FDR) across proteins within cancer_abbr",
    "Model residual SD (sigma)",
    "Model-standardized effect size: es_adj = effect / sigma",
    "Small-sample bias-corrected standardized effect (Hedges-type): g_adj = J(df) * es_adj, where J(v) = Gamma(v/2) / (sqrt(v/2) * Gamma((v-1)/2))",
    "Model formula actually used for this protein/cancer subset",
    "Direction label derived from sign of effect (Up if effect>0, Down if effect<0)",
    "Mapped gene name(s) corresponding to the UniProt ID (if available in UniProt Knowledgebase)"
  ),
  stringsAsFactors = FALSE
)

# output
res.lm.tbls <- list(
  Instruction = instru_cn,
  Diff.report = res_main,
  Diff.report.filter = res_main_filter
)
rio::export(res.lm.tbls, 'output/tumor_nontumor_analysis_reports_v1202.xlsx')
saveRDS(res.lm.tbls, 'output/tumor_nontumor_analysis_reports_v1202.rds')
# res.lm.tbls <- loadRDS('output/tumor_nontumor_analysis_reports_v1202.rds')

rio::export(impute.data, 'output/T_NT_analysis_imputated_data_v1202.xlsx')
saveRDS(impute.data, 'output/T_NT_analysis_imputated_data_v1202.rds')
# impute.data <- readRDS('output/T_NT_analysis_imputated_data_v1202.rds')

plyr::ldply(res.all.list, function(D) D$impute.data) %>%
  pivot_wider(id_cols = c(patient_ID, Dataset, sample_type,
                          cancer_abbr, cancer_subtype, Gender, Age),
              names_from = 'protein') %>%
  mutate(ID = str_c(patient_ID, Dataset, '_', sample_type), .before = 1) %>%
  pull(patient_ID) %>% setequal(df1b$patient_ID) # TRUE
t.nt.input <- df1b
saveRDS(t.nt.input, 'output/tumor_nontumor_input_v1202.rds')




View(res_main %>% filter(p_adj_BH < 0.05, #abs(effect) >= 0.5,
                         abs(g_adj) >= 0.5) %>% 
       count(cancer_abbr))

pdf('T_NT_datamatrix_paired_cancer_split_V1202.pdf', width = 15, height = 5)
# Symmetric color scale around 0
vmax <- 3
hm_cols <- colorRampPalette(brewer.pal(11, "RdYlBu")[c(10:6, 5:1)])(101)  # red→blue diverging
hm_breaks <- seq(-vmax, vmax, length.out = 101)

for(ca in unique(impute.data$cancer_abbr)){
  cat('Drawing ', ca, '...\n')
  tmp <- impute.data %>%
    filter(cancer_abbr == ca) %>%
    column_to_rownames('ID') %>% 
    select(-(patient_ID:Age)) %>% 
    t() %>% removeRowsAllNa()
  # prot_cv05 <- apply(tmp, 1, cv) %>% .[. > quantile(apply(tmp, 1, cv), 0.5)] %>% names()
  ann_col <- impute.data %>%
    filter(cancer_abbr == ca) %>%
    column_to_rownames('ID') %>% 
    select(patient_ID:Age) %>% 
    arrange(cancer_abbr, cancer_subtype, sample_type, Age) %>% 
    select(-cancer_subtype)
  ann_colors <- list(
    # cancer_subtype = cancer_color,
    cancer_abbr = cancer_color,
    sample_type = ggsci::pal_d3()(length(unique(ann_col$sample_type))) %>% setNames(sort(unique(ann_col$sample_type))),
    Gender = ggsci::pal_aaas()(length(unique(ann_col$Gender))) %>% setNames(sort(unique(ann_col$Gender))),
    Dataset = ggsci::pal_jama()(length(unique(ann_col$Dataset))) %>% setNames(sort(unique(ann_col$Dataset)))
  )
  
  pheatmap::pheatmap(
    # tmp[prot_cv05, ],
    tmp,
    scale = 'row',
    color = hm_cols, breaks = hm_breaks,
    cluster_rows = T, cluster_cols = T,
    clustering_method = 'complete',
    show_rownames = F, show_colnames = T,
    fontsize_row = 6, fontsize_col = 10,
    border_color = NA, na_col = "#CCCCCC",
    # annotation_row = ann_row,
    # annotation_row = prm_speprot %>% column_to_rownames('protein') %>% select(cancer_abbr),
    annotation_col = ann_col %>% select(-cancer_abbr),
    annotation_colors = ann_colors,
    fontsize = 7,
    main = str_c("T-NT matrix ", ca, ' - ', str_c(dim(tmp), c('prots', 'runs'), collapse = ' ')),
  ) %>% print()
}
graphics.off()

# post-correlation
pdf('T_NT_datamatrix_paired_cancer_cor_split_V1202.pdf', width = 15, height = 5)
# Symmetric color scale around 0
vmax <- 1
hm_cols <- colorRampPalette(brewer.pal(11, "RdYlBu"))(101)  # red→blue diverging
hm_breaks <- seq(0, vmax, length.out = 101)

post.cor.list <- list()
for(ca in unique(impute.data$cancer_abbr)[1:25]){
  cat('Correlation of ', ca, '...\n')
  tmp <- impute.data %>%
    filter(cancer_abbr == ca) %>%
    column_to_rownames('ID') %>% 
    select(-(patient_ID:Age)) %>% 
    t() %>% removeRowsAllNa()
  tmp.cor <- cor(tmp, method = 'pearson', use = "pairwise.complete.obs")
  post.cor.list[[ca]] <- tmp.cor
  
  ann_col <- impute.data %>%
    filter(cancer_abbr == ca) %>%
    column_to_rownames('ID') %>% 
    select(patient_ID:Age) %>% 
    arrange(cancer_abbr, cancer_subtype, sample_type, Age) %>% 
    select(-cancer_subtype)
  ann_colors <- list(
    # cancer_subtype = cancer_color,
    cancer_abbr = cancer_color,
    sample_type = ggsci::pal_d3()(length(unique(ann_col$sample_type))) %>% setNames(sort(unique(ann_col$sample_type))),
    Gender = ggsci::pal_aaas()(length(unique(ann_col$Gender))) %>% setNames(sort(unique(ann_col$Gender))),
    Dataset = ggsci::pal_jama()(length(unique(ann_col$Dataset))) %>% setNames(sort(unique(ann_col$Dataset)))
  )
  
  pheatmap::pheatmap(
    # tmp[prot_cv05, ],
    tmp.cor,
    scale = 'none',
    color = hm_cols, #breaks = hm_breaks,
    cluster_rows = T, cluster_cols = T,
    clustering_method = 'complete',
    show_rownames = F, show_colnames = T,
    fontsize_row = 6, fontsize_col = 10,
    border_color = NA, na_col = "#CCCCCC",
    # annotation_row = ann_row,
    # annotation_row = prm_speprot %>% column_to_rownames('protein') %>% select(cancer_abbr),
    annotation_col = ann_col %>% select(-cancer_abbr),
    annotation_colors = ann_colors,
    fontsize = 7,
    main = str_c("T-NT matrix ", ca, ' - ', str_c(dim(tmp), c('prots', 'runs'), collapse = ' ')),
  ) %>% print()
}
graphics.off()



# 3.DEP numbers -----
## all DEPs -----
expr <- readRDS('output/expr.rds')
meta <- readRDS('output/T_NT_info_check.rds')
impute.data <- readRDS('output/T_NT_analysis_imputated_data_v1202.rds')
res.lm.tbls <- readRDS('output/tumor_nontumor_analysis_reports_v1202.rds')
res_main <- res.lm.tbls$Diff.report
res_main_filter <- rio::import('output/tumor_nontumor_analysis_reports_v1202.xlsx', sheet = 'Diff.report.filter')
# res_main_filter %<>% mutate(direction = if_else(effect >= 0, "Up", "Down"))

### volcano plot -------
p1 <- ggplot(res_main, aes(x = effect, y = -log10(p_adj_BH))) +
  facet_wrap(~ cancer_abbr, scale = 'free') +
  geom_point(aes(color = (p_adj_BH < 0.05) & (abs(effect) >= 0.5))) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic(base_size = 10)
p2 <- ggplot(res_main, aes(x = g_adj, y = -log10(p_adj_BH))) +
  facet_wrap(~ cancer_abbr, scale = 'free') +
  geom_point(aes(color = (p_adj_BH < 0.05) & (abs(g_adj) >= 0.5))) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic(base_size = 10)

pdf('output/T_NT_volcano_paired_v1202.pdf', width = 12, height = 10)
print(p1)
print(p2)
graphics.off()

### DEPs bar -------
dep_stat <- res_main_filter %>%
  count(cancer_abbr, direction) %>%
  mutate(direction = factor(direction, levels = c("Up", "Down"))) %>%
  complete(cancer_abbr, direction, fill = list(n = 0)) %>%
  mutate(
    signed_n = if_else(direction == "Down", -n, n),
    label_y  = if_else(n == 0, NA_real_, signed_n / 2),
    label_txt = if_else(n == 0, "", as.character(n))
  )  %>%
  inner_join(meta %>% distinct(cancer, cancer_abbr), .)

p_dys <- ggplot(dep_stat, aes(x = cancer_abbr, y = signed_n, fill = direction)) +
  geom_col(width = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(y = label_y, label = label_txt), size = 3.5, hjust = 0.5, vjust = 0.5) +
  scale_fill_manual(values = c(Up = "#AB2947", Down = "#17B544")) +
  scale_y_continuous(labels = function(x) abs(x)) +
  labs(
    y = "# tumor dysregulated proteins"
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(colour = "black"),
    axis.line.x = element_blank(),
    plot.subtitle = element_text(size = 10, hjust = 0, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, vjust = 1, hjust = 1, color = "black", angle = 45),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 10, hjust = 0.5, color = "black", angle = 90),
    axis.text.y = element_text(size = 10, color = "black", angle = 0),
    axis.ticks.y = element_blank()
  )
ggsave('output/TPHP_tumor_dysregulated_proteins_v1202.pdf', p_dys, width = 10, height = 4.6)

p_dys <- ggplot(dep_stat, aes(x = cancer, y = signed_n, fill = direction)) +
  geom_col(width = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(y = label_y, label = label_txt), size = 3.5, hjust = 0.5, vjust = 0.5) +
  scale_fill_manual(values = c(Up = "#AB2947", Down = "#17B544")) +
  scale_y_continuous(labels = function(x) abs(x)) +
  labs(
    y = "# tumor dysregulated proteins"
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(colour = "black"),
    axis.line.x = element_blank(),
    plot.subtitle = element_text(size = 10, hjust = 0, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, vjust = 1, hjust = 1, color = "black", angle = 45),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 10, hjust = 0.5, color = "black", angle = 90),
    axis.text.y = element_text(size = 10, color = "black", angle = 0),
    axis.ticks.y = element_blank()
  )
ggsave('output/TPHP_tumor_dysregulated_proteins_fullname_v1202.pdf', p_dys, width = 10, height = 4.6)


### upset -----
protls_upset <- split(res_main_filter$protein, res_main_filter$cancer_abbr)
pdf('output/TPHP_CA_ADJ_dys_overlap_upset_top5_v1202.pdf', width = 10, height = 5)
print(
  UpSetR::upset(
    UpSetR::fromList(protls_upset),
    order.by = 'freq',# order.by = c('freq', 'degree'),
    mainbar.y.label = 'Protein number',
    # queries = list(
    #   list(query = UpSetR::intersects, params = list('TE', 'UTE', 'MAG', 'CU', 'MU'), color = "orange", active = T, query.name = 'The greatest DEPs intersection of top 5 organs')),
    # query.legend = "top",
    # group.by = "sets",
    # empty.intersections = "on",
    set_size.show = T, number.angles = 0, text.scale = c(1, 1, 1, 1, 1, 1.5)
  )
)
graphics.off()

pdf('TPHP_CA_ADJ_dys_overlap_upset_25sets_95intersects_v1202.pdf', width = 21.0, height = 29.7, paper = 'a4r')
print(UpSetR::upset(UpSetR::fromList(protls_upset), order.by = 'freq', nsets = length(protls_upset), nintersects = 95, set_size.show = T))
graphics.off()


## tumor-specific proteins -----
# A protein is "tumor-specific" if it appears in exactly one cancer in the filtered table
# up/down distinct
up <- res_main_filter %>% filter(direction == 'Up')
down <- res_main_filter %>% filter(direction == 'Down')
up_counts <- up %>% count(protein, name = "n_cancers") %>% mutate(direction = 'Up')
down_counts <- down %>% count(protein, name = "n_cancers") %>% mutate(direction = 'Down')
prot_counts <- rbind(up_counts, down_counts)

# prot_counts <- res_main_filter %>% count(protein, name = "n_cancers")
cs_map <- res_main_filter %>%
  inner_join(prot_counts %>% filter(n_cancers == 1), by = c("protein", "direction"))
length(unique(cs_map$cancer_abbr)) # 25
length(unique(cs_map$protein)) # 2878
length(unique(cs_map$protein[cs_map$direction == 'Up'])) # 1547
length(unique(cs_map$protein[cs_map$direction == 'Down'])) # 1431
cs_counts <- cs_map %>% count(cancer_abbr, direction) %>%
  group_by(cancer_abbr) %>%
  arrange(cancer_abbr, direction) %>%
  ungroup()
rio::export(cs_map, 'T_NT_DEP_tumor_specific_v1202.xlsx')
split(cs_map$protein, cs_map$cancer_abbr) %>% 
  create_dataframe_from_list() %>% 
  rio::export('TPHP_CA_targetlist_3.csv')

# Order cancers by total counts (ascending). Plot Down as negative for symmetric bars.
cs_counts_plot <- cs_counts %>%
  group_by(cancer_abbr) %>% mutate(total = sum(n)) %>% ungroup() %>%
  mutate(
    cancer_abbr = forcats::fct_reorder(cancer_abbr, total, .fun = sum, .desc = FALSE),
    n_signed    = if_else(direction == "Down", -n, n)
  )
p_bar <- ggplot(cs_counts_plot, aes(x = cancer_abbr, y = n_signed, fill = direction)) +
  geom_col(width = 0.9, color = "black") +
  scale_fill_manual(values = c(Up = "#AB2947", Down = "#17B544")) +
  coord_flip() +
  labs(title = "Tumor-specific proteins", y = "# proteins", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    text = element_text(size = 12)
  ) +
  geom_text(aes(label = abs(n_signed)), hjust = if_else(cs_counts_plot$direction=="Down", 1.1, -0.1),
            size = 3)
ggsave("output/tumor_specific_proteins_bar_v1202.pdf", p_bar, width = 8, height = 5)

# 4.Tumor-specific analysis -----
## 4.1 use intensity matrix -----
## Visualization B: heatmap of top-|g_adj| tumor-specific proteins by cancer
# Choose top-K per cancer for readability (change K as needed)
top_k_per_cancer <- 3

cs_topk <- cs_map %>%
  group_by(cancer_abbr) %>%
  slice_max(order_by = abs(g_adj),
            #order_by = g_adj,
            n = top_k_per_cancer, with_ties = FALSE) %>%
  ungroup() %>% 
  left_join(dfprot %>% rename(protein = Protein.Group)) %>% 
  mutate(label = str_glue("{Genes} ({protein})"),
         specific.type = str_c('specific.', tolower(direction)))

# Build a matrix of g_adj: rows = proteins, cols = cancers (NA if not that cancer)
hm_df <- cs_topk %>%
  select(protein, cancer_abbr, g_adj) %>%
  pivot_wider(names_from = cancer_abbr, values_from = g_adj)

# Keep an informative rowname; here we use protein as-is (you can append gene if available)
hm_mat <- hm_df %>%
  column_to_rownames("protein") %>%
  as.matrix()

# Symmetric color scale around 0
vmax <- 5
hm_cols <- colorRampPalette(brewer.pal(11, "RdYlBu")[c(10:6, 5:2)])(101)  # red→blue diverging
# hm_cols <- colorRampPalette(c("#D73027", "#FFFFBF", "#4575B4"))(101)
hm_breaks <- seq(-vmax, vmax, length.out = 101)
hm_cols[abs(hm_breaks) < 0.5] <- 'white'

# Optional row annotation by direction in its cancer (Up/Down).
# Build a vector mapping protein→direction from cs_topk
# row_dir <- cs_topk %>% select(protein, direction) %>% distinct()
# ann_row <- data.frame(Direction = row_dir$direction)
# row.names(ann_row) <- row_dir$protein
for(direct in unique(sort(cs_topk$direction))){
  cs_topk_di <- cs_topk %>% filter(direction == direct)
  ann_row <- cs_topk_di %>%
    arrange(cancer_abbr, desc(direction)) %>% 
    select(label, cancer_abbr, direction) %>% 
    column_to_rownames('label')
  # ann_colors <- list(Direction = c(Up = "#AB2947", Down = "#17B544"))
  ann_col <- data.frame(FileName = meta %>%
                          filter(patient_ID %in% impute.data$patient_ID) %>%
                          pull(FileName)) %>% 
    inner_join(meta %>% select(FileName, cancer_abbr, sample_type, Gender, Age), .) %>% 
    filter(cancer_abbr != 'DLBCL') %>% 
    mutate(#cancer_subtype = cancer_abbr, .after = cancer_abbr,
      cancer_abbr = ifelse(str_detect(cancer_abbr, 'BRCA'), 'BRCA', cancer_abbr)) %>% 
    column_to_rownames('FileName') %>% 
    arrange(cancer_abbr,# cancer_subtype,
            desc(sample_type))
  ann_colors <- list(
    # cancer_subtype = cancer_color,
    # cancer_abbr = cancer_color %>% append(c('OC' = '#E58989')),
    cancer_abbr = cancer_color,
    sample_type = ggsci::pal_d3()(length(unique(ann_col$sample_type))) %>% setNames(sort(unique(ann_col$sample_type))),
    Gender = ggsci::pal_aaas()(length(unique(ann_col$Gender))) %>% setNames(sort(unique(ann_col$Gender))),
    direction = c(Up = "#D73027", Down = "#4575B4")
  )
  
  hm_mat_expr <- expr[rownames(hm_mat), rownames(ann_col)]
  rownames(hm_mat_expr) <- cs_topk_di %>% pull(label, protein) %>% .[rownames(hm_mat)]
  
  # compute gap positions by cancer_abbr runs (in the current display order)
  # rows
  rle_row <- rle(as.character(ann_row$cancer_abbr))
  gaps_row_idx <- cumsum(rle_row$lengths)
  gaps_row_idx <- gaps_row_idx[-length(gaps_row_idx)]  # remove the last break
  
  # cols
  rle_col <- rle(as.character(ann_col$cancer_abbr))
  gaps_col_idx <- cumsum(rle_col$lengths)
  gaps_col_idx <- gaps_col_idx[-length(gaps_col_idx)]
  
  pheatmap::pheatmap(
    hm_mat_expr[rownames(ann_row), ], scale = 'row',
    color = hm_cols, breaks = hm_breaks,
    cluster_rows = F, cluster_cols = F,
    show_rownames = T, show_colnames = F,
    fontsize_row = 6, fontsize_col = 10,
    border_color = NA, na_col = "#CCCCCC",
    annotation_row = ann_row,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    gaps_row = gaps_row_idx,
    gaps_col = gaps_col_idx,
    fontsize = 7,
    main = paste0("Top-", top_k_per_cancer, " tumor-specific proteins per cancer (g_adj)"),
    filename = str_c("tumor_specific_proteins_heatmap_v1202", direct, ".pdf"), width = 10, height = 4
  )
}
ann_row <- cs_topk %>%
  arrange(cancer_abbr, desc(specific.type)) %>% 
  select(label, cancer_abbr, specific.type) %>% 
  pivot_wider(id_cols = label, names_from = specific.type, values_from = cancer_abbr) %>% 
  column_to_rownames('label')
# ann_colors <- list(Direction = c(Up = "#AB2947", Down = "#17B544"))
ann_col <- data.frame(FileName = meta %>%
                        filter(patient_ID %in% impute.data$patient_ID) %>%
                        pull(FileName)) %>% 
  inner_join(meta %>% select(FileName, cancer_abbr, sample_type, Gender, Age), .) %>% 
  filter(cancer_abbr != 'DLBCL') %>% 
  mutate(#cancer_subtype = cancer_abbr, .after = cancer_abbr,
    cancer_abbr = ifelse(str_detect(cancer_abbr, 'BRCA'), 'BRCA', cancer_abbr)) %>% 
  column_to_rownames('FileName') %>% 
  arrange(cancer_abbr,# cancer_subtype,
          desc(sample_type))
ann_colors <- list(
  # cancer_subtype = cancer_color,
  cancer_abbr = cancer_color,
  specific.up = cancer_color,
  specific.down = cancer_color,
  sample_type = ggsci::pal_d3()(length(unique(ann_col$sample_type))) %>% setNames(sort(unique(ann_col$sample_type))),
  Gender = ggsci::pal_aaas()(length(unique(ann_col$Gender))) %>% setNames(sort(unique(ann_col$Gender)))#,
  # direction = c(Up = "#D73027", Down = "#4575B4")
)

hm_mat_expr <- expr[rownames(hm_mat), rownames(ann_col)]
rownames(hm_mat_expr) <- cs_topk %>% pull(label, protein) %>% .[rownames(hm_mat)]

# compute gap positions by cancer_abbr runs (in the current display order)
# rows
rle_row <- rle(as.character(ann_row$cancer_abbr))
gaps_row_idx <- cumsum(rle_row$lengths)
gaps_row_idx <- gaps_row_idx[-length(gaps_row_idx)]  # remove the last break

# cols
rle_col <- rle(as.character(ann_col$cancer_abbr))
gaps_col_idx <- cumsum(rle_col$lengths)
gaps_col_idx <- gaps_col_idx[-length(gaps_col_idx)]

pheatmap::pheatmap(
  hm_mat_expr[rownames(ann_row), ], scale = 'row',
  color = hm_cols, breaks = hm_breaks,
  cluster_rows = F, cluster_cols = F,
  show_rownames = T, show_colnames = F,
  fontsize_row = 6, fontsize_col = 10,
  border_color = NA, na_col = "#CCCCCC",
  annotation_row = ann_row,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  gaps_row = gaps_row_idx,
  gaps_col = gaps_col_idx,
  fontsize = 7,
  main = paste0("Top-", top_k_per_cancer, " tumor-specific proteins per tumor (g_adj)"),
  filename = "output/tumor_specific_proteins_heatmap_v1202.pdf", width = 10, height = 4
)


## 4.2 use logFC matrix -----
data_spec_ana <- res_main %>%
  inner_join(dfprot %>% rename(protein = Protein.Group)) %>% 
  filter(protein %in% cs_topk$protein) %>% 
  mutate(label = str_glue("{Genes} ({protein})"))

hm_mat_adjg <- data_spec_ana %>%
  arrange(cancer_abbr) %>% 
  pivot_wider(id_cols = label, names_from = cancer_abbr, values_from = g_adj) %>%
  column_to_rownames('label') %>% 
  .[rownames(ann_row), ]
hm_mat_adjp <- data_spec_ana %>%
  arrange(cancer_abbr) %>% 
  mutate(p_adj_BH = ifelse(p_adj_BH < 0.001, '***',
                           ifelse(p_adj_BH < 0.01, '**',
                                  ifelse(p_adj_BH < 0.05, '*', '')))) %>% 
  pivot_wider(id_cols = label, names_from = cancer_abbr, values_from = p_adj_BH,
              values_fill = '') %>%
  column_to_rownames('label') %>% 
  .[rownames(ann_row), ]

pheatmap::pheatmap(
  t(hm_mat_adjg), scale = 'none',
  color = hm_cols, breaks = hm_breaks,
  cluster_rows = F, cluster_cols = F,
  show_rownames = T, show_colnames = T,
  font_size = 8, angle_col = 45,
  display_numbers = t(hm_mat_adjp), fontsize_number = 6,
  border_color = NA, na_col = "#CCCCCC",
  annotation_col = ann_row,
  annotation_row = data.frame(cancer_abbr = colnames(hm_mat_adjg),
                              row.names = colnames(hm_mat_adjg)),
  annotation_colors = ann_colors,
  gaps_row = seq_along(colnames(hm_mat_adjg)),
  gaps_col = gaps_row_idx,
  main = paste0("Top-", top_k_per_cancer, " tumor-specific proteins per tumor (g_adj)"),
  filename = "output/tumor_specific_proteins_adjg_heatmap_v1202.pdf", width = 15, height = 5
)


# data0[, 'Q12816', drop=F] %>% rownames_to_column('FileName') %>% inner_join(meta, .) %>% rio::export('check_Q12816_raw.xlsx')

# 5.PRM validated -----
prm.data.raw <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/PRM/20230313TPHP_PRM_prot_matrix_48.csv')
prm.matrix.filterna <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/PRM/20230318_PUH_PRM_submatrix_combine_5_missing.csv')
df_diff_all.prm <- readxl::read_excel('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/PRM/20230318_TPHP_PRM_dysregulated_filter50NAByOrgan.xlsx')
# prm_prot <- c('O15079', 'O15194', 'O95858', 'O95990', 'P17947', 'Q8TCG5', 'Q92521', 'Q9H1C7', 'Q9HCJ6', 'Q9Y692', 'Q9UH73', 'O75554', 'O15204', 'P56749', 'Q6UXB2', 'Q9UIK4', 'O14757', 'Q02548', 'Q15910', 'Q7L513', 'Q96RG2', 'Q9H967', 'Q9NS87', 'Q9H0Q3', 'P07202', 'Q969H4')
prm_dysprot <- df_diff_all.prm %>%
  filter(pAdj_t < 0.05, abs(log2FC) > log2(2)) %>%
  rename(anatomical_classification = organ) %>%
  inner_join(info.all %>%
               filter(sample_type %in% c('T', 'NT')) %>%
               mutate(cancer_abbr = ifelse(str_detect(cancer_abbr, 'BRCA'), 'BRCA', cancer_abbr),
                      cancer_abbr = ifelse(str_detect(cancer_abbr, '(HGS|CC)OC'), 'OC', cancer_abbr)) %>%
               distinct(anatomical_classification, cancer_abbr) %>% drop_na())

up.prm <- prm_dysprot %>%
  filter(log2FC > 0) %>%
  filter(protein %in% res_main_filter$protein) %>% 
  # inner_join(res_main_filter) %>%
  inner_join(dfprot %>% rename(protein = Protein.Group), .)
up.prm %<>% count(protein) %>% filter(n == 1) %>% semi_join(up.prm, .)
down.prm <- prm_dysprot %>%
  filter(log2FC < 0) %>%
  filter(protein %in% res_main_filter$protein) %>% 
  # inner_join(res_main_filter) %>%
  inner_join(dfprot %>% rename(protein = Protein.Group), .)
down.prm %<>% count(protein) %>% filter(n == 1) %>% semi_join(down.prm, .)
prm_speprot <- rbind(up.prm, down.prm)
hm_mat_expr <- expr[unique(prm_speprot$protein), ]


pheatmap::pheatmap(
  hm_mat_expr, scale = 'row',
  color = hm_cols, breaks = hm_breaks,
  cluster_rows = F, cluster_cols = F,
  show_rownames = T, show_colnames = F,
  fontsize_row = 6, fontsize_col = 10,
  border_color = NA, na_col = "#CCCCCC",
  # annotation_row = ann_row,
  # annotation_row = prm_speprot %>% column_to_rownames('protein') %>% select(cancer_abbr),
  # annotation_col = ann_col,
  annotation_colors = ann_colors,
  fontsize = 7,
  main ="PRM validated proteins",
  filename = "output/PRM_validated_proteins_heatmap_v1202.pdf", width = 10, height = 3
)


## V2 -----
data <- impute.data %>% column_to_rownames('ID') %>% 
  select(-(patient_ID:Age)) %>% t()
prm_speprot_wide <- prm_speprot %>% 
  mutate(specific.type = case_when(log2FC > 0 ~ 'specific.up',
                                   log2FC < 0 ~ 'specific.down')) %>% 
  arrange(cancer_abbr, desc(specific.type)) %>% 
  select(protein, cancer_abbr, specific.type) %>% 
  pivot_wider(id_cols = protein, names_from = specific.type, values_from = cancer_abbr,
              values_fn = function(x) 1, values_fill = 0)
ann_col <- impute.data %>% column_to_rownames('ID') %>% 
  select(patient_ID:Age) %>% 
  arrange(cancer_abbr, cancer_subtype, sample_type, Age) %>% 
  select(-cancer_subtype)
ann_colors <- list(
  # cancer_subtype = cancer_color,
  cancer_abbr = cancer_color,
  sample_type = ggsci::pal_d3()(length(unique(ann_col$sample_type))) %>% setNames(sort(unique(ann_col$sample_type))),
  Gender = ggsci::pal_aaas()(length(unique(ann_col$Gender))) %>% setNames(sort(unique(ann_col$Gender))),
  Dataset = ggsci::pal_jama()(length(unique(ann_col$Dataset))) %>% setNames(sort(unique(ann_col$Dataset)))
)

hm_mat_expr <- data[intersect(prm_speprot$protein, rownames(data)), rownames(ann_col)]

# Symmetric color scale around 0
vmax <- 3
hm_cols <- colorRampPalette(brewer.pal(11, "RdYlBu")[c(10:6, 5:1)])(101)  # red→blue diverging
hm_breaks <- seq(-vmax, vmax, length.out = 101)

pheatmap::pheatmap(
  hm_mat_expr, scale = 'row',
  color = hm_cols, breaks = hm_breaks,
  cluster_rows = F, cluster_cols = F,
  show_rownames = T, show_colnames = F,
  fontsize_row = 6, fontsize_col = 10,
  border_color = NA, na_col = "#CCCCCC",
  # annotation_row = ann_row,
  annotation_row = prm_speprot_wide %>% column_to_rownames('protein'),
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  fontsize = 7,
  main ="PRM validated proteins",
  filename = "PRM_validated_proteins_heatmap_v1202_V2.pdf", width = 10, height = 3
)
graphics.off()

bp_mat_expr <- impute.data %>% 
  select(ID:Age, all_of(intersect(prm_speprot$protein, rownames(data)))) %>% 
  pivot_longer(cols = intersect(prm_speprot$protein, rownames(data)),
               names_to = 'protein', values_to = 'Log2Intensity', values_drop_na = T) %>% 
  inner_join(
    prm_speprot %>%
      distinct(cancer_abbr, protein) %>% 
      rename(cancer.valid = cancer_abbr) %>%
      group_by(protein) %>%
      reframe(cancer.valid = str_c(cancer.valid, collapse = ';')), .
  ) %>% 
  inner_join(dfprot %>% rename(protein = Protein.Group), .) %>% 
  mutate(sample_type = factor(sample_type, c('T', 'NT')))

# ggplot(bp_mat_expr) +
#   facet_wrap(~ str_c(protein, '-', cancer.valid), ncol = 1) +
#   aes(x = cancer_abbr, y = Log2Intensity, color = sample_type) +
#   geom_boxplot() +
#   geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
#   theme_classic()

# Build p-values per (protein, cancer_abbr, cancer.valid), comparing sample_type
pvals <- bp_mat_expr %>%
  group_by(Genes, protein, cancer_abbr, cancer.valid) %>%
  rstatix::wilcox_test(Log2Intensity ~ sample_type, paired = T) %>%
  rstatix::add_xy_position(x = "cancer_abbr", dodge = 0.75) %>%
  mutate(
    p.label = sprintf('%.2e', p),
    # p.label = rstatix::p_format(p, digits = 3),
    panel = str_c(protein, "-", cancer.valid)  # match your facet
  )

# p_box <- ggplot(bp_mat_expr, aes(x = cancer_abbr, y = Log2Intensity, color = sample_type)) +
#   facet_wrap(~ str_c(Genes, ".", protein, "-", cancer.valid), ncol = 1, scales = 'free_y') +
#   geom_boxplot(position = position_dodge(width = 0.75)) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
#   labs(x = 'Cancer', y = 'Log2Intensity') +
#   stat_pvalue_manual(
#     data = pvals,
#     label = "p.label",
#     y.position = "y.position", xmin = "xmin", xmax = "xmax",
#     hide.ns = F, tip.length = 0,
#     inherit.aes = FALSE        # <-- key to avoid length mismatch
#   ) +
#   theme_classic() +
#   theme(text = element_text(size = 10),
#         strip.background = element_blank())
p_box <- ggplot(bp_mat_expr, aes(x = cancer_abbr, y = Log2Intensity, color = sample_type)) +
  facet_wrap(~ str_c(cancer.valid, "-", Genes, ".", protein),
             ncol = 1, scales = "free_y", strip.position = "top") +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  labs(x = "Cancer", y = "Log2Intensity") +
  stat_pvalue_manual(
    data = pvals, label = "p.label", size = 3,
    y.position = "y.position", xmin = "xmin", xmax = "xmax",
    hide.ns = FALSE, tip.length = 0,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    strip.text.y.right = element_text(size = 8, angle = 0, hjust = 0, margin = margin(l = 2)),
    # strip.background = element_blank(),
  )
ggsave('output/PRM_validated_proteins_box_v1202.pdf', p_box, width = 20, height = 1*9)




impute.data %>% 
  select(ID:Age, O00481, Q53RD9) %>% 
  # filter(cancer_abbr == 'CCOC') %>% 
  pivot_longer(cols = c(O00481, Q53RD9),
               names_to = 'protein', values_to = 'Log2Intensity', values_drop_na = T) %>% 
  inner_join(dfprot %>% rename(protein = Protein.Group), .) %>% 
  mutate(sample_type = factor(sample_type, c('T', 'NT'))) %>% 
  ggplot() +
  facet_wrap(~ str_c(Genes, ".", protein),
             ncol = 1, scales = "free_y", strip.position = "top") +
  aes(x = cancer_abbr, y = Log2Intensity, color = sample_type) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
  theme_classic()


impute.data %>% 
  filter(cancer_abbr == 'CCOC', sample_type == 'NT') %>% 
  pull(O00481) -> x1
impute.data %>% 
  filter(cancer_abbr == 'CCOC', sample_type == 'NT') %>% 
  pull(Q53RD9) -> x2
plot(impute.data %>% 
       filter(cancer_abbr == 'CCOC', sample_type == 'NT') %>% 
       pivot_longer(cols = -(ID:Age)) %>% pull(value) %>% density(na.rm = T))
abline(v = x1, col = 'red1')
abline(v = x2, col = 'blue3')


# 6.New analysis --------
## Small helpers -------------
.scale01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

# pick top-N by two keys without overfitting number of params
.pick_top <- function(df, n = 200) {
  df %>% arrange(desc(sharing_n), desc(abs_g)) %>% slice_head(n = n)
}

## 6.1 Load data ------
# meta <- rio::import('//172.16.13.136/tphp/code/6_3paired_tumor_nontumor/output/T_NT_info_check.xlsx')
meta <- rio::import('output/T_NT_info_check.xlsx')
t.nt.input <- readRDS('output/tumor_nontumor_input_v1202.rds')
res.lm.tbls <- readRDS('output/tumor_nontumor_analysis_reports_v1202.rds')
impute.data <- readRDS('output/T_NT_analysis_imputated_data_v1202.rds')

res_main_filter <- res.lm.tbls$Diff.report.filter
filtered12 <- res_main_filter %>% select(cancer_abbr, protein, Genes, g_adj, p_adj_BH, direction)


## 6.2 Pie-arc plot -----
# wide matrix (proteins × cancers)
mat_wide <- filtered12 %>%
  dplyr::select(protein, cancer_abbr, g_adj) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(names_from = cancer_abbr, values_from = g_adj) %>%
  tibble::column_to_rownames("protein") %>%
  as.matrix()

# hygiene: drop rows/cols with <2 finite values
rw <- rowSums(is.finite(mat_wide)) >= 2
cw <- colSums(is.finite(mat_wide[rw, , drop = FALSE])) >= 2
M  <- mat_wide[rw, cw, drop = FALSE]

# Euclidean (NA→0) so Ward.D2 is valid
Me <- M; Me[!is.finite(Me)] <- 0
D_eu <- stats::dist(t(Me), method = "euclidean")
HC_eu_ward <- stats::hclust(D_eu, method = "ward.D2")
HC_eu_avg  <- stats::hclust(D_eu, method = "average")
HC_eu_comp <- stats::hclust(D_eu, method = "complete")

# Cosine via normalized columns (NA→0), then D = 1 - cos
Mc <- M; Mc[!is.finite(Mc)] <- 0
nrm <- sqrt(colSums(Mc^2)); keep_cos <- nrm > 0
Mc  <- sweep(Mc[, keep_cos, drop = FALSE], 2, nrm[keep_cos], "/")
S_cos <- crossprod(Mc); S_cos[S_cos > 1] <- 1; S_cos[S_cos < -1] <- -1
D_cos <- stats::as.dist(1 - S_cos)
HC_cos_avg  <- stats::hclust(D_cos, method = "average")
HC_cos_comp <- stats::hclust(D_cos, method = "complete")

# Pearson(pairwise complete); NA pairs -> rho=0 -> D=1
S_sp <- suppressWarnings(stats::cor(M, method = "pearson", use = "pairwise.complete.obs"))
S_sp[!is.finite(S_sp)] <- 0
D_sp <- stats::as.dist(1 - S_sp)
HC_sp_avg  <- stats::hclust(D_sp, method = "average")
HC_sp_comp <- stats::hclust(D_sp, method = "complete")

# Pearson + Ward.D2 (valid): Ward on z-scored rank vectors
rank_norm <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  if (all(is.na(r))) return(r)
  r[is.na(r)] <- median(r, na.rm = TRUE)   # impute missing ranks by median
  r <- r - mean(r)                         # center
  s <- sqrt(sum(r^2)); if (s > 0) r/s else r
}
Mz <- apply(M, 2, rank_norm)
D_spW <- stats::dist(t(Mz), method = "euclidean")
HC_sp_ward <- stats::hclust(D_spW, method = "ward.D2")

# choose one dendrogram for downstream ordering
hc_use <- HC_sp_comp  # <-- you can change to HC_sp_ward / HC_eu_ward / etc.
pdf('HC_sp_comp_v1202.pdf', height = 5, width = 4)
op <- par(mar = c(2, 4, 2, 8))
plot(as.dendrogram(HC_sp_comp), main = "Clustered by g_adj", xlab = "1-r", sub = "", horiz = T)
par(op)
graphics.off()

pdf('HC_sp_ward2_v1202.pdf', height = 5, width = 4)
op <- par(mar = c(2, 4, 2, 8))
plot(as.dendrogram(HC_sp_ward), main = "Clustered by g_adj", xlab = "1-r", sub = "", horiz = T)
par(op)
graphics.off()


### Pairwise shared proteins + arc diagram (pie nodes: Up/Down) ----
prot_sets <- filtered12 %>%
  dplyr::distinct(cancer_abbr, protein) %>%
  dplyr::group_by(cancer_abbr) %>%
  dplyr::summarise(proteins = list(unique(protein)), .groups = "drop")

C <- prot_sets$cancer_abbr
pairs <- if (length(C) >= 2) t(utils::combn(C, 2)) else matrix(character(0), ncol = 2)

shared_tbl <- purrr::map_dfr(seq_len(nrow(pairs)), function(i) {
  a <- pairs[i, 1]; b <- pairs[i, 2]
  pa <- prot_sets$proteins[[match(a, C)]]
  pb <- prot_sets$proteins[[match(b, C)]]
  inter <- base::intersect(pa, pb)
  tibble::tibble(from = a, to = b, n_shared = length(inter), proteins = list(inter))
}) %>% dplyr::filter(n_shared > 0)

# ggraph-style
# Node order & Up/Down counts
edgelist <- as.matrix(shared_tbl %>% dplyr::select(from, to))
edgelist <- apply(edgelist, 2, as.character)
nodes_all <- sort(unique(c(edgelist)))
ord <- hc_use$labels[hc_use$order]
nodes_order <- c(ord[ord %in% nodes_all], setdiff(nodes_all, ord))

updown <- filtered12 %>%
  dplyr::filter(cancer_abbr %in% nodes_order) %>%
  dplyr::count(cancer_abbr, direction, name = "n") %>%
  tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  dplyr::rename(Up = `Up`, Down = `Down`) %>% 
  arrange(Up + Down)
nodes_order <- updown$cancer_abbr

updown <- dplyr::full_join(tibble::tibble(cancer_abbr = nodes_order), updown, by = "cancer_abbr") %>%
  dplyr::mutate(dplyr::across(c(Up, Down), ~ tidyr::replace_na(.x, 0L)))

nodes_df <- tibble::tibble(name = nodes_order, cancer = nodes_order)
edges_df <- shared_tbl %>%
  dplyr::mutate(from = match(from, nodes_order),
                to   = match(to,   nodes_order),
                weight = n_shared)

g <- tidygraph::tbl_graph(nodes = nodes_df, edges = edges_df, directed = FALSE)

edge_protein_tbl <- filtered12 %>%
  filter(cancer_abbr %in% nodes_order) %>%
  distinct(cancer_abbr, protein, direction) %>%
  group_by(direction, protein) %>%
  summarise({
    cs <- sort(unique(cancer_abbr))
    if (length(cs) >= 2) {
      pm <- utils::combn(cs, 2)
      tibble::tibble(from = pm[1,], to = pm[2,])
    } else tibble::tibble(from=character(0), to=character(0))
  }, .groups = "drop")

edges_df2 <- edge_protein_tbl %>%
  mutate(from = match(from, nodes_order),
         to   = match(to,   nodes_order)) %>%    # protein, direction retained
  dplyr::count(from, to, direction, name = "weight")

# # sanity check: should see values >1 somewhere
# edges_df2 %>% count(from, to, direction, name = "n_proteins") %>%
#   arrange(desc(n_proteins)) %>% print(n = 20)

g2 <- tidygraph::tbl_graph(nodes = nodes_df, edges = edges_df2, directed = FALSE)
g  <- g2  # override for downstream layout/plot

# Layout & pie data
lay <- ggraph::create_layout(g, layout = "linear")

# Enlarge horizontal spacing
SPACING <- 2.0
lay$x <- lay$x * SPACING

lay_df <- as.data.frame(lay)
if ("label" %in% names(lay_df)) {
  lay_df <- dplyr::rename(lay_df, name = label)
}

pie_df <- lay_df %>%
  dplyr::select(name, x, y) %>%
  dplyr::left_join(updown, by = c("name" = "cancer_abbr"))

# ---- Adaptive radius based on min neighbor spacing (pre-scaling)
dx <- diff(sort(unique(lay_df$x)))
min_dx <- if (length(dx)) min(dx) else 1
RADIUS <- 0.9 * (min_dx / SPACING)
pie_df$pie_r <- RADIUS

# --- Central labels (at exact pie center)
# Show Up and Down stacked at the center; no reliance on slice angles.
center_label_df <- pie_df %>%
  dplyr::mutate(
    label = dplyr::case_when(
      Up > 0 & Down > 0 ~ paste0(Up, "\n", Down),
      Up > 0            ~ as.character(Up),
      Down > 0          ~ as.character(Down),
      TRUE              ~ "0"
    )
  )

# arcdiagram
### Pies only ----
pie_df2 <- pie_df %>%
  rename(x2 = y, y2 = x) %>% 
  rename(x = x2, y = y2) %>% 
  mutate(size = Up + Down,
         size.rel = (size / max(size)) ^ 0.2,
         .RATIO_Up = Up / size,
         .RATIO_Down = Down / size)
center_label_df2 <- center_label_df %>%
  rename(x2 = y, y2 = x) %>% 
  rename(x = x2, y = y2) %>% 
  mutate(size = Up + Down) %>% 
  inner_join(meta %>% distinct(cancer, cancer_abbr) %>% mutate(name = cancer_abbr)) %>% 
  mutate(x = x + 6,
         # label = str_glue("{cancer_abbr} - {cancer}\nn={size} ({Up} up, {Down} down)"),
         label = str_glue("{cancer_abbr}\n({Up} up, {Down} down)"))

rho <- 0.6
phi0 <- pi/2
s    <- -1     # clockwise like a clock

ratio_lab_df2 <- pie_df2 %>%
  transmute(name, x, y, r = size.rel, Up, Down) %>%
  pivot_longer(c(Up, Down), names_to = "dir", values_to = "value") %>%
  group_by(name) %>%
  mutate(
    dir   = factor(dir, levels = c("Up", "Down")),
    value = as.numeric(value),
    total = sum(value),
    prop  = value / total,
    end   = cumsum(prop),
    start = end - prop,
    mid   = (start + end) / 2, 
    theta = phi0 + s * 2 * pi * mid,
    x_lab = x + rho * r * cos(theta),
    y_lab = y + rho * r * sin(theta),
    label = scales::percent(prop, accuracy = 1)
  ) %>%
  ungroup()

pie.chart <- ggplot() +
  scatterpie::geom_scatterpie(
    data = pie_df2,
    aes(x = x, y = y, r = size.rel),
    cols = c("Up", "Down"),
    color = NA,
    alpha = 0.95
  ) +
  # geom_text(
  #   data = ratio_lab_df2,
  #   aes(x = x_lab, y = y_lab, label = label),
  #   size = 3, fontface = "bold", show.legend = FALSE
  # ) +
  geom_text(
    data = center_label_df2,
    aes(x = x, y = y, label = label),
    colour = "black", fontface = "bold", size = 5, lineheight = 0.95,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = sample_color[1:2] %>% setNames(c('Up', 'Down')),
                    name = "Dysregulated proteins") +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  coord_fixed(clip = "off") +
  theme(plot.margin = margin(10, 5, 10, 5), legend.position = "top")
ggsave('output/T_NT_DEP_pie_v1202.pdf', pie.chart, width = 7, height = 15)




## 6.2B Pie-arc plot tumor-specific -----
cs_map <- res_main_filter %>%
  inner_join(prot_counts %>% filter(n_cancers == 1), by = "protein") %>% 
  arrange(cancer_abbr)

# wide matrix (proteins × cancers)
mat_wide <- filtered12 %>%
  dplyr::filter(protein %in% cs_map$protein) %>% 
  dplyr::select(protein, cancer_abbr, g_adj) %>%
  dplyr::distinct() %>%
  # dplyr::arrange(cancer_abbr) %>% 
  tidyr::pivot_wider(names_from = cancer_abbr, values_from = g_adj) %>%
  tibble::column_to_rownames("protein") %>%
  as.matrix()

### Pairwise shared proteins + arc diagram (pie nodes: Up/Down) ----
prot_sets <- filtered12 %>%
  dplyr::filter(protein %in% cs_map$protein) %>% 
  dplyr::distinct(cancer_abbr, protein) %>%
  dplyr::group_by(cancer_abbr) %>%
  dplyr::summarise(proteins = list(unique(protein)), .groups = "drop")

C <- prot_sets$cancer_abbr
pairs <- if (length(C) >= 2) t(utils::combn(C, 2)) else matrix(character(0), ncol = 2)

shared_tbl <- purrr::map_dfr(seq_len(nrow(pairs)), function(i) {
  a <- pairs[i, 1]; b <- pairs[i, 2]
  pa <- prot_sets$proteins[[match(a, C)]]
  pb <- prot_sets$proteins[[match(b, C)]]
  inter <- base::intersect(pa, pb)
  tibble::tibble(from = a, to = b, n_shared = length(inter), proteins = list(inter))
}) %>% dplyr::filter(n_shared > 0)

if (nrow(shared_tbl) > 0) {
  # ggraph-style
  # Node order & Up/Down counts
  edgelist <- as.matrix(shared_tbl %>% dplyr::select(from, to))
  edgelist <- apply(edgelist, 2, as.character)
  nodes_all <- sort(unique(c(edgelist)))
  nodes_order <- cs_map %>% count(cancer_abbr) %>% arrange(desc(n)) %>% pull(cancer_abbr)
  
  updown <- filtered12 %>%
    dplyr::filter(protein %in% cs_map$protein) %>% 
    dplyr::filter(cancer_abbr %in% nodes_order) %>%
    dplyr::count(cancer_abbr, direction, name = "n") %>%
    tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
    dplyr::rename(Up = `Up`, Down = `Down`)
  
  updown <- dplyr::full_join(tibble::tibble(cancer_abbr = nodes_order), updown, by = "cancer_abbr") %>%
    dplyr::mutate(dplyr::across(c(Up, Down), ~ tidyr::replace_na(.x, 0L)))
  
  nodes_df <- tibble::tibble(name = nodes_order, cancer = nodes_order)
  edges_df <- shared_tbl %>%
    dplyr::mutate(from = match(from, nodes_order),
                  to   = match(to,   nodes_order),
                  weight = n_shared)
  
  g <- tidygraph::tbl_graph(nodes = nodes_df, edges = edges_df, directed = FALSE)
  
  edge_protein_tbl <- filtered12 %>%
    filter(cancer_abbr %in% nodes_order) %>%
    distinct(cancer_abbr, protein, direction) %>%
    group_by(direction, protein) %>%
    summarise({
      cs <- sort(unique(cancer_abbr))
      if (length(cs) >= 2) {
        pm <- utils::combn(cs, 2)
        tibble::tibble(from = pm[1,], to = pm[2,])
      } else tibble::tibble(from=character(0), to=character(0))
    }, .groups = "drop")
  
  edges_df2 <- edge_protein_tbl %>%
    mutate(from = match(from, nodes_order),
           to   = match(to,   nodes_order)) %>%    # protein, direction retained
    dplyr::count(from, to, direction, name = "weight")
  
  # # sanity check: should see values >1 somewhere
  # edges_df2 %>% count(from, to, direction, name = "n_proteins") %>%
  #   arrange(desc(n_proteins)) %>% print(n = 20)
  
  
  # g2 <- tidygraph::tbl_graph(nodes = nodes_df, edges = edges_df2, directed = FALSE)
  # g  <- g2  # override for downstream layout/plot
  
  # Layout & pie data
  lay <- ggraph::create_layout(g, layout = "linear")
  
  # Enlarge horizontal spacing
  SPACING <- 2.0
  lay$x <- lay$x * SPACING
  
  lay_df <- as.data.frame(lay)
  if ("label" %in% names(lay_df)) {
    lay_df <- dplyr::rename(lay_df, name = label)
  }
  
  pie_df <- lay_df %>%
    dplyr::select(name, x, y) %>%
    dplyr::left_join(updown, by = c("name" = "cancer_abbr"))
  
  # ---- Adaptive radius based on min neighbor spacing (pre-scaling)
  dx <- diff(sort(unique(lay_df$x)))
  min_dx <- if (length(dx)) min(dx) else 1
  RADIUS <- 0.9 * (min_dx / SPACING)
  pie_df$pie_r <- RADIUS
  
  # --- Central labels (at exact pie center)
  # Show Up and Down stacked at the center; no reliance on slice angles.
  center_label_df <- pie_df %>%
    dplyr::mutate(
      label = dplyr::case_when(
        Up > 0 & Down > 0 ~ paste0(Up, "\n", Down),
        Up > 0            ~ as.character(Up),
        Down > 0          ~ as.character(Down),
        TRUE              ~ "0"
      )
    )
}

### Pies only ----
pie_df2 <- pie_df %>%
  # rename(x2 = y, y2 = x) %>% 
  # rename(x = x2, y = y2) %>% 
  mutate(size = Up + Down,
         size.rel = (size / max(size)) ^ 0.2,
         .RATIO_Up = Up / size,
         .RATIO_Down = Down / size)
center_label_df2 <- center_label_df %>%
  # rename(x2 = y, y2 = x) %>% 
  # rename(x = x2, y = y2) %>% 
  mutate(size = Up + Down) %>% 
  inner_join(meta %>% distinct(cancer, cancer_abbr) %>% mutate(name = cancer_abbr)) %>% 
  mutate(y = y - 2,
         label = str_glue("{cancer_abbr}\n({Up} up, {Down} down)"))

rho <- 0.6
phi0 <- pi/2
s    <- -1     # clockwise like a clock

ratio_lab_df2 <- pie_df2 %>%
  transmute(name, x, y, r = size.rel, Up, Down) %>%
  pivot_longer(c(Up, Down), names_to = "dir", values_to = "value") %>%
  group_by(name) %>%
  mutate(
    dir   = factor(dir, levels = c("Up", "Down")),
    value = as.numeric(value),
    total = sum(value),
    prop  = value / total,
    end   = cumsum(prop),
    start = end - prop,
    mid   = (start + end) / 2, 
    theta = phi0 + s * 2 * pi * mid,
    x_lab = x + rho * r * cos(theta),
    y_lab = y + rho * r * sin(theta),
    label = scales::percent(prop, accuracy = 1)
  ) %>%
  ungroup()

pie.chart <- ggplot() +
  scatterpie::geom_scatterpie(
    data = pie_df2,
    aes(x = x, y = y, r = size.rel),
    cols = c("Up", "Down"),
    color = NA,
    alpha = 0.95
  ) +
  # geom_text(
  #   data = ratio_lab_df2,
  #   aes(x = x_lab, y = y_lab, label = label),
  #   size = 3, fontface = "bold", show.legend = FALSE
  # ) +
  geom_text(
    data = center_label_df2,
    aes(x = x, y = y, label = label),
    colour = "black", fontface = "bold", size = 2, lineheight = 0.95,
    show.legend = FALSE
  ) +
  # scale_fill_manual(values = c(Up = "#D73027", Down = "#4575B4"),
  #                   name = "Cancer-unique DEPs") +
  scale_fill_manual(values = sample_color[1:2] %>% setNames(c('Up', 'Down')),
                    name = "Tumor-specific DEPs") +
  theme_void() +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  coord_fixed(clip = "off") +
  theme(plot.margin = margin(2, 5, 2, 5), legend.position = "top")
ggsave('output/T_NT_DEP_specific_pie_v1202.pdf', pie.chart, width = 20, height = 7)


## 6.3 Sankey (skipped) -----
### Sharing-number groups → g:Profiler → Sankey =====
# prot_sharing.up <- filtered12 %>%
#   dplyr::filter(direction == 'Up') %>% 
#   dplyr::distinct(protein, Genes, direction, cancer_abbr) %>%
#   dplyr::count(protein, Genes, direction, name = "sharing_n")
# prot_sharing.down <- filtered12 %>%
#   dplyr::filter(direction == 'Down') %>% 
#   dplyr::distinct(protein, Genes, direction, cancer_abbr) %>%
#   dplyr::count(protein, Genes, direction, name = "sharing_n")
# prot_sharing <- rbind(prot_sharing.up, prot_sharing.down)

# map_uniprot_to_symbol <- function(up) {
#   m <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(up), keytype = "UNIPROT", columns = c("SYMBOL"))
#   m <- dplyr::filter(m, !is.na(SYMBOL)) %>% dplyr::distinct(UNIPROT, SYMBOL)
#   tibble::tibble(protein = m$UNIPROT, SYMBOL = m$SYMBOL)
# }
# map_up2sym <- map_uniprot_to_symbol(filtered12$protein)

prot_meta <- filtered12 %>%
  dplyr::group_by(protein, Genes, direction) %>%
  dplyr::summarise(g_adj_med = median(g_adj, na.rm = TRUE),
                   p_adj_med = median(p_adj_BH, na.rm = TRUE),
                   sharing_n = dplyr::n(), .groups = "drop") %>% 
  arrange(sharing_n, direction, p_adj_med)
rio::export(prot_meta, 'T_NT_shared_DEP_1_22.xlsx')
# prot_meta <- rio::import('T_NT_shared_DEP_1_22.xlsx')

prot_meta %>% count(sharing_n)
prot_meta %>% filter(sharing_n >= 20)
prot_meta %>% filter(sharing_n >= 20) %>% pull(Genes) %>% unique() %>% writeClipboard()

# Use metascape instead of local database
res.metascape <- rio::import('metascape/v4_20_of_25_cancer_DEP/metascape_result.xlsx', sheet = 2)
top_terms <- res.metascape %>% filter(str_detect(GroupID, '_Summary$'))
top_terms_separate <- top_terms %>% separate_rows(Symbols) %>% 
  dplyr::select(-Genes) %>% dplyr::rename(Genes = Symbols) %>% 
  inner_join(prot_meta)

# Sankey
## ---- pick enriched proteins that are direction-consistent (1 direction only)
prot_pick <- unique(top_terms_separate$protein)
length(prot_pick) # 16
consensus_dir <- filtered12 %>%
  filter(protein %in% prot_pick) %>%
  distinct(protein, direction)

links_pc <- filtered12 %>% dplyr::filter(protein %in% prot_pick) %>%
  dplyr::distinct(protein, cancer_abbr) %>%
  dplyr::count(source = protein, target = cancer_abbr, name = "value") %>%
  dplyr::mutate(group = target)

pc <- links_pc %>%
  filter(source %in% consensus_dir$protein) %>%
  transmute(protein = source, cancer = target) %>%
  distinct()

pt <- top_terms_separate %>%
  filter(protein %in% consensus_dir$protein) %>%
  transmute(protein, term = Description)

triples <- pc %>%
  inner_join(pt) %>%
  inner_join(consensus_dir) %>%
  left_join(top_terms_separate %>% distinct(protein, Genes)) %>%
  mutate(
    protein_lab = str_c(Genes, ' (', protein, ')'),
    term_lab    = term, 42,
    alluvium_id = paste(protein, cancer, term_lab, sep = "|"),
    weight      = 1L # one unit per protein-cancer-term occurrence
  )

## ---- build links (Protein→Cancer, Cancer→Term) with LinkGroup = direction
L1 <- triples %>% count(source = protein_lab, target = cancer,  group = direction, name = "value")
L2 <- triples %>% count(source = cancer,      target = term_lab, group = direction, name = "value")

## ---- nodes & indices
nodes <- tibble::tibble(name = unique(c(L1$source, L1$target, L2$target))) %>%
  mutate(NodeGroup = case_when(
    name %in% unique(triples$protein_lab) ~ "protein",
    name %in% unique(triples$cancer)      ~ "cancer",
    TRUE                                  ~ "term"
  ))

ndx <- setNames(seq_len(nrow(nodes)) - 1L, nodes$name)
links_all <- dplyr::bind_rows(
  L1 %>% mutate(source = ndx[source], target = ndx[target]),
  L2 %>% mutate(source = ndx[source], target = ndx[target])
)

colourScale <- sprintf(
  "d3.scaleOrdinal().domain(%s).range(%s)",
  jsonlite::toJSON(c("up","down","protein","term","cancer","other")),
  jsonlite::toJSON(c("#d73027","#4575b4","#8da0cb","#66c2a5","#666666","#B0B0B0"))
)

## ---- Before plotting: optional minimum flow threshold (keep as 1 if no filtering)
# Recommended 2 or 3 to significantly reduce noise
min_link_value <- 1L
links_all <- links_all[links_all$value >= min_link_value, , drop = FALSE]

## ---- Plot dimensions: adapt height automatically based on node count
# node_pad ↑ original 10 → 22, increases vertical spacing
# font_size smaller labels
# each node ~28 px tall + margin
node_pad    <- 22
font_size   <- 10
nP <- sum(nodes$NodeGroup == "protein")
nC <- sum(nodes$NodeGroup == "cancer")
nT <- sum(nodes$NodeGroup == "term")
nMax <- max(nP, nC, nT)
height_px <- 800 #max(800, 28 * nMax + 160)
width_px  <- 1500

## ---- Sankey diagram
# sinksRight = TRUE: fix sinks on right to reduce crossings
# margin: extra space for long term labels
sankey <- networkD3::sankeyNetwork(
  Links = as.data.frame(links_all),
  Nodes = as.data.frame(nodes),
  Source = "source", Target = "target", Value = "value",
  NodeID = "name", NodeGroup = "NodeGroup", LinkGroup = "group",
  sinksRight = FALSE, fontSize = 12, nodeWidth = 18, nodePadding = 10,
  width = width_px, height = height_px,                                  # ↑ 仅加宽，保持高度
  margin = list(top = 20, right = 380, bottom = 20, left = 300),# ↑ 给两侧标签更多空间
  colourScale = colourScale
)

## —— After plotting: front-end visual refinement only (no data/ID changes)
sankey <- htmlwidgets::onRender(
  sankey,
  "
  function(el,x){
    //var maxLen = 28;  // label truncation length (adjust as needed)
    //var nodes = d3.select(el).selectAll('.node');
    // 1) Semi-transparent blended links to reduce clutter
    d3.select(el).selectAll('.link')
      .style('stroke-opacity', 0.22)
      .style('mix-blend-mode', 'multiply')
      .append('title')
      .text(function(d){
        var s = (d.source && d.source.name) ? d.source.name : '';
        var t = (d.target && d.target.name) ? d.target.name : '';
        var g = (d.group ? (' ['+d.group+']') : '');
        return s + ' → ' + t + ': ' + d.value + g;
      });
  }"
)

## ---- vector PDF (and SVG) via Chrome print
tmp_html <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(sankey, tmp_html, selfcontained = TRUE)
webshot2::webshot(tmp_html, file = "T_NT_shared_DEP_over_16cancers_sankey.pdf", vwidth = width_px, vheight = height_px, zoom = 1)
webshot2::webshot(tmp_html, file = "T_NT_shared_DEP_over_16cancers_sankey.webp", vwidth = width_px, vheight = height_px, zoom = 1)

# not use ggalluvial -style here
if(FALSE){
  # ggalluvial -style
  # ## ---- pick enriched proteins that are direction-consistent (1 direction only)
  # prot_pick <- unique(top_terms_separate$protein)
  # consensus_dir <- filtered12 %>%
  #   filter(protein %in% prot_pick) %>%
  #   distinct(protein, direction)
  # 
  # pc <- links_pc %>%
  #   filter(source %in% consensus_dir$protein) %>%
  #   transmute(protein = source, cancer = target) %>%
  #   distinct()
  # 
  # pt <- top_terms_separate %>%
  #   filter(protein %in% consensus_dir$protein) %>%
  #   transmute(protein, term = Description)
  # 
  # triples <- pc %>%
  #   inner_join(pt) %>%
  #   inner_join(consensus_dir) %>%
  #   left_join(top_terms_separate %>% distinct(protein, Genes)) %>%
  #   mutate(
  #     protein_lab = str_c(Genes, ' (', protein, ')'),
  #     term_lab    = term, 42,
  #     alluvium_id = paste(protein, cancer, term_lab, sep = "|"),
  #     weight      = 1L # one unit per protein-cancer-term occurrence
  #   )
  # 
  # wide_df <- triples %>%
  #   select(alluvium_id, protein_lab, cancer, term_lab, direction, weight)
  # 
  # lodes <- ggalluvial::to_lodes_form(
  #   as.data.frame(wide_df),
  #   axes    = c("protein_lab", "cancer", "term_lab"),
  #   id      = "alluvium_id",
  #   weights = "weight",
  #   key     = "axis",
  #   value   = "stratum"
  # )
  # 
  # lodes$axis <- factor(lodes$axis, levels = c("protein_lab","cancer","term_lab"))
  # lodes$node_type <- ifelse(lodes$axis == "protein_lab", "protein",
  #                           ifelse(lodes$axis == "cancer", "cancer", "term"))
  # 
  # p <- ggplot(lodes,
  #             aes(x = axis, stratum = stratum, alluvium = alluvium_id, y = weight)) +
  #   # flows by direction (vector paths)
  #   ggalluvial::geom_flow(aes(color = direction), alpha = 0.35, show.legend = TRUE) +
  #   # node rectangles by type
  #   ggalluvial::geom_stratum(aes(fill = node_type), width = 0.1, color = "#666666") +
  #   # node labels: proteins (gene symbols), cancers, terms (truncated)
  #   ggalluvial::geom_text(stat = "stratum", aes(label = stratum), size = 3, color = "#222222") +
  #   scale_color_manual(
  #     name   = "Direction",
  #     values = c("up" = "#d73027", "down" = "#4575b4")
  #   ) +
  #   scale_fill_manual(
  #     name   = "Node type",
  #     values = c("protein" = "#8da0cb", "cancer" = "#666666", "term" = "#66c2a5")
  #   ) +
  #   scale_x_discrete(labels = c("Protein", "Cancer", "Term"), expand = c(0.12, 0.02)) +
  #   guides(fill = guide_legend(order = 2), color = guide_legend(order = 1)) +
  #   theme_minimal(base_size = 11) +
  #   theme(
  #     panel.grid = element_blank(),
  #     axis.title = element_blank(),
  #     axis.text.x = element_text(size = 11),
  #     legend.position = "right",
  #     plot.margin = margin(10, 20, 10, 20)
  #   )
  
}

tbl.list <- list(
  DEP_list = filtered12,
  DEP_share = pie_df %>% select(name, Up, Down),
  DEP_share_N = prot_meta,
  Top_terms = top_terms_separate
)
rio::export(tbl.list, 'T_NT_shared_DEP_source.xlsx')


## 6.4 Jaccard ----

### --- helper functions ---------------------------------------------------------

#' Jaccard similarity for two atomic vectors
jaccard <- function(A, B, empty_value = NA_real_) {
  # Returns |A ∩ B| / |A ∪ B|.
  # - Removes duplicates and NA.
  # - If both sets are empty, returns `empty_value` (default NA_real_).
  A <- unique(stats::na.omit(A))
  B <- unique(stats::na.omit(B))
  inter <- length(intersect(A, B))
  uni   <- length(union(A, B))
  if (uni == 0) return(empty_value)
  inter / uni
}

#' Pairwise Jaccard matrix for a list of sets
jaccard_matrix <- function(sets, empty_value = NA_real_) {
  # `sets` is a named list; names used as row/col names in the result.
  stopifnot(is.list(sets), length(sets) > 0)
  nm <- names(sets)
  if (is.null(nm)) nm <- paste0("S", seq_along(sets))
  n  <- length(sets)
  M  <- matrix(1.0, n, n, dimnames = list(nm, nm))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      M[i, j] <- M[j, i] <- jaccard(sets[[i]], sets[[j]], empty_value)
    }
  }
  diag(M) <- 1.0
  M
}

#' Jaccard from an incidence (membership) matrix
jaccard_from_incidence <- function(X, empty_value = NA_real_) {
  # Columns = sets (e.g., cancer types), rows = universe elements (e.g., proteins).
  # Entries: TRUE/1 if row-element belongs to the column-set.
  X <- X != 0
  # Intersections (c_ij) via cross-product of logical matrix
  C <- crossprod(X)                           # integer matrix
  s <- colSums(X)                             # set sizes
  U <- outer(s, s, "+") - C                   # unions (u_ij)
  J <- C / U
  diag(J) <- 1.0
  J[U == 0] <- empty_value                    # handle columns with all FALSE
  dimnames(J) <- list(colnames(X), colnames(X))
  J
}

#' Symmetrize, clip to [0,1], preserve matrix shape and dimnames
.sym_clip01 <- function(M) {
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  stopifnot(nrow(M) == ncol(M))
  d  <- dim(M)
  dn <- dimnames(M)
  M  <- (M + t(M)) / 2
  # Clip diagonal first
  diag(M) <- pmin(1, pmax(0, diag(M)))
  # Elementwise clip while preserving matrix structure
  M[] <- pmin(1, pmax(0, M))
  dim(M) <- d
  dimnames(M) <- dn
  M
}

#' Build a "dist" object from a similarity matrix S (values in [0,1])
#' Uses d_ij = 1 - S_ij; NA in S is set to 0 similarity after symmetrization.
.dist_from_similarity <- function(S) {
  S <- as.matrix(S)
  storage.mode(S) <- "double"
  stopifnot(nrow(S) == ncol(S))
  d  <- dim(S)
  dn <- dimnames(S)
  # sanitize similarity
  S[!is.finite(S)] <- NA_real_
  diag(S) <- 1
  S <- (S + t(S)) / 2
  # Preserve matrix on clipping
  S[] <- pmin(1, pmax(0, S))
  S[is.na(S)] <- 0
  dim(S) <- d
  dimnames(S) <- dn
  
  n <- nrow(S)
  Dmat <- 1 - S
  diag(Dmat) <- 0
  Dmat <- (Dmat + t(Dmat)) / 2
  dv <- Dmat[lower.tri(Dmat, diag = FALSE)]
  structure(dv, Size = n, Diag = FALSE, Upper = FALSE, Labels = rownames(S), class = "dist")
}

#' Convert a Jaccard similarity matrix to a Jaccard distance object
.as_dist_from_jaccard <- function(J) {
  .dist_from_similarity(.sym_clip01(J))
}

#' Build a dist for rows/cols from a similarity matrix J under a chosen metric
#' metric: "jaccard" | "euclidean" | "maximum" | "manhattan" | "canberra" | "minkowski" |
#'         "correlation" | "cosine" | dist-object | function(mat, axis) -> dist
#' axis: "rows" or "cols"
.as_dist_generic <- function(J, metric = "jaccard", axis = c("rows","cols")) {
  axis <- match.arg(axis)
  
  # Local sanitizer; avoids dependency on hidden helpers when sourcing subsets
  sanitize_sim <- function(M) {
    M <- as.matrix(M)
    storage.mode(M) <- "double"
    stopifnot(nrow(M) == ncol(M))
    d  <- dim(M)
    dn <- dimnames(M)
    M[!is.finite(M)] <- NA_real_
    diag(M) <- 1
    M <- (M + t(M)) / 2
    M[] <- pmin(1, pmax(0, M))
    M[is.na(M)] <- 0
    dim(M) <- d
    dimnames(M) <- dn
    M
  }
  
  # Keep a sanitized similarity for distance construction (no NA, diag=1)
  Jsim <- sanitize_sim(J)
  Jdis <- 1 - Jsim
  Xr   <- if (axis == "rows") Jdis else t(Jdis)  # features for Lp distances
  Xsim <- if (axis == "rows") Jsim else t(Jsim)  # similarity for jaccard/corr/cosine
  
  if (inherits(metric, "dist")) {
    n <- if (axis == "rows") nrow(J) else ncol(J)
    stopifnot(isTRUE(attr(metric, "Size") == n))
    return(metric)
  }
  if (is.function(metric)) {
    out <- metric(if (axis == "rows") Jsim else t(Jsim), axis)
    stopifnot(inherits(out, "dist"))
    return(out)
  }
  stopifnot(is.character(metric), length(metric) == 1)
  metric <- tolower(metric)
  
  if (metric == "jaccard") {
    out <- .dist_from_similarity(Xsim)
  } else if (metric %in% c("euclidean","maximum","manhattan","canberra","minkowski")) {
    # Lp distances on dissimilarity features
    out <- stats::dist(Xr, method = metric)
  } else if (metric == "correlation") {
    # 1 - Pearson correlation between row (or column) similarity profiles
    C <- stats::cor(t(if (axis == "rows") Xsim else t(Xsim)), use = "pairwise.complete.obs")
    out <- .dist_from_similarity(.sym_clip01(C))
  } else if (metric == "cosine") {
    # 1 - cosine similarity between row (or column) similarity profiles
    Z <- if (axis == "rows") Xsim else t(Xsim)
    nr <- sqrt(rowSums(Z * Z))
    Z  <- Z / pmax(nr, .Machine$double.eps)
    S  <- tcrossprod(Z)
    out <- .dist_from_similarity(.sym_clip01(S))
  } else {
    stop("Unsupported metric: ", metric)
  }
  
  # Defensive checks to prevent hclust() errors
  n <- if (axis == "rows") nrow(J) else ncol(J)
  stopifnot(inherits(out, "dist"), isTRUE(attr(out, "Size") == n))
  if (any(!is.finite(out))) stop("Distance contains non-finite values; check input J or metric.")
  out
}

#' Label a metric for facet titles
.label_metric <- function(m) {
  if (inherits(m, "dist")) return("dist")
  if (is.function(m))      return("fun")
  as.character(m)
}

#' Heatmaps of a Jaccard similarity matrix with synchronized clustering
#'
#' Rows and columns share the SAME distance and the SAME agglomeration method.
#' Supports trying multiple methods in one call and returns a named list of ggplot objects.
#'
#' J: square Jaccard similarity matrix to display
plot_jaccard_heatmap <- function(
    J,
    distance              = c("jaccard"),      # e.g., "jaccard", "euclidean"
    method                = c("average"),      # one or more of: "ward.D","ward.D2","single","complete","average","mcquitty","median","centroid"
    ann_col               = NULL,
    ann_colors            = NULL,
    cluster_rows          = TRUE,
    cluster_cols          = TRUE,
    show_colnames         = TRUE,
    show_rownames         = TRUE,
    main                  = "Jaccard similarity",
    ...                                        # passed to pheatmap (color, breaks, fontsize, etc.)
) {
  if (!requireNamespace("pheatmap", quietly = TRUE))
    stop("Install 'pheatmap' first: install.packages('pheatmap')")
  if (!requireNamespace("ggplotify", quietly = TRUE))
    stop("Install 'ggplotify' first: install.packages('ggplotify')")
  
  # Matrix for plotting (preserve NA to show gaps if any)
  J_plot <- .sym_clip01(J)
  # Matrix for distance construction (no NA); prevents hclust() length errors
  J_safe <- J_plot
  J_safe[!is.finite(J_safe)] <- 0
  diag(J_safe) <- 1
  
  allowed_methods <- c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
  if (length(setdiff(method, allowed_methods)))
    stop("Unsupported 'method'. Allowed: ", paste(allowed_methods, collapse = ", "))
  
  combos <- expand.grid(
    dist = I(as.list(distance)),
    meth = method,
    stringsAsFactors = FALSE
  )
  
  out <- vector("list", nrow(combos))
  hclust.list <- vector("list", nrow(combos))
  for (i in seq_len(nrow(combos))) {
    dmet <- combos$dist[[i]]
    mth  <- combos$meth[[i]]
    
    # For Ward/centroid/median, use Euclidean distances as recommended by stats::hclust
    use_metric_rows <- dmet
    use_metric_cols <- dmet
    # if (mth %in% c("ward.D","ward.D2","centroid","median")) {
    #   use_metric_rows <- "euclidean"
    #   use_metric_cols <- "euclidean"
    # }
    
    D_rows <- if (isTRUE(cluster_rows)) .as_dist_generic(J_safe, use_metric_rows, "rows") else NULL
    D_cols <- if (isTRUE(cluster_cols)) .as_dist_generic(J_safe, use_metric_cols, "cols") else NULL
    
    HC_rows <- if (isTRUE(cluster_rows)) stats::hclust(D_rows, method = mth) else FALSE
    HC_cols <- if (isTRUE(cluster_cols)) stats::hclust(D_cols, method = mth) else FALSE
    
    ph <- pheatmap::pheatmap(
      mat               = J_plot,
      cluster_rows      = HC_rows,            # hclust objects override clustering_method
      cluster_cols      = HC_cols,
      annotation_col    = ann_col,
      annotation_colors = ann_colors,
      show_colnames     = show_colnames,
      show_rownames     = show_rownames,
      # display_numbers = T, number_color = '#333333',fontsize_number = 8, 
      fontsize = 8,
      main              = sprintf("%s\nrows=cols: %s-%s", main, .label_metric(dmet), mth),
      silent            = TRUE,
      ...
    )
    gg <- ggplotify::as.ggplot(ph$gtable)
    key <- sprintf("dist:%s|method:%s", .label_metric(dmet), mth)
    out[[i]] <- gg
    names(out)[i] <- key
    hclust.list[[i]] <- HC_cols
    names(hclust.list)[i] <- key
  }
  
  class(out) <- c("ggplot.list", class(out))
  list(plots = out, hclust.list = hclust.list)
}

### main -------
dep_list <- split(filtered12$protein, filtered12$cancer_abbr)
J <- jaccard_matrix(dep_list)
saveRDS(J, 'Jaccard.matrix.rds')
# J <- readRDS('Jaccard.matrix.rds')
print(round(J, 3))
ann_col <- filtered12 %>% dplyr::distinct(cancer_abbr) %>% set_rownames(.$cancer_abbr)
ann_colors <- list(
  cancer_abbr = cancer_color#[-which(names(cancer_color) %in% c('CCOC', 'HGSOC'))] %>%
    # append(cancer_color['HGSOC'] %>% setNames('OC'))
)

# Rows and columns will use the SAME distance and SAME hclust method in each panel.
res.jaccard <- plot_jaccard_heatmap(
  J,
  distance   = c("jaccard"),    # try one or more distance definitions
  method     = c("complete", "average"),  # try multiple agglomeration methods
  # ann_col    = ann_col,
  # ann_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE
)
saveRDS(res.jaccard, 'res.jaccard.rds')
# res.jaccard <- readRDS('res.jaccard.rds')
plots <- res.jaccard$plots
length(plots) # 2
plot_jac_hclut <- ggpubr::ggarrange(plotlist = plots, ncol = 2,
                                    nrow = length(plots) / 2, common.legend = T)
ggsave('output/T_NT_shared_DEP_jaccard_heatmap_v1202.pdf', plot_jac_hclut,
       width = 6 * 2, height = 6 * length(plots) / 2)


## 6.5 Selected proteins boxplot (old) -------
## ---- pick enriched proteins that are direction-consistent (1 direction only)
prot_box_show <- c('P09758', 'P11387', 'P09874')

bp_mat_expr <- impute.data %>% 
  select(ID:Age, all_of(prot_box_show)) %>% 
  pivot_longer(cols = prot_box_show,
               names_to = 'protein', values_to = 'Log2Intensity', values_drop_na = T) %>% 
  # inner_join(dfprot %>% rename(protein = Protein.Group), .) %>% 
  inner_join(filtered12 %>%
               select(canbcer_abbr, protein, Genes, g_adj, p_adj_BH) %>%
               mutate(label = str_c(protein, '_', Genes))) %>% 
  mutate(sample_type = factor(sample_type, c('T', 'NT')))

# Build p-values per (protein, cancer_abbr), comparing sample_type
pvals <- bp_mat_expr %>%
  group_by(Genes, protein, cancer_abbr) %>%
  rstatix::wilcox_test(Log2Intensity ~ sample_type, paired = T) %>%
  rstatix::add_xy_position(x = "cancer_abbr", dodge = 0.75) %>%
  mutate(
    p.label = sprintf('%.2e', p),
    # p.label = rstatix::p_format(p, digits = 3),
    panel = str_c(protein)  # match your facet
  )

p_box <- ggplot(bp_mat_expr, aes(x = cancer_abbr, y = Log2Intensity, color = sample_type)) +
  facet_wrap(~ str_c(Genes, ".", protein),
             ncol = 1, scales = "free_y", strip.position = "top") +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  labs(x = "Cancer", y = "Log2Intensity") +
  stat_pvalue_manual(
    data = pvals, label = "p.label", size = 3,
    y.position = "y.position", xmin = "xmin", xmax = "xmax",
    hide.ns = FALSE, tip.length = 0,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    strip.text.y.right = element_text(size = 8, angle = 0, hjust = 0, margin = margin(l = 2)),
    # strip.background = element_blank(),
  )
ggsave('T_NT_DEP_selected_boxplot_v1202.pdf', p_box, width = 15, height = 2*3)


# 7. Cancer-common DEPs -----
## 7.1 Read data -----
res.lm.tbls <- readRDS('output/tumor_nontumor_analysis_reports_v1202.rds')
res_main_filter <- res.lm.tbls$Diff.report.filter
res_main_filter %<>% #mutate(direction = if_else(effect >= 0, "Up", "Down")) %>% 
  inner_join(dfprot %>% rename(protein = Protein.Group), .)
filtered12 <- res_main_filter %>% select(cancer_abbr, protein, Genes, g_adj, p_adj_BH, direction)
dep_list <- split(filtered12$protein, filtered12$cancer_abbr)

# ## 7.2 Skipped: compare enrich -----
library(clusterProfiler)
library(org.Hs.eg.db)

### 7.2.1 compareCluster with DEP-list ------
#### gobp -------
ego_dep <- compareCluster(
  geneClusters = dep_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = FALSE
)
# saveRDS(ego_dep, 'T_NT_DEP_ego_compare.rds')
# ego_dep <- readRDS('T_NT_DEP_ego_compare.rds')
# dotplot(ego_dep, color = "p.adjust", showCategory = 2, font.size = 10)

#### kegg -----
tmp <- clusterProfiler::bitr_kegg(unique(sort(filtered12$protein)), fromType='uniprot', toType='ncbi-proteinid', organism='hsa') %>% 
  setNames(c('protein', 'proteinid'))
filtered12.pid <- tmp1 <- tmp %>%
  count(protein) %>% filter(n == 1) %>% semi_join(tmp, .) %>% 
  count(proteinid) %>% filter(n == 1) %>% semi_join(tmp, .) %>% 
  inner_join(filtered12, .)
# saveRDS(filtered12.pid, 'filtered12.pid.rds')
# filtered12.pid <- readRDS('filtered12.pid.rds')

# KEGG <- clusterProfiler::enrichKEGG(tmp$`ncbi-proteinid`, organism = 'hsa',
#                                     keyType = 'ncbi-proteinid',
#                                     minGSSize = 10, maxGSSize = 500,
#                                     pvalueCutoff = 0.05, pAdjustMethod = 'BH')
# enrichplot::dotplot(KEGG)

dep_pid_list <- split(tmp1$proteinid, tmp1$cancer_abbr)
# saveRDS(dep_pid_list, 'dep_pid_list.rds')
# dep_pid_list <- readRDS(('dep_pid_list.rds'))
ekegg_dep <- compareCluster(dep_pid_list, fun = enrichKEGG,
                            pAdjustMethod = "BH",
                            keyType = 'ncbi-proteinid',
                            minGSSize = 10, maxGSSize = 500,
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)
# saveRDS(ekegg_dep, 'T_NT_DEP_ekegg_compare.rds')
# ekegg_dep <- readRDS('T_NT_DEP_ekegg_compare.rds')
# dotplot(ekegg_dep)


pacman::p_unload(clusterProfiler, enrichplot, ReactomePA, AnnotationDbi, BiocGenerics, IRanges, S4Vectors, Biobase)

### 7.2.2 compute generality & directionality per pathway term ----------------------
#### gobp ----
ego_dep_df <- ego_dep %>%
  filter(pvalue < 0.05, qvalue < 0.05) %>% # in fact redundant
  as.data.frame() %>%
  separate_rows(geneID, sep = '/') %>% 
  rename(cancer_abbr = Cluster,
         protein = geneID) %>% 
  mutate(Pathway.Size = as.numeric(str_extract(BgRatio, '^\\d+'))) %>% 
  inner_join(filtered12 %>% distinct(cancer_abbr, protein, direction), .)

go_dep_summary <- ego_dep_df %>%
  group_by(ID, Description, Pathway.Size) %>%
  summarise(
    generality = n_distinct(cancer_abbr),
    # n_up = n_distinct(cancer_abbr[direction == "Up"]), # cannot be distinct
    # n_down = n_distinct(cancer_abbr[direction == "Down"]), # follow the reference
    n_up = sum(direction == "Up", na.rm = TRUE),
    n_down = sum(direction == "Down", na.rm = TRUE),
    # geneID in compareCluster is "/"-delimited character of SYMBOLs per cancer
    protein.list = list(unique(protein)),
    # max_count = max(Count, na.rm = TRUE),
    med_padj = median(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    directionality = n_up - n_down,
    n_protein = lengths(protein.list),
    label = if_else(generality > 20, Description, NA_character_)
  )

if(FALSE){ # old styles
  # p_go_dep <- ggplot(go_dep_summary) +
  #   aes(x = directionality, y = generality,
  #       colour = -log10(med_padj),
  #       size = Pathway.Size) +
  #   geom_point(alpha = 0.85) +
  #   geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  #   ggrepel::geom_text_repel(aes(label = label), color = 'black', min.segment.length = 1,
  #                            max.overlaps = Inf, size = 1.5) +
  #   scale_size_continuous(name = "Proteins in GO") +
  #   scale_colour_viridis_c(name = "-log10(adj. p)", end = 0.8, direction = -1) +
  #   scale_x_continuous(breaks = scales::pretty_breaks()) +
  #   labs(
  #     x = "Directionality (Up-enriched cancers - Down-enriched cancers)",
  #     y = "Generality (number of cancers enriched)",
  #     title = "Common GO functions among 17 cancers",
  #     subtitle = "clusterProfiler::compareCluster on DEPs"
  #   ) +
  #   theme_bw() +
  #   theme(plot.title = element_text(face = "bold"))
  # p_go_dep <- ggplot(go_dep_summary) +
  #   aes(
  #     x = directionality,
  #     y = generality,
  #     colour = -log10(med_padj),
  #     size   = Pathway.Size
  #   ) +
  #   ggbeeswarm::geom_quasirandom(
  #     alpha = 0.85,
  #     shape = 21,
  #     stroke = 0.3,
  #     fill = NA,
  #     width = 0.15,   # horizontal spread
  #     varwidth = FALSE
  #   ) +
  #   # ggplot(go_dep_summary) +
  #   aes(
  #     x = directionality,
  #     y = generality,
  #     colour = -log10(med_padj),
  #     size   = Pathway.Size
  #   ) +
  #   geom_point(
  #     alpha = 0.85,
  #     shape = 21,        # hollow-ish point
  #     stroke = 0.3,
  #     fill = NA,
  #     position = position_jitter(width = 0.08, height = 0.08)
  #   ) +
  #   geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  #   ggrepel::geom_text_repel(
  #     aes(label = label),
  #     color = "black",
  #     min.segment.length = 1,
  #     max.overlaps = Inf,
  #     size = 1.5
  #   ) +
  #   scale_size_continuous(name = "Proteins in GO") +
  #   scale_colour_viridis_c(name = "-log10(adj. p)", end = 0.8, direction = -1) +
  #   scale_x_continuous(breaks = scales::pretty_breaks()) +
  #   labs(
  #     x = "Directionality (Up-enriched cancers - Down-enriched cancers)",
  #     y = "Generality (number of cancers enriched)",
  #     title = "Common GO functions among 17 cancers",
  #     subtitle = "clusterProfiler::compareCluster on DEPs"
  #   ) +
  #   theme_bw() +
  #   theme(plot.title = element_text(face = "bold"))
  # p_go_dep <- p_go_dep <- ggplot(go_dep_summary) +
  #   aes(
  #     x    = directionality,
  #     y    = generality,
  #     fill = -log10(med_padj),   # 用 fill 而不是 colour
  #     size = Pathway.Size
  #   ) +
  #   # 点：±0.4 抖动
  #   geom_point(
  #     shape    = 21,                     # 有描边、可用 fill
  #     colour   = "grey25",               # 固定描边颜色，便于区分
  #     stroke   = 0.25,                   # 细描边
  #     alpha    = 0.85,
  #     position = position_jitter(width = 0.4, height = 0.4)
  #   ) +
  #   # 竖虚线
  #   geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  #   # 标签也跟着抖动后的位置走
  #   ggrepel::geom_text_repel(
  #     aes(label = label),
  #     color              = "black",
  #     size               = 1.5,
  #     min.segment.length = 1,
  #     max.overlaps       = Inf,
  #     position           = position_jitter(width = 0.4, height = 0.4)
  #   ) +
  #   # size 拉开：小的更小，大的更大
  #   scale_size_continuous(
  #     name  = "Proteins in GO",
  #     range = c(0.5, 6)   # 可根据实际分布再调
  #   ) +
  #   # 用 fill 的 viridis
  #   scale_fill_viridis_c(
  #     name      = "-log10(adj. p)",
  #     end       = 0.8,
  #     direction = -1
  #   ) +
  #   scale_x_continuous(breaks = scales::pretty_breaks()) +
  #   labs(
  #     x        = "Directionality (Up-enriched cancers - Down-enriched cancers)",
  #     y        = "Generality (number of cancers enriched)",
  #     title    = "Common GO functions among 17 cancers",
  #     subtitle = "clusterProfiler::compareCluster on DEPs"
  #   ) +
  #   theme_bw() +
  #   theme(
  #     plot.title = element_text(face = "bold")
  #   )
} # old style figures

# 统一的抖动位置（点和标签都用这个）
pos_jit <- position_jitter(width = 0.4, height = 0.4, seed = 42)

p_go_dep <- ggplot(go_dep_summary) +
  aes(
    x     = directionality,
    y     = generality,
    fill  = -log10(med_padj),  # 用 fill 映射显著性
    size  = Pathway.Size#,
    # alpha = Pathway.Size       # 用 alpha 映射大小，后面反转
  ) +
  # 抖动后的散点
  geom_point(
    alpha = 0.5,
    shape    = 21,             # 可填充、有描边
    colour   = "grey25",       # 固定描边颜色，避免和 fill 冲突
    stroke   = 0.25,
    position = pos_jit
  ) +
  # 中心竖虚线
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    colour     = "grey60"
  ) +
  # 标签：与散点同一个抖动 + 强制画短线
  ggrepel::geom_text_repel(
    aes(label = label),
    position           = pos_jit,   # 与点对齐
    color              = "black",
    size               = 2,
    min.segment.length = 0,         # 一定画线
    segment.size       = 0.25,
    segment.color      = "grey30",
    segment.alpha      = 0.8,
    max.overlaps       = Inf,
    seed               = 42         # ggrepel 自身也设个 seed
  ) +
  # size: 拉开差异
  scale_size_continuous(
    name  = "Proteins in GO",
    range = c(0.5, 6)               # 小的更小
  ) +
  # alpha: 大size → 更透明
  # scale_alpha_continuous(
  #   range = c(0.45, 0.9),           # 0.35(大点) ~ 0.9(小点更实)
  #   trans = "reverse",
  #   guide = "none"                  # 不需要单独的alpha图例
  # ) +
  # fill: viridis
  scale_fill_viridis_c(
    name      = "-Log10(p.adjust.median)",
    end       = 0.8,
    direction = -1
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    x        = "Directionality (Up-enriched cancers - Down-enriched cancers)",
    y        = "Generality (number of cancers enriched)",
    title    = "Common GO functions among 25 cancers",
    subtitle = "Labeled at least in 21 cancers"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )
ggsave(filename = "common_GO_bubble_plot.pdf", plot = p_go_dep, width = 9, height = 6)


#### kegg ----
ekegg_dep_df <- ekegg_dep %>%
  filter(pvalue < 0.05, qvalue < 0.05) %>% # in fact redundant
  as.data.frame() %>%
  separate_rows(geneID, sep = '/') %>% 
  rename(cancer_abbr = Cluster, proteinid = geneID) %>% 
  mutate(Pathway.Size = as.numeric(str_extract(BgRatio, '^\\d+'))) %>% 
  inner_join(filtered12.pid %>% distinct(cancer_abbr, protein, proteinid, direction), .)

kegg_dep_summary <- ekegg_dep_df %>%
  group_by(ID, Description, Pathway.Size) %>%
  summarise(
    generality = n_distinct(cancer_abbr),
    n_up = sum(direction == "Up", na.rm = TRUE),
    n_down = sum(direction == "Down", na.rm = TRUE),
    # geneID in compareCluster is "/"-delimited character of SYMBOLs per cancer
    protein.list = list(unique(protein)),
    # max_count = max(Count, na.rm = TRUE),
    med_padj = median(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    directionality = n_up - n_down,
    n_protein = lengths(protein.list),
    label = if_else(generality >= 15, Description, NA_character_)
  )

# 统一的抖动位置（点和标签都用这个）
pos_jit <- position_jitter(width = 0.4, height = 0.4, seed = 42)

p_kegg_dep <- ggplot(kegg_dep_summary) +
  aes(
    x     = directionality,
    y     = generality,
    fill  = -log10(med_padj),  # 用 fill 映射显著性
    size  = Pathway.Size#,
    # alpha = Pathway.Size       # 用 alpha 映射大小，后面反转
  ) +
  # 抖动后的散点
  geom_point(
    alpha = 0.5,
    shape    = 21,             # 可填充、有描边
    colour   = "grey25",       # 固定描边颜色，避免和 fill 冲突
    stroke   = 0.25,
    position = pos_jit
  ) +
  # 中心竖虚线
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    colour     = "grey60"
  ) +
  # 标签：与散点同一个抖动 + 强制画短线
  ggrepel::geom_text_repel(
    aes(label = label),
    position           = pos_jit,   # 与点对齐
    color              = "black",
    size               = 2,
    min.segment.length = 0,         # 一定画线
    segment.size       = 0.25,
    segment.color      = "grey30",
    segment.alpha      = 0.8,
    max.overlaps       = Inf,
    seed               = 42         # ggrepel 自身也设个 seed
  ) +
  # size: 拉开差异
  scale_size_continuous(
    name  = "Proteins in KEGG",
    range = c(0.5, 6)               # 小的更小
  ) +
  # alpha: 大size → 更透明
  # scale_alpha_continuous(
  #   range = c(0.45, 0.9),           # 0.35(大点) ~ 0.9(小点更实)
  #   trans = "reverse",
  #   guide = "none"                  # 不需要单独的alpha图例
  # ) +
  # fill: viridis
  scale_fill_viridis_c(
    name      = "-Log10(p.adjust.median)",
    end       = 0.8,
    direction = -1
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    x        = "Directionality (Up-enriched cancers - Down-enriched cancers)",
    y        = "Generality (number of cancers enriched)",
    title    = "Common KEGG pathways among 25 cancers",
    subtitle = "Labeled at least in 15 cancers"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )
ggsave(filename = "common_KEGG_bubble_plot.pdf", plot = p_kegg_dep, width = 9, height = 6)


#### output ------
list(GOBP = go_dep_summary %>% select(-protein.list),
     KEGG = kegg_dep_summary %>% select(-protein.list)) %>% 
  rio::export('common_DEP_bubble_plot.xlsx')

list(GOBP = go_dep_summary,
     KEGG = kegg_dep_summary) %>% 
  saveRDS('common_DEP_bubble_plot.rds')
# DEP.common.pathway <- readRDS('common_DEP_bubble_plot.rds')

## 7.2x compare enrich - KEGG -----
library(clusterProfiler)
library(org.Hs.eg.db)

# Build direction-specific clusters and a common universe
filtered12.pid <- readRDS('filtered12.pid.rds')
dep_pid_list <- filtered12.pid %>%
  dplyr::filter(direction %in% c("Up", "Down"),
                !is.na(cancer_abbr), !is.na(proteinid)) %>%
  dplyr::distinct(cancer_abbr, direction, proteinid) %>%
  dplyr::group_by(cancer_abbr, direction) %>%
  dplyr::summarise(genes = list(unique(proteinid)), .groups = "drop") %>%
  dplyr::mutate(Cluster = paste(cancer_abbr, direction, sep = "|")) %>%
  { rlang::set_names(.$genes, .$Cluster) }
# saveRDS(dep_pid_list, 'dep_pid_list_up_down.rds')
# dep_pid_list <- readRDS('dep_pid_list_up_down.rds')

# bg <- filtered12.pid %>%
#   dplyr::filter(!is.na(proteinid)) %>%
#   dplyr::pull(proteinid) %>%
#   unique()

# ORA with the same universe across all clusters (valid for direction)
ekegg_dep <- compareCluster(
  dep_pid_list, fun = enrichKEGG,
  pAdjustMethod = "BH",
  keyType = "ncbi-proteinid",
  minGSSize = 10, maxGSSize = 500,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2#,
  # universe = bg
)
# saveRDS(ekegg_dep, 'T_NT_DEP_ekegg_compare_up_down.rds')
# ekegg_dep <- readRDS('T_NT_DEP_ekegg_compare_up_down.rds')
# dotplot(ekegg_dep)
pacman::p_unload(clusterProfiler, enrichplot, ReactomePA, AnnotationDbi, BiocGenerics, IRanges, S4Vectors, Biobase)


ekegg_dep_df <- ekegg_dep %>%
  clusterProfiler::filter(pvalue < 0.05, qvalue < 0.05, # in fact redundant
         !category %in% c('Human Diseases', 'Organismal Systems')) %>%
  as.data.frame() %>%
  tidyr::separate(Cluster, into = c("cancer_abbr", "direction"), sep = "\\|") %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  rename(proteinid = geneID) %>% 
  mutate(Pathway.Size = as.numeric(str_extract(BgRatio, '^\\d+'))) %>% 
  inner_join(filtered12.pid %>% distinct(cancer_abbr, protein, proteinid, direction), .)

# Summarize by term; direction counts are by DISTINCT cancers (not proteins)
kegg_dep_summary <- ekegg_dep_df %>%
  # dplyr::distinct(ID, Description, Pathway.Size, cancer_abbr, direction, protein, proteinid, .keep_all = TRUE) %>%
  dplyr::group_by(category, subcategory, ID, Description, Pathway.Size) %>%
  dplyr::summarise(
    generality = dplyr::n_distinct(cancer_abbr),
    n_up   = dplyr::n_distinct(cancer_abbr[direction == "Up"]),
    n_down = dplyr::n_distinct(cancer_abbr[direction == "Down"]),
    protein.list = list(unique(protein[!is.na(protein)])),
    RichFactor = stats::median(RichFactor, na.rm = TRUE),
    FoldEnrichment = stats::median(FoldEnrichment, na.rm = TRUE),
    zScore = stats::median(zScore, na.rm = TRUE),
    pvalue = stats::median(pvalue, na.rm = TRUE),
    p.adjust = stats::median(p.adjust, na.rm = TRUE),
    qvalue = stats::median(qvalue, na.rm = TRUE),
    Count = stats::median(Count, na.rm = TRUE),
    geneID = str_c(sort(unique(proteinid)), collapse = '/'), # actually proteinid
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    directionality = n_up - n_down,
    n_protein = lengths(protein.list),
    label = dplyr::if_else(generality >= 15, Description, NA_character_)
  ) %>% 
  as.data.frame()

# Unified jitter position (used for both points and labels)
pos_jit <- position_jitter(width = 0.5, height = 0.5, seed = 42)

p_kegg_dep <- ggplot(kegg_dep_summary) +
  aes(
    x = directionality,
    y = generality,
    fill = -log10(p.adjust),    # Fill mapped to significance
    size = Pathway.Size
  ) +
  geom_point(                   # Jittered points
    alpha = 0.5,
    shape = 21,                 # Fillable circle with border
    colour = "grey25",
    stroke = 0.25,
    position = pos_jit
  ) +
  geom_vline(                   # Center dashed line
    xintercept = 0,
    linetype = "dashed",
    colour = "grey60"
  ) +
  ggrepel::geom_text_repel(     # Labels, aligned with jittered points
    aes(label = label),
    position = pos_jit,
    color = "black",
    size = 2,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.color = "grey30",
    segment.alpha = 0.8,
    max.overlaps = Inf,
    seed = 42
  ) +
  scale_size_continuous(        # Adjust size scaling
    name = "Proteins in KEGG",
    range = c(0.5, 6)
  ) +
  scale_fill_viridis_c(         # Fill color mapping
    name = "-Log10(p.adjust.median)",
    end = 0.8,
    direction = -1
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    x = "Directionality (Up-enriched cancers - Down-enriched cancers)",
    y = "Generality (number of enriched cancers)",
    title = "Common KEGG pathways among 25 cancers",
    subtitle = "Labelled in at least 15 cancers"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

ggsave("common_KEGG_bubble_plot.pdf", plot = p_kegg_dep, width = 9, height = 6)
rio::export(kegg_dep_summary, 'kegg_dep_summary.xlsx')

# devtools::install_github('ievaKer/aPEAR')
library(aPEAR)
kegg_id <- kegg_dep_summary %>% filter(generality >= 10, directionality >= 10) %>% pull(ID)
set.seed(1)
p1 <- enrichmentNetwork(
  kegg_dep_summary %>% filter(ID %in% kegg_id),
  clustMethod = 'hier', # method for detecting pathway clusters. Available methods: 'markov', 'hier' and 'spectral'. Using 'spectral' method requires that you have Spectrum package installed.
  simMethod = 'cosine', # method for calculating similarity matrix between the pathways. Available methods: 'jaccard', 'cosine' and 'cor'
  clustNameMethod = 'pagerank', # method for selecting cluster names. Available methods: 'pagerank', 'hits' and 'none'
  drawEllipses = TRUE, fontSize = 4, nodeSize = 'Pathway.Size',
  colorBy = 'p.adjust', colorType = 'pval',# pCutoff = -5,
  minClusterSize = 2
) +
  scale_color_viridis_c(name = '-Log10(p.adjust)', option = 'C', end = 0.8, limits = c(-10, -2)) +
  scale_size_continuous(name = "Pathway size", range = c(2, 6), limits = c(20, 240))

kegg_id <- kegg_dep_summary %>% filter(generality >= 10, directionality <= -10) %>% pull(ID)
set.seed(1)
p2 <- enrichmentNetwork(
  kegg_dep_summary %>% filter(ID %in% kegg_id),
  clustMethod = 'hier', # method for detecting pathway clusters. Available methods: 'markov', 'hier' and 'spectral'. Using 'spectral' method requires that you have Spectrum package installed.
  simMethod = 'cosine', # method for calculating similarity matrix between the pathways. Available methods: 'jaccard', 'cosine' and 'cor'
  clustNameMethod = 'pagerank', # method for selecting cluster names. Available methods: 'pagerank', 'hits' and 'none'
  drawEllipses = TRUE, fontSize = 4, nodeSize = 'Pathway.Size',
  colorBy = 'p.adjust', colorType = 'pval',# pCutoff = -5,
  minClusterSize = 2
) +
  scale_color_viridis_c(name = '-Log10(p.adjust)', option = 'C', end = 0.8, limits = c(-10, -2)) +
  scale_size_continuous(name = "Pathway size", range = c(2, 6), limits = c(20, 240))

plot.network <- ggpubr::ggarrange(p1, p2, labels = c('Up', 'Down'), common.legend = TRUE)
ggsave('common_KEGG_network_plot.pdf', plot.network, width = 10, height = 5)


## 7.2x compare enrich - GOBP -----
library(clusterProfiler)
library(org.Hs.eg.db)

# Build direction-specific clusters and a common universe
dep_list <- filtered12 %>%
  dplyr::filter(direction %in% c("Up", "Down"),
                !is.na(cancer_abbr), !is.na(protein)) %>%
  dplyr::distinct(cancer_abbr, direction, protein) %>%
  dplyr::group_by(cancer_abbr, direction) %>%
  dplyr::summarise(genes = list(unique(protein)), .groups = "drop") %>%
  dplyr::mutate(Cluster = paste(cancer_abbr, direction, sep = "|")) %>%
  { rlang::set_names(.$genes, .$Cluster) }

# bg <- filtered12 %>%
#   dplyr::filter(!is.na(protein)) %>%
#   dplyr::pull(protein) %>%
#   unique()

# ORA with the same universe across all clusters (valid for direction)
ego_dep <- compareCluster(
  geneClusters = dep_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = FALSE#,
  # universe = bg
)
# saveRDS(ego_dep, 'T_NT_DEP_ego_compare_up_down.rds')
# ego_dep <- readRDS('T_NT_DEP_ego_compare_up_down.rds')
pacman::p_unload(clusterProfiler, enrichplot, ReactomePA, AnnotationDbi, BiocGenerics, IRanges, S4Vectors, Biobase)

ego_dep_df <- ego_dep %>%
  clusterProfiler::filter(pvalue < 0.05, qvalue < 0.05) %>% # in fact redundant
  as.data.frame() %>%
  tidyr::separate(Cluster, into = c("cancer_abbr", "direction"), sep = "\\|") %>%
  separate_rows(geneID, sep = '/') %>% 
  rename(protein = geneID) %>% 
  mutate(Pathway.Size = as.numeric(str_extract(BgRatio, '^\\d+'))) %>% 
  inner_join(filtered12 %>% distinct(cancer_abbr, protein, direction), .)


go_dep_summary <- ego_dep_df %>%
  group_by(ID, Description, Pathway.Size) %>%
  summarise(
    generality = n_distinct(cancer_abbr),
    n_up = n_distinct(cancer_abbr[direction == "Up"]), # cannot be distinct
    n_down = n_distinct(cancer_abbr[direction == "Down"]), # follow the reference
    # n_up = sum(direction == "Up", na.rm = TRUE),
    # n_down = sum(direction == "Down", na.rm = TRUE),
    # geneID in compareCluster is "/"-delimited character of SYMBOLs per cancer
    protein.list = list(unique(protein)),
    # max_count = max(Count, na.rm = TRUE),
    med_padj = median(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    directionality = n_up - n_down,
    n_protein = lengths(protein.list),
    label = if_else(generality > 15, Description, NA_character_)
  )


# Unified jitter position (used for both points and labels)
pos_jit <- position_jitter(width = 0.4, height = 0.4, seed = 42)

p_go_dep <- ggplot(go_dep_summary) +
  aes(
    x = directionality,
    y = generality,
    fill = -log10(med_padj),    # Fill mapped to significance
    size = Pathway.Size
  ) +
  geom_point(                   # Jittered points
    alpha = 0.5,
    shape = 21,                 # Fillable circle with border
    colour = "grey25",
    stroke = 0.25,
    position = pos_jit
  ) +
  geom_vline(                   # Center dashed line
    xintercept = 0,
    linetype = "dashed",
    colour = "grey60"
  ) +
  ggrepel::geom_text_repel(     # Labels, aligned with jittered points
    aes(label = label),
    position = pos_jit,
    color = "black",
    size = 2,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.color = "grey30",
    segment.alpha = 0.8,
    max.overlaps = Inf,
    seed = 42
  ) +
  scale_size_continuous(        # Adjust size scaling
    name = "Proteins in GO",
    range = c(0.5, 6)
  ) +
  scale_fill_viridis_c(         # Fill color mapping
    name = "-Log10(p.adjust.median)",
    end = 0.8,
    direction = -1
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    x = "Directionality (Up-enriched cancers - Down-enriched cancers)",
    y = "Generality (number of enriched cancers)",
    title = "Common GO pathways among 25 cancers",
    subtitle = "Labelled in at least 15 cancers"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

ggsave("common_GO_bubble_plot.pdf", plot = p_go_dep, width = 9, height = 6)
rio::export(go_dep_summary, 'go_dep_summary.xlsx')


## 7.3 high-related cancers, common enrich ------
# dist_mat <- as.dist(1 - J)
# hc <- hclust(dist_mat, method = "average")
cl <- cutree(res.jaccard$hclust.list$`dist:jaccard|method:average`, h = 0.7)
J_cls <- split(colnames(J), cl)
cancer.core <- J_cls$`1`
# [1] "BRCA" "CESC" "COCA" "ENCA" "GIST" "LUCA" "MUT"  "READ" "TGCT" "TOCA"

# Reduce(intersect, dep_list[cancer.core])
dep_list1 <- dep_list[cancer.core]

library(clusterProfiler)
library(org.Hs.eg.db)

### 7.3.1 compareCluster with DEP-list ------
#### gobp -------
ego_dep1 <- compareCluster(
  geneClusters = dep_list1,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = FALSE
)
# saveRDS(ego_dep1, 'T_NT_DEP_ego1_compare.rds')
# ego_dep1 <- readRDS('T_NT_DEP_ego1_compare.rds')
# dotplot(ego_dep1)

#### kegg -----
# tmp <- clusterProfiler::bitr_kegg(unique(sort(filtered12$protein)), fromType='uniprot', toType='ncbi-proteinid', organism='hsa') %>% 
#   setNames(c('protein', 'proteinid'))
# filtered12.pid <- tmp1 <- tmp %>% count(protein) %>% filter(n == 1) %>% semi_join(tmp, .) %>% 
#   inner_join(filtered12, .)
# saveRDS(filtered12.pid, 'filtered12.pid.rds')
filtered12.pid <- readRDS('filtered12.pid.rds')

# KEGG <- clusterProfiler::enrichKEGG(tmp$`ncbi-proteinid`, organism = 'hsa',
#                                     keyType = 'ncbi-proteinid',
#                                     minGSSize = 10, maxGSSize = 500,
#                                     pvalueCutoff = 0.05, pAdjustMethod = 'BH')
# enrichplot::dotplot(KEGG)

dep_pid_list <- readRDS('dep_pid_list.rds')
dep_pid_list1 <- dep_pid_list[cancer.core]
ekegg_dep1 <- compareCluster(dep_pid_list1, fun = enrichKEGG,
                             pAdjustMethod = "BH",
                             keyType = 'ncbi-proteinid',
                             minGSSize = 10, maxGSSize = 500,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2)
# saveRDS(ekegg_dep1, 'T_NT_DEP_ekegg1_compare.rds')
# ekegg_dep1 <- readRDS('T_NT_DEP_ekegg1_compare.rds')
# dotplot(ekegg_dep1)


pacman::p_unload(clusterProfiler, enrichplot, ReactomePA, AnnotationDbi, BiocGenerics, IRanges, S4Vectors, Biobase)

### 7.3.2 compute generality & directionality per pathway term ----------------------
#### gobp ----
ego_dep1_df <- ego_dep1 %>%
  filter(pvalue < 0.05, qvalue < 0.05) %>% # in fact redundant
  as.data.frame() %>%
  separate_rows(geneID, sep = '/') %>% 
  rename(cancer_abbr = Cluster,
         protein = geneID) %>% 
  mutate(Pathway.Size = as.numeric(str_extract(BgRatio, '^\\d+'))) %>% 
  inner_join(filtered12 %>% distinct(cancer_abbr, protein, direction), .)

go_dep_summary <- ego_dep1_df %>%
  group_by(ID, Description, Pathway.Size) %>%
  summarise(
    generality = n_distinct(cancer_abbr),
    n_up = sum(direction == "Up", na.rm = TRUE),
    n_down = sum(direction == "Down", na.rm = TRUE),
    # geneID in compareCluster is "/"-delimited character of SYMBOLs per cancer
    protein.list = list(unique(protein)),
    # max_count = max(Count, na.rm = TRUE),
    med_padj = median(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    directionality = n_up - n_down,
    n_protein = lengths(protein.list),
    label = if_else(generality >= 5, Description, NA_character_)
  )

# 统一的抖动位置（点和标签都用这个）
pos_jit <- position_jitter(width = 0.4, height = 0.4, seed = 42)

p_go_dep <- ggplot(go_dep_summary) +
  aes(
    x     = directionality,
    y     = generality,
    fill  = -log10(med_padj),  # 用 fill 映射显著性
    size  = Pathway.Size#,
    # alpha = Pathway.Size       # 用 alpha 映射大小，后面反转
  ) +
  # 抖动后的散点
  geom_point(
    alpha = 0.5,
    shape    = 21,             # 可填充、有描边
    colour   = "grey25",       # 固定描边颜色，避免和 fill 冲突
    stroke   = 0.25,
    position = pos_jit
  ) +
  # 中心竖虚线
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    colour     = "grey60"
  ) +
  # 标签：与散点同一个抖动 + 强制画短线
  ggrepel::geom_text_repel(
    aes(label = label),
    position           = pos_jit,   # 与点对齐
    color              = "black",
    size               = 2,
    min.segment.length = 0,         # 一定画线
    segment.size       = 0.25,
    segment.color      = "grey30",
    segment.alpha      = 0.8,
    max.overlaps       = Inf,
    seed               = 42         # ggrepel 自身也设个 seed
  ) +
  # size: 拉开差异
  scale_size_continuous(
    name  = "Proteins in GO",
    range = c(0.5, 6)               # 小的更小
  ) +
  # alpha: 大size → 更透明
  # scale_alpha_continuous(
  #   range = c(0.45, 0.9),           # 0.35(大点) ~ 0.9(小点更实)
  #   trans = "reverse",
  #   guide = "none"                  # 不需要单独的alpha图例
  # ) +
  # fill: viridis
  scale_fill_viridis_c(
    name      = "-Log10(p.adjust.median)",
    end       = 0.8,
    direction = -1
  ) +
  scale_x_continuous(limits = c(-5, 15),
                     breaks = scales::pretty_breaks()) +
  scale_y_continuous(limits = c(0, 12)) +
  labs(
    x        = "Directionality (Up-enriched cancers - Down-enriched cancers)",
    y        = "Generality (number of cancers enriched)",
    title    = "Common GO functions among 10 cancers",
    subtitle = "Labeled genarality == 10, directionality >= 5"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )
ggsave(filename = "Cancer10_common_GO_bubble_plot.pdf", plot = p_go_dep, width = 9, height = 6)


#### kegg ----
ekegg_dep1_df <- ekegg_dep1 %>%
  filter(pvalue < 0.05, qvalue < 0.05) %>% # in fact redundant
  as.data.frame() %>%
  separate_rows(geneID, sep = '/') %>% 
  rename(cancer_abbr = Cluster,
         proteinid = geneID) %>% 
  mutate(Pathway.Size = as.numeric(str_extract(BgRatio, '^\\d+'))) %>% 
  inner_join(filtered12.pid %>% distinct(cancer_abbr, protein, proteinid, direction), .)

kegg_dep_summary <- ekegg_dep1_df %>%
  group_by(ID, Description, Pathway.Size) %>%
  summarise(
    generality = n_distinct(cancer_abbr),
    n_up = sum(direction == "Up", na.rm = TRUE),
    n_down = sum(direction == "Down", na.rm = TRUE),
    # geneID in compareCluster is "/"-delimited character of SYMBOLs per cancer
    protein.list = list(unique(protein)),
    # max_count = max(Count, na.rm = TRUE),
    med_padj = median(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    directionality = n_up - n_down,
    n_protein = lengths(protein.list),
    label = if_else(generality >= 8, Description, NA_character_)
  )

# 统一的抖动位置（点和标签都用这个）
pos_jit <- position_jitter(width = 0.4, height = 0.4, seed = 42)

p_kegg_dep <- ggplot(kegg_dep_summary) +
  aes(
    x     = directionality,
    y     = generality,
    fill  = -log10(med_padj),  # 用 fill 映射显著性
    size  = Pathway.Size#,
    # alpha = Pathway.Size       # 用 alpha 映射大小，后面反转
  ) +
  # 抖动后的散点
  geom_point(
    alpha = 0.5,
    shape    = 21,             # 可填充、有描边
    colour   = "grey25",       # 固定描边颜色，避免和 fill 冲突
    stroke   = 0.25,
    position = pos_jit
  ) +
  # 中心竖虚线
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    colour     = "grey60"
  ) +
  # 标签：与散点同一个抖动 + 强制画短线
  ggrepel::geom_text_repel(
    aes(label = label),
    position           = pos_jit,   # 与点对齐
    color              = "black",
    size               = 2,
    min.segment.length = 0,         # 一定画线
    segment.size       = 0.25,
    segment.color      = "grey30",
    segment.alpha      = 0.8,
    max.overlaps       = Inf,
    seed               = 42         # ggrepel 自身也设个 seed
  ) +
  # size: 拉开差异
  scale_size_continuous(
    name  = "Proteins in KEGG",
    range = c(0.5, 6)               # 小的更小
  ) +
  # alpha: 大size → 更透明
  # scale_alpha_continuous(
  #   range = c(0.45, 0.9),           # 0.35(大点) ~ 0.9(小点更实)
  #   trans = "reverse",
  #   guide = "none"                  # 不需要单独的alpha图例
  # ) +
  # fill: viridis
  scale_fill_viridis_c(
    name      = "-Log10(p.adjust.median)",
    end       = 0.8,
    direction = -1
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    x        = "Directionality (Up-enriched cancers - Down-enriched cancers)",
    y        = "Generality (number of cancers enriched)",
    title    = "Common KEGG pathways among 10 cancers",
    subtitle = "Labeled at least in 8 cancers"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )
ggsave(filename = "cancer10_common_KEGG_bubble_plot.pdf", plot = p_kegg_dep, width = 9, height = 6)


#### output ------
list(GOBP = go_dep_summary %>% select(-protein.list),
     KEGG = kegg_dep_summary %>% select(-protein.list)) %>% 
  rio::export('cancer10_common_DEP_bubble_plot.xlsx')

list(GOBP = go_dep_summary,
     KEGG = kegg_dep_summary) %>% 
  saveRDS('cancer10_common_DEP_bubble_plot.rds')

## 7.4 commonly up/down --------
prot_meta <- rio::import('T_NT_shared_DEP_1_22.xlsx')
25*0.6 == 15
25*0.8 == 20
25*0.88 == 22

### 7.4.1 20 up/down DEPs -----
# DEPs dysregulated in at least 20 kinds of cancer
sharen.cutoff <- 22
com.up.down <- prot_meta %>% filter(sharing_n >= sharen.cutoff)
dim(com.up.down) # 179 6

matfc.common <- filtered12 %>% filter(protein %in% com.up.down$protein) %>% 
  mutate(label = str_c(protein, '_', Genes)) %>% 
  mutate(size = abs(g_adj), color = -log10(p_adj_BH),
         direction = ifelse(g_adj > 0, 'up', 'down')) %>% 
  mutate(color = ifelse(color < 8, color, 8))
p <- ggplot(matfc.common, aes(x = cancer_abbr, y = label))+
  geom_point(aes(size = size, fill = color, shape = direction), color = '#FFFFFF', alpha=1) +
  labs(x = '', y = '', title = '') +
  # scale_color_continuous(
  #   name = 'padj', low = brewer.pal(9, 'OrRd')[4], high = brewer.pal(9, 'OrRd')[9]
  # ) +
  # scale_fill_gradientn(
  #   name = '-Log10(p.adj)',
  #   colors = c(brewer.pal(9, 'OrRd')[2],
  #              colorRampPalette(c(brewer.pal(9, 'OrRd')[4], brewer.pal(9, 'OrRd')[9]))(8)),
  #   limits = c(1, 8),
  # ) +
  scale_fill_viridis_c(name = '-Log10(p.adj)', end = 0.8, direction = -1) +
  # scale_color_aaas() +
  # scale_radius(range = c(1, 4.8), name = 'log2FC') +
  scale_shape_manual(values = c(up = 24, down = 25)) +
  scale_size(range = c(.5, 2.5), name = '|g_adj|',
             limits=c(floor(quantile(matfc.common$size)[1]), ceiling(quantile(matfc.common$size)[5]))) +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2, override.aes = list(shape = 24, fill = "white", color = "black")),
    shape = guide_legend(order = 3, override.aes = list(shape = c(25, 24), fill = "white", color = "black"))
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, vjust = 0.5 ,color = 'black', angle = 45),
        axis.text.y = element_text(size = 6, vjust = 0.5 ,color = 'black', angle = 0))
ggsave('output/TPHP_DEP_shared_atleast22_cancer.pdf', p, width = 5, height = 4)

matfc.common %>% 
  rio::export('TPHP_DEP_shared_atleast22_cancer.xlsx')

# ### 7.4.2 Skipped: 15 up/down DEP pathways -----
sharen.cutoff <- 15
com.up.down <- prot_meta %>% filter(sharing_n >= sharen.cutoff)
matfc.common <- filtered12 %>% filter(protein %in% com.up.down$protein) %>% 
  mutate(label = str_c(protein, '_', Genes)) %>% 
  pivot_wider(id_cols = label, names_from = cancer_abbr, values_from = g_adj)

# heatmap
matfc.hm <- matfc.common %>% column_to_rownames('label') %>% as.matrix()
vmax <- 15 # set suitable breaks based on g_adj
hm_cols <- colorRampPalette(brewer.pal(11, "RdYlBu")[c(10:6, 5:1)])(101)  # red→blue diverging
hm_breaks <- seq(-vmax, vmax, length.out = 101)
pheatmap::pheatmap(
  matfc.hm, scale = 'none',
  color = hm_cols, breaks = hm_breaks,
  cluster_rows = T, cluster_cols = T,
  show_rownames = F, show_colnames = T,
  # fontsize_row = 6, fontsize_col = 10,
  fontsize = 10,
  border_color = NA, na_col = "#CCCCCC",
  # annotation_row = ann_row,
  # annotation_col = ann_col,
  # annotation_colors = ann_colors,
  # gaps_row = gaps_row_idx,
  # gaps_col = gaps_col_idx,
  main = paste0('Cancer-common (at least ', sharen.cutoff, ') up/down DEPs', ' (g_adj)'),
  filename = "cancer_common15_DEP_heatmap.pdf", width = 6, height = 4
)

#### 1) pathway analysis --------------------------------------------
# com.up <- com.up.down %>% filter(direction == 'Up')
# com.down <- com.up.down %>% filter(direction == 'Down')
writeClipboard(com.up.down$protein) # for metascape input

metascape_file <- "metascape/15_common_up_down/metascape_result.xlsx"

res.metascape <- rio::import(metascape_file, sheet = 2)

## 1) Only _Summary
metascape_sum <- res.metascape %>%
  filter(str_detect(GroupID, "_Summary$")) %>% 
  slice(1:5)

## 2) Split Symbols, then join with your table
metascape_long <- metascape_sum %>%
  separate_rows(Symbols) %>%
  select(-Genes) %>%
  rename(Genes = Symbols)

## 2a) Some people's com.up.down uses gene names, others use protein
if ("Genes" %in% names(com.up.down)) {
  metascape_long <- metascape_long %>%
    inner_join(com.up.down, by = "Genes")
} else if ("protein" %in% names(com.up.down)) {
  ## Treat gene names in metascape as protein to join
  metascape_long <- metascape_long %>%
    rename(protein = Genes) %>%
    inner_join(com.up.down, by = "protein")
} else {
  stop("com.up.down has neither Genes nor protein, cannot join with metascape.")
}

## —— Check here to avoid discovering 0 rows after long run
if (nrow(metascape_long) == 0L) {
  stop("inner_join result is empty: genes/proteins in metascape do not match those in com.up.down. Fix keys first.")
}

## 3) direction takes filtered12 as reference
if ("g_adj" %in% names(filtered12)) {
  dir_tbl <- filtered12 %>%
    distinct(protein, g_adj) %>%
    mutate(direction = if_else(g_adj > 0, "up", "down"))
} else if ("direction" %in% names(filtered12)) {
  dir_tbl <- filtered12 %>%
    distinct(protein, direction)
} else {
  stop("filtered12 has neither g_adj nor direction, cannot determine up/down regulation.")
}

## 4) Protein → Term
prot_pick <- unique(metascape_long$protein)

pt <- metascape_long %>%
  filter(protein %in% prot_pick) %>%
  transmute(protein, term = Description)

## 5) triples
triples <- pt %>%
  inner_join(dir_tbl, by = "protein") %>%
  left_join(
    metascape_long %>% distinct(protein, Genes),
    by = "protein"
  ) %>%
  mutate(
    protein_lab = str_c(Genes, " (", protein, ")"),
    term_lab    = term,
    alluvium_id = paste(protein, term_lab, sep = "|"),
    weight      = 1L
  )

if (nrow(triples) == 0L) {
  stop("triples is empty: proteins in filtered12 do not appear in metascape.")
}

#### 2) Aggregate first → reduce lines --------------------------------------------
flows <- triples %>%
  count(protein_lab, term_lab, direction, name = "weight")

## If too many, increase threshold
min_flow <- 1L     # change to 2L or 3L to make it cleaner
flows <- flows %>% filter(weight >= min_flow)

#### 3) Order protein by heatmap ------------------------------------------
lbl2pro <- com.up.down %>% mutate(label = str_c(protein, '_', Genes)) %>%
  pull(protein, label)
pro_ordered <- lbl2pro[plot.hm$tree_row$labels[plot.hm$tree_row$order]]

protein_nodes <- flows %>%
  distinct(protein_lab) %>%
  mutate(ord = match(protein_lab, pro_ordered)) %>%
  arrange(is.na(ord), ord, protein_lab) %>%
  pull(protein_lab)

term_nodes <- flows %>%
  distinct(term_lab) %>%
  arrange(term_lab) %>%
  pull(term_lab)

#### 4) Wide → lodes (note id = "flow_id") ------------------------
flows_long <- flows %>%
  mutate(
    protein_lab = factor(protein_lab, levels = protein_nodes),
    term_lab    = factor(term_lab,    levels = term_nodes),
    flow_id     = paste(protein_lab, term_lab, direction, sep = "|")
  )

lodes <- ggalluvial::to_lodes_form(
  flows_long,
  axes    = c("protein_lab", "term_lab"),
  id      = "flow_id",     # ← later use it as alluvium
  weights = "weight",
  key     = "axis",
  value   = "stratum"
) %>%
  mutate(
    axis      = factor(axis, levels = c("protein_lab", "term_lab")),
    node_type = if_else(axis == "protein_lab", "protein", "term")
  )

#### 5) Colors: protein=direction, term=each its own color --------------------------
get_term_palette <- function(n){
  base_n <- min(9L, n)
  base   <- RColorBrewer::brewer.pal(base_n, "Set3")
  if (n <= base_n) return(base[seq_len(n)])
  colorRampPalette(base)(n)
}

term_cols <- get_term_palette(length(term_nodes))
names(term_cols) <- term_nodes

fill_values <- c(
  "dir__up"   = "#d73027",
  "dir__down" = "#4575b4",
  term_cols
)

#### 6) Theme ----------------------------------------------------------
sankey_theme_sci <- function(base_size = 10L){
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid   = element_blank(),
      axis.title   = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks   = element_blank(),
      legend.position = "right",
      plot.margin  = margin(6, 26, 6, 6)
    )
}

#### 7) Actual plot ------------------------------------------------------
p <- ggplot(
  lodes %>%
    mutate(
      fill_key = case_when(
        node_type == "protein" & direction == "up"   ~ "dir__up",
        node_type == "protein" & direction == "down" ~ "dir__down",
        node_type == "term"                          ~ as.character(stratum),
        TRUE                                         ~ "dir__down"
      )
    ),
  aes(
    x        = axis,
    stratum  = stratum,
    alluvium = flow_id,     # ← use the id above, so protein_lab won't be missing
    y        = weight
  )
) +
  ggalluvial::geom_flow(
    aes(color = direction),
    stat          = "alluvium",
    lode.guidance = "frontback",
    knot.pos      = 0.35,
    alpha         = 0.15,   # ← more transparent
    size          = 0.25,   # ← thinner
    lineend       = "round"
  ) +
  ggalluvial::geom_stratum(
    aes(fill = fill_key),
    width = 0.10,
    color = "#666666",
    size  = 0.20,
    alpha = 0.95
  ) +
  geom_text(
    data = dplyr::filter(lodes, axis == "term_lab"),
    aes(label = after_stat(stratum)),
    stat  = ggalluvial::StatStratum,
    size  = 2.7,
    color = "#1f1f1f",
    hjust = -0.05
  ) +
  scale_x_discrete(
    labels = c("protein_lab" = "Protein", "term_lab" = "Pathway / Term"),
    expand = c(0.10, 0.02)
  ) +
  scale_color_manual(
    name   = "Direction",
    values = c("up" = "#d73027", "down" = "#4575b4")
  ) +
  scale_fill_manual(
    values = fill_values,
    guide  = "none"
  ) +
  sankey_theme_sci(base_size = 10L) +
  coord_cartesian(clip = "off")

p

#### 8) Export and keep transparency --------------------------------------------
ggplot2::ggsave(
  "cancer_common15_DEP_heatmap_sankey.pdf",
  p,
  width  = 9,
  height = 5.5,
  device = cairo_pdf   # ← use cairo, do not use default pdf :contentReference[oaicite:2]{index=2}
)


# ## 7.5 Skipped: common DEPs which druggable ----
DEP.common.pathway <- readRDS('common_DEP_bubble_plot.rds')
DEP.common.pathway.rna <- DEP.common.pathway %>% plyr::ldply(.id = 'Dataset') %>% 
  filter(str_detect(Description, 'transcription|RNA'),
         generality >= 10)
DEP.common.pathway.rna.split <- split(unique(unlist(DEP.common.pathway.rna$protein.list)), DEP.common.pathway.rna$Description)

canditar.rna <- unique(unlist(DEP.common.pathway.rna$protein.list)) # target candidates
df.canditar.rna <- plyr::ldply(DEP.common.pathway.rna.split, .id = 'Description', function(x) data.frame(Uniprot = x)) %>%
  inner_join(DEP.common.pathway.rna %>% select(-protein.list, -label))

## 7.6 Selected proteins boxplot (new) -------
# ## ---- pick enriched proteins that are direction-consistent (1 direction only)
# # After load data
# prot_box_show <- c('P09874', 'P11387', 'P09758')
# df_expr <- impute.data %>% 
#   select(patient_ID, cancer_abbr, sample_type, all_of(prot_box_show)) %>% 
#   mutate(sample_type = factor(sample_type, c('NT', 'T')))
# 
# df_log10 <- df_expr %>%
#   mutate_if(is.numeric, function(x) log10(2^x))
# 
# # na_filling <- log10(min(df_prot %>% select(-(cancer_abbr:sample_type)), na.rm = T)) + 
# #   log10(0.5)
# 
# # Get annotation info for all proteins
# protinfo1 <- filtered12 %>%
#   filter(protein %in% prot_box_show, direction == 'Up') %>% # only up-DEPs
#   rename(UniprotID = protein)
# 
# # 'P11387', 'P09758' show be dysregulated in the same cancer types
# ca_selector <- protinfo1 %>% filter(UniprotID %in% c('P11387', 'P09758')) %>% 
#   count(cancer_abbr) %>% filter(n > 1) %>% 
#   pull(cancer_abbr)
# 
# protinfo1 <- protinfo1 %>% filter(cancer_abbr %in% ca_selector | !(UniprotID %in% c('P11387', 'P09758')))
# 
# 
# # Create plots for each protein (all cancer types in one plot)
# plots <- list()
# for(protein in prot_box_show){
#   # Get cancer types and p-values for this protein
#   protein_info2 <- protinfo1 %>% filter(UniprotID == protein)
#   cancer_list <- protein_info2$cancer_abbr
#   anno_list <- lapply(protein_info2$p_adj_BH, function(x) sprintf("%1.2e", x))
#   
#   # Create plot with all relevant cancer types
#   plots[[protein]] <- df_log10 %>%
#     select(patient_ID, cancer_abbr, sample_type, all_of(protein)) %>%
#     plot_CA_boxviolin_multi(cancer_abbrs = cancer_list, 
#                             na_cutoff = 1, 
#                             refer = 'local',
#                             anno_txt = anno_list,
#                             protein_id = protein)
# }
# 
# # Combine plots
# p <- ggpubr::ggarrange(plotlist = lapply(plots, function(x) x$p), 
#                        nrow = length(plots))
# ggsave('TPHP_target_drugs_paired_box_3prots_20251202.pdf', p, width = 4, height = 6)
# ggsave('TPHP_target_drugs_paired_box_3prots_20251202_narrow.pdf', p, width = 2.5, height = 6)
# ggsave('TPHP_target_drugs_paired_box_3prots_20251202_wide.pdf', p, width = 10, height = 6)
# 
# list(protinfo = protinfo1,
#      expr = df_log10) %>% 
#   rio::export('TPHP_target_drugs_paired_box_3prots_20251202.xlsx')

## ---- review version style
# After load data
prot_box_show <- c('P11387', 'P09758', 'P09874')
df_expr <- impute.data %>% 
  select(patient_ID, cancer_abbr, sample_type, all_of(prot_box_show)) %>% 
  mutate(sample_type = factor(sample_type, c('NT', 'T')))

df_log10 <- df_expr %>%
  mutate_if(is.numeric, function(x) log10(2^x))

# Get annotation info for all proteins
protinfo1 <- filtered12 %>%
  filter(protein %in% prot_box_show, direction == 'Up') %>% # only up-DEPs
  rename(UniprotID = protein)

# 'P11387', 'P09758' show be dysregulated in the same cancer types
ca_selector <- protinfo1 %>% filter(UniprotID %in% c('P11387', 'P09758')) %>% 
  count(cancer_abbr) %>% filter(n > 1) %>% 
  pull(cancer_abbr)

protinfo1 <- protinfo1 %>% filter(cancer_abbr %in% ca_selector | !(UniprotID %in% c('P11387', 'P09758')))
protinfo1 %<>%
  filter(UniprotID %in% c('P09874') & cancer_abbr %in% c('BRCA', 'OC', 'ENCA', 'CESC') |
           UniprotID %in% c('P09758', 'P11387') & cancer_abbr %in% c('BRCA', 'ENCA')) %>% 
  mutate(cancer_abbr = factor(cancer_abbr, c('BRCA', 'OC', 'ENCA', 'CESC')))

# Create plots for each protein (all cancer types in one plot)
plots <- list()
for(protein in prot_box_show){
  # Get cancer types and p-values for this protein
  protein_info2 <- protinfo1 %>% filter(UniprotID == protein)
  cancer_list <- protein_info2$cancer_abbr
  anno_list <- lapply(protein_info2$p_adj_BH, function(x) sprintf("%1.2e", x))
  
  # Create plot with all relevant cancer types
  plots[[protein]] <- df_log10 %>%
    select(patient_ID, cancer_abbr, sample_type, all_of(protein)) %>%
    plot_CA_boxviolin_multi(cancer_abbrs = cancer_list, 
                            na_cutoff = 1, 
                            refer = 'local',
                            anno_txt = anno_list,
                            protein_id = protein)
}

# Combine plots
p <- ggpubr::ggarrange(plotlist = lapply(plots, function(x) x$p), nrow = 1, widths = c(1, 1, 1.7))
ggsave('TPHP_target_drugs_paired_box_3prots_20251202.pdf', p, width = 7, height = 2.5)

list(protinfo = protinfo1,
     expr = df_log10) %>% 
  rio::export('TPHP_target_drugs_paired_box_3prots_20251202.xlsx')


## 7.7 Up-DEPs (Re-new) -------
# up-DEPs
df_dep <- rio::import('../6_3paired_tumor_nontumor/output/tumor_nontumor_analysis_reports_v1202.xlsx', sheet = 3)
up <- df_dep %>% filter(direction == 'Up')
protls_upset <- plyr::dlply(up, 'cancer_abbr', function(dfsub) dfsub$protein)
pdf('TPHP_CA_ADJ_dys_overlap_upset_top5_up_v1202.pdf', width = 10, height = 5)
print(
  UpSetR::upset(
    UpSetR::fromList(protls_upset),
    order.by = 'freq',# order.by = c('freq', 'degree'),
    mainbar.y.label = 'Protein number',
    # queries = list(
    #   list(query = UpSetR::intersects, params = list('TE', 'UTE', 'MAG', 'CU', 'MU'), color = "orange", active = T, query.name = 'The greatest DEPs intersection of top 5 organs')),
    # query.legend = "top",
    # group.by = "sets",
    # empty.intersections = "on",
    set_size.show = T, number.angles = 0, text.scale = c(1, 1, 1, 1, 1, 1.5)
  )
)
graphics.off()

pdf('TPHP_CA_ADJ_dys_overlap_upset_20sets_95intersects_up_v1202.pdf', width = 21.0, height = 29.7, paper = 'a4r')
print(
  UpSetR::upset(
    UpSetR::fromList(protls_upset),
    order.by = 'freq',
    mainbar.y.label = 'Protein number',
    nsets = length(protls_upset), nintersects = 95,
    set_size.show = T, number.angles = 0, text.scale = c(1, 1, 1, 1, 1, 1)
  )
)
graphics.off()


up_wide <- up %>% pivot_wider(id_cols = 'protein', names_from = 'cancer_abbr', values_from = 'direction')
up_n <- apply(up_wide, 1, function(x) sum(!is.na(x[-1])))
table(up_n)
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22 
# 1547  772  624  483  425  355  329  298  328  255  301  280  284  242  200  200  157  144  106   63   26    7 
plot(table(up_n), type = 'b', xlab = '# shared tumor types', ylab = '# upregulated proteins')
# ggplot style
df_line <- data.frame(n_organ = as.numeric(names(table(up_n))), n_protein = as.vector(table(up_n)))

p <- ggplot(df_line) +
  aes(x = n_organ, y = n_protein) +
  geom_line(size = 1, color = "#EF562D") +
  geom_point(size = 2, color = "#000000") +
  geom_text(aes(label = n_protein), vjust = -0.5, hjust = 0, size = 2.4) +
  labs(x = '# shared tumor types', y = '# upregulated proteins') +
  theme_minimal() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # axis.line = element_line(colour = 'black'),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = 10, vjust = 0.5 ,color = 'black', angle = 0),
        axis.text.x = element_text(size = 8, vjust = 0.5 ,color = 'black', angle = 0),
        axis.text.y = element_text(size = 8, vjust = 0.5 ,color = 'black', angle = 0),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10, vjust = 0.5 ,color = 'black', angle = 90),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )
ggsave('sF7B_TPHP_DEP_up_overlap_v1202.pdf', p, width = 5, height = 2)
rio::export(df_line, 'sF7B_TPHP_DEP_up_overlap_v1202.xlsx')

quantile(up_n)
# 0%  25%  50%  75% 100% 
# 1    2    5   11   22 
up21 <- up_wide[up_n >= 21, ]
up21_long <- up21 %>%
  pivot_longer(cols = -protein, names_to = 'cancer_abbr', values_to = 'direction', values_drop_na = T) %>% 
  semi_join(up, .)

organ_order <- up21_long %>% count(cancer_abbr) %>% arrange(desc(n)) %>% pull(cancer_abbr) %>% unname()
up21 %<>% select(protein, all_of(organ_order))

up21_gadj <- up21 %>% column_to_rownames('protein') %>% as.matrix()
up21_gadj[!is.na(up21_gadj)] <- 1
up21_gadj <- apply(up21_gadj, 1:2, as.numeric)
for(x in rownames(up21_gadj)){
  for(y in colnames(up21_gadj)){
    if(!is.na(up21_gadj[x, y])){
      up21_gadj[x, y] <- up21_long %>% filter(protein == all_of(x), cancer_abbr == all_of(y)) %>% pull(g_adj)
    }
  }
}
# colnames(up21_gadj) %>% writeClipboard()
# up21_gadj %<>% arrange(desc(BRCA), desc(CESC), desc(ENCA), desc(MUT), desc(OC), desc(TOCA), desc(TGCT), desc(READ), desc(PECA), desc(LARCA), desc(GBCA), desc(DLBCL), desc(THCA), desc(LUCA), desc(HCC), desc(GC), desc(COCA), desc(PACA), desc(GBM), desc(GIST))
up21_gadj %<>% as.data.frame() %>% arrange(across(everything(), desc))


up21_long$protein <- factor(up21_long$protein, levels = rev(rownames(up21_gadj)), ordered = T)
up21_long$cancer_abbr <- factor(up21_long$cancer_abbr, levels = organ_order, ordered = T)
up21_long %<>% arrange(protein, cancer_abbr)



#uniprot -> uniprot_genename
# df_uniprot <- rio::import('Y:/members/jiangwenhao/TPHP/20220908/drug/uniprot/uniprot-download_true_fields_accession_2Cid_2Cprotein_name_2Cgene_na-2022.09.21-14.18.23.22.xlsx')
# df_uniprot12 <- df_uniprot %>% filter(Entry %in% up21_long$protein)
up21_long <- up21_long %>% 
  mutate(label = str_c(protein, Genes, sep = '_'))

up21_long$label <- factor(up21_long$label, levels = unique(up21_long$label), ordered = T)


library(RColorBrewer)
up21_long$color <- -log10(up21_long$p_adj_BH)
up21_long$size <- up21_long$g_adj
p <- ggplot(up21_long, aes(x = cancer_abbr, y = label))+
  geom_point(aes(size = size, color = color), alpha=1) +
  labs(x = 'Protein', y = 'Cancer name') +
  # scale_color_continuous(
  #   name = 'p_adj_BH', low = brewer.pal(9, 'OrRd')[4], high = brewer.pal(9, 'OrRd')[9]
  # ) +
  scale_colour_gradientn(
    name = 'p_adj_BH',
    colors = c(brewer.pal(9, 'OrRd')[2],
               colorRampPalette(c(brewer.pal(9, 'OrRd')[4], brewer.pal(9, 'OrRd')[9]))(8)),
    limits = c(1, NA),
  ) + 
  # scale_radius(range = c(1, 4.8), name = 'g_adj') +
  scale_size(range = c(.1, 2), limits=c(floor(quantile(up21_long$size)[1]), ceiling(quantile(up21_long$size)[5])), name = 'g_adj') +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  ) +
  ggtitle("DEPs upregulated in at least 12 kinds of cancer") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # axis.line = element_line(colour = 'black'),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, vjust = 0.5 ,color = 'black', angle = 45),
        axis.text.y = element_text(size = 6, vjust = 0.5 ,color = 'black', angle = 0),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )
ggsave('TPHP_DEP_up_atleast21_tumors_v1202.pdf', p, width = 21.0, height = 29.7, units = 'cm')

up21_gadj %>%
  rownames_to_column('protein') %>%
  # left_join(dfprot %>% select(`Protein ID`, Gene, `Entry Name`, `Protein Description`), by = c(prot = 'Protein ID')) %>% 
  inner_join(up21_long %>% select(protein, Genes, label) %>% distinct(), .) %>%
  slice(nrow(.):1) %>% 
  rio::export('TPHP_DEP_up_atleast21_tumors_v1202.xlsx')

# 8.PRM -----
## 8.1 DIA-specific PRM-validated -----
### read dia specific proteins -----
# dia_target <- read.csv('TPHP_CA_targetlist_3.csv', check.names = F)
# dia_specific <- dia_target %>% pivot_longer(cols = everything(), names_to = 'cancer_abbr', values_to = 'prot')
cs_map <- rio::import('T_NT_DEP_tumor_specific_v1202.xlsx')
dia_specific <- cs_map %>%
  mutate(type = str_c(direction, 'regulation')) %>% 
  select(cancer_abbr, protein, type) %>%
  rename(prot = protein)

### read prm data ------
df_pair_prm <- readxl::read_excel('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/PRM/PUH_PRM_prot_info_914pair_48prot.xlsx', col_types = c(rep('guess', 10), rep('numeric', 48)), na = '') %>%
  as.data.frame()
length(unique(df_pair_prm$DIA_ID)) # 914
length(unique(df_pair_prm$patient_ID)) # 457

prm <- readxl::read_excel('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/PRM/20230318_TPHP_PRM_dysregulated_filter50NAByOrgan_DEP.xlsx', sheet = 'all')
prm %<>% rename(class_abbr = organ) %>% inner_join(meta %>% distinct(class_abbr, cancer_abbr))

prm_specific_prot <- prm %>% count(prot) %>% filter(n == 1) %>% pull(prot)
# prm_specific_up <- prm %>% filter(prot %in% prm_specific_prot) %>% filter(type == 'Upregulation') %>% select(cancer_abbr, prot)
# prm_specific_down <- prm %>% filter(prot %in% prm_specific_prot) %>% filter(type == 'Downregulation') %>% select(cancer_abbr, prot)

# stat
# length(intersect(dia_specific$prot, prm.data.raw$prot)) # 20
# length(intersect(dia_specific$prot, prm_specific_prot)) # 6
# length(intersect(dia_specific$prot, prm_specific_up$prot)) # 6


prm_specific <- prm %>% count(prot, type) %>% filter(n == 1) %>% semi_join(prm, .)
length(unique(prm_specific$prot)) # 11

prm_specific_up <- prm %>%
  filter(prot %in% prm_specific$prot) %>%
  filter(type == 'Upregulation')
prm_specific_down <- prm %>%
  filter(prot %in% prm_specific$prot) %>%
  filter(type == 'Downregulation')


dia_specific_up <- dia_specific %>% filter(type == 'Upregulation')
dia_specific_down <- dia_specific %>% filter(type == 'Downregulation')

# stat
length(intersect(dia_specific$prot, prm.data.raw$prot)) # 20
length(intersect(dia_specific_up$prot, prm_specific_up$prot))
dia_specific_up %>% inner_join(prm_specific_up) %>% nrow() # 6
length(intersect(dia_specific_down$prot, prm_specific_down$prot)) # 0
dia_specific_down %>% inner_join(prm_specific_down) %>% nrow() # 0

prot_validated <- intersect(dia_specific$prot, prm_specific_up$prot)
tbl1 <- prm_specific_up %>% filter(prot %in% prot_validated) %>% mutate(data = 'PRM')
tbl2 <- dia_specific %>% filter(prot %in% prot_validated) %>% mutate(data = 'DIA')
tbl <- rbind(tbl1, tbl2)
tbl_validated <- tbl %>% select(cancer_abbr, prot, data) %>%
  pivot_wider(id_cols = c('cancer_abbr', 'prot'), names_from = 'data', values_from = 'data')

### boxplot -----
protinfo <- tbl_validated %>% select(cancer_abbr, prot)
protinfo %<>% left_join(prm, by = c('cancer_abbr', 'prot')) %>% rename(pAdj = p_value)
protinfo$sample_type <- 'T'
protinfo %<>% arrange(cancer_abbr, desc(sample_type))


df_prot <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/PRM/20230318_PUH_PRM_submatrix_combine_5_missing.csv') # use submatrix combine matrix
df_prot %<>% left_join(meta %>% distinct(tissue_name, cancer_abbr)) %>% relocate(cancer_abbr, .after = tissue_name) %>% 
  mutate(cancer_abbr = ifelse(organ == 'ovary', 'OC', cancer_abbr),
         cancer_abbr = ifelse(organ == 'mammary gland', 'BRCA', cancer_abbr),
         cancer_abbr = ifelse(organ == 'prostate', 'PRCA', cancer_abbr)) %>% 
  filter(!is.na(cancer_abbr)) %>% 
  mutate(sample_type = ifelse(cancer_type == 'adjacent', 'NT', 'T'), .after = cancer_type)

colnames(df_prot)[1:8]


# add gene name
protinfo %<>% inner_join(dfprot %>% rename(prot = Protein.Group))

# add peptide sequence
pep <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/PRM/20230313TPHP_PRM_pep_matrix_57.csv')
protinfo %<>% inner_join(pep %>% select(Peptide, UniprotID) %>% rename(prot = UniprotID) %>% distinct())

# filter missing ratio < 50%
na_ratio <- c()
for(i in 1:nrow(protinfo)){
  this_prot <- protinfo$prot[i]
  this_cancer_abbr <- protinfo$cancer_abbr[i]
  this_cancer <- protinfo$sample_type[i]
  
  x <- df_prot %>% filter(cancer_abbr == this_cancer_abbr, sample_type == this_cancer) %>%
    pull(all_of(this_prot))
  na_ratio[i] <- sum(is.na(x)) / length(x)
}
protinfo %<>% filter(all_of(na_ratio) < 0.5)
df_prot %<>% dplyr::select(1:patient_ID, all_of(protinfo$prot))


df_info <- df_prot %>% dplyr::select(1:patient_ID)
rownames(df_info) <- df_info$DIA_ID


# remark
df_prot %<>% mutate(sample_type = factor(sample_type, levels = c('T', 'NT')))

# data prepare
df_expr <- df_prot %>%
  select(-DIA_ID) %>%
  select(cancer_abbr, sample_type, everything())


# source('TPHP_source_QC.R')
df_log10 <- df_expr %>% mutate_if(is.numeric, log10)
na_filling <- log10(min(df_prot %>% select(-(1:patient_ID)), na.rm = T)) + log10(0.5)

plots <- list()
for(i in 1:nrow(protinfo)){
  # anno_txt <- str_glue("log2FC={round(protinfo$log2FC[i], 2)},\npAdj={signif(protinfo$pAdj[i], 3)}")
  # anno_txt <- sprintf("%1.2f (%1.2e)", 2 ^ protinfo$log2FC[i], protinfo$pAdj[i])
  anno_txt <- sprintf("%1.2e", protinfo$pAdj[i])
  plots[[i]] <- df_log10 %>%
    select(cancer_abbr, sample_type, patient_ID, all_of(protinfo$prot[i])) %>%
    plot_CA_boxviolin_new(cancer_abbr = protinfo$cancer_abbr[i], pepseq = protinfo$Peptide[i], na_cutoff = 1, refer = 'local',
                          anno_txt = anno_txt, drawline = T)
}

p <- eval(parse(text = stringr::str_c('ggpubr::ggarrange(', stringr::str_c('plots[[', 1:length(plots), ']]$p', collapse = ', '), ', nrow = length(plots))')))
ggsave('TPHP_PRM_CA_specific_expr_box_20251202.pdf', p, width = 297, height = 210, units = "mm")

# write tables
df_log10_output <- df_log10
pm_output <- df_log10_output[, -(1:6)]
pm_output[is.na(pm_output)] <- na_filling + log10(rnorm(sum(is.na(pm_output)), mean = 1, sd = 0.05))
df_log10_output[, -(1:6)] <- pm_output
list(log10 = df_log10, log10_fillna = df_log10_output) %>%
  rio::export('TPHP_PRM_CA_specific_expr_box_20251202.xlsx')

## 8.2 DIA-specific PRM-detected -----
### read prm data ------
prm <- readxl::read_excel('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/PRM/20230318_TPHP_PRM_dysregulated_filter50NAByOrgan.xlsx')
prm %<>% rename(anatomical_classification = organ, prot = protein) %>% inner_join(meta %>% distinct(anatomical_classification, cancer_abbr))

prot_validating <- intersect(prm$prot, dia_specific$prot)
tbl1 <- prm %>% filter(prot %in% prot_validating) %>% mutate(data = 'PRM') %>%
  select(cancer_abbr, prot) %>% mutate(data = 'PRM')
tbl2 <- dia_specific %>% filter(prot %in% prot_validating) %>% mutate(data = 'DIA')
tbl <- rbind(tbl1, tbl2)
tbl_validating <- tbl %>% select(cancer_abbr, prot, data) %>%
  pivot_wider(id_cols = c('cancer_abbr', 'prot'), names_from = 'data', values_from = 'data')

### boxplot -----
protinfo <- tbl_validating %>% select(cancer_abbr, prot)
protinfo %<>% left_join(prm, by = c('cancer_abbr', 'prot')) %>% rename(pAdj = pAdj_t)
protinfo$sample_type <- 'T'
protinfo %<>% arrange(cancer_abbr, desc(sample_type))


df_prot <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/PRM/20230318_PUH_PRM_submatrix_combine_5_missing.csv') # use submatrix combine matrix
df_prot %<>% left_join(meta %>% distinct(tissue_name, cancer_abbr)) %>% relocate(cancer_abbr, .after = tissue_name) %>% 
  mutate(cancer_abbr = ifelse(organ == 'ovary', 'OC', cancer_abbr),
         cancer_abbr = ifelse(organ == 'mammary gland', 'BRCA', cancer_abbr),
         cancer_abbr = ifelse(organ == 'prostate', 'PRCA', cancer_abbr)) %>% 
  filter(!is.na(cancer_abbr)) %>% 
  mutate(sample_type = ifelse(cancer_type == 'adjacent', 'NT', 'T'), .after = cancer_type)

colnames(df_prot)[1:8]

# add gene name
protinfo %<>% inner_join(dfprot %>% rename(prot = Protein.Group))

# add peptide sequence
pep <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/20220908/PRM/20230313TPHP_PRM_pep_matrix_57.csv')
protinfo %<>% inner_join(pep %>% select(Peptide, UniprotID) %>% rename(prot = UniprotID) %>% distinct())


# filter missing ratio < 50%
na_ratio <- c()
for(i in 1:nrow(protinfo)){
  this_prot <- protinfo$prot[i]
  this_cancer_abbr <- protinfo$cancer_abbr[i]
  this_cancer <- protinfo$sample_type[i]
  
  x <- df_prot %>% filter(cancer_abbr == this_cancer_abbr, sample_type == this_cancer) %>%
    pull(all_of(this_prot))
  na_ratio[i] <- sum(is.na(x)) / length(x)
}
protinfo %<>% filter(all_of(na_ratio) < 0.5) %>% arrange(prot, cancer_abbr)
df_prot %<>% dplyr::select(1:patient_ID, all_of(protinfo$prot))


df_info <- df_prot %>% dplyr::select(1:patient_ID)
rownames(df_info) <- df_info$DIA_ID


# remark
df_prot %<>% mutate(sample_type = factor(sample_type, levels = c('T', 'NT')))

# data prepare
df_expr <- df_prot %>%
  select(-DIA_ID) %>%
  select(cancer_abbr, sample_type, everything())


# source('TPHP_source_QC.R')
df_log10 <- df_expr %>% mutate_if(is.numeric, log10)
na_filling <- log10(min(df_prot %>% select(-(1:patient_ID)), na.rm = T)) + log10(0.5)

plots <- list()
for(i in 1:nrow(protinfo)){
  # anno_txt <- str_glue("log2FC={round(protinfo$log2FC[i], 2)},\npAdj={signif(protinfo$pAdj[i], 3)}")
  # anno_txt <- sprintf("%1.2f (%1.2e)", 2 ^ protinfo$log2FC[i], protinfo$pAdj[i])
  anno_txt <- sprintf("%1.2e", protinfo$pAdj[i])
  plots[[i]] <- df_log10 %>%
    select(cancer_abbr, sample_type, patient_ID, all_of(protinfo$prot[i])) %>%
    plot_CA_boxviolin_new(cancer_abbr = protinfo$cancer_abbr[i], na_cutoff = 1, refer = 'local',
                          anno_txt = anno_txt, drawline = T)
}

p <- eval(parse(text = stringr::str_c('ggpubr::ggarrange(', stringr::str_c('plots[[', 1:length(plots), ']]$p', collapse = ', '), ', nrow = length(plots))')))
ggsave('TPHP_PRM_CA_specific_testing_expr_box_20251202.pdf', p, width = 297, height = 210*10, units = "mm", limitsize = FALSE)

# write tables
df_log10_output <- df_log10
pm_output <- df_log10_output[, -(1:6)]
pm_output[is.na(pm_output)] <- na_filling + log10(rnorm(sum(is.na(pm_output)), mean = 1, sd = 0.05))
df_log10_output[, -(1:6)] <- pm_output
list(protinfo = protinfo %>% select(-(p_wilcox:sample_type)),
     log10 = df_log10, log10_fillna = df_log10_output) %>%
  rio::export('TPHP_PRM_CA_specific_testing_expr_box_2025202.xlsx')



# 9.TCGA Curated Pathways ------
# pathway 10 (Table S3)
TCGA_pathway10 <- read_excel_allsheets('//192.168.99.100/share/members/jiangwenhao/TPHP/input/1-s2.0-S0092867418303593-mmc3.xlsx')
TCGA_pathway10_ <- head(TCGA_pathway10, 10)
colnames(TCGA_pathway10_[["Cell Cycle"]]) <- TCGA_pathway10_[["Cell Cycle"]][2, ]
TCGA_pathway10_[["Cell Cycle"]] %<>% slice(-(1:2))


df_pathway <- plyr::ldply(TCGA_pathway10_, .id = 'Pathway')
df_pathway[!(df_pathway$`OG/TSG` %in% c('OG', 'TSG')), 'OG/TSG'] <- NA


# # check
# df_pathway %>% count(`OG/TSG`)
# TCGA_mat <- read_excel_allsheets('//172.16.13.136/share/members/jiangwenhao/TPHP/input/1-s2.0-S0092867418303593-mmc4.xlsx')
# ptn_na <- df_pathway %>% filter(is.na(`OG/TSG`)) %>% pull(Gene) %>% str_c(collapse = '|')
# str_subset(TCGA_mat$`Alteration level`[2, ], ptn_na)
# ptn <- df_pathway %>% filter(!is.na(`OG/TSG`)) %>% pull(Gene) %>% str_c(collapse = '|')
# str_subset(TCGA_mat$`Alteration level`[2, ], ptn)
# 
# df_check <- TCGA_mat$`Alteration level`
# colnames(df_check) <- df_check[2, ]
# df_check %<>% slice(-(1:2))
# colnames(df_check)[-1] %>% str_split('\\.') %>% sapply(tail, 1) %>% unique() %>% intersect((df_pathway %>% filter(is.na(`OG/TSG`)) %>% pull(Gene)))


# # read carcinoma-adjacent difference matrix
# df1 <- read.csv('//172.16.13.136/share/members/yuel/2022/tables/20230226_PUH_submatrix_combine_05NA_paired_diff_matrix_v2.csv', check.names = F) %>% select(-1) # raw scale
# df1[df1$organ == 'cervix_uteri', 'organ'] <- 'cervix uteri'
# colnames(df1)[1:7]
# # [1] ""DIA_ID"                    "sample_type"               "tissue_name"
# # [4] "patient_ID"                "anatomical_classification" "organ"                     "A0A075B6H7"
# 
# 
# # organ x protein matrix (raw scale)
# df_organ <- df1 %>% select(-(DIA_ID:anatomical_classification)) %>%
#   group_by(organ) %>%
#   summarise_all(median, na.rm = T) # median
# 
# df_organ[2, "organ"] == "cervix uteri" # TRUE
# df_organ[2, "organ"] <- "cervix uteri"


## DEP -------
df_dep <- rio::import('../6_3paired_tumor_nontumor/output/tumor_nontumor_analysis_reports_v1202.xlsx', sheet = 3)
# gene_uni_map <- df_dep %>% distinct(protein, Genes) %>% setNames(c('UNIPROT', 'SYMBOL'))

# map uniprotid to symbol
gene_uni_map <- clusterProfiler::bitr(df_dep$protein, fromType = "UNIPROT",toType = "SYMBOL",OrgDb = "org.Hs.eg.db",drop = T)
gene_uni_map %<>% filter(SYMBOL %in% df_pathway$Gene | SYMBOL %in% df_pathway$`Aliases Used`)
gene_uni_map %>% count(UNIPROT) %>% filter(n > 1) # zero rows
gene_uni_map %>% filter(SYMBOL %in% df_pathway$`Aliases Used`) %>% pull(SYMBOL)
# "MST1"

# df_pathway %>% filter(Gene %in% gene_uni_map$SYMBOL)
# df_pathway %>% filter(Gene %in% gene_uni_map$SYMBOL) %>% count(`OG/TSG`)
# intersect(gene_uni_map$SYMBOL, df_pathway$Gene)
# df_pathway %>% filter(Gene %in% gene_uni_map$SYMBOL) %>% count(`OG/TSG`)

# # Swap `Aliases Used` and `Gene` to map `Gene`
# gene_map_alias <- gene_uni_map %>% filter(SYMBOL %in% df_pathway$`Aliases Used`) %>% pull(SYMBOL)
# pos <- which(df_pathway$`Aliases Used` %in% gene_map_alias)
# swap <- df_pathway$Gene[pos]
# df_pathway[pos, 'Gene'] <- df_pathway$`Aliases Used`[pos]
# df_pathway[pos, 'Aliases Used'] <- swap

# # Curated DEPs; remain all UniprotID in PHU dataset
# df_dep1 <- df_dep %>%
#   right_join(., gene_uni_map, by = c('protein' = 'UNIPROT')) # match gene name; remain all UniprotID in PHU dataset
# df_pathway %>%
#   left_join(., gene_uni_map, by = c('SYMBOL' = 'Gene')) %>% # match 10 TCGA curated pathways
#   rename(Gene = SYMBOL) # rename column
# 
# colnames(df_dep1)
# df_dep1 %>% count(organ, Pathway)



# Curated DEPs; remain all UniprotID in PHU dataset
df_dep1 <- df_dep %>%
  right_join(., gene_uni_map, by = c('protein' = 'UNIPROT')) %>% # match gene name; remain all UniprotID in PHU dataset
  inner_join(., df_pathway, by = c('SYMBOL' = 'Gene')) %>% # match 10 TCGA curated pathways
  rename(Gene = SYMBOL) # rename column

df_dep2 <- df_dep %>%
  right_join(., gene_uni_map, by = c('protein' = 'UNIPROT')) %>% # match gene name; remain all UniprotID in PHU dataset
  inner_join(., df_pathway, by = c('SYMBOL' = 'Aliases Used')) %>% # match 10 TCGA curated pathways; with Aliases Used
  rename(`Aliases Used` = SYMBOL) # rename column

df <- bind_rows(df_dep1, df_dep2)
length(unique(df$protein)) == length(unique(gene_uni_map$UNIPROT)) # TRUE


# statistics
df_n <- df %>% count(cancer_abbr, Pathway)
df_n_wide <- df_n %>% pivot_wider(Pathway, names_from = 'cancer_abbr', values_from = 'n')
df_n_wide <- df %>%
  distinct(protein, Pathway) %>%
  count(Pathway, name = 'Total DEPs') %>%
  inner_join(df_n_wide, ., by = 'Pathway')


# add total proteins instead of DEPs
impute.data <- rio::import('output/T_NT_analysis_imputated_data_v1202.xlsx')

gene_uni_map_all <- clusterProfiler::bitr(na.omit(union(df_pathway$Gene, df_pathway$`Aliases Used`)), fromType = "SYMBOL",toType = "UNIPROT",OrgDb = "org.Hs.eg.db",drop = T)
gene_uni_map_all %<>% filter(UNIPROT %in% colnames(impute.data)) # detected in PHU CA dataset
df_map1 <- gene_uni_map_all %>% inner_join(df_pathway, by = c(SYMBOL = 'Gene'))
df_map2 <- gene_uni_map_all %>% inner_join(df_pathway, by = c(SYMBOL = 'Aliases Used'))
df_map <- bind_rows(df_map1, df_map2)
df_n_wide <- df_map %>%
  distinct(UNIPROT, Pathway) %>%
  count(Pathway, name = 'Total detected') %>%
  inner_join(df_n_wide, ., by = 'Pathway')

# map to abbr
# df_abbr <- read.delim('//172.16.13.136/share/members/jiangwenhao/TPHP/input/sample_types_abbr_20230113.txt', stringsAsFactors = F, check.names = F, na.strings = '')
# h <- hash::hash(keys = df_abbr$'Entire', values = df_abbr$'Abbr')
# colnames(df_n_wide) %<>% sapply(function(e){
#   ifelse(!is.null(h[[tolower(e)]]), h[[tolower(e)]], e)
# })
rio::export(df, 'output/TPHP_CA_TCGA_curatedPathway_20251202.xlsx')
# rio::export(df_n_wide, 'TPHP_CA_TCGA_curatedPathway_stat_20251101.xlsx')

# heatmap
df_n_long <- df_n_wide %>% pivot_longer(-Pathway, names_to = 'cancer_abbr', values_to = 'ProteinNumber')
df_n_long$Pathway %<>% factor(levels = unique(.), ordered = T)
df_n_long[is.na(df_n_long)] <- 0

#set order
# organs_ordered <- c('UTE', 'REC', 'PR', 'LU', 'KI', 'CO', 'SI', 'MAG', 'FT', 'LTH', 'OV', 'CE', 'PA', 'MU', 'THR', 'CU', 'THY', 'ST', 'GB', 'TE', 'LI', 'LN', 'ESO', 'TON', 'PEN')
heat_cor <- df_n_wide %>% column_to_rownames('Pathway') %>% select(-`Total DEPs`, -`Total detected`) %>%
  apply(1:2, as.numeric)
heat_cor[is.na(heat_cor)] <- 0
heat_cor <- cor(heat_cor, method = 'spearman')
heat_dist <- as.dist(1 - heat_cor) # 1-Spearman's rho
heat_tree <- hclust(heat_dist, method="complete")
organs_ordered <- heat_tree$labels[heat_tree$order]
organs_ordered %<>% rev() # for mapping to ggplot
df_n_long$cancer_abbr %<>% factor(levels = union(organs_ordered, .), ordered = T)

df_n_wide %<>% select(Pathway, all_of(organs_ordered), `Total DEPs`, `Total detected`)
rio::export(df_n_wide, 'TPHP_CA_TCGA_curatedPathway_stat_20251202.xlsx')

X <- df_pathway %>% select(Gene, Pathway) %>% count(Pathway, name = 'Total in TCGA') %>% 
  inner_join(df_n_wide, ., by = 'Pathway')
rio::export(X, 'TPHP_CA_TCGA_curatedPathway_stat_20251202_v2.xlsx')

# heatmap
df_heat <- df_n_long %>% filter(cancer_abbr %in% organs_ordered)

p <- ggplot(df_heat) + 
  geom_tile(aes(cancer_abbr, Pathway, fill= ProteinNumber)) +
  geom_text(aes(cancer_abbr, Pathway, label = ifelse(ProteinNumber > 0, ProteinNumber, NA)), size = 5, color = 'black') +
  scale_fill_gradient(low="white", high="#3388AA") +
  labs(
    x = 'Cancer type',
    y = 'TCGA curated pathways'
  )+
  scale_x_discrete(expand = c(0, 0))+
  coord_cartesian(xlim = c(0, length(unique(df_n_long$cancer_abbr))) * 1.05,
                  ylim = c(0, length(unique(df_n_long$Pathway)) * 1.3))+
  #theme_bw()+
  theme(panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12,color = 'black'),
        axis.line = element_line(color = 'black'),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        plot.subtitle = element_text(size = 30, hjust = 0, color = 'black')
  )+
  theme(strip.text.x = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 12, color = 'black'), legend.position = 'right',
        legend.title = element_text(size = 15, color = 'black'))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, vjust = 0.5 ,color = 'black', angle = 90),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )+
  theme(axis.title.y = element_text(size = 12, hjust = 0.5, color = 'black', angle = 90),
        axis.text.y = element_text(size = 10,color = 'black', angle = 0),
        axis.ticks.y = element_blank()
  )

p1 <- p +
  annotate("text", x = length(unique(df_heat$cancer_abbr)) + 1, y = 1:length(unique(df_n_long$Pathway)), parse = T, label = filter(df_n_long, cancer_abbr == 'Total DEPs') %>% pull(ProteinNumber), size = 5, color = "black") +
  annotate("text", x = length(unique(df_heat$cancer_abbr)) + 2, y = 1:length(unique(df_n_long$Pathway)), parse = T, label = filter(df_n_long, cancer_abbr == 'Total detected') %>% pull(ProteinNumber), size = 5, color = "black") +
  annotate("text", x = length(unique(df_heat$cancer_abbr)) + 1, y = length(unique(df_n_long$Pathway)) * 1.05, label = 'Total DEPs', size = 5, color = "black", vjust = 0.5, hjust = 0, angle = 90) +
  annotate("text", x = length(unique(df_heat$cancer_abbr)) + 2, y = length(unique(df_n_long$Pathway)) * 1.05, label = 'Total detected', size = 5, color = "black", vjust = 0.5, hjust = 0, angle = 90)

ggsave('TPHP_CA_TCGA_curatedPathway_stat_20251202.pdf', p1, width = 12, height = 6.5)

# 10.RTKs ------
df_dea <- rio::import('../6_3paired_tumor_nontumor/output/tumor_nontumor_analysis_reports_v1202.xlsx', sheet = 2)

df <- rio::import('input/RTKs.xlsx')
rtk_id <- df$From
df_label <- df_dea %>%
  filter(protein %in% rtk_id) %>%
  select(cancer_abbr, protein, Genes, g_adj, p_adj_BH)

df <- df_label %>%
  select(cancer_abbr, protein, g_adj) %>%
  pivot_wider(names_from = 'cancer_abbr', values_from = 'g_adj') %>%
  column_to_rownames('protein')


df_label_p <- df_label %>%
  select(protein, p_adj_BH, cancer_abbr)

df_label_p$p_label[df_label_p$p_adj_BH < 0.001] <- '***'
df_label_p$p_label[df_label_p$p_adj_BH >= 0.001 & df_label_p$p_adj_BH < 0.01] <- '**'
df_label_p$p_label[df_label_p$p_adj_BH >= 0.01 & df_label_p$p_adj_BH < 0.05] <- '*'
df_label_p %<>%
  pivot_wider(protein, names_from = cancer_abbr, values_from = p_label) %>%
  column_to_rownames('protein')
df_label_p[is.na(df_label_p)] <- ''


clusterProfiler::bitr(rownames(df), fromType = "UNIPROT",toType = "SYMBOL",OrgDb = "org.Hs.eg.db",drop = T) %>% 
  rename(protein = UNIPROT) %>% inner_join(df_label %>% distinct(protein, Genes))
rownames(df) <- clusterProfiler::bitr(rownames(df), fromType = "UNIPROT",toType = "SYMBOL",OrgDb = "org.Hs.eg.db",drop = T) %>% unite('label', UNIPROT, SYMBOL) %>% pull(label)
rownames(df_label_p) <- clusterProfiler::bitr(rownames(df_label_p), fromType = "UNIPROT",toType = "SYMBOL",OrgDb = "org.Hs.eg.db",drop = T) %>% unite('label', UNIPROT, SYMBOL) %>% pull(label)

identical(rownames(df), rownames(df_label_p)) # TRUE


dys_cancer_abbr_num <- apply(df_label_p, 1, function(e) sum(str_detect(e, '\\*'))) %>% unname
df <- df[which(dys_cancer_abbr_num != 0), ]
df_label_p <- df_label_p[which(dys_cancer_abbr_num != 0), ]

quantile(df, na.rm = T)
#         0%        25%        50%        75%       100% 
# -1.8073838 -0.1615253  0.3263752  0.7077229  2.1868337 
hist(as.matrix(df))


#breaks
my_breaks <- seq(-2, 2, by = 0.01)

#colors
my_colors <- c(colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[1:3]))(length(my_breaks)/2/8*6),
               colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[3:5], 'white'))(length(my_breaks)/2/8*2),
               colorRampPalette(colors = c('white', RColorBrewer::brewer.pal(11, "PiYG")[7:9]))(length(my_breaks)/2/8*2),
               colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[9:11]))(length(my_breaks)/2/8*6)
)
my_colors %<>% rev() # high intensity should be hot while low to be cold

df_imp <- df; df_imp[is.na(df_imp)] <- 0
a <- pheatmap::pheatmap(df_imp, fontsize = 16,
                        # color = colorRampPalette(colors = c("#74ADD1", "#ABD9E9", "#E0F3F8", "white", "#FDAE61", "#F46D43", "#D73027"))(100),
                        color = my_colors,
                        breaks = my_breaks,
                        display_numbers = df_label_p, fontsize_number = 25, number_color = '#000000',angle_col = 90,
                        cluster_rows = T, cluster_cols = T,na_col = "#DDDDDD",
                        filename = 'TPHP_rtks_heatmap_v1202.pdf',
                        width = 16, height = 9)

list(g_adj = df[a$tree_row$order, a$tree_col$order],
     p_adj_BH = df_label_p[a$tree_row$order, a$tree_col$order],
     long_data = df_label) %>%
  rio::export('TPHP_rtks_heatmap_v1202.xlsx')

# 11. imputated data proteins ------
t.nt.input <- readRDS('output/tumor_nontumor_input_v1202.rds')
t.nt.input.list <- split(t.nt.input, t.nt.input$cancer_abbr)
t.nt.input.identity <- plyr::ldply(t.nt.input.list, function(sub.input){
  sub.input_imp <- preproc_fill_na_by_group(
    df = sub.input, group_col = "cancer_abbr", na_cutoff = 0.50, noise_sd = 0, seed = 2022, min_log2 = min_log2_tnt
  )
  sample_type_labels <- sub.input_imp$sample_type
  sub.input_imp %<>% select(-(1:Dataset))
  sub.input_mask <- sub.input_imp < min_log2_tnt
  data.frame(sample_type = sample_type_labels, proteins = rowSums(!sub.input_mask))
}, .id = 'cancer_abbr')
t.nt.input.identity %<>% rename(`# proteins` = proteins)

plot_identity <- ggplot(t.nt.input.identity) +
  aes(x = `# proteins`, y = cancer_abbr, color = sample_type) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.1, position = position_jitterdodge(jitter.width = 0.3, seed = 0)) +
  scale_color_manual(values = sample_color) +
  labs(y = '') +
  theme_classic(base_size = 10)
ggsave('T_NT_filterNA_50_proteins.pdf', plot_identity, width = 4, height = 5.5)
rio::export(t.nt.input.identity, 'T_NT_filterNA_50_proteins.xlsx')

# EOF ----
pacman::p_unload(clusterProfiler, enrichplot, ReactomePA, AnnotationDbi, BiocGenerics, IRanges, S4Vectors, Biobase)

sapply(res.all.list, function(X){
  length(unique(X$impute.data$protein))
})
# BRCA  CESC  COCA DLBCL  ENCA  ESCA  FTCA  GBCA   GBM    GC  GIST   HCC LARCA  LUCA   MUT    OC  PACA  PECA  PRCA    RC  READ  TGCT  THCA 
# 6220  6268  7242  7585  6777  6173  6451  6660  7144  7107  7287  6862  6876  7111  6368  6619  7289  6235  6126  7258  6791  6813  6011 
# THYM  TOCA 
# 6807  6171 



