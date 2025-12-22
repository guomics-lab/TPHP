pkgs <- c("tidyverse","arrow","dplyr","tidyr","lmerTest","lme4","plyr")
sapply(pkgs, function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
  return(TRUE)
})

library(tidyverse)

# helper functions ------
safe_div <- function(a, b) {
  ifelse(is.finite(a) & is.finite(b) & b != 0, a / b, NA_real_)
}

J <- function(df) {
  ifelse(
    is.finite(df) && df > 1,
    exp(lgamma(df/2) - 0.5*log(df/2) - lgamma((df-1)/2)),
    NA_real_
  )
}

# Prepare data ------
paired.data <- arrow::read_parquet('data/tumor_compare_data.parquet')
compare.data <- paired.data %>% 
  select(patient_ID, sample_type, cancer_abbr, cancer_subtype,
         Gender, Age, Dataset, everything()) %>% 
  mutate(
    patient_ID = as.character(patient_ID),
    sample_type = factor(as.character(sample_type), levels = c("NT", "T")),
    cancer_abbr = as.character(cancer_abbr),
    cancer_subtype = as.character(cancer_subtype),
    Gender = factor(Gender),
    Age = as.integer(Age),
    Dataset = factor(Dataset)
  )

# Linear mixed model --------
vec_ca <- sort(unique(compare.data$cancer_abbr))

res.all.list <- lapply(vec_ca, function(nm) {
  dfsub <- compare.data %>% filter(cancer_abbr == nm)
  
  long_dfsub <- dfsub %>%
    tidyr::pivot_longer(
      cols = -(patient_ID:Dataset),
      names_to = "protein",
      values_to = "value"
    )
  
  wide_dfsub <- long_dfsub %>%
    tidyr::pivot_wider(
      id_cols = c(cancer_abbr, cancer_subtype, patient_ID,
                  Gender, Age, Dataset, protein),
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
    st_ok <- length(unique(na.omit(d$sample_type))) >= 2
    if (!st_ok) return(data.frame())
    
    g2         <- all(c("F", "M") %in% unique(na.omit(d$Gender)))
    subtype_ok <- length(unique(na.omit(d$cancer_subtype))) >= 2
    ds_ok      <- length(unique(na.omit(d$Dataset))) >= 2
    age_ok <- {
      x <- d$Age_c
      x <- x[is.finite(x)]
      length(unique(na.omit(x))) >= 2
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
    # --- Fit Linear Mixed-Effects Models ---
    
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
res.df <- dplyr::bind_rows(lapply(res.all.list, function(D) D$res))
res_main <- res.df %>% mutate(direction = ifelse(g_adj > 0, "Up", "Down"))
res_main_filter <- res_main %>% filter(abs(g_adj) >= 0.5, p_adj_BH < 0.05)

DEA <- list(Diff.report = res_main,
            Diff.report.filter = res_main_filter)

saveRDS(DEA, 'output/compare_report_output.rds')


# Version check ----
pkgs_ver <- sapply(pkgs, function(pkg) as.character(packageVersion(pkg)))

txt <- c(
  paste("R: ", R.version$major, ".", R.version$minor, sep = ""),
  "",
  paste0(pkgs, ": ", pkgs_ver)
)

writeLines(txt, "output/package_versions.txt")

