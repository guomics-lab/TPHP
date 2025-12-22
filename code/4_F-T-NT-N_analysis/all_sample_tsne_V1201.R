# rm(list = ls())
# pacman::p_unload(pacman::p_loaded(), character.only = T)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source('//172.16.13.136/share/members/jiangwenhao/code/myQC.R')
source('../source/source_code.R')


dfprot <- rio::import('../../TPL/libs/20220616_fraglib/protein.tsv')
reported_protinfo <- rio::import('../0_process_DIA-NN_data/output/protein_info_from_reports_V1009.csv')


## helper functions -----
my_uniprot2symbol <- function(vec_uni, refer = 'uniprot'){
  if(refer == 'uniprot'){
    # crawl protein name from uniprotid on uniprot.org
    vec_prot <- lapply(vec_uni, function(uniprotid){
      # my_html <- rvest::read_html(stringr::str_glue('https://www.uniprot.org/uniprotkb/{uniprotid}/entry'))
      my_html <- rvest::read_html(stringr::str_glue('https://rest.uniprot.org/genecentric/{uniprotid}')) # 2022-08-08 update; get json
      my_json <- my_html %>%
        rvest::html_text() %>%
        jsonlite::fromJSON()
      
      symbol <- my_json$canonicalProtein$geneName
      return(symbol)
    }) %>% unlist
  }
  if(refer == 'local'){
    vec_prot <- clusterProfiler::bitr(my_uni, fromType = "UNIPROT",toType = "SYMBOL",OrgDb = "org.Hs.eg.db",drop = T) %>% pull(SYMBOL)
  }
  return(vec_prot)
}



# 1. all samples t-SNE -------------------------

## 1.1 t-SNE ---------------------------------------------------------------
info <- rio::import("../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V6.xlsx")
info %<>% filter(!Tissue_heterogeneity_abnormal, !Low_protein_IDs, !FileName %in% pancancer_low_pearson)
# DIA-NN protein-group matrices
# Original (log2 expected in file); if raw, we transform below
file_list <- list.files(path = 'input/matrices_input', full.names = T)
names(file_list) <- c(
  'original', 'qimp.cmb.mo', 'qimp.cmb', 'qimp.lma', 'qimp.lma.s2',
  'filtered', 'imp', 'qimp'
)
pm_list <- lapply(file_list, function(X){
  message('Reading file ', X, '...\n')
  rio::import(X) %>% column_to_rownames(colnames(.)[1])
})
names(pm_list) <- names(file_list)

# transpose
sapply(pm_list, dim)
# qimp.cmb.mo qimp.cmb qimp.lma qimp.lma.s2   imp  qimp
# [1,]       12754    12754    12754       12754 12754 12754
# [2,]        2997     2997     2997        2997  2997  2997
pm_list$filtered %<>% t()

# log2-transform
sapply(pm_list, quantile, na.rm = T)
#          original qimp.cmb.mo  qimp.cmb  qimp.lma qimp.lma.s2     filtered       imp     qimp
# 0%   1.277158e+01    5.921994 -1.310841  3.949413    3.161037 1.277158e+01  2.674865  6.38579
# 25%  1.362713e+03    6.407567  6.388891  6.386089    6.401242 1.367858e+03  2.674865  6.38579
# 50%  3.397383e+03    6.490974  6.971537  6.892888    7.064453 3.406983e+03  2.674865  6.38579
# 75%  1.076838e+04   12.365733 12.357216 12.365974   12.446083 1.079430e+04 11.688395 12.35507
# 100% 4.662632e+08   22.872299 33.539784 24.994231   24.906838 4.662632e+08 28.796569 22.61057
pm_list$original %<>% log2()
pm_list$filtered %<>% log2()
sapply(pm_list, quantile, na.rm = T)
#       original qimp.cmb.mo  qimp.cmb  qimp.lma qimp.lma.s2  filtered       imp     qimp
# 0%    3.674865    5.921994 -1.310841  3.949413    3.161037  3.674865  2.674865  6.38579
# 25%  10.412266    6.407567  6.388891  6.386089    6.401242 10.417703  2.674865  6.38579
# 50%  11.730208    6.490974  6.971537  6.892888    7.064453 11.734279  2.674865  6.38579
# 75%  13.394513   12.365733 12.357216 12.365974   12.446083 13.397983 11.688395 12.35507
# 100% 28.796569   22.872299 33.539784 24.994231   24.906838 28.796569 28.796569 22.61057


# Filter to samples present in metadata
info1 <- info %>% dplyr::filter(sample_type != 'p')
mat_list <- lapply(pm_list, function(X){
  X[, info1$FileName, drop = FALSE]
})
mat_list1 <- mat_list[-which(names(mat_list) %in% c('original', 'filtered'))]
# saveRDS(mat_list1, 'output/mat_list1.rds')
# saveRDS(info1, 'output/info1.rds')

# colors of major tissue types
# df_color <- rio::import('../0_process_DIA-NN_data/input/labels/PUH_tissue_colorset_20251009.xlsx')
# organ_color <- str_c('#', df_color$color) %>%
#   setNames(df_color$class_abbr) %>% .[!is.na(names(.))]
class_abbr2full <- info1 %>% pull(anatomical_classification, class_abbr)
major_color <- organ_color %>% setNames(class_abbr2full[names(.)])



## primary tSNE ------
ge.plot.tsne<-function(data,type,title="",label=NULL,seed=2023){
  col2<-rep(brewer.pal(9,"Set1")[1:9],3)[1:4]
  cl <- data.frame(col2,row.names =unique(type),stringsAsFactors = F)
  cl2 <- cl[match(type,row.names(cl)),1]
  set2<-c(rep(0,times=4))
  st <- data.frame(set2,row.names =unique(type),stringsAsFactors = F)
  st2<-st[match(type,row.names(cl)),1]
  df10 <- data
  df10[is.na(df10)] <- 0 #
  df10 <- t(apply(df10, 1, scale))  #
  colnames(df10) <- type
  set.seed(seed)
  tsne <- Rtsne(
    t(df10), dims = 2,
    # perplexity = 20, # not more than (ncol(data)-1)/3-1 # for V1 params
    perplexity = 10, eta = 200, exaggeration_factor = 16, theta = 0.5, max_iter = 2000, # tested params (V2)
    verbose = T , check_duplicates = FALSE
  )
  pdf(paste0(title,"_TSNE.pdf"),height = 7,width = 10)
  if(is.null(label)){
    par(mai=c(0.4,0.5,0.2,2.4))
    plot(tsne$Y,col=cl2, main = "tsne", pch =st2,cex=0.6,cex.axis=2,cex.lab=2)
    legend(x=5.1,y=5,legend=row.names(cl),pch=set2, cex=0.8,pt.cex=1,col=cl$col2,ncol=1,xpd=T,
           x.intersp=0.5,bty='n')
  }else{
    plot(tsne$Y,col=cl2, main = "tsne", pch = 20,cex=2,cex.axis=2,cex.lab=2,xlab='t-SNE 1',ylab='t-SNE 2')
    text(tsne$Y,pos = 1, labels = label, col= "DimGrey",cex = 0.8)
    legend("topright",legend=row.names(cl), fill=cl$col2, lty=1,lwd=1)
  }
  dev.off()
  
  
  return(tsne$Y)
}

tsne_list <- lapply(names(mat_list1), function(nm){
  message('Primary tSNE of matrix ', nm, '...\n')
  dat1 <- mat_list[[nm]]
  tsne<-ge.plot.tsne(dat1, type = info1$sample_type,
                     title = str_c('output/20251009_PUH_tsneV2_2772samples_sample_type_', nm))
  tsne_df<-as.data.frame(cbind(info1,tsne))
  rio::export(tsne_df, str_c('output/20251009_PUH_tsneV2_2772samples_sample_type_', nm, '.xlsx'))
  return(tsne_df)
})
names(tsne_list) <- names(mat_list1)
saveRDS(tsne_list, 'all_tissues_tsneV2_list.rds')

# New plot with stat_ellipse
# change sample type labels
major_type_labeled <- c('brain', 'kidney', 'liver', 'lung', 'muscle', 'pancreas',
                        'thyroid', 'small intestine', 'stomach', 'thymus', 'other tissues')
tmp_lbl_unmapped <- setdiff(major_type_labeled, info1$anatomical_classification) # "other tissues
ann_tsne <- info1 %>% mutate(major_type = ifelse(anatomical_classification %in% major_type_labeled, anatomical_classification, 'other tissues'))

tmp_lbl_mapped <- intersect(major_type_labeled, names(major_color))
major_type_labeled_colors <- major_color[tmp_lbl_mapped] %>%
  append(#major_color %>% .[-which(names(.) %in% tmp_lbl_mapped)] %>% head(length(tmp_lbl_unmapped)) %>% 
    c('#DDDDDD') %>%
      setNames(tmp_lbl_unmapped))

for(nm in names(tsne_list)){
  message('Scatter plot generating of ', nm, '...\n')
  
  tsne_df <- tsne_list[[nm]] %>%
    inner_join(ann_tsne) %>%
    mutate(sample_type = factor(sample_type, c('F', 'T', 'NT', 'N')),
           major_type = factor(major_type, major_type_labeled))
  tsne_df %<>% dplyr::rename(X = `1`, Y = `2`)
  
  p <- ggplot(tsne_df) +
    aes(x = X, y = Y, color = sample_type) +
    geom_point(shape = "circle", size = 1.5) +
    labs(x = "t-SNE 1", y = "t-SNE 2") +
    guides(color = guide_legend(title = "Sample type")) +
    scale_color_manual(values = c(`T` = '#E41A1C', `NT` = '#377EB8', `N` =  '#4DAF4A', `F` = '#984EA3')) +
    theme_minimal() +
    theme(legend.position = 'right',
          legend.text = element_text(size = 20, color = "black"),
          legend.title = element_text(size=20,color="black"),
          panel.grid = element_blank(),
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 20,color = "black"),
          axis.text.x = element_text(size = 20,color = "black", angle = 0),
          axis.title.x=element_text(size = 22, hjust = 0.5, color = "black"),
          axis.title.y=element_text(size = 22, hjust = 0.5, color = "black"),
          plot.subtitle=element_text(size = 25, hjust = 0, color = "black")
    )

  p <- p +
    stat_ellipse(aes(fill = sample_type), level = 0.8, geom = 'polygon', size = 0.5, alpha = 0.2) +
    scale_fill_manual(values = c(T = '#E41A1C', NT = '#377EB8', N =  '#4DAF4A', F = '#984EA3'))

  ggsave(str_c('output/20251009_PUH_2772samples_12797proteins_', nm, '_tSNEV2_sampleType_ellipse.pdf'),
         p, height = 8, width = 11)
  
  p <- ggplot(tsne_df) +
    aes(x = X, y = Y, color = major_type) +
    stat_ellipse(aes(fill = major_type), type = 't', level = 0.9, geom = 'polygon', size = 0.5, alpha = 0.2,
                 data = tsne_df %>% filter(major_type %in% c('brain', 'liver'))) +
    geom_point(shape = "circle", size = 1.5) +
    labs(x = "t-SNE 1", y = "t-SNE 2", subtitle = "Type of ellipses: multivariate t-distribution; level = 0.9") +
    guides(color = guide_legend(title = "Tissue type", override.aes = list(size = 5)),
           fill = guide_legend(title = "Ellipses")) +
    scale_color_manual(values = major_type_labeled_colors) +
    scale_fill_manual(values = major_type_labeled_colors) +
    theme_minimal() +
    theme(legend.position = 'right',
          legend.text = element_text(size = 20, color = "black"),
          legend.title = element_text(size=20,color="black"),
          panel.grid = element_blank(),
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 20,color = "black"),
          axis.text.x = element_text(size = 20,color = "black", angle = 0),
          axis.title.x=element_text(size = 22, hjust = 0.5, color = "black"),
          axis.title.y=element_text(size = 22, hjust = 0.5, color = "black"),
          plot.subtitle=element_text(size = 8, hjust = 0, color = "black")
    )
  
  ggsave(str_c('output/20251009_PUH_2772samples_12797proteins_', nm, '_tSNEV2_someMajorType_ellipse.pdf'),
         p, height = 8, width = 11)
}


## 1.2 DR algorithms --------
# 0) Source your script (adjust path as needed)
source('drpipeline_pro_20250829.R')

# # 1) Minimal run (PCA + t-SNE + UMAP) with defaults
# res <- ge.plot(
#   data   = dat1,                 # matrix: features x samples
#   group  = df1$sample_type,      # factor/character vector, length = ncol(dat1)
#   title  = "output/20250822_PUH_2847samples_cancer_type",
#   seed   = 2023,
#   methods = c("PCA","tSNE","UMAP")
# )
# res$files   # Named vector of PDF paths

# 2) Define grids explicitly
tsne_grid_custom <- list(
  perplexity = c(5, 10, 20, 30, 40, 50),   # will be auto-capped at floor((n-1)/3)
  eta = c(200, 800, 1200),
  exaggeration_factor = c(12, 16),
  theta = c(0.5),                   # Barnes-Hut accuracy; 0.5 is a common balance
  max_iter = c(1000, 2000)
)

umap_grid_custom <- list(
  n_neighbors = c(15, 30, 50, 80),  # consider ~log(n) to ~n/20 for large n
  min_dist    = c(0.01, 0.1, 0.3),
  metric      = c("euclidean", "cosine"),
  n_epochs    = c(200, 500)
)

res.dr_list <- lapply(names(mat_list1) %>% setNames(names(mat_list1)), function(nm){
  message('A series of tSNE: matrix ', nm, '...\n')
  dat1 <- mat_list1[[nm]]
  res.grid <- ge.plot(
    data   = as.matrix(dat1),
    group  = info1$sample_type,
    label  = info1$sample_id,
    title  = str_c("output/20251009_PUH_advanced_grid_", nm),
    seed   = 2023,
    methods   = c("PCA","tSNE","UMAP"),
    n_pcs     = 50, pca_center = TRUE, pca_scale = TRUE,
    tsne_grid = tsne_grid_custom,
    umap_grid = umap_grid_custom,
    legends_mode = "all_last", unique_legend = TRUE, legends_ncol = 1,
    width = 7, height = 6
  )
  return(res.grid)
})
saveRDS(res.dr_list, 'all_samples_DR_grid_list.rds')


res2$files


# res_eval <- ge.plot(
#   data   = dat1,
#   group  = df1$sample_type,
#   title  = "output/20250822_PUH_eval",
#   seed   = 2023,
#   methods = c("PCA","tSNE","UMAP"),
#   tsne_grid = list(perplexity=c(20,40), eta=c(200,1000),
#                    exaggeration_factor=c(12), theta=c(0.5), max_iter=c(1000)),
#   umap_grid = list(n_neighbors=c(15,30,50),
#                    min_dist=c(0.01,0.1), metric=c("euclidean","cosine"), n_epochs=c(200)),
#   eval = list(
#     run = TRUE,
#     k_values = c(5,10,15,30),
#     labels = df1$sample_type,                      # or another ground-truth vector if available
#     save_csv = "output/20250822_PUH_eval_metrics.csv"
#   )
# )
# 
# # Inspect the evaluation table in memory
# head(res_eval$evaluation)



### DR tuning test ------
# cfg <- list(
#   log1p          = FALSE,
#   center         = TRUE,
#   scale_features = TRUE,
#   n_pcs          = 50,
#   tsne_grid = expand.grid(
#     perplexity   = c(5, 10, 30, 50),
#     eta          = c(200, 500, 1000),
#     exaggeration = c(12, 16),
#     theta        = c(0.5),
#     max_iter     = c(1000, 2000),
#     stringsAsFactors = FALSE
#   ),
#   outdir         = "dimred_tuning_outputs"
# )

# res <- ge.plot(
#   data   = dat1,     # features x samples
#   type   = df1$sample_type,   # sample labels
#   title  = "my_project",
#   label  = NULL,
#   seed   = 2023,
#   grid   = cfg$tsne_grid,
#   log1p  = cfg$log1p,
#   center = cfg$center,
#   scale_features = cfg$scale_features,
#   n_pcs  = cfg$n_pcs,       # if want to pre-PCA; set NULL to let Rtsne do internal PCA
#   outdir = cfg$outdir
# )
# res$params$pdf


## 1.3 subplots ----
# DR_Minimal_Extract_And_Highlight.R
# Purpose: (1) From *grid_custom* variables, select the 46th UMAP and 4th t-SNE
#          and return a single combined data.frame with columns:
#          method ("UMAP"/"t_SNE"), FileName, Dim1, Dim2.
#          (2) Build highlight-style subplots per anatomical_classification that
#          show ALL points (fixed method-wide axes); non-target classes in gray;
#          target class colored by a single GLOBAL tissue_name palette shared by UMAP & t_SNE.
#
# Notes:
# - This script does not rerun DR; it consumes ge.plot() outputs in `res$coords`.
# - `res$coords$umap` must contain columns: sample_id, UM1, UM2, n_neighbors, min_dist, metric, n_epochs
#   `res$coords$tsne` must contain columns: sample_id, Dim1, Dim2, perplexity, eta, exaggeration_factor, theta, max_iter
# - If `sample_id` is not the desired `FileName`, pass a mapping via `sample_id_to_FileName`.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  # Optional but used for broader palettes when available
  # library(ggsci); library(RColorBrewer)
})

### Section 1 — Minimal extraction from grid indices -----
.param_cols <- list(
  UMAP = c("n_neighbors","min_dist","metric","n_epochs"),
  tSNE = c("perplexity","eta","exaggeration_factor","theta","max_iter")
)

.expand_grid_from_list <- function(grid_list){
  stopifnot(is.list(grid_list))
  g <- do.call(tidyr::expand_grid, grid_list)
  g %>% dplyr::mutate(across(where(is.factor), as.character))
}

.filter_coords_by_params <- function(coords_df, method, param_row){
  pcols <- .param_cols[[method]]
  stopifnot(all(pcols %in% names(coords_df)))
  stopifnot(all(pcols %in% names(param_row)))
  out <- coords_df
  for (nm in pcols){
    v <- param_row[[nm]][[1]]
    if (is.numeric(out[[nm]])) {
      out <- dplyr::filter(out, dplyr::near(.data[[nm]], as.numeric(v), tol = 1e-12))
    } else {
      out <- dplyr::filter(out, .data[[nm]] == as.character(v))
    }
  }
  out
}

# Optional mapping: sample_id -> FileName
.map_FileName <- function(df, sample_id_to_FileName){
  if (is.null(sample_id_to_FileName)) {
    df$FileName <- df$sample_id
    return(df)
  }
  if (is.data.frame(sample_id_to_FileName)) {
    stopifnot(all(c("sample_id","FileName") %in% names(sample_id_to_FileName)))
    return(dplyr::left_join(df, sample_id_to_FileName, by = "sample_id"))
  }
  if (is.character(sample_id_to_FileName) && !is.null(names(sample_id_to_FileName))) {
    df$FileName <- unname(sample_id_to_FileName[df$sample_id])
    return(df)
  }
  stop("sample_id_to_FileName must be NULL, a named character vector, or a data.frame with sample_id & FileName")
}

extract_minimal_coords_from_grids <- function(res,
                                              tsne_grid_custom,
                                              umap_grid_custom,
                                              umap_idx = 46L,
                                              tsne_idx = 4L,
                                              sample_id_to_FileName = NULL){
  stopifnot(is.list(res), !is.null(res$coords))
  if (is.null(res$coords$umap) || is.null(res$coords$tsne))
    stop("res$coords$umap and res$coords$tsne must both be present")
  
  g_umap <- .expand_grid_from_list(umap_grid_custom)
  g_tsne <- .expand_grid_from_list(tsne_grid_custom)
  if (umap_idx < 1L || umap_idx > nrow(g_umap)) stop(sprintf("UMAP index %d out of range [1,%d]", umap_idx, nrow(g_umap)))
  if (tsne_idx < 1L || tsne_idx > nrow(g_tsne)) stop(sprintf("t-SNE index %d out of range [1,%d]", tsne_idx, nrow(g_tsne)))
  
  row_umap <- g_umap[umap_idx, , drop = FALSE]
  row_tsne <- g_tsne[tsne_idx, , drop = FALSE]
  
  sel_umap <- .filter_coords_by_params(res$coords$umap, method = "UMAP", param_row = row_umap)
  sel_tsne <- .filter_coords_by_params(res$coords$tsne, method = "tSNE", param_row = row_tsne)
  
  if (!nrow(sel_umap)) stop("No UMAP rows matched grid index 46. Ensure res was computed with the same umap_grid_custom.")
  if (!nrow(sel_tsne)) stop("No t-SNE rows matched grid index 4. Ensure res was computed with the same tsne_grid_custom.")
  
  sel_umap <- .map_FileName(sel_umap, sample_id_to_FileName) %>%
    transmute(method = "UMAP", FileName, Dim1 = UM1, Dim2 = UM2)
  sel_tsne <- .map_FileName(sel_tsne, sample_id_to_FileName) %>%
    transmute(method = "t_SNE", FileName, Dim1, Dim2)
  
  bind_rows(sel_umap, sel_tsne)
}

### Section 2 — Global high-contrast palette for tissue_name (shared across methods) ----
# Length-safe palette builder: always returns exactly one color per unique value.
# Tries wide qualitative pools (ggsci/RColorBrewer) and falls back to a 360-color wheel; 
# guarantees names(pal) aligns with unique levels.
build_distinct_palette <- function(values) {
  lv <- sort(unique(ifelse(is.na(values), "NA", as.character(values))))
  K  <- length(lv)
  
  pool <- character(0)
  
  # ggsci pools (optional)
  if (requireNamespace("ggsci", quietly = TRUE)) {
    pool <- c(
      tryCatch(ggsci::pal_igv("default")(51),     error = function(e) NULL),
      tryCatch(ggsci::pal_d3("category20")(20),    error = function(e) NULL),
      tryCatch(ggsci::pal_npg("nrc")(10),          error = function(e) NULL),
      tryCatch(ggsci::pal_uchicago("default")(9),  error = function(e) NULL),
      tryCatch(ggsci::pal_lancet()(9),              error = function(e) NULL),
      tryCatch(ggsci::pal_jama()(7),                error = function(e) NULL),
      tryCatch(ggsci::pal_tron("legacy")(7),       error = function(e) NULL),
      tryCatch(ggsci::pal_nejm()(8),                error = function(e) NULL)
    )
    pool <- unlist(pool, use.names = FALSE)
  }
  
  # RColorBrewer qualitative (optional fallback)
  if (!length(pool) && requireNamespace("RColorBrewer", quietly = TRUE)) {
    qual <- c("Dark2","Set1","Set2","Paired","Accent","Set3","Pastel1","Pastel2")
    pool <- unlist(lapply(qual, function(p) {
      nmax <- RColorBrewer::brewer.pal.info[p, "maxcolors"]
      RColorBrewer::brewer.pal(nmax, p)
    }), use.names = FALSE)
  }
  
  # Last resort: large HSV wheel
  if (!length(pool)) {
    pool <- grDevices::rainbow(360, s = 0.9, v = 0.9)
  }
  
  # Ensure at least K colors, then pick exactly K evenly across the pool
  if (length(pool) < K) pool <- rep(pool, length.out = K)
  idx <- floor(seq(1, length(pool), length.out = K))
  pal <- pool[idx]
  
  stopifnot(length(pal) == K)
  names(pal) <- lv
  pal
}

### ---- Section 3 — Highlight pages per anatomical_classification ----
# For a single method ("UMAP" or "t_SNE"):
# - show ALL points as background in gray;
# - overlay only the target anatomical_classification in color (fill=tissue_name) using a GLOBAL palette;
# - keep fixed method-wide axes across its pages.
highlight_anatomical_pages <- function(df, method_label, pal_tissue,
                                       point_size = 2.4, stroke = 0.6,
                                       bg_alpha = 0.35, fg_alpha = 0.98) {
  stopifnot(all(c("method","Dim1","Dim2","tissue_name","anatomical_classification") %in% names(df)))
  
  d0 <- df %>%
    filter(method == !!method_label) %>%
    dplyr::mutate(
      tissue_name = ifelse(is.na(tissue_name), "NA", as.character(tissue_name)),
      anatomical_classification = as.character(anatomical_classification)
    )
  if (!nrow(d0)) return(list())
  
  xr <- range(d0$Dim1, na.rm = TRUE)
  yr <- range(d0$Dim2, na.rm = TRUE)
  anat_levels <- sort(unique(d0$anatomical_classification))
  
  pages <- vector("list", length(anat_levels))
  names(pages) <- anat_levels
  
  for (i in seq_along(anat_levels)) {
    anat <- anat_levels[i]
    d0$.__is_target <- (d0$anatomical_classification == anat)
    
    p <- ggplot() +
      # Background: all points in gray
      geom_point(
        data = d0, aes(Dim1, Dim2),
        color = "grey78", alpha = bg_alpha, size = point_size
      ) +
      # Overlay: target class as FILLED circles with a dark stroke for contrast
      geom_point(
        data = d0 %>% filter(.__is_target),
        aes(Dim1, Dim2, fill = tissue_name),
        shape = 21, color = "grey15", stroke = stroke,
        alpha = fg_alpha, size = point_size + 0.8
      ) +
      scale_fill_manual(values = pal_tissue, drop = FALSE, name = "tissue_name") +
      coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
      labs(title = sprintf("%s | highlight: %s", method_label, anat),
           x = "Dim1", y = "Dim2") +
      theme_bw(base_size = 11) +
      theme(legend.position = "right", plot.title = element_text(face = "bold"))
    
    pages[[i]] <- p
  }
  pages
}

### ---- Section 4 — Export helper (optional) ----
export_highlight_pdf <- function(clu_sub, pal_tissue_global,
                                 file = "Highlight_by_AnatomicalClassification.pdf",
                                 width = 7, height = 6) {
  umap_pages <- highlight_anatomical_pages(clu_sub, "UMAP", pal_tissue_global)
  tsne_pages <- highlight_anatomical_pages(clu_sub, "t_SNE", pal_tissue_global)
  pdf(file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  for (p in c(umap_pages, tsne_pages)) print(p)
  invisible(file)
}

### ---- Section 5 — Example (commented) ----
# 1) Extract combined coordinates by selecting 46th UMAP & 4th t-SNE
df_combined <- extract_minimal_coords_from_grids(
  res = res2,
  tsne_grid_custom = tsne_grid_custom,
  umap_grid_custom = umap_grid_custom
)
rio::export(df_combined, 'Highlight_by_AnatomicalClassification.xlsx')

# 2) Join metadata (must include FileName, tissue_name, anatomical_classification, ...)
clu_sub <- df_combined %>% dplyr::inner_join(info, by = "FileName") %>% 
  dplyr::mutate(tissue_name = str_c(tissue_name, '_', sample_type))

# 3) Build a single GLOBAL palette for tissue_name (shared by UMAP & t_SNE and all subplots)
pal_tissue_global <- build_distinct_palette(clu_sub$tissue_name)

# 4) Export highlight pages to a single multi-page PDF
export_highlight_pdf(clu_sub, pal_tissue_global,
                     file = "Highlight_by_AnatomicalClassification.pdf",
                     width = 7, height = 6)

# 2.For DEA ------
# tsne_df <- rio::import('Highlight_by_AnatomicalClassification.xlsx')
nm <- 'qimp.lma.s2'
tsne_df <- rio::import('output/20251009_PUH_tsneV2_2772samples_sample_type_qimp.lma.s2.xlsx')
df <- mat_list1[[nm]]

# remove brain and 13 plasma samples
# tsne_df <- rio::import('20230710_PUH_tsne_1731samples_cancer_type_v2.xlsx')
# tsne_df[tsne_df$sample_type == 'carcinoma', 'sample_type'] <- 'T'
# tsne_df[tsne_df$sample_type == 'adjacent', 'sample_type'] <- 'NT'
# tsne_df$sample_type %<>% factor(levels = c('T', 'NT', 'N', 'F'))
# tsne_df %<>% rename(X = `1`, Y = `2`)

rm_index <- which(colnames(df) %in% (tsne_df %>% filter(anatomical_classification %in% c('brain', 'plasma')) %>% pull(FileName)))
mat <- 2 ^ df[, -rm_index] %>% as.matrix()
class(mat) # "matrix" "array"
dim(mat) # 12797  2664

# 按 F/C/NC/N 分别控制 NA ratio
f_index <- which(colnames(mat) %in% (tsne_df %>% filter(sample_type == 'F') %>% pull(FileName)))
c_index <- which(colnames(mat) %in% (tsne_df %>% filter(sample_type == 'T') %>% pull(FileName)))
nc_index <- which(colnames(mat) %in% (tsne_df %>% filter(sample_type == 'NT') %>% pull(FileName)))
n_index <- which(colnames(mat) %in% (tsne_df %>% filter(sample_type == 'N') %>% pull(FileName)))

condition1 <- apply(mat[, f_index], 1, function(x) any(!is.na(x))) # any F not NA
condition2 <- apply(mat[, c_index], 1, function(x) any(!is.na(x))) # any T not NA
condition3 <- apply(mat[, nc_index], 1, function(x) any(!is.na(x))) # any NT not NA
condition4 <- apply(mat[, n_index], 1, function(x) any(!is.na(x))) # any N not NA
mat <- mat[condition1 & condition2 & condition3 & condition4, ]
dim(mat) # 12797  2664

# Wilcox test
df_compare_mean <- as.data.frame(t(mat)) %>% rownames_to_column('FileName')
df_compare_mean %<>% add_column(Group = NA, .before = 2)
df_compare_mean[f_index, 'Group'] <- 'F'
df_compare_mean[c_index, 'Group'] <- 'T'
df_compare_mean[nc_index, 'Group'] <- 'NT'
df_compare_mean[n_index, 'Group'] <- 'N'

df_p_ls <- list()
for(j in 3:ncol(df_compare_mean)){
  cat(j, '...\r')
  dfsub <- df_compare_mean %>% dplyr::select(2, all_of(j))
  prot_name <- colnames(dfsub)[2]
  colnames(dfsub)[2] <- 'Protein'
  my_comparisons <- list(c('F', 'T'), c('T', 'NT'), c('NT', 'N'))
  df_p <- plyr::ldply(my_comparisons, function(comp){
    x1 <- na.omit(dfsub %>% filter(Group == comp[1]) %>% pull(Protein))
    x2 <- na.omit(dfsub %>% filter(Group == comp[2]) %>% pull(Protein))
    p <- wilcox.test(x1, x2, paired = F)$p.value
    data.frame(Protein = prot_name, Group1 = comp[1], Group2 = comp[2], p = p) %>% return()
  })
  df_p$p.adj <- p.adjust(df_p$p, method = 'BH')
  df_p_ls[[j-2]] <- df_p
  # ggpubr::compare_means(formula = Protein ~ Group, data = dfsub, method = 'wilcox.test', p.adjust.method = 'BH', comparisons = my_comparisons)
}
names(df_p_ls) <- colnames(df_compare_mean)[-(1:2)]
df_test <- plyr::ldply(df_p_ls)
df_test %<>% dplyr::mutate(Protein = .id, .id = NULL)

df_test$Group1 <- sapply(df_test$Group1, function(x) switch(x, F = 'f.', T = 'c.', NT = 'nc.'))
df_test$Group2 <- sapply(df_test$Group2, function(x) switch(x, T = 'c', NT = 'nc', N = 'n'))
df_test$Pair_P <- str_c('P_', df_test$Group1, df_test$Group2)
df_test$Pair_adjP <- str_c('adjP_', df_test$Group1, df_test$Group2)
df_test %<>% dplyr::select(-Group1, -Group2)
df_test1 <- df_test %>% pivot_wider(id_cols = 'Protein', names_from = 'Pair_P', values_from = 'p')
df_test2 <- df_test %>% pivot_wider(id_cols = 'Protein', names_from = 'Pair_adjP', values_from = 'p.adj')
df_wilcox_test <- df_test1 %>% left_join(df_test2, by = 'Protein')

log2fc_f.c <- log2fc_c.nc <- log2fc_nc.n <- c()
for(i in 1:nrow(mat)){
  cat('Processing protein ', i, '...\r')
  x1 <- na.omit(mat[i, f_index])
  x2 <- na.omit(mat[i, c_index])
  x3 <- na.omit(mat[i, nc_index])
  x4 <- na.omit(mat[i, n_index])
  
  log2fc_f.c[i] <- log2(mean(x1) / mean(x2))
  log2fc_c.nc[i] <- log2(mean(x2) / mean(x3))
  log2fc_nc.n[i] <- log2(mean(x3) / mean(x4))
}
df_wilcox_test$log2FC_f.c <- log2fc_f.c
df_wilcox_test$log2FC_c.nc <- log2fc_c.nc
df_wilcox_test$log2FC_nc.n <- log2fc_nc.n

p_significant <- df_wilcox_test %>% dplyr::select(contains('adjP')) %>% apply(1, function(x) all(x < 0.05)) # wheather adjp < 0.05
fc_significant <- df_wilcox_test %>% dplyr::select(contains('log2FC')) %>% apply(1, function(x) {(all(abs(x) >= log2(1.2))) & (length(unique(sign(x))) == 1)}) # wheather |FC| >= 1.5 and with the same signature
df_wilcox_test$isSignificant <- ifelse(p_significant & fc_significant, T, F)
df_wilcox_test$Regulation <- ifelse(df_wilcox_test$log2FC_f.c > 0, 'Up', 'Down')
df_wilcox_test[!df_wilcox_test$isSignificant, 'Regulation'] <- NA # cannot describe insignificance
table(df_wilcox_test$isSignificant)
# FALSE  TRUE 
# 11491  1306 

rio::export(df_wilcox_test, 'PUH_allSample_differentiation_rmBrainPlasma_20251201.xlsx')

# volcano plot
df_dea <- df_wilcox_test
df_dea %<>% dplyr::mutate(
  isSignificant_f.c = abs(log2FC_f.c) >= log2(1.2) & adjP_f.c < 0.05,
  isSignificant_c.nc = abs(log2FC_c.nc) >= log2(1.2) & adjP_c.nc < 0.05,
  isSignificant_nc.n = abs(log2FC_nc.n) >= log2(1.2) & adjP_nc.n < 0.05
)

p1 <- ggplot() +
  geom_point(data = df_dea %>% filter(!isSignificant_f.c), aes(x = log2FC_f.c, y = -log10(adjP_f.c), size = abs(log2FC_f.c)), color = 'gray', alpha = 0.5) +
  geom_point(data = df_dea %>% filter(isSignificant_f.c & log2FC_f.c > 0), aes(x = log2FC_f.c, y = -log10(adjP_f.c), size = abs(log2FC_f.c)), color = '#BC3C29', alpha = 0.8) +
  geom_point(data = df_dea %>% filter(isSignificant_f.c & log2FC_f.c < 0), aes(x = log2FC_f.c, y = -log10(adjP_f.c), size = abs(log2FC_f.c)), color = '#0072B5', alpha = 0.8) +
  scale_size(name = '|log2 (Fold change)|', limits = c(0, 10), breaks = c(2, 4, 6), labels = c(2, 4, 6), range = c(2, 6)) +
  geom_hline(yintercept = -log10(0.05), lty = 4, lwd = 0.6, alpha = 0.8) +
  geom_vline(xintercept = log2(1.2) * c(1, -1),lty = 4,lwd = 0.6, alpha = 0.8) +
  theme_classic() +
  theme(legend.position = 'right',
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.45),
  )+
  labs(x = "log2 (Fold change)", y = "-log10 (Adj. P value)") +
  guides(size = guide_legend(title = '|log2 (Fold change)|'),) +
  annotate("text", x = -4, y = 20, label = str_c("Down-regulated:\nn = ", sum(df_dea$isSignificant_f.c & df_dea$log2FC_f.c < 0)), size = 5) +
  annotate("text", x = 4, y = 20, label = str_c("Up-regulated:\nn = ", sum(df_dea$isSignificant_f.c & df_dea$log2FC_f.c > 0)), size = 5)

ggsave('PUH_allSample_differentiation_rmBrainPlasma_volcano_fetal-tumor_20251201.pdf', p1, width = 8, height = 7)
ggsave('PUH_allSample_differentiation_rmBrainPlasma_volcano_fetal-tumor_wolegend_20251201.pdf', p1+theme(legend.position = 'none'), width = 7, height = 7)



p2 <- ggplot() +
  geom_point(data = df_dea %>% filter(!isSignificant_c.nc), aes(x = log2FC_c.nc, y = -log10(adjP_c.nc), size = abs(log2FC_c.nc)), color = 'gray', alpha = 0.5) +
  geom_point(data = df_dea %>% filter(isSignificant_c.nc & log2FC_c.nc > 0), aes(x = log2FC_c.nc, y = -log10(adjP_c.nc), size = abs(log2FC_c.nc)), color = '#BC3C29', alpha = 0.8) +
  geom_point(data = df_dea %>% filter(isSignificant_c.nc & log2FC_c.nc < 0), aes(x = log2FC_c.nc, y = -log10(adjP_c.nc), size = abs(log2FC_c.nc)), color = '#0072B5', alpha = 0.8) +
  scale_size(name = '|log2 (Fold change)|', limits = c(0, 10), breaks = c(2, 4, 6), labels = c(2, 4, 6), range = c(2, 6)) +
  geom_hline(yintercept = -log10(0.05), lty = 4, lwd = 0.6, alpha = 0.8) +
  geom_vline(xintercept = log2(1.2) * c(1, -1),lty = 4,lwd = 0.6, alpha = 0.8) +
  theme_classic() +
  theme(legend.position = 'right',
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.45),
  )+
  labs(x = "log2 (Fold change)", y = "-log10 (Adj. P value)") +
  guides(size = guide_legend(title = '|log2 (Fold change)|'),) +
  annotate("text", x = -4, y = 20, label = str_c("Down-regulated:\nn = ", sum(df_dea$isSignificant_c.nc & df_dea$log2FC_c.nc < 0)), size = 5) +
  annotate("text", x = 4, y = 20, label = str_c("Up-regulated:\nn = ", sum(df_dea$isSignificant_c.nc & df_dea$log2FC_c.nc > 0)), size = 5)

ggsave('PUH_allSample_differentiation_rmBrainPlasma_volcano_tumor-nontumor_20251201.pdf', p2, width = 8, height = 7)
ggsave('PUH_allSample_differentiation_rmBrainPlasma_volcano_tumor-nontumor_wolegend_20251201.pdf', p2+theme(legend.position = 'none'), width = 7, height = 7)



p3 <- ggplot() +
  geom_point(data = df_dea %>% filter(!isSignificant_nc.n), aes(x = log2FC_nc.n, y = -log10(adjP_nc.n), size = abs(log2FC_nc.n)), color = 'gray', alpha = 0.5) +
  geom_point(data = df_dea %>% filter(isSignificant_nc.n & log2FC_nc.n > 0), aes(x = log2FC_nc.n, y = -log10(adjP_nc.n), size = abs(log2FC_nc.n)), color = '#BC3C29', alpha = 0.8) +
  geom_point(data = df_dea %>% filter(isSignificant_nc.n & log2FC_nc.n < 0), aes(x = log2FC_nc.n, y = -log10(adjP_nc.n), size = abs(log2FC_nc.n)), color = '#0072B5', alpha = 0.8) +
  scale_size(name = '|log2 (Fold change)|', limits = c(0, 10), breaks = c(2, 4, 6), labels = c(2, 4, 6), range = c(2, 6)) +
  geom_hline(yintercept = -log10(0.05), lty = 4, lwd = 0.6, alpha = 0.8) +
  geom_vline(xintercept = log2(1.2) * c(1, -1),lty = 4,lwd = 0.6, alpha = 0.8) +
  theme_classic() +
  theme(legend.position = 'right',
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.45),
  )+
  labs(x = "log2 (Fold change)", y = "-log10 (Adj. P value)") +
  guides(size = guide_legend(title = '|log2 (Fold change)|'),) +
  annotate("text", x = -4, y = 20, label = str_c("Down-regulated:\nn = ", sum(df_dea$isSignificant_nc.n & df_dea$log2FC_nc.n < 0)), size = 5) +
  annotate("text", x = 4, y = 20, label = str_c("Up-regulated:\nn = ", sum(df_dea$isSignificant_nc.n & df_dea$log2FC_nc.n > 0)), size = 5)

ggsave('PUH_allSample_differentiation_rmBrainPlasma_volcano_nontumor-N_20251201.pdf', p3, width = 8, height = 7)
ggsave('PUH_allSample_differentiation_rmBrainPlasma_volcano_nontumor-N_wolegend_20251201.pdf', p3+theme(legend.position = 'none'), width = 7, height = 7)

df_dea %>% rio::export('PUH_allSample_differentiation_rmBrainPlasma_addLabels_20251201.xlsx')

up_f.c <- df_dea$Protein[df_dea$isSignificant_f.c & df_dea$log2FC_f.c > 0]
up_c.nc <- df_dea$Protein[df_dea$isSignificant_c.nc & df_dea$log2FC_c.nc > 0]
up_nc.n <- df_dea$Protein[df_dea$isSignificant_nc.n & df_dea$log2FC_nc.n > 0]
down_f.c <- df_dea$Protein[df_dea$isSignificant_f.c & df_dea$log2FC_f.c < 0]
down_c.nc <- df_dea$Protein[df_dea$isSignificant_c.nc & df_dea$log2FC_c.nc < 0]
down_nc.n <- df_dea$Protein[df_dea$isSignificant_nc.n & df_dea$log2FC_nc.n < 0]
up_list <- list(`F vs. T` = up_f.c,
                `T vs. NT` = up_c.nc,
                `NT vs. N` = up_nc.n)
down_list <- list(`F vs. T` = down_f.c,
                  `T vs. NT` = down_c.nc,
                  `NT vs. N` = down_nc.n)



p1 <- VennDiagram::venn.diagram(x = up_list,
                                resolution = 300,
                                col='#A81B3C', 
                                fill='#A81B3C',alpha=0.8,
                                # fill=c("#EB495C", "#2DCC62FD", "#3C65D6"), 
                                main="Upregulation",
                                #sub = rep,
                                main.cex = 4, 
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff", 
                                filename = NULL,disable.logging = T,
                                width = 5000, height = 5000
)
p2 <- VennDiagram::venn.diagram(x = down_list,
                                resolution = 300,
                                col='#19199E', 
                                fill='#19199E',alpha=0.8,
                                # fill=c("#EB495C", "#2DCC62FD", "#3C65D6"), 
                                main="Downregulation",
                                #sub = rep,
                                main.cex = 4, 
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff", 
                                filename = NULL,disable.logging = T,
                                width = 5000, height = 5000
)
pdf('PUH_allSample_differentiation_rmBrainPlasma_VennDiagram_20251201.pdf', width = 20, height = 20)
grid::grid.newpage(); grid::grid.draw(p1)
grid::grid.newpage(); grid::grid.draw(p2)
graphics.off()


## 2.2 heatmap --------------------------------------------------------------
df_wilcox_test <- rio::import('PUH_allSample_differentiation_rmBrainPlasma_20251201.xlsx')

prot_diff_up <- df_wilcox_test %>% filter(isSignificant) %>% filter(Regulation == 'Up') %>% pull(Protein) # 1244
prot_diff_down <- df_wilcox_test %>% filter(isSignificant) %>% filter(Regulation == 'Down') %>% pull(Protein) # 62

file_info <- info1 %>%
  filter(anatomical_classification != 'plasma') %>% 
  dplyr::select(FileName, tissue_name, sample_type, anatomical_classification) %>% 
  rename(organ = anatomical_classification) %>% 
  dplyr::mutate(sample_type = sapply(sample_type, function(x) {switch(x, adjacent = 'NT', carcinoma = 'T', x)}))
# file_info[file_info$tissue_name %in% c('fetal telencephalon', 'fetal mesencephalon', 'fetal myelencephalon', 'fetal metencephalon'), 'organ'] <- 'brain' # set fetal brain to organ brain

# file_info[which(file_info$organ == 'brain'), 'sample_type'] <- 'Brain'
rownames(file_info) <- file_info$FileName
file_info %<>% dplyr::select(-FileName) %>%
  dplyr::mutate(sample_type = factor(sample_type, levels = c('F', 'T', 'NT', 'N'), ordered = T)) %>% 
  arrange(sample_type, organ)

file_info1 <- file_info %>% slice(-which(organ == 'brain'))
file_info2 <- file_info %>% slice(which(organ == 'brain'))
file_info <- rbind(file_info1, file_info2)

heatmat <- cbind(
  mat[c(prot_diff_up, prot_diff_down), ],
  df[c(prot_diff_up, prot_diff_down), rownames(file_info %>% filter(organ == 'brain', sample_type == 'F'))], # brain_F
  df[c(prot_diff_up, prot_diff_down), rownames(file_info %>% filter(organ == 'brain', sample_type == 'T'))], # brain_C
  df[c(prot_diff_up, prot_diff_down), rownames(file_info %>% filter(organ == 'brain', sample_type == 'NT'))], # brain_NC
  df[c(prot_diff_up, prot_diff_down), rownames(file_info %>% filter(organ == 'brain', sample_type == 'N'))] # brain_N
)
dim(heatmat) # 1306 2772


heatmat1 <- heatmat[apply(heatmat[, f_index], 1, function(x) sum(is.na(x)) / length(x)) <= 0.5, ] # missing ratio <= 50% in fetal
heatmat2 <- heatmat1[apply(heatmat1[, c_index], 1, function(x) sum(is.na(x)) / length(x)) <= 0.5, ] # missing ratio <= 50% in cancerous
heatmat3 <- heatmat2[apply(heatmat2[, nc_index], 1, function(x) sum(is.na(x)) / length(x)) <= 0.5, ] # missing ratio <= 50% in NT
heatmat4 <- heatmat3[apply(heatmat3[, n_index], 1, function(x) sum(is.na(x)) / length(x)) <= 0.5, ] # missing ratio <= 50% in N
heatmat4 <- log2(heatmat4)
mat_scale <- apply(heatmat4, 1, scale) %>% t() %>% data.frame()
quantile(mat_scale, na.rm = T)
colnames(mat_scale) <- colnames(heatmat)

#breaks
my_breaks <- seq(-4, 4, by = 0.01)

#colors
my_colors <- c(
  colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[1:3]))(length(my_breaks)/2/8*6),
  colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[3:5], 'white'))(length(my_breaks)/2/8*2),
  colorRampPalette(colors = c('white', RColorBrewer::brewer.pal(11, "PiYG")[7:9]))(length(my_breaks)/2/8*2),
  colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[9:11]))(length(my_breaks)/2/8*6)
)
my_colors %<>% rev() # high intensity should be hot while low to be cold




# file_info %>% dplyr::count(organ) # 80 organs
heat <- mat_scale[, rownames(file_info)]

df_color <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/input/PUH_tissue_colorset_20230210.xlsx')
organ_colors <- str_c('#', df_color$color)
set.seed(2023)
organ_colors %<>% append(sample(organ_colors, length(unique(file_info$organ)) - length(organ_colors)))
names(organ_colors) <- sort(unique(file_info$organ))
sample_type_colors <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#F39800') %>% setNames(c('T', 'NT', 'N', 'F', 'Brain'))
ann_colors <- list(sample_type = sample_type_colors, organ = organ_colors)
a <- pheatmap::pheatmap(
  heat,
  color = my_colors,
  legend_breaks = seq(-4,4,2),
  breaks = my_breaks,
  fontsize_col = 8,
  border_color = F,
  annotation_col = file_info,
  # annotation_row = ann_row,
  gaps_col = c(which(!duplicated(file_info1$sample_type))[-1], (which(!duplicated(file_info2$sample_type)) + nrow(file_info1))) - 1,
  # gaps_col = cumsum(table(ann_col$tissue))[-length(cumsum(table(ann_col$tissue)))],
  # gaps_row = cumsum(table(row.tissue))[-length(cumsum(table(row.tissue)))],
  fontsize_row=10,
  annotation_colors = ann_colors,
  na_col = "grey70",#scale = "row",
  cluster_rows = T, cluster_cols = F,
  show_rownames = T, show_colnames = F,
  filename = 'PUH_allSample_differentiation_rmBrainPlasma_FCNcN_DEP_50NA_addBrain_addGlio_heatmap_20251201.pdf',
  width = 15, height = 17
)

list(matrix = heat[a$tree_row$order, ], label = file_info) %>%
  rio::export('PUH_allSample_differentiation_rmBrainPlasma_FCNcN_DEP_50NA_addBrain_addGlio_heatmap_20251201.xlsx')




## get median matrix
df_lbl <- file_info %>% rownames_to_column('FileName')
df_lbl[df_lbl$organ != 'brain', 'organ'] <- 'notBrain'
names(df_lbl$sample_type) <- NULL

filename_split <- plyr::dlply(df_lbl, c('sample_type', 'organ'), function(dfsub){
  dfsub$FileName
})


mat <- df[rownames(mat), ]
dim(mat) # 12797  2772

mat_median <- sapply(filename_split, function(files){
  apply(mat[, files], 1, function(x) median(x, na.rm = T))
})
mat_median <- mat_median[, c('F.notBrain', 'T.notBrain', 'NT.notBrain', 'N.notBrain', 'T.brain', 'NT.brain', 'N.brain')] # columns arrange
mat_median_scale <- apply(mat_median, 1, scale) %>% t() %>% data.frame()
quantile(mat_median_scale, na.rm = T)
colnames(mat_median_scale) <- colnames(mat_median)

#breaks
my_breaks <- seq(-2, 2, by = 0.01)

#colors
my_colors <- c(
  colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[1:3]))(length(my_breaks)/2/8*6),
  colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[3:5], 'white'))(length(my_breaks)/2/8*2),
  colorRampPalette(colors = c('white', RColorBrewer::brewer.pal(11, "PiYG")[7:9]))(length(my_breaks)/2/8*2),
  colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[9:11]))(length(my_breaks)/2/8*6)
)
my_colors %<>% rev() # high intensity should be hot while low to be cold


# file_info %>% dplyr::count(organ) # 80 organs
heat <- mat_median_scale[c(prot_diff_up, prot_diff_down), ]
ann_col <- data.frame(Sample_type = c('F', 'T', 'NT', 'N', 'T', 'NT', 'N'),
                      row.names = colnames(heat))

df_color <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/input/PUH_tissue_colorset_20230210.xlsx')
organ_colors <- str_c('#', df_color$color)
set.seed(2023)
organ_colors %<>% append(sample(organ_colors, length(unique(file_info$organ)) - length(organ_colors)))
names(organ_colors) <- sort(unique(file_info$organ))
Sample_type_colors <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#F39800') %>% setNames(c('T', 'NT', 'N', 'F', 'Brain'))
ann_colors <- list(Sample_type = Sample_type_colors)


b <- pheatmap::pheatmap(
  heat,
  color = my_colors,
  legend_breaks = seq(-2,2,2),
  breaks = my_breaks,
  fontsize_col = 6,
  border_color = F,
  annotation_col = ann_col,
  # annotation_row = ann_row,
  gaps_col = 4,
  # gaps_col = cumsum(table(ann_col$tissue))[-length(cumsum(table(ann_col$tissue)))],
  # gaps_row = cumsum(table(row.tissue))[-length(cumsum(table(row.tissue)))],
  fontsize_row=10,
  annotation_colors = ann_colors,
  na_col = "grey70",#scale = "row",
  cluster_rows = T, cluster_cols = F,
  show_rownames = T, show_colnames = F,
  filename = 'PUH_allSample_differentiation_rmBrainPlasma_FCNcN_DEP_50NA_addBrain_addGlio_median_heatmap_20251201.pdf',
  width = 7, height = 20
)

b.same <- pheatmap::pheatmap(
  heat[b$tree_row$order, ],
  color = my_colors,
  legend_breaks = seq(-2,2,2),
  breaks = my_breaks,
  fontsize_col = 6,
  border_color = F,
  annotation_col = ann_col,
  # annotation_row = ann_row,
  gaps_col = 4,
  # gaps_col = cumsum(table(ann_col$tissue))[-length(cumsum(table(ann_col$tissue)))],
  # gaps_row = cumsum(table(row.tissue))[-length(cumsum(table(row.tissue)))],
  fontsize_row=10,
  annotation_colors = ann_colors,
  na_col = "grey70",#scale = "row",
  cluster_rows = F, cluster_cols = F,
  show_rownames = T, show_colnames = F,
  filename = 'PUH_allSample_differentiation_rmBrainPlasma_FCNcN_DEP_50NA_addBrain_addGlio_median_heatmap_20251201.pdf',
  width = 5, height = 25
)


protinfo <- heat[b$tree_row$order, ] %>%
  rownames_to_column('Protein') %>% 
  left_join(dfprot %>% dplyr::select(`Protein ID`, Gene, `Entry Name`, `Protein Description`), by = c(Protein = 'Protein ID'))%>%
  left_join(df_wilcox_test)


list(matrix = protinfo, label = ann_col) %>% rio::export('PUH_allSample_differentiation_rmBrainPlasma_FCNcN_DEP_50NA_addBrain_addGlio_median_heatmap_20251201.xlsx')



## 2.3 GO enrichment --------------------------------------------------------
# save.image('all_sample_tsne.RData')
df_wilcox_test <- rio::import('PUH_allSample_differentiation_rmBrainPlasma_20251201.xlsx')
libprot <- rio::import('../../TPL/libs/20220616_fraglib/protein.tsv')
prot2gene <- libprot$Gene
names(prot2gene) <- libprot$`Protein ID`

prot_diff_up <- df_wilcox_test %>% filter(isSignificant) %>% filter(Regulation == 'Up') %>% pull(Protein) # 1244
prot_diff_down <- df_wilcox_test %>% filter(isSignificant) %>% filter(Regulation == 'Down') %>% pull(Protein) # 62

library(clusterProfiler)
library(org.Hs.eg.db)

#BP
goall_BP_up <- enrichGO(prot_diff_up, OrgDb = org.Hs.eg.db,
                        ont = "BP", keyType = "UNIPROT", pvalueCutoff = 0.05,
                        readable = T, minGSSize = 1, maxGSSize = 30000)
goall_BP_up2 <- clusterProfiler::simplify(goall_BP_up, cutoff=0.7, by="p.adjust", select_fun = min)
df_GOBP_up1 <- DOSE::gsfilter(goall_BP_up, min = 10, max = 5000)[]
df_GOBP_up2 <- DOSE::gsfilter(goall_BP_up2, min = 10, max = 5000)[]

goall_BP_down <- enrichGO(prot_diff_down, OrgDb = org.Hs.eg.db,
                          ont = "BP", keyType = "UNIPROT", pvalueCutoff = 0.05,
                          readable = T, minGSSize = 1, maxGSSize = 30000)
goall_BP_down2 <- clusterProfiler::simplify(goall_BP_down, cutoff=0.7, by="p.adjust", select_fun = min)
df_GOBP_down1 <- DOSE::gsfilter(goall_BP_down, min = 10, max = 5000)[]
df_GOBP_down2 <- DOSE::gsfilter(goall_BP_down2, min = 10, max = 5000)[]

df_GOBP <- plyr::ldply(list(Up = df_GOBP_up2, Down = df_GOBP_down2), .id = 'Regulation')
df_GOBP$pvalue <- -log10(df_GOBP$pvalue)
df_GOBP$pvalue <- ifelse(df_GOBP$Regulation == 'Up', df_GOBP$pvalue, -df_GOBP$pvalue)
df_GOBP %<>% arrange(desc(pvalue))
# df_GOBP$Count <- as.factor(df_GOBP$Count)
# df_GOBP$Description %<>% str_to_sentence()

# library(ggcharts)
quantile(df_GOBP$pvalue)
# 0%       25%       50%       75%      100% 
# 2.565715  3.272864  4.992356  9.480528 29.822978 
GO_BP <- ggplot(df_GOBP,aes(pvalue,Description))+
  geom_point(aes(size=Count,color=Regulation))+
  labs(x="-log10 (p value)",y="Pathway name",title="GO Enrichment Biological Process")+
  scale_x_continuous(breaks=seq(-18, 8, 2), labels = abs(seq(-18, 8, 2))) +
  scale_y_discrete(limits=rev(df_GOBP$Description))+
  scale_color_manual(values=c(Up = "#BC3C29", Down = "#293CBC"))+
  # scale_size_manual(values = c('9' = 5, '11' = 7)) +
  theme_classic()+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black"))+
  theme(legend.text=element_text(size=15),legend.title = element_text(size=15))+
  theme(plot.title = element_text(size = 15, face = "bold"))
ggsave('PUH_allSample_differentiation_rmBrainPlasma_GOBP_v2_20251201.pdf', GO_BP, width = 20, height = 30)
rio::export(df_GOBP, 'PUH_allSample_differentiation_rmBrainPlasma_GOBP_v2_20251201.xlsx')

GO_BP <- ggplot((df_GOBP %>% dplyr::slice(1:5, (nrow(df_GOBP)-4):nrow(df_GOBP))),aes(pvalue,Description))+
  geom_point(aes(size=Count,color=Regulation))+
  labs(x="-log10 (p value)",y="Pathway name",title="GO Enrichment Biological Process")+
  # scale_x_continuous(breaks=seq(-10, 14, 4), labels = abs(seq(-10, 14, 4))) +
  scale_y_discrete(limits=rev(df_GOBP %>% dplyr::slice(1:5, (nrow(df_GOBP)-4):nrow(df_GOBP)) %>% pull(Description)))+
  scale_color_manual(values=c(Up = "#BC3C29", Down = "#293CBC"))+
  # scale_size_manual(values = c('9' = 5, '11' = 7)) +
  theme_classic()+
  theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black"))+
  theme(legend.text=element_text(size=15),legend.title = element_text(size=15))+
  theme(plot.title = element_text(size = 15, face = "bold"))
ggsave('PUH_allSample_differentiation_rmBrainPlasma_GOBP_top20UpDown_v2_20251201.pdf', GO_BP, width = 15, height = 5)

# save(GO_BP, df_GOBP, df_GOBP_up1, df_GOBP_up2, df_GOBP_down1, df_GOBP_down2, goall_BP_up, goall_BP_up2, goall_BP_down, goall_BP_down2, file = 'all_sample_tsne_GOBP.RData')
# load('all_sample_tsne_GOBP.RData')


### boxplot ----
prot_selected <- c('P08294', 'P11216', 'P01034')
df_box <- heatmat[prot_selected, ] %>% t() %>% as.data.frame() %>% 
  rownames_to_column('FileName') %>% 
  left_join(tsne_df %>% dplyr::select(FileName, sample_type) %>% distinct()) %>% 
  pivot_longer(cols = prot_selected, names_to = 'Protein', values_to = 'Intensity')

df_box$sample_type %<>% as.character()
# df_box <- df_compare_mean %>% dplyr::select(1:2, all_of(prot_selected)) %>% 
#   pivot_longer(cols = prot_selected, names_to = 'Protein', values_to = 'Intensity')
# df_box[df_box$sample_type == 'T', 'sample_type'] <- 'T'
# df_box[df_box$sample_type == 'NT', 'sample_type'] <- 'NT'

brains <- tsne_df %>% filter(anatomical_classification %in% c('brain')) %>% pull(FileName)
df_box$sample_type[(df_box$FileName %in% brains)] %<>% str_c('.b')
df_box$sample_type %<>% factor(levels = c('F', 'T', 'NT', 'N', 'T.b', 'NT.b', 'N.b'), ordered = T)
df_box$color <- str_replace(df_box$sample_type, '.b', '')
df_box$fill <- ifelse(str_detect(df_box$sample_type, '.b$'), 'Brain', 'Others')

df_box_log2 <- df_box
df_box_log2$Intensity %<>% log2()

sample_type_colors <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#AAAAAA', 'white') %>% setNames(c('T', 'NT', 'N', 'F', 'Brain', 'Others'))

plot_ls <- plyr::dlply(df_box_log2, 'Protein', function(dfsub){
  p <- ggplot(dfsub) +
    aes(x = sample_type, y = Intensity) +
    geom_jitter(alpha = 0.4, size = 0.1) +
    geom_boxplot(aes(color = color, fill = fill), outlier.shape = NA, width = 0.4) +
    scale_color_manual(values = sample_type_colors) +
    scale_fill_manual(values = sample_type_colors) +
    labs(
      x = '', y = 'log2 Intensity', title = dfsub$Protein[1]
    )+
    theme(#panel.grid = element_blank(),
      #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text = element_text(size = 12,color = 'black'),
      axis.line = element_line(color = 'black'),
      axis.line.x = element_blank(),
      plot.subtitle = element_text(size = 30, hjust = 0, color = 'black')
    )+
    theme(legend.text = element_text(size = 12, color = 'black'), legend.position = 'right',
          legend.title = element_text(size = 15, color = 'black'))+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, vjust = 0.5 ,color = 'black', angle = 0),
          #axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    )+
    theme(axis.title.y = element_text(size = 12, hjust = 0.5, color = 'black', angle = 90),
          axis.text.y = element_text(size = 10,color = 'black', angle = 0),
          # axis.ticks.y = element_blank()
    )
  
  p + ggpubr::geom_signif	(#aes(x = sample_type, y = Intensity),
    # data = dfsub,
    comparisons = list(c('F', 'T'),
                       c('T', 'NT'),
                       c('NT', 'N'),
                       c('T.b', 'NT.b'),
                       c('NT.b', 'N.b')
    ),
    step_increase = 0.05,
    #map_signif_level = T,
    # test = function(x1, x2) t.test(x1, x2, alternative = 'less')) # 单尾
    test = wilcox.test) # 双尾
})

p <- ggpubr::ggarrange(plotlist = plot_ls, ncol = 1)
ggsave('PUH_P08294_P11216_P01034_expr_20251201.pdf', p, width = 10, height = 10)

df_p <- plyr::ddply(df_box_log2, 'Protein', function(dfsub){
  # dfsub <- df_box_log2 %>% filter(Protein == 'P01034')
  my_comparisons <- list(c('F', 'T'), c('T', 'NT'), c('NT', 'N'), c('T.b', 'NT.b'), c('NT.b', 'N.b'))
  
  df_p <- plyr::ldply(my_comparisons, function(comp){
    x1 <- na.omit(dfsub %>% filter(sample_type == comp[1]) %>% pull(Intensity))
    x2 <- na.omit(dfsub %>% filter(sample_type == comp[2]) %>% pull(Intensity))
    p <- wilcox.test(x1, x2, paired = F, var.equal = F)$p.value
    data.frame(Protein = dfsub$Protein[1], Group1 = comp[1], Group2 = comp[2], p = p) %>% return()
  })
  df_p$p.adj <- p.adjust(df_p$p, method = 'BH')
  return(df_p)
})
rio::export(df_p, 'PUH_P08294_P11216_P01034_expr_20251201.xlsx')


### heatmap -----
# 126/222 DEPs with missing ratio <= 50% in fetal, tumor, non-tumor, and normal
# (95/187 down + 31/35 up)
prot_diff_up <- df_wilcox_test %>% filter(isSignificant) %>% filter(Regulation == 'Up') %>% pull(Protein) # 35
prot_diff_down <- df_wilcox_test %>% filter(isSignificant) %>% filter(Regulation == 'Down') %>% pull(Protein) # 187


pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)

file_info <- df1 %>%
  filter(anatomical_classification != 'plasma') %>% 
  dplyr::select(FileName, tissue_name, sample_type, anatomical_classification) %>% 
  rename(organ = anatomical_classification) %>% 
  dplyr::mutate(sample_type = sapply(sample_type, function(x) {switch(x, adjacent = 'NT', carcinoma = 'T', x)}))
# file_info[file_info$tissue_name %in% c('fetal telencephalon', 'fetal mesencephalon', 'fetal myelencephalon', 'fetal metencephalon'), 'organ'] <- 'brain' # set fetal brain to organ brain

# file_info[which(file_info$organ == 'brain'), 'sample_type'] <- 'Brain'
rownames(file_info) <- file_info$FileName
file_info %<>% dplyr::select(-FileName) %>%
  dplyr::mutate(sample_type = factor(sample_type, levels = c('F', 'T', 'NT', 'N'), ordered = T)) %>% 
  arrange(sample_type, organ)

file_info1 <- file_info %>% slice(-which(organ == 'brain'))
file_info2 <- file_info %>% slice(which(organ == 'brain'))
file_info <- rbind(file_info1, file_info2)

heatmat <- cbind(
  mat[c(prot_diff_up, prot_diff_down), ],
  # df[c(prot_diff_up, prot_diff_down), rownames(file_info %>% filter(organ == 'brain', sample_type == 'F'))], # brain_F
  df[c(prot_diff_up, prot_diff_down), rownames(file_info %>% filter(organ == 'brain', sample_type == 'T'))], # brain_C
  df[c(prot_diff_up, prot_diff_down), rownames(file_info %>% filter(organ == 'brain', sample_type == 'NT'))], # brain_NC
  df[c(prot_diff_up, prot_diff_down), rownames(file_info %>% filter(organ == 'brain', sample_type == 'N'))] # brain_N
)
dim(heatmat) # 537 1786


heatmat1 <- heatmat[apply(heatmat[, f_index], 1, function(x) sum(is.na(x)) / length(x)) <= 0.5, ] # missing ratio <= 50% in fetal
heatmat2 <- heatmat1[apply(heatmat1[, c_index], 1, function(x) sum(is.na(x)) / length(x)) <= 0.5, ] # missing ratio <= 50% in cancerous
heatmat3 <- heatmat2[apply(heatmat2[, nc_index], 1, function(x) sum(is.na(x)) / length(x)) <= 0.5, ] # missing ratio <= 50% in noncancerous
heatmat4 <- heatmat3[apply(heatmat3[, n_index], 1, function(x) sum(is.na(x)) / length(x)) <= 0.5, ] # missing ratio <= 50% in normal
heatmat4 <- log2(heatmat4)
mat_scale <- apply(heatmat4, 1, scale) %>% t() %>% data.frame()
quantile(mat_scale, na.rm = T)
colnames(mat_scale) <- colnames(heatmat)

#breaks
my_breaks <- seq(-4, 4, by = 0.01)

#colors
my_colors <- c(
  colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[1:3]))(length(my_breaks)/2/8*6),
  colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[3:5], 'white'))(length(my_breaks)/2/8*2),
  colorRampPalette(colors = c('white', RColorBrewer::brewer.pal(11, "PiYG")[7:9]))(length(my_breaks)/2/8*2),
  colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[9:11]))(length(my_breaks)/2/8*6)
)
my_colors %<>% rev() # high intensity should be hot while low to be cold




# file_info %>% dplyr::count(organ) # 65 organs
heat <- mat_scale[, rownames(file_info)]

df_color <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/input/PUH_tissue_colorset_20230210.xlsx')
organ_colors <- str_c('#', df_color$color)
set.seed(2023)
organ_colors %<>% append(sample(organ_colors, length(unique(file_info$organ)) - length(organ_colors)))
names(organ_colors) <- sort(unique(file_info$organ))
sample_type_colors <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#F39800') %>% setNames(c('T', 'NT', 'N', 'F', 'Brain'))
ann_colors <- list(sample_type = sample_type_colors, organ = organ_colors)

pheatmap::pheatmap(
  heat,
  color = my_colors,
  legend_breaks = seq(-4,4,2),
  breaks = my_breaks,
  fontsize_col = 8,
  border_color = F,
  annotation_col = file_info,
  # annotation_row = ann_row,
  gaps_col = c(which(!duplicated(file_info1$sample_type))[-1], (which(!duplicated(file_info2$sample_type)) + nrow(file_info1))) - 1,
  # gaps_col = cumsum(table(ann_col$tissue))[-length(cumsum(table(ann_col$tissue)))],
  # gaps_row = cumsum(table(row.tissue))[-length(cumsum(table(row.tissue)))],
  fontsize_row=10,
  annotation_colors = ann_colors,
  na_col = "grey70",#scale = "row",
  cluster_rows = T, cluster_cols = F,
  show_rownames = T, show_colnames = F,
  filename = 'PUH_allSample_differentiation_FCNcN_DEP_50NA_addBrain_addGlio_heatmap_20251201.pdf',
  width = 15, height = 17
)
dim(heat) # 347 2847

### ssGSEA ------
# BiocManager::install('GSVA')
# BiocManager::install('GSEABase')
# remotes::install_github("HenrikBengtsson/future.apply@develop")


# BiocManager::install("GSVA")
# install.packages(c("matrixStats","irlba","rsvd"))
library(GSVA)
library(GSEABase)

# transform prot x file matrix to prot x sampleType matrix
df_lbl <- file_info %>% rownames_to_column('FileName')
df_lbl[df_lbl$organ != 'brain', 'organ'] <- 'notBrain'
names(df_lbl$sample_type) <- NULL

filename_split <- plyr::dlply(df_lbl, c('sample_type', 'organ'), function(dfsub){
  dfsub$FileName
})

mat_median <- sapply(filename_split, function(files){
  rowMedians(heatmat4[, files], na.rm = T)
})
rownames(mat_median) <- rownames(heatmat4)
# mat_median <- mat_median[, c('F.notBrain', 'T.notBrain', 'NT.notBrain', 'N.notBrain', 'F.brain', 'T.brain', 'NT.brain', 'N.brain')] # arrange
mat_median <- mat_median[, c('F.notBrain', 'T.notBrain', 'NT.notBrain', 'N.notBrain',  'T.brain', 'NT.brain', 'N.brain')] # arrange
mat_median[is.na(mat_median)] <- min(mat_median, na.rm = T) + log2(0.8)

# # uniprotid to symbol
# writeClipboard(rownames(mat_median))
# df_query <- rio::import('idmapping_2025_10_06.xlsx')
# # P0DOX7, P0DOX5, P0DOX2, P0DOX8 not found
# pid2gnm.app <- c('P0DOX2' = 'IGA2',
#                  'P0DOX5' = 'IGG1',
#                  'P0DOX7' = 'IGK',
#                  'P0DOX8' = 'IGL1')
# symbol_vec <- c()
# for(i in 1:nrow(mat_median)){
#   message('Search the protein ', i)
#   if(!is.na(pid2gnm.app[rownames(mat_median)[i]])){
#     symbol_vec[i] <- pid2gnm.app[rownames(mat_median)[i]] %>% unname()
#   } else{
#     symbol_vec[i] <- my_uniprot2symbol(rownames(mat_median)[i])
#   }
# }
df_query_map1 <- data.frame(From = rownames(mat_median), symbol = symbol_vec) %>% full_join(df_query)

df_query_map2 <- data.frame(Protein.Group = rownames(mat_median)) %>% 
  left_join(reported_protinfo)

df_query_map <- df_query_map1 %>% dplyr::mutate(Protein.Group = From) %>% full_join(df_query_map2)


# 
# df_query_map %>% filter(symbol != To)
# # From symbol     To
# # 1 P0C0L5    C4B  C4B_2
# # 2 P12532 CKMT1A CKMT1B
df_query_map %>% filter(symbol != Genes)
df_query_map %<>% dplyr::distinct(From, symbol)

# Verify that all UniprotIDs are uniquely mapped to symbols.
df_query_map %>% dplyr::count(From) %>% dplyr::count(n)
#   n  nn
# 1 1 347
setdiff(rownames(mat_median), df_query_map$From)
# character(0)

rownames(mat_median) <- symbol_vec


#####  download  http://www.gsea-msigdb.org/gsea/msigdb/index.jsp
#  human： http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
#  mouse： http://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp
# h1gmt <- getGmt("h.all.v7.4.symbols.gmt") # hallmark
gogmt <- getGmt("c5.go.bp.v2025.1.Hs.symbols.gmt") # GO

ssgsea <- gsva(param = ssgseaParam(
  exprData = mat_median, # a matrix of expression values where rows correspond to genes and columns correspond to samples.
  geneSets = gogmt#, minSize = 15, maxSize = 500
)
# expr = mat_median,
# gset.idx.list = gogmt,# min.sz=15, max.sz=500,
# method = 'ssgsea',
# kcdf = 'Gaussian', abs.ranking = T,
)

# for pathway rename
library(msigdbr)
C5.meta <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gs_description, gs_exact_source) %>%     # description=Canonical term；exact_source=GO:ID
  dplyr::distinct()
gs_names <- sapply(gogmt, setName)
gs_sizes <- sapply(gogmt, function(gs) length(geneIds(gs)))
C5.meta2 <- C5.meta %>%
  inner_join(tibble(gs_name = gs_names, size = gs_sizes), by = "gs_name") %>%
  dplyr::mutate(term  = str_remove(gs_name, '^GOBP_'),
                label = paste0(term, " (", gs_exact_source, ", n=", size, ")"))
gs2label <- setNames(C5.meta2$label, C5.meta2$gs_name)


# min-max normalization
ssgsea.1 <- ssgsea
for (i in colnames(ssgsea)) {
  #i <- colnames(ssgsea)[1]
  ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
}
rownames(ssgsea.1) <- gs2label[rownames(ssgsea.1)]

pheatmap::pheatmap(
  ssgsea.1,
  scale = "none",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  color =colorRampPalette(c("blue", "white","red"))(100),
  cellwidth = 10, cellheight = 15,
  fontsize = 10
)




# limma deps pathway 
library(limma)

# set groups
group <- factor(c('F', 'T', 'NT', 'N', rep('Anyway', 4)))
group <- factor(c('F', 'T', 'NT', 'N', 'F', rep('Brain', 3)))
group <- factor(c('F', 'T', 'NT', 'N', rep('Brain', 3)))
# group <- factor(c('F', 'T', 'NT', 'N'))

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
rownames(design) <- colnames(ssgsea.1)[1:7]

compare <- makeContrasts(contrasts = c('F - T', 'T - NT', 'NT - N'), levels = design)
fit <- lmFit(ssgsea.1[, 1:7], design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, number=Inf) %>% dplyr::filter(P.Value < 0.01)

Diff %>% rownames_to_column('Pathway') %>% rio::export('PUH_allSample_differentiation_FCNcN_DEP_50NA_addBrain_addGlio_rmPlasma_GSVA_limmaFilter_20251201.xlsx')


heat_ssgsea.1 <- ssgsea.1[rownames(Diff), ]
heat_ssgsea.1 <- heat_ssgsea.1[apply(heat_ssgsea.1, 1, function(x) identical(sort(x[1:4]), x[1:4])|identical(sort(x[1:4], decreasing = T), x[1:4])), ]
a <- pheatmap::pheatmap( # to perform a "hc" tree
  heat_ssgsea.1[, 1:4],
  scale = "none",
  cluster_rows = T,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  color =colorRampPalette(c("blue", "white","red"))(100),
  cellwidth = 10, cellheight = 15,
  fontsize = 10,
  filename = 'PUH_allSample_differentiation_FCNcN_DEP_50NA_addBrain_addGlio_rmPlasma_GSVA_heatmap_20251201.pdf',
  # width = 10, height = 15
)

heat_ssgsea.1 <- heat_ssgsea.1[a$tree_row$order, ]
ann_col <- data.frame(
  sample_type = c('F', 'T', 'NT', 'N', 'T', 'NT', 'N'),
  row.names = colnames(heat_ssgsea.1)
)
b <- pheatmap::pheatmap(
  heat_ssgsea.1,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  scale = "none",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  gaps_col = 4,
  color = colorRampPalette(c("blue", "white","red"))(100),
  cellwidth = 10, cellheight = 15,
  fontsize = 10,
  filename = 'PUH_allSample_differentiation_FCNcN_DEP_50NA_addBrain_addGlio_rmPlasma_GSVA_heatmap_v2_20251201.pdf',
  # width = 10, height = 15
)

heat_ssgsea.1 %>%
  as.data.frame() %>%
  rownames_to_column('Pathway') %>%
  rio::export('PUH_allSample_differentiation_FCNcN_DEP_50NA_addBrain_addGlio_rmPlasma_GSVA_heatmap_v2_20251201.xlsx')




# not finished -----
### GOBP terms clusterProfiler::simplify
# BiocManager::install('RamiGO')
library(rrvgo)
# x <- rownames(Diff)[1] %>% str_replace('^GOBP_', '') %>% str_replace_all('_', ' ')

Diff1 <- Diff[rownames(heat_ssgsea.1), ]

msigdbC52GO <- function(x){
  # QuickGO version 2022-03-14
  x %<>% str_replace_all(' ', '%20')
  my_html <- rvest::read_html(stringr::str_glue('https://www.ebi.ac.uk/QuickGO/services/internal/search/ontology?query={x}')) # get json
  # for geneproduct:  https://www.ebi.ac.uk/QuickGO/services/geneproduct/search?query=
  my_json <- my_html %>%
    rvest::html_text() %>%
    jsonlite::fromJSON()
  
  return(my_json$results[1, 1:2])
}

GO_query_rlt <- list()
for(i in 1:nrow(Diff1)){
  cat(i, '...\r')
  x <- rownames(Diff1)[i] %>% str_replace('^GOBP_', '') %>% str_replace_all('_', ' ')
  GO_query_rlt[[i]] <- msigdbC52GO(x)
}
names(GO_query_rlt) <- rownames(Diff1)
df_GO <- plyr::ldply(GO_query_rlt, .id = 'msigdb_name')

# go_analysis <- Diff1
simMatrix <- calculateSimMatrix(df_GO$id,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
# scores <- setNames(-log10(Diff1$adj.P.Val), df_GO$id)
scores <- setNames(-log10(Diff1$P.Value), df_GO$id)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

df_reduce <- reducedTerms %>%
  left_join(df_GO, by = c(go = 'id')) %>%
  dplyr::rename(score_parent = parent,
                score_parentTerm = parentTerm)
df_reduce <- plyr::ddply(df_reduce, 'cluster', function(dfsub){
  dfsub %<>% arrange(desc(size))
  dfsub$parent <- dfsub$go[1]
  dfsub$parentTerm <- dfsub$term[1]
  return(dfsub)
})

# rio::export(df_reduce, 'PUH_allSample_differentiation_FCNcN_DEP_50NA_addBrain_addGlio_rmPlasma_GSVA_limmaFilter_termReduce_20251201.xlsx')

# pdf('PUH_allSample_differentiation_FCNcN_DEP_50NA_addBrain_addGlio_rmPlasma_GSVA_limmaFilter_termReduce_20251201.pdf', width = 10, height = 10)
# par(mfrow = c(2, 2))
heatmapPlot(simMatrix,
            df_reduce ,row_cluster = F,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=10)
scatterPlot(simMatrix, reducedTerms)
treemapPlot(reducedTerms)
wordcloudPlot(reducedTerms, min.freq=1, colors="black")
# graphics.off()


msigdb2go <- df_GO$id
names(msigdb2go) <- df_GO$msigdb_name
rownames(heat_ssgsea.1) <- msigdb2go[rownames(heat_ssgsea.1)]

#
pheatmap::pheatmap(
  heat_ssgsea.1,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  scale = "none",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  gaps_col = 4,
  color = colorRampPalette(c("blue", "white","red"))(100),
  cellwidth = 10, cellheight = 15,
  fontsize = 10,
  filename = 'PUH_allSample_differentiation_FCNcN_DEP_50NA_addBrain_addGlio_rmPlasma_GSVA_heatmap_v3_20251201.pdf',
  width = 10, height = 15
)


#### 2.4.3.x GOBP terms clusterProfiler::simplify
# 1. First, deduplicate the matrix by row.
# 2. Calculate up-regulation and down-regulation separately.
##### 

# Deduplicate the matrix by row and filter for GO terms with completely identical expression patterns.
same_terms <- heat_ssgsea.1 %>%
  as.data.frame %>% 
  rownames_to_column('go') %>% 
  dplyr::mutate(go = str_replace(go, '^GO\\.', 'GO:')) %>% 
  plyr::dlply(., colnames(ssgsea.1), function(dfsub){
    dfsub$go
  })
length(same_terms) # 7
same_terms1 <- same_terms[sapply(same_terms, length) > 1] # 筛选GO terms > 1个的结果
names(same_terms1) <- seq_along(same_terms1)

Diff2 <- Diff1 %>% # combine
  rownames_to_column('msigdb_name') %>% 
  left_join(df_GO)
rownames(Diff2) <- Diff2$msigdb_name

df_same_term <- plyr::ldply(same_terms1, function(term_vec){
  Diff2 %>% filter(id %in% term_vec)
}, .id = 'Group')
df_same_term$Group %<>% as.character()

term_relations <- as.list(GO.db::GOBPANCESTOR)[df_same_term$id]
df_term_relations <- plyr::ldply(term_relations, function(x) {
  data.frame(from = x)
}, .id = 'to')
df_term_relations$to %<>% as.character()
df_term_relations %<>%
  dplyr::select(from, to) %>% 
  dplyr::filter(from %in% unlist(same_terms1), to %in% unlist(same_terms1))

# add root
df_term_relations <-
  data.frame(from = 'all', to = unique(c(as.character(df_term_relations$to), df_term_relations$from))) %>%
  rbind(df_term_relations)

# add Group
vertices <- data.frame(
  name = df_same_term$id,
  group = df_same_term$Group
  # cluster = sample(letters[1:4], length(name), replace=T),
  # value = 1
) %>% dplyr::filter(name %in% unlist(df_term_relations))
vertices %<>% rbind(data.frame(name = 'all', group = 0))

v_colors <- RColorBrewer::brewer.pal(9, 'Set1') %>% setNames(c(1:8, 0))
# 1         2         3         4         5         6         7         8         0 
"#E41A1C"; "#377EB8"; "#4DAF4A"; "#984EA3"; "#FF7F00"; "#FFFF33"; "#A65628"; "#F781BF"; "#999999" 

mygraph <- igraph::graph_from_data_frame(df_term_relations, vertices = vertices)
# pdf('identical_expression_GO_network_20251201.pdf', width = 10, height = 10)
set.seed(2023)
plot(mygraph, vertex.color = v_colors[vertices$group], # Node color
     vertex.frame.color = "Forestgreen",            # Node border color
     # vertex.shape=c("circle","square"),             # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
     vertex.size=15,                          # Size of the node (default is 15)
     vertex.size2=NA,                               # The second size of the node (e.g. for a rectangle))
     vertex.label.dist=0,                           # Distance between the label and the vertex
     vertex.label.degree=0 ,                        # The position of the label in relation to the vertex (use pi)
     edge.color='black',           # Edge color
     edge.width=1,                        # Edge width, defaults to 1
     edge.arrow.size=0.5,                           # Arrow size, defaults to 1
     edge.arrow.width=0.5,                          # Arrow width, defaults to 1
     edge.lty=c("solid")                           # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
     #edge.curved=c(rep(0,5), rep(1,5))            # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
)
# graphics.off()






list(
  v = vertices %>% left_join(dplyr::select(df_same_term, -Group), by = c(name = 'id')),
  g = df_term_relations
) %>% rio::export('identical_expression_GO_network_20251201.xlsx')

terms_rm <- vertices$name %>% setdiff(
  c('GO:0000041', 'GO:1902904', 'GO:1902742', 'GO:0071478', 'GO:0043254', 'GO:0060267', 'GO:0003014', 'GO:0021761', 'GO:0051962', 'GO:0010720', 'GO:1904406', 'GO:0050779', 'GO:0010985', 'all')
)

df_same_term_filtered <- df_same_term %>%
  dplyr::filter(!(id %in% terms_rm))





# Calculate up-regulation and down-regulation separately.
Diff_down <- Diff2 %>% filter(logFC > 0, id %in% df_same_term_filtered$id)
Diff_up <- Diff2 %>% filter(logFC < 0, id %in% df_same_term_filtered$id)
df_GO_down <- df_GO %>% filter(msigdb_name %in% rownames(Diff_down))
df_GO_up <- df_GO %>% filter(msigdb_name %in% rownames(Diff_up))

#down
simMatrix_down <- calculateSimMatrix(df_GO_down$id,
                                     orgdb="org.Hs.eg.db",
                                     ont="BP",
                                     method="Rel")
scores_down <- setNames(-log10(Diff_down$P.Value), df_GO_down$id)
reducedTerms_down <- reduceSimMatrix(simMatrix_down,
                                     scores_down,
                                     threshold=0.7,
                                     orgdb="org.Hs.eg.db")

#up
simMatrix_up <- calculateSimMatrix(df_GO_up$id,
                                   orgdb="org.Hs.eg.db",
                                   ont="BP",
                                   method="Rel")
scores_up <- setNames(-log10(Diff_up$P.Value), df_GO_up$id)
reducedTerms_up <- reduceSimMatrix(simMatrix_up,
                                   scores_up,
                                   threshold=0.7,
                                   orgdb="org.Hs.eg.db")

list(down = reducedTerms_down %>% arrange(cluster), up = reducedTerms_up %>% arrange(cluster)) %>%
  rio::export('PUH_allSample_differentiation_FCNcN_DEP_50NA_addBrain_addGlio_rmPlasma_GSVA_limmaFilter_termReduce_up-down_20251201.xlsx')


pdf('PUH_allSample_differentiation_FCNcN_DEP_50NA_addBrain_addGlio_rmPlasma_GSVA_limmaFilter_termReduce_20251201.pdf', width = 10, height = 10)
# par(mfrow = c(2, 2))
heatmapPlot(simMatrix_up,
            df_reduce ,row_cluster = F,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=10)
scatterPlot(simMatrix_up, reducedTerms_up)
treemapPlot(reducedTerms_up)
wordcloudPlot(reducedTerms_up, min.freq=1, colors="black")
heatmapPlot(simMatrix_down,
            df_reduce ,row_cluster = F,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=10)
scatterPlot(simMatrix_down, reducedTerms_down)
treemapPlot(reducedTerms_down)
wordcloudPlot(reducedTerms_down, min.freq=1, colors="black")
graphics.off()



#
heat_sorted <- heat_ssgsea.1[c(rownames(reducedTerms_up), rownames(reducedTerms_down)), ]
go2lbl <- str_glue("{str_remove(df_GO$msigdb_name, '^GOBP_') %>% str_replace_all('_', ' ') %>% str_to_sentence()} ({df_GO$id})")
names(go2lbl) <- df_GO$id
rownames(heat_sorted) <- go2lbl[rownames(heat_sorted)]

pheatmap::pheatmap(
  heat_sorted[c(2, 9, 7, 1, 3, 8, 11, 4:6, 10, 12, 14:16, 19:20, 13, 17:18), ],
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  scale = "none",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  gaps_col = 4,
  color = colorRampPalette(c("blue", "white","red"))(100),
  cellwidth = 10, cellheight = 15,
  fontsize = 10,
  filename = 'PUH_allSample_differentiation_FCNcN_DEP_50NA_addBrain_addGlio_rmPlasma_GSVA_heatmap_20251201.pdf',
  width = 8, height = 8
)


















