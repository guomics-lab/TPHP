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
source('../source/my_fun.R')

dfprot <- rio::import('../0_process_DIA-NN_data/output/protein_info_from_reports.csv')
# df_abbr <- read.delim('//172.16.13.136/share/members/jiangwenhao/TPHP/input/sample_types_abbr_20230113.txt', stringsAsFactors = F, check.names = F, na.strings = '')
# h <- hash::hash(keys = df_abbr$'Entire', values = df_abbr$'Abbr')


## helper functions -----
pca_plot <- function(M, title, info,
                     color_var = "Dataset",
                     key = 'FileName',
                     label_var = "anatomical_classification",
                     n_top = 5000L,
                     abbrev_min = 3L) {
  stopifnot(is.matrix(M) || is.data.frame(M))
  X <- as.matrix(M)
  
  # 1) Use top-variance features to stabilize PCs
  if (nrow(X) > n_top) {
    v <- matrixStats::rowVars(X)
    keep <- order(v, decreasing = TRUE)[seq_len(n_top)]
    X <- X[keep, , drop = FALSE]
  }
  
  # 2) Samples are columns (features × samples) → transpose for prcomp
  #    Remove any columns with non-finite values as a guard
  samp <- colnames(X)
  good_cols <- colSums(is.finite(X)) == nrow(X)
  if (!all(good_cols)) X <- X[, good_cols, drop = FALSE]
  samp <- colnames(X)
  
  # 3) Align metadata rows to samples, trying common keys in order
  meta <- tibble::as_tibble(info)
  meta$.rown <- rownames(info)
  if(is.null(key)){
    key <- c("sample_id", "file_id", ".rown")[c("sample_id","file_id",".rown") %in% names(meta)][1]
  }
  if (!is.na(key) && all(samp %in% meta[[key]])) {
    meta <- meta[match(samp, meta[[key]]), , drop = FALSE]
  } else if (nrow(meta) == length(samp)) {
    # assume already aligned
  } else {
    stop("Cannot align `info` to samples: provide a column matching colnames(M).")
  }
  
  # 4) PCA
  pc <- stats::prcomp(t(X), center = TRUE, scale. = TRUE)
  var_exp <- (pc$sdev^2) / sum(pc$sdev^2)
  subtitle <- sprintf("Variance: PC1=%.1f%%, PC2=%.1f%%", 100 * var_exp[1], 100 * var_exp[2])
  
  # 5) Build plotting/data frame
  df <- cbind(
    data.frame(sample = samp, PC1 = unclass(pc$x[, 1]), PC2 = unclass(pc$x[, 2])),
    meta
  )
  
  # Color mapping (factor to ensure consistent legend order across panels)
  if (!is.null(color_var) && color_var %in% names(df)) {
    df[[color_var]] <- as.factor(df[[color_var]])
  } else {
    color_var <- NULL
  }
  
  # Abbreviated labels per class (drawn only if abbrev_min > 0)
  draw_labels <- !is.null(label_var) && label_var %in% names(df) && abbrev_min > 0L
  if (draw_labels) {
    df$.lbl <- base::abbreviate(as.character(df[[label_var]]),
                                minlength = abbrev_min, strict = TRUE)
  }
  
  # 6) Plot
  df <- tibble::as_tibble(df, .name_repair = "unique")  # ensure unique column names
  
  p <- ggplot2::ggplot(
    df,
    if (!is.null(color_var))
      ggplot2::aes(x = PC1, y = PC2, color = .data[[color_var]])
    else
      ggplot2::aes(x = PC1, y = PC2)
  ) +
    # {
    #   if (!is.null(label_var) && label_var %in% names(df) && abbrev_min > 0L) {
    #     # draw abbreviated class labels
    #     ggplot2::geom_text(ggplot2::aes(label = .data$.lbl), alpha = 0.85, size = 3)
    #   } else {
    #     # no labels requested -> emulate points using a small text glyph
    #     ggplot2::geom_text(label = ".", alpha = 0.9, size = 3)
    #   }
    # } +
    {
      ggplot2::geom_point(alpha = 0.85, size = 1.5)
    } +
    ggplot2::labs(title = title, subtitle = subtitle, x = "PC1", y = "PC2",
                  color = color_var) +
    ggplot2::theme_bw(base_size = 12)
  
  # 7) Optional: enforce a stable color palette across panels if available
  if (!is.null(color_var)) {
    lv <- levels(df[[color_var]])
    if (exists(".Dataset_palette", mode = "function")) {
      pal <- .Dataset_palette(length(lv))
      p <- p + ggplot2::scale_color_manual(drop = FALSE, values = stats::setNames(pal, lv))
    } else {
      p <- p + ggplot2::scale_color_discrete(drop = FALSE)
    }
  }
  
  list(pca_plot = p, pca_df = df)
}

pca_plot_from_df <- function(df,
                             title,
                             color_var = "Dataset",
                             label_var = NULL,
                             abbrev_min = 3L,
                             subtitle = NULL) {
  stopifnot(is.data.frame(df))
  if (!all(c("PC1", "PC2") %in% names(df))) {
    stop("Input data frame must contain columns 'PC1' and 'PC2'.")
  }
  
  
  df <- tibble::as_tibble(df, .name_repair = "unique")
  
  
  # Color mapping
  if (!is.null(color_var) && color_var %in% names(df)) {
    df[[color_var]] <- as.factor(df[[color_var]])
  } else {
    color_var <- NULL
  }
  
  
  # Labels: use precomputed '.lbl' if present; otherwise create if possible
  draw_labels <- !is.null(label_var) && label_var %in% names(df) && abbrev_min > 0L
  if (draw_labels && !(".lbl" %in% names(df))) {
    df$.lbl <- base::abbreviate(as.character(df[[label_var]]),
                                minlength = abbrev_min, strict = TRUE)
  }
  
  
  aes_map <- if (!is.null(color_var)) {
    ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data[[color_var]])
  } else {
    ggplot2::aes(x = .data$PC1, y = .data$PC2)
  }
  
  
  p <- ggplot2::ggplot(df, aes_map) +
    {
      if (draw_labels && ".lbl" %in% names(df)) {
        ggplot2::geom_text(ggplot2::aes(label = .data$.lbl), alpha = 0.85, size = 3)
      } else {
        ggplot2::geom_point(alpha = 0.85, size = 1.5)
      }
    } +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = "PC1", y = "PC2", color = color_var) +
    ggplot2::theme_bw(base_size = 12)
  
  
  # Optional: stable color palette if available
  if (!is.null(color_var)) {
    lv <- levels(df[[color_var]])
    if (exists(".Dataset_palette", mode = "function")) {
      pal <- .Dataset_palette(length(lv))
      p <- p + ggplot2::scale_color_manual(drop = FALSE,
                                           values = stats::setNames(pal, lv))
    } else {
      p <- p + ggplot2::scale_color_discrete(drop = FALSE)
    }
  }
  
  
  return(p)
}

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
info_pan <- rio::import('../6_pancancer_harmonizR_QC/output/PUH_pancancer_sample_information_2datesets_1146files.xlsx')
info.all <- rio::import('input/20251009_PUH_sample_information_2997files.xlsx')
# datall <- rio::import('input/imputed_raw_12754proteins_2997samples.csv') %>% column_to_rownames('V1')
datall <- rio::import('../2_batch_effect_evaluation_correction/output/batch_corrected_data_combat3.csv') %>% column_to_rownames('V1')

out_dir   <- "subtype_cluster_V1009_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# 2. Load data ----
## 2.1 Pan-cancer data -----
meta <- info.all %>% filter(sample_type == 'T')
expr <- datall[, meta$FileName]


dim(meta) # 1105   18
dim(expr) # 12754  1105


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
dim(mat_consensus) # 6377 1105
quantile(mat_consensus)
# 0%       25%       50%       75%      100% 
# 3.368983  9.455099 12.306143 14.488626 25.897229 

ann_df <- meta %>%
  select(FileName, tissue_name, class_abbr, cancer_abbr, instrument, new) %>%
  mutate(across(-FileName, \(x) factor(x))) %>%
  column_to_rownames("FileName")

ann_colors <- lapply(ann_df, function(f){
  lv <- levels(f)
  setNames(get_scicol(length(lv)), lv)
})

a <- pheatmap::pheatmap(
  mat_consensus,
  color = c(RColorBrewer::brewer.pal(11,"RdYlBu")[9:7],
            RColorBrewer::brewer.pal(11,"RdYlBu")[4:2]),
  scale = "row", annotation_col = ann_df, annotation_colors = ann_colors,
  cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, show_colnames = FALSE,
  fontsize_col = 8, filename = "consensus_input_heatmap.pdf", width = 10, height = 10
)

## 4.2 Consensus -------
saveRDS(mat_consensus, 'mat_consensus.1105tumors_6377proteins.cmb.rds')
mat_consensus <- readRDS('mat_consensus.1105tumors_6377proteins.cmb.rds')


# ** columns=items/samples; and rows are features/proteins **
library(ConsensusClusterPlus)
res.consensus <- ConsensusClusterPlus(mat_consensus, maxK = 20, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg = 'km', distance = 'euclidean', seed = 8, plot = 'pdf', writeTable = T, title = 'pancancer_1105s_6377p_cmb_consen', verbose = T)
# save(res.consensus, file = 'pancancer_1105s_6377p_cmb_consen.RData') # top50% CV

# load('pancancer_1105s_6377p_cmb_consen.RData')
# 

### clustering stats -----
# REFERNCE:
# DESeq2 RNAseq_unsupervised.R [line 170-226]
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
# https://github.com/IARCbioinfo/RNAseq_analysis_scripts/blob/master/RNAseq_unsupervised.R
library(optparse)
require(ConsensusClusterPlus)
library(ade4)
library(DESeq2)
library(fpc)
library(cluster)

maxK = 20
d = mat_consensus
clusters = res.consensus

# prettycolors = c(1,2,rgb(0,152/255,219/255),rgb(233/255,131/255,0),rgb(0,155/255,118/255),rgb(0.5,0.5,0.5),rgb(0,124/255,146/255),rgb(178/255,111/255,22/255),rgb(234/255,171/255,0),rgb(83/255,40/255,79/255))
prettycolors <- colorRampPalette(c('red', 'green', 'blue'))(maxK)


# compute clustering stats
stcl   = lapply(2:maxK, function(i) cluster.stats(dist(t(d)),clusters[[i]]$consensusClass) )
ldunn = sapply(1:(maxK-1), function(i) stcl[[i]]$dunn )
lwbr  = sapply(1:(maxK-1), function(i) stcl[[i]]$wb.ratio ) #c(stcl.B$wb.ratio,stcl.hc$wb.ratio,stcl.Whc$wb.ratio,stcl.W2hc$wb.ratio,stcl.km$wb.ratio)
lch   = sapply(1:(maxK-1), function(i) stcl[[i]]$ch ) #c(stcl.B$ch,stcl.hc$ch,stcl.Whc$ch,stcl.W2hc$ch,stcl.km$ch)
lsil = vector("list",(maxK-1))

pdf("pancancer_1105s_6377p_cmb_consen_silhouette.pdf",h=4*ceiling(maxK / ceiling(sqrt(maxK))),w=4*ceiling(sqrt(maxK)))
par(mfrow=c(ceiling(maxK / ceiling(sqrt(maxK))),ceiling(sqrt(maxK))))
for(i in 2:maxK){
  cat(i, '...\r')
  sil = silhouette(clusters[[i]]$consensusClass,dist(t(d),method = "euclidean"))
  sizes = table(clusters[[i]]$consensusClass)
  plot( sil ,col=rep( clusters[[i]]$clrs[[3]],rep=sizes) ,main=paste("K=",i))
  lsil[[i-1]]=sil
}
dev.off()

msil = sapply(1:(maxK-1), function(i) mean( lsil[[i]][,3] ) )
cdl = lapply(2:maxK, function(i) as.dist(1-clusters[[i]]$consensusMatrix ) )
md = dist( t(d),method = "euclidean")
corl =sapply(cdl, cor,md)

# plot clustering stats
pdf("pancancer_1105s_6377p_cmb_consen_cluster_separation_stats.pdf",h=3.5,w=3.5*5)
par(mfrow=c(1,5),family="Times")
co = rep(1,(maxK-1))
co[which.max(ldunn)]=2
barplot(ldunn ,names.arg = paste("K =",2:maxK) ,las=2,ylab="Dunn index",col= co)
co = rep(1,(maxK-1))
co[which.min(lwbr)]=2
barplot(lwbr ,names.arg = paste("K =",2:maxK) ,las=2,ylab="Within-between SS ratio",col= co)
co = rep(1,(maxK-1))
co[which.max(lch)]=2
barplot(lch ,names.arg = paste("K =",2:maxK) ,las=2,ylab="Calinski-Harabasz index",col= co)
co = rep(1,(maxK-1))
co[which.max(msil)]=2
barplot(msil,names.arg = paste("K =",2:maxK) ,las=2,ylab="Mean Silhouette",col= co)
co = rep(1,(maxK-1))
co[which.max(corl)]=2
barplot(corl,names.arg = paste("K =",2:maxK) ,las=2,ylab="Mean cophenetic distance",col= co)
dev.off()


# plot PCA with clusters
set.seed(2023)
pca <- ade4::dudi.pca(t(d),scannf = F,nf = 10,center = T, scale = F)

pdf("pancancer_1105s_6377p_cmb_consen_PCA.pdf",h=4*ceiling(maxK / ceiling(sqrt(maxK))),w=4*ceiling(sqrt(maxK)))
par(mfrow=c(ceiling(maxK / ceiling(sqrt(maxK))),ceiling(sqrt(maxK))),family="Times")
for(i in 2:(maxK)) s.class(pca$li,as.factor(clusters[[i]]$consensusClass),col=prettycolors,xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
dev.off()

# plot t-SNE with clusters
set.seed(2023)
dtsne = sweep(d, 1, apply(d,1,median,na.rm=T))
tsne <- Rtsne::Rtsne(t(dtsne), dims = 2, perplexity = floor((nrow(t(dtsne))-1)/3-1), verbose=T, max_iter = 1000,check_duplicates = F,theta = 0)

pdf("pancancer_1105s_6377p_cmb_consen_tSNE.pdf",h=4*ceiling(maxK / ceiling(sqrt(maxK))),w=4*ceiling(sqrt(maxK)))
par(mfrow=c(ceiling(maxK / ceiling(sqrt(maxK))),ceiling(sqrt(maxK))),family="Times")
for(i in 2:(maxK)) s.class(tsne$Y,as.factor(clusters[[i]]$consensusClass),col=prettycolors,xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
dev.off()



for(k in 2:20){
  message("k = ", k)
  sil <- lsil[[k-1]] # silhouette results of k=2 stored in lsil[[k-1]]
  df_sil <- as.data.frame(sil) %>%
    mutate(FileName = names(res.consensus[[k]]$consensusClass)) %>%
    inner_join(info.all)
  
  ann_df <- df_sil %>%
    select(FileName, cluster, tissue_name, anatomical_classification, instrument) %>%
    mutate(across(-FileName, \(x) factor(x))) %>%
    column_to_rownames("FileName")
  
  ann_colors <- lapply(ann_df, function(f){
    lv <- levels(f)
    setNames(get_scicol(length(lv)), lv)
  })
  
  a <- pheatmap::pheatmap(
    mat_consensus,
    color = c(RColorBrewer::brewer.pal(11,"RdYlBu")[9:7],
              RColorBrewer::brewer.pal(11,"RdYlBu")[4:2]),
    scale = "row", annotation_col = ann_df, annotation_colors = ann_colors,
    cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, show_colnames = FALSE, clustering_method = 'ward.D2',
    fontsize_col = 8, filename = str_glue("pancancer_1105s_6377p_cmb_consen_k{k}_heatmap.pdf"), width = 10, height = 5
  )
}




# ===============

## 4.3 select core samples -------
k <- 5
sil <- lsil[[k-1]] # silhouette results of k=20 stored in lsil[[k-1]]
df_sil <- as.data.frame(sil) %>% mutate(FileName = names(res.consensus[[k]]$consensusClass))

dim(df_sil) # 1105    4
df_sil %>% dplyr::count(cluster)
# cluster   n
# 1       1 178
# 2       2 364
# 3       3  74
# 4       4 398
# 5       5  91

df_sil1 <- df_sil %>% group_by(cluster) %>%
  arrange(desc(sil_width), .by_group = T) %>% 
  ungroup() %>% 
  filter(sil_width > 0)

dim(df_sil1) # 1014    4
df_sil1 %>% dplyr::count(cluster)
# cluster     n
# 1       1 178
# 2       2 312
# 3       3  68
# 4       4 365
# 5       5  91

df_sil1 %<>% left_join(info.all)
df_sil1$cluster %<>% factor()
df_sil1$neighbor %<>% factor()
df_sil1 %<>% as.data.frame()
rownames(df_sil1) <- df_sil1$FileName

clustering <- res.consensus[[k]]$consensusClass
consensus_mat <- res.consensus[[k]]$consensusMatrix
colnames(consensus_mat) <- rownames(consensus_mat) <- names(clustering)


# heatmap
ann_colors <- list(cancer_abbr = cancer_color)

# df_sil1 %>% distinct(cancer_abbr, tissue_name) %>% count(cancer_abbr) %>% filter(n > 1) %>% semi_join(df_sil1, .) %>% distinct(cancer_abbr, tissue_name) %>% arrange(cancer_abbr) %>% dplyr::mutate(cancer_abbr = factor(cancer_abbr))


a <- pheatmap::pheatmap( # simple heatmap
  consensus_mat,
  # color = c(brewer.pal(11,"RdYlBu")[9:7],brewer.pal(11,"RdYlBu")[4:2]), fontsize_col = 8,
  scale = 'none',
  annotation_col = as.data.frame(as.factor(clustering)),
  cluster_rows = T, cluster_cols = T,
  show_rownames = F, show_colnames = F
)

a_ <- pheatmap::pheatmap( # simple heatmap
  consensus_mat[df_sil1$FileName, df_sil1$FileName],
  # color = c(brewer.pal(11,"RdYlBu")[9:7],brewer.pal(11,"RdYlBu")[4:2]), fontsize_col = 8,
  scale = 'none',
  annotation_col = df_sil1 %>% select(cluster, cancer_abbr, tissue_name),
  annotation_colors = ann_colors,
  cluster_rows = T, cluster_cols = T,
  show_rownames = F, show_colnames = F,
  fontsize = 10,
  filename = 'PTDB_pan-cancer_6377_proteins_1014_coreSamples_reproduce_5clusters.pdf'
)


#breaks
quantile(scale(mat_consensus[, df_sil1$FileName]))
# my_breaks <- unique(c(seq(-1.5, 0, length.out = 251), 0, seq(0, 6, length.out = 251)))
my_breaks <- seq(-1.5, 1.5, by = 0.01)

#colors
my_colors <- rev(c(colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[1:5], 'white'))(length(my_breaks)/2),
                   colorRampPalette(colors = c('white', RColorBrewer::brewer.pal(11, "PiYG")[7:11]))(length(my_breaks)/2)))

dim(mat_consensus[, df_sil1$FileName]) # 6377 1014
b_ <- pheatmap::pheatmap( # simple heatmap
  mat_consensus[, df_sil1$FileName],
  # color = c(brewer.pal(11,"RdYlBu")[9:7],brewer.pal(11,"RdYlBu")[4:2]),
  color = my_colors,
  legend_breaks = quantile(my_breaks),
  breaks = my_breaks,
  fontsize_col = 8,
  scale = 'row',
  annotation_col = df_sil1 %>% dplyr::select(cluster, cancer_abbr, tissue_name) %>% 
    dplyr::mutate(cancer_abbr = factor(cancer_abbr)),
  annotation_colors = ann_colors,
  cutree_rows = 4, cutree_cols = 4,
  cluster_rows = T, cluster_cols = T,
  clustering_method = 'ward.D2',
  show_rownames = F, show_colnames = F,
  fontsize = 10,
  filename = 'PTDB_pan-cancer_6377_proteins_1014_coreSamples.pdf'
)


df_sil1 %>% rio::export('PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters.xlsx')


## 4.4 select signatures -------
mat_sig <- mat_consensus[, df_sil1$FileName]
cls <- unique(df_sil1$cluster)

### DEA --------
df_dea_ls <- list()
for(i in 1:nrow(mat_sig)){ # for each protein
  cat(i, '...\r')
  log2FC <- P <- c()
  for(ii in seq_along(cls)){ # for each cluster
    cl <- cls[ii]
    fnames1 <- df_sil1 %>% filter(cluster == cl) %>% pull(FileName)
    fnames2 <- setdiff(df_sil1$FileName, fnames1)
    x1 <- mat_sig[i, fnames1]
    x2 <- mat_sig[i, fnames2]
    
    # fold-change
    log2FC[ii] <- log2(mean(2^x1) / mean(2^x2))
    
    # wilcoxon rank test
    P[ii] <- wilcox.test(x1, x2, alternative = 'two.sided', paired = F)$p.value
  }
  
  adj.P <- p.adjust(P, 'BH')
  names(log2FC) <- str_c('log2FC_', cls, 'vs.')
  names(P) <- str_c('P_', cls, 'vs.')
  names(adj.P) <- str_c('adj.P_', cls, 'vs.')
  df_tmp <- c(log2FC, P, adj.P) %>%
    as.data.frame() %>%
    setNames(rownames(mat_sig)[i]) %>%
    t() %>% as.data.frame()
  
  df_dea_ls[[i]] <- df_tmp
}
names(df_dea_ls) <- rownames(mat_sig)
df_dea <- plyr::ldply(df_dea_ls, .id = 'Protein') %>% mutate(Protein = as.character(Protein))


df_dea_ <- dfprot %>% rename(Protein = Protein.Group) %>% right_join(df_dea)
rio::export(df_dea_, 'PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_DEA.xlsx')

df_dea_signif3 <- df_dea_ %>% mutate_if(is.numeric, signif, digits = 3)
rio::export(df_dea_signif3, 'PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_DEA_signif3.xlsx')


# DEP cutoff: |FC| >= 1.5, adj.P < 0.05
#DEP_up
prot_up_ls <- list()
for(i in seq_along(cls)){
  cl <- cls[i]
  df_tmp <- df_dea %>% select(Protein, matches(str_c(as.character(cl), 'vs.')))
  # isSignificant <- abs(df_tmp$log2FC) >= 2 & df_tmp$adj.P < 0.05
  isUp <- df_tmp$log2FC >= 1.5 & df_tmp$adj.P < 0.05 # select only up-regulated proteins
  prot_up_ls[[i]] <- df_tmp$Protein[isUp]
}

prot_specific_up_ls <- lapply(seq_along(prot_up_ls), function(i){
  setdiff(prot_up_ls[[i]], unlist(prot_up_ls[-i]))
})

names(prot_up_ls) <- names(prot_specific_up_ls) <- cls
dep_up <- plyr::ldply(prot_up_ls, .id = 'cluster', function(x){
  x %>% as.data.frame() %>% setNames('Protein')
})
dep_up %<>% left_join(df_dea_, by = 'Protein')

sapply(prot_up_ls, length)
# 1    2    3    4    5 
# 293  782 1003  224  606 

sapply(prot_specific_up_ls, length)
# 1   2   3   4   5 
# 280 770 939 222 561 

dep_specific_up <- plyr::ldply(prot_specific_up_ls, .id = 'cluster', function(x){
  x %>% as.data.frame() %>% setNames('Protein')
})
dep_specific_up %<>% left_join(df_dea_, by = 'Protein')
rownames(dep_specific_up) <- dep_specific_up$Protein
dim(dep_specific_up) # 2772   18


#DEP_down
prot_down_ls <- list()
for(i in seq_along(cls)){
  cl <- cls[i]
  df_tmp <- df_dea %>% select(Protein, matches(str_c(as.character(cl), 'vs.')))
  # isSignificant <- abs(df_tmp$log2FC) >= 2 & df_tmp$adj.P < 0.05
  isDown <- df_tmp$log2FC <= -1.5 & df_tmp$adj.P < 0.05 # select only down-regulated proteins
  prot_down_ls[[i]] <- df_tmp$Protein[isDown]
}

prot_specific_down_ls <- lapply(seq_along(prot_down_ls), function(i){
  setdiff(prot_down_ls[[i]], unlist(prot_down_ls[-i]))
})

names(prot_down_ls) <- names(prot_specific_down_ls) <- cls
dep_down <- plyr::ldply(prot_down_ls, .id = 'cluster', function(x){
  x %>% as.data.frame() %>% setNames('Protein')
})
dep_down %<>% left_join(df_dea_, by = 'Protein')

sapply(prot_down_ls, length)
# 1    2    3    4    5 
# 975  579 1373 1050 1406 

sapply(prot_specific_down_ls, length)
# 1   2   3   4   5 
# 469 322 485 355 493

dep_specific_down <- plyr::ldply(prot_specific_down_ls, .id = 'cluster', function(x){
  x %>% as.data.frame() %>% setNames('Protein')
})
dep_specific_down %<>% left_join(df_dea_, by = 'Protein')
rownames(dep_specific_down) <- dep_specific_down$Protein
dim(dep_specific_down) # 2124   18



# heatmap -- only upregulated proteins
#breaks
quantile(scale(mat_sig[dep_specific_up$Protein, ]))
# my_breaks <- unique(c(seq(-1.5, 0, length.out = 251), 0, seq(0, 3, length.out = 251)))
my_breaks <- seq(-1.5, 1.5, by = 0.01)

#colors
my_colors <- rev(c(colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[1:5], 'white'))(length(my_breaks)/2),
                   colorRampPalette(colors = c('white', RColorBrewer::brewer.pal(11, "PiYG")[7:11]))(length(my_breaks)/2)))

ann_colors <- list(
  cancer_abbr = cancer_color,
  # cluster = RColorBrewer::brewer.pal(nlevels(df_sil1$cluster), 'Spectral') %>% setNames(levels(df_sil1$cluster)),
  cluster = c('#CA1D1F', '#F4AA63', '#006835', '#7E2F8E', '#2E7CAF') %>% setNames(levels(df_sil1$cluster))#,
  #Type = c(up = 'red3', down = 'blue3')
)

ph_spe_prot <- pheatmap::pheatmap(
  mat_sig[dep_specific_up$Protein, ],
  color = my_colors,
  legend_breaks = quantile(my_breaks),
  breaks = my_breaks,
  scale = 'row',
  annotation_row = dep_specific_up %>%
    # as.data.frame() %>% column_to_rownames('Protein') %>% 
    select(cluster) %>%
    mutate(cluster = factor(as.numeric(cluster), levels = 1:max(as.numeric(cluster)))),
  annotation_col = df_sil1 %>% dplyr::select(cluster, cancer_abbr, tissue_name),
  annotation_colors = ann_colors,
  cluster_rows = F, cluster_cols = F,
  show_rownames = F, show_colnames = F,
  fontsize = 10,
  filename = 'PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_2772specific-up.pdf', width = 10, height = 9
)


rio::export(list(DEP_up = dep_up, DEP_up_specific = dep_specific_up,
                 DEP_down = dep_down, DEP_down_specific = dep_specific_down),
            'PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_DEP_up_down_v2.xlsx')




# enriched proteins
#median matrix
df_median <- 2 ^ mat_sig %>%
  t() %>% as.data.frame() %>% 
  rownames_to_column('FileName') %>% 
  inner_join(df_sil1 %>% select(FileName, cluster), ., by = 'FileName')

mat_median <- df_median %>% select(-FileName) %>% 
  group_by(cluster) %>% 
  summarise_all(median) %>% 
  ungroup()

colnames(mat_median)[1] <- 'Group.1'


clas.all <- ulhen_class(mat_median, fct = 3)
colnames(mat_median)[1] <- colnames(clas.all)[1] <- 'cluster'

df_class <- clas.all %>% pivot_longer(-cluster, names_to = 'Protein', values_to = 'specificity')
df_class$specificity %<>% str_to_sentence()
df_class %>% count(specificity)
#   specificity                  n
# 1 Group enriched  30501
# 2 Not detected        1
# 3 Tissue enriched  1383

df_class %>% filter(specificity == 'Tissue enriched') %>% count(cluster)
#   cluster     n
# 1 1         170
# 2 2         260
# 3 3         504
# 4 4          51
# 5 5         398

df_class_enrich <- df_class %>% filter(specificity == 'Tissue enriched') %>% arrange(cluster)

rio::export(
  list(matrix_median = mat_median, specificity = df_class, cluster_enriched = df_class_enrich),
  'PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_ulhen_v2.xlsx')


list(cluster_enriched = df_class_enrich$Protein, DEP_spec = dep_specific_up$Protein) %>%
  VennDiagram::venn.diagram(category.names = names(.),
                            filename = NULL,disable.logging = T,
                            output=TRUE,
                            cat.cex = 2,
                            cex = 1.5,
                            cat.default.pos = "outer") %>% grid::grid.draw()





dep_specific_up <- readxl::read_excel('PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_DEP_up_down_v2.xlsx', sheet = 'DEP_up_specific')
df_class_enrich <- readxl::read_excel('PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_ulhen_v2.xlsx', sheet = 'cluster_enriched')

# df_spec <- df_class_enrich %>% inner_join(dep_specific_up)
# rio::export(df_spec, 'PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_334signature.xlsx')
# df_spec %>% count(cluster)
# #   cluster     n
# # 1 1          12
# # 2 2          60
# # 3 3           3
# # 4 4          30
# # 5 5          62
# # 6 6           7
# # 7 7          80
# # 8 8          80
# 
# prot_signitures <- df_spec$Protein # 334

df_spec <- df_class_enrich %>% select(cluster, Protein) %>% mutate(Type = 'up-regulated') %>%
  rbind(dep_specific_up %>% select(cluster, Protein) %>% mutate(Type = 'cluster-enriched')) %>%
  distinct(cluster, Protein, .keep_all = T) %>% as.data.frame()
df_spec %>% count(Protein) %>% filter(n > 1) %>% semi_join(df_spec, .) %>% arrange(Protein, cluster)

df_spec <- df_spec %>% count(Protein) %>% filter(n > 1) %>% anti_join(df_spec, .) # delete enriched and DEP differences
rownames(df_spec) <- df_spec$Protein
prot_signitures <- df_spec$Protein # 3014


ann_row <- df_spec %>% 
  arrange(cluster, Type) %>%
  dplyr::select(cluster) %>% 
  mutate(cluster = factor(as.numeric(cluster), levels = 1:max(as.numeric(cluster))))
ann_colors <- list(
  cancer_abbr=cancer_color,
  # cluster = RColorBrewer::brewer.pal(nlevels(df_sil1$cluster), 'Spectral') %>% setNames(levels(df_sil1$cluster)),
  cluster = c('#CA1D1F', '#F4AA63', '#006835', '#7E2F8E', '#2E7CAF') %>% setNames(levels(df_sil1$cluster))#,
  # Type = c(up = 'red3', down = 'blue3')
)



ph_spe_prot <- pheatmap::pheatmap(
  mat_sig[prot_signitures, ],
  color = my_colors,
  legend_breaks = quantile(my_breaks),
  breaks = my_breaks,
  scale = 'row',
  # annotation_col = df_sil1 %>% select(cluster, tissue_name),
  # annotation_row = df_spec %>% select(cluster, Type),
  # annotation_colors = ann_colors,
  annotation_row = dep_specific_up %>%
    as.data.frame() %>% column_to_rownames('Protein') %>% 
    select(cluster) %>%
    mutate(cluster = factor(as.numeric(cluster), levels = 1:max(as.numeric(cluster)))),
  annotation_col = df_sil1 %>% dplyr::select(cluster, cancer_abbr, tissue_name),
  annotation_colors = ann_colors,
  cluster_rows = F, cluster_cols = F,
  show_rownames = F, show_colnames = F,
  fontsize = 10,
  filename = 'PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_3014signatures.pdf'
)

rio::export(df_spec, 'PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_3014signatures.xlsx')
mat_sig %>% as.data.frame() %>% rownames_to_column() %>% rio::export('PTDB_pan-cancer_6377_proteins_1014_coreSamples_matrix.tsv')





pheatmap::pheatmap(
  mat_sig[prot_signitures, ],
  color = my_colors,
  legend_breaks = quantile(my_breaks),
  breaks = my_breaks,
  scale = 'row',
  # annotation_col = df_sil1 %>% select(anatomical_classification, cancer_abbr, cluster) %>% rename(Cancer = cancer_abbr, Cluster = cluster),
  # annotation_row = df_spec %>% select(cluster, Type) %>% rename(Cluster = cluster),
  # annotation_colors = ann_colors,
  annotation_row = ann_row,
  annotation_col = df_sil1 %>% dplyr::select(cluster, cancer_abbr, tissue_name),
  annotation_colors = ann_colors,
  cluster_rows = F, cluster_cols = F,
  show_rownames = F, show_colnames = F,
  fontsize = 10,
  filename = 'PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_3014signatures_v2.pdf'
)



# 5. Post-analyse of molecular subtypes ------------
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)

df_sil1 <- rio::import('PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters.xlsx')

## 5.1 Cancer Type Distribution of Core-samples 
### 5.1.1 signatures distribution ---------
# pie chart
df_pie <- df_sil1 %>%
  mutate(cluster = factor(as.numeric(cluster)),
         cancer_abbr = factor(cancer_abbr)) %>%
  count(cancer_abbr, cluster)


plot_ls <- list()
for(i in seq_along(levels(df_pie$cancer_abbr))){
  x <- levels(df_pie$cancer_abbr)[i]
  tmp <- df_pie %>% filter(cancer_abbr == x)
  tmp2 <- tmp %>% 
    mutate(csum = rev(cumsum(rev(n))), 
           pos = n/2 + lead(csum, 1),
           pos = if_else(is.na(pos), n/2, pos))
  
  plot_ls[[i]] <- ggplot(tmp, aes(x = "" , y = n, fill = cluster)) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y", start = 0) +
    # scale_fill_brewer(palette = "Set2") +
    scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(nlevels(df_pie$cluster)) %>% setNames(levels(df_pie$cluster))) +
    ggrepel::geom_label_repel(data = tmp2,
                              aes(y = pos, label = paste0(round(100*n/sum(n), 1), "%")),
                              size = 4.5, nudge_x = 1, show.legend = FALSE) +
    guides(fill = guide_legend(title = "Group")) +
    theme_void()+
    labs(title = x)
}
length(plot_ls) # 30
p <- ggpubr::ggarrange(plotlist = plot_ls, nrow = 6, ncol = 5)
ggsave('PTDB_pancancer_6377_proteins_1014_coreSamples_5clusters_ratio_pie.pdf', p, width = 3*5, height = 3*6)

#20251024改成等长bar plot

# 计算百分比数据
df_bar <- df_sil1 %>%
  mutate(cluster = factor(as.numeric(cluster)),
         cancer_abbr = factor(cancer_abbr)) %>%
  count(cancer_abbr, cluster) %>%
  group_by(cancer_abbr) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()

# 创建颜色调色板

df_bar$cluster<-paste0("P",df_bar$cluster)

df_bar$cancer_abbr<-factor(df_bar$cancer_abbr,levels = c("ENCA","CESC","TOCA","THCA","LARCA","TGCT","GC","GIST","PRCA","PECA","PACA","CCOC",
                                                         "HGSOC","MUT","BRCA-HER2+","BRCA-LumB-HER2+","BRCA-LumA","BRCA-LumB-HER2-","BRCA-TNBC",
                                                         "LUCA","THYM","DLBCL","HCC","READ","COCA","RC","GBCA","FTCA","ESCA","GBM"))





cluster_colors <- c('#CA1D1F', '#F4AA63', '#006835', '#7E2F8E', '#2E7CAF') %>% 
  setNames(levels(df_bar$cluster))

# 创建堆叠百分比条形图
p <- ggplot(df_bar, aes(x = percentage, y = cancer_abbr, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "Cluster Distribution by Cancer Type",
    x = "Percentage",
    y = "Cancer Type",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

# 显示图形
print(p)

ggsave("20251028_PTDB_pancancer_6377_proteins_1014_coreSamples_5clusters_ratio_bar_V2.pdf",p,width = 6,height = 6)




#######################################

### 5.1.2 cancer distribution --------
#20251024 颜色改成指定颜色
color_mat<-rio::import('./input/labels/PUH_cancer_colorset_20251009.xlsx')
color_mat$color <- str_c('#', color_mat$color) %>% setNames(color_mat$cancer_abbr)
df_pie <- df_sil1 %>%
  mutate(cluster = factor(as.numeric(cluster)),
         cancer_abbr = factor(cancer_abbr)) %>%
  count(cancer_abbr, cluster)
df_pie$color<-color_mat$color[match(df_pie$cancer_abbr,color_mat$cancer_abbr)]

plot_ls <- list()
for(i in seq_along(levels(df_pie$cluster))){
  x <- levels(df_pie$cluster)[i]
  tmp <- df_pie %>% filter(cluster == x)
  tmp2 <- tmp %>% 
    mutate(csum = rev(cumsum(rev(n))), 
           pos = n/2 + lead(csum, 1),
           pos = if_else(is.na(pos), n/2, pos))
  
  plot_ls[[i]] <- ggplot(tmp, aes(x = "" , y = n, fill = cancer_abbr)) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y", start = 0) +
    # scale_fill_brewer(palette = "Set2") +
    scale_fill_manual(values = df_pie$color[df_pie$cluster==x]) +
    # ggrepel::geom_label_repel(data = tmp2,
    #                           aes(y = pos, label = paste0(round(100*n/sum(n), 1), "%")),
    #                           size = 4.5, nudge_x = 1, show.legend = FALSE,max.overlaps = Inf ) +
    guides(fill = guide_legend(title = "Group")) +
    theme_void()+
    labs(title = str_c('P', x))
}

length(plot_ls) # 5
p <- ggpubr::ggarrange(plotlist = plot_ls, nrow = 2, ncol = 3)
ggsave('20251024_PTDB_pancancer_6377_proteins_1014_coreSamples_5clusters_cancer_ratio_pie_del_label_v2.pdf', p, width = 5*3, height = 5*2)

df_pie %>% rio::export('PTDB_pancancer_6377_proteins_1014_coreSamples_5clusters_cancer_ratio_pie.xlsx')


sum_n <- df_pie %>% group_by(cluster) %>% summarise(x = sum(n)) %>% pull(x)
df_pie %>%
  group_by(cluster) %>%
  arrange(desc(n), .by_group = T) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(`ratio (%)` = 100 * n / sum_n) %>% 
  arrange(desc(`ratio (%)`))
# # A tibble: 4 × 4
# cancer_abbr cluster     n `ratio (%)`
# <fct>       <fct>   <int>       <dbl>
# 1 RC          5          47        51.6
# 2 GIST        3          30        44.1
# 3 CESC        1          58        32.6
# 4 COCA        2          56        17.9
# 5 GC          4          58        15.9

### 5.1.3 Cancer Hallmark --------
#### a. select signatures --------
mat_sig <- mat_consensus[, df_sil1$FileName]
df_sil1$cluster <- factor(str_c('P', df_sil1$cluster), levels = str_c('P', sort(as.numeric(unique(df_sil1$cluster)))))
df_sil1 %<>% mutate_at(vars(neighbor), factor)
df_sil1 %<>% mutate(cancer. = str_c(cancer_abbr, " - ", cancer), .after = cancer_abbr)
rownames(df_sil1) <- df_sil1$FileName
# saveRDS(df_sil1, 'df_sil1.rds')
df_sil1 <- readRDS('df_sil1.rds')
cls <- unique(df_sil1$cluster)

df_dea_ls <- list()
for(i in 1:nrow(mat_sig)){ # for each protein
  cat(i, '...\r')
  log2FC <- P <- c()
  for(ii in seq_along(cls)){ # for each cluster
    cl <- cls[ii]
    fnames1 <- df_sil1 %>% filter(cluster == cl) %>% pull(FileName)
    fnames2 <- setdiff(df_sil1$FileName, fnames1)
    x1 <- mat_sig[i, fnames1]
    x2 <- mat_sig[i, fnames2]
    
    # fold-change
    log2FC[ii] <- log2(mean(2^x1) / mean(2^x2))
    
    # wilcoxon rank test
    P[ii] <- wilcox.test(x1, x2, alternative = 'two.sided', paired = F)$p.value
  }
  
  adj.P <- p.adjust(P, 'BH')
  names(log2FC) <- str_c('log2FC_', cls, 'vs.')
  names(P) <- str_c('P_', cls, 'vs.')
  names(adj.P) <- str_c('adj.P_', cls, 'vs.')
  df_tmp <- c(log2FC, P, adj.P) %>%
    as.data.frame() %>%
    setNames(rownames(mat_sig)[i]) %>%
    t() %>% as.data.frame()
  
  df_dea_ls[[i]] <- df_tmp
}
names(df_dea_ls) <- rownames(mat_sig)
df_dea <- plyr::ldply(df_dea_ls, .id = 'Protein') %>% mutate(Protein = as.character(Protein))

df_dea_ <- dfprot %>% rename(Protein = Protein.Group) %>% right_join(df_dea)
rio::export(df_dea_, 'PTDB_pan-cancer_6377_proteins_1014_core-samples_5clusters_DEA.xlsx')

df_dea_signif3 <- df_dea_ %>% mutate_if(is.numeric, signif, digits = 3)
rio::export(df_dea_signif3, 'PTDB_pan-cancer_6377_proteins_1014_core-samples_5clusters_DEA_signif3.xlsx')


# DEP cutoff: |FC| >= 1.5, adj.P < 0.05
#DEP_up
prot_up_ls <- list()
for(i in seq_along(cls)){
  cl <- cls[i]
  df_tmp <- df_dea %>% select(Protein, matches(str_c(as.character(cl), 'vs.')))
  # isSignificant <- abs(df_tmp$log2FC) >= 2 & df_tmp$adj.P < 0.05
  isUp <- df_tmp$log2FC >= 1.5 & df_tmp$adj.P < 0.05 # select only up-regulated proteins
  prot_up_ls[[i]] <- df_tmp$Protein[isUp]
}

prot_specific_up_ls <- lapply(seq_along(prot_up_ls), function(i){
  setdiff(prot_up_ls[[i]], unlist(prot_up_ls[-i]))
})

names(prot_up_ls) <- names(prot_specific_up_ls) <- cls
dep_up <- plyr::ldply(prot_up_ls, .id = 'cluster', function(x){
  x %>% as.data.frame() %>% setNames('Protein')
})
dep_up %<>% left_join(df_dea_, by = 'Protein')

sapply(prot_up_ls, length)
# P1   P2   P3   P4   P5 
# 293  782 1003  224  606 

sapply(prot_specific_up_ls, length)
# P1  P2  P3  P4  P5 
# 280 770 939 222 561 

dep_specific_up <- plyr::ldply(prot_specific_up_ls, .id = 'cluster', function(x){
  x %>% as.data.frame() %>% setNames('Protein')
})
dep_specific_up %<>% left_join(df_dea_, by = 'Protein')
rownames(dep_specific_up) <- dep_specific_up$Protein
dim(dep_specific_up) # 2772  18


#DEP_down
prot_down_ls <- list()
for(i in seq_along(cls)){
  cl <- cls[i]
  df_tmp <- df_dea %>% select(Protein, matches(str_c(as.character(cl), 'vs.')))
  # isSignificant <- abs(df_tmp$log2FC) >= 2 & df_tmp$adj.P < 0.05
  isDown <- df_tmp$log2FC <= -1.5 & df_tmp$adj.P < 0.05 # select only down-regulated proteins
  prot_down_ls[[i]] <- df_tmp$Protein[isDown]
}

prot_specific_down_ls <- lapply(seq_along(prot_down_ls), function(i){
  setdiff(prot_down_ls[[i]], unlist(prot_down_ls[-i]))
})

names(prot_down_ls) <- names(prot_specific_down_ls) <- cls
dep_down <- plyr::ldply(prot_down_ls, .id = 'cluster', function(x){
  x %>% as.data.frame() %>% setNames('Protein')
})
dep_down %<>% left_join(df_dea_, by = 'Protein')

sapply(prot_down_ls, length)
# P1   P2   P3   P4   P5 
# 975  579 1373 1050 1406 

sapply(prot_specific_down_ls, length)
# P1  P2  P3  P4  P5 
# 469 322 485 355 493 

dep_specific_down <- plyr::ldply(prot_specific_down_ls, .id = 'cluster', function(x){
  x %>% as.data.frame() %>% setNames('Protein')
})
dep_specific_down %<>% left_join(df_dea_, by = 'Protein')
rownames(dep_specific_down) <- dep_specific_down$Protein
dim(dep_specific_down) # 2124   18



# heatmap -- only upregulated proteins
#breaks
quantile(scale(mat_sig[dep_specific_up$Protein, ]))
# my_breaks <- unique(c(seq(-1.5, 0, length.out = 251), 0, seq(0, 3, length.out = 251)))
my_breaks <- seq(-1.5, 1.5, by = 0.01)

#colors
my_colors <- rev(c(colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[1:5], 'white'))(length(my_breaks)/2),
                   colorRampPalette(colors = c('white', RColorBrewer::brewer.pal(11, "PiYG")[7:11]))(length(my_breaks)/2)))

ph_spe_prot <- pheatmap::pheatmap(
  mat_sig[dep_specific_up$Protein, ],
  color = my_colors,
  legend_breaks = quantile(my_breaks),
  breaks = my_breaks,
  scale = 'row',
  annotation_col = df_sil1 %>% column_to_rownames('FileName') %>% select(cluster, cancer_abbr),
  annotation_row = dep_specific_up %>% select(cluster),
  # annotation_colors = ann_colors,
  cluster_rows = F, cluster_cols = F,
  show_rownames = F, show_colnames = F,
  fontsize = 10,
  filename = 'PTDB_pan-cancer_6377_proteins_1014_core-samples_5clusters_2772specific-up.pdf'
)


rio::export(list(DEP_up = dep_up, DEP_up_specific = dep_specific_up,
                 DEP_down = dep_down, DEP_down_specific = dep_specific_down),
            'PTDB_pan-cancer_6377_proteins_1014_core-samples_5clusters_DEP_up_down_v2.xlsx')


# enriched proteins
#median matrix
df_median <- 2 ^ mat_sig %>%
  t() %>% as.data.frame() %>% 
  rownames_to_column('FileName') %>% 
  inner_join(df_sil1 %>% select(FileName, cluster), ., by = 'FileName')

mat_median <- df_median %>% select(-FileName) %>% 
  group_by(cluster) %>% 
  summarise_all(median) %>% 
  ungroup()

colnames(mat_median)[1] <- 'Group.1'


clas.all <- ulhen_class(mat_median, fct = 3)
colnames(mat_median)[1] <- colnames(clas.all)[1] <- 'cluster'

df_class <- clas.all %>% pivot_longer(-cluster, names_to = 'Protein', values_to = 'specificity')
df_class$specificity %<>% str_to_sentence()
df_class %>% count(specificity)
#   specificity                  n
# 1 Group enriched  30501
# 2 Not detected        1
# 3 Tissue enriched  1383

df_class %>% filter(specificity == 'Tissue enriched') %>% count(cluster)
#   cluster     n
# 1 P1        170
# 2 P2        260
# 3 P3        504
# 4 P4         51
# 5 P5        398

df_class_enrich <- df_class %>% filter(specificity == 'Tissue enriched') %>% arrange(cluster)

rio::export(
  list(matrix_median = mat_median, specificity = df_class, cluster_enriched = df_class_enrich),
  'PTDB_pan-cancer_6377_proteins_1014_core-samples_5clusters_ulhen_v2.xlsx')


list(cluster_enriched = df_class_enrich$Protein, DEP_spec = dep_specific_up$Protein) %>%
  VennDiagram::venn.diagram(category.names = names(.),
                            filename = NULL,disable.logging = T,
                            output=TRUE,
                            cat.cex = 2,
                            cex = 1.5,
                            cat.default.pos = "outer") %>% grid::grid.draw()

#### b. signatures expression --------
dep_specific_up <- readxl::read_excel('PTDB_pan-cancer_6377_proteins_1014_core-samples_5clusters_DEP_up_down_v2.xlsx', sheet = 'DEP_up_specific')
df_class_enrich <- readxl::read_excel('PTDB_pan-cancer_6377_proteins_1014_core-samples_5clusters_ulhen_v2.xlsx', sheet = 'cluster_enriched')

df_spec <- df_class_enrich %>% select(cluster, Protein) %>% mutate(Type = 'up-regulated') %>%
  rbind(dep_specific_up %>% select(cluster, Protein) %>% mutate(Type = 'cluster-enriched')) %>%
  distinct(cluster, Protein, .keep_all = T) %>% as.data.frame()
df_spec %>% count(Protein) %>% filter(n > 1) %>% semi_join(df_spec, .) %>% arrange(Protein, cluster)

df_spec <- df_spec %>% count(Protein) %>% filter(n > 1) %>% anti_join(df_spec, .) # delete enriched and DEP differences
rownames(df_spec) <- df_spec$Protein
prot_signitures <- df_spec$Protein # 3014


ann_colors <- list(
  cancer_abbr=cancer_color,
  # cluster = RColorBrewer::brewer.pal(nlevels(df_sil1$cluster), 'Spectral') %>% setNames(levels(df_sil1$cluster)),
  cluster = c('#CA1D1F', '#F4AA63', '#006835', '#7E2F8E', '#2E7CAF') %>% setNames(levels(df_sil1$cluster))#,
  # Type = c(up = 'red3', down = 'blue3')
)
ph_spe_prot <- pheatmap::pheatmap(
  mat_sig[prot_signitures, ],
  color = my_colors,
  legend_breaks = quantile(my_breaks),
  breaks = my_breaks,
  scale = 'row',
  annotation_col = df_sil1 %>% column_to_rownames('FileName') %>% select(cluster, cancer_abbr),
  annotation_row = df_spec %>% select(cluster, Type),
  annotation_colors = ann_colors,
  cluster_rows = F, cluster_cols = F,
  show_rownames = F, show_colnames = F,
  fontsize = 10,
  filename = 'PTDB_pan-cancer_6377_proteins_1014_core-samples_5clusters_3014signatures.pdf'
)

rio::export(df_spec, 'PTDB_pan-cancer_6377_proteins_1014_core-samples_5clusters_3014signatures.xlsx')
mat_sig %>% as.data.frame() %>% rownames_to_column() %>% rio::export('PTDB_pan-cancer_6377_proteins_1014_core-samples_matrix.tsv')


# df_color <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/input/PUH_tissue_colorset_20230210.xlsx')
# cancer_color <- df_sil1 %>%
#   distinct(class_abbr, cancer_abbr) %>%
#   rename(tissue = class_abbr) %>% 
#   full_join(df_color)
# rio::export(cancer_color, 'PUH_cancer_colorset.xlsx')
# dfc_color <- rio::import('PUH_cancer_colorset_20251009.xlsx')
# cancer_color <- str_c('#', dfc_color$color) %>% setNames(dfc_color$cancer_abbr)
# df_sil1$cancer_abbr %>% setdiff(names(cancer_color)) # character(0)

ann_colors <- list(
  cancer_abbr=cancer_color,
  # cluster = RColorBrewer::brewer.pal(nlevels(df_sil1$cluster), 'Spectral') %>% setNames(levels(df_sil1$cluster)),
  cluster = c('#CA1D1F', '#F4AA63', '#006835', '#7E2F8E', '#2E7CAF') %>% setNames(levels(df_sil1$cluster))#,
  # Type = c(up = 'red3', down = 'blue3')
)

pheatmap::pheatmap(
  mat_sig[prot_signitures, ],
  color = my_colors,
  legend_breaks = quantile(my_breaks),
  breaks = my_breaks,
  scale = 'row',
  annotation_col = df_sil1 %>% column_to_rownames('FileName') %>% select(cancer_abbr, cluster),
  annotation_row = df_spec %>% select(cluster, Type),
  annotation_colors = ann_colors,
  cluster_rows = F, cluster_cols = F,
  show_rownames = F, show_colnames = F,
  fontsize = 10,
  filename = 'PTDB_pan-cancer_6377_proteins_1014_core-samples_5clusters_3014signatures_v2.pdf'
)

#### c. ssGSEA with MSigDB **H** datasets -----
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)

mat_sig <- rio::import('PTDB_pan-cancer_6377_proteins_1014_core-samples_matrix.tsv') %>% column_to_rownames() %>% as.matrix()
df_sil1 <- readRDS('df_sil1.rds')
rownames(df_sil1) <- df_sil1$FileName

dep_specific_up <- readxl::read_excel('PTDB_pan-cancer_6377_proteins_1014_core-samples_5clusters_DEP_up_down_v2.xlsx', sheet = 'DEP_up_specific')
dep_specific_down <- readxl::read_excel('PTDB_pan-cancer_6377_proteins_1014_core-samples_5clusters_DEP_up_down_v2.xlsx', sheet = 'DEP_down_specific')

dep_specific <- list(up = dep_specific_up, down = dep_specific_down) %>%
  plyr::ldply(.id = 'Type')

dim(dep_specific) # 4896   19
length(unique(dep_specific$Protein)) # 3873
dep_specific1 <- dep_specific %>% distinct(Protein, .keep_all = T) %>% arrange(cluster)
dim(dep_specific1) # 3873   19
rio::export(dep_specific1, 'PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_3873signatures_up_down.xlsx')

dep_specific2 <- dep_specific1 %>%
  mutate(log2FC_major = str_glue("log2FC_{cluster}vs."),
         adj.P_major = str_glue("adj.P_{cluster}vs."),
         .after = Genes)
dep_specific2$log2FC_major <- sapply(seq_along(dep_specific2$log2FC_major), function(x){
  y <- dep_specific2$log2FC_major[x]
  dep_specific2[[y]][x]
})
dep_specific2$adj.P_major <- sapply(seq_along(dep_specific2$adj.P_major), function(x){
  y <- dep_specific2$adj.P_major[x]
  dep_specific2[[y]][x]
})

mat_specg <- mat_specp <- mat_sig[dep_specific1$Protein, ] # g: gene matrix; p: protein matrix
rownames(mat_specg) <- dep_specific1$Gene
dim(mat_specp) # 3873  1014
dim(df_sil1) # 1014 21

library(GSVA)
library(GSEABase)


##### Load GSEA files. ###
## Downloaded from http://www.gsea-msigdb.org/gsea/msigdb/index.jsp ##
#  human： http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
#  mouse： http://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp
# h1gmt <- getGmt("h.all.v7.4.symbols.gmt") # hallmark
h1gmt <- getGmt("h.all.v2025.1.Hs.symbols.gmt") # gmt dataset
tmp <- sapply(seq_along(h1gmt@.Data), function(i) h1gmt@.Data[[i]]@geneIds)
names(tmp) <- names(h1gmt)
unlist(tmp) %>% unique() %>% length() # 4384
unlist(tmp) %>% intersect(rownames(mat_specg)) %>% length() # 1032

ssgsea <- GSVA::gsva(
  param = ssgseaParam(exprData = mat_specg, geneSets = h1gmt),#, minSize = 15, maxSize = 500
  # expr = mat_specg,
  # gset.idx.list = h1gmt,# min.sz=15, max.sz=500,
  # method = 'ssgsea',
  # kcdf = 'Gaussian', abs.ranking = T
)

# save(ssgsea, file = 'ssgsea_pancancer_hallmark.RData')
load('ssgsea_pancancer_hallmark.RData')


# min-max normalization
ssgsea.1 <- ssgsea
for (i in colnames(ssgsea)) {
  #i <- colnames(ssgsea)[1]
  ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
  
}

pheatmap::pheatmap(
  ssgsea.1,
  annotation_col = df_sil1 %>% dplyr::select(cluster, cancer_abbr),
  annotation_colors = ann_colors,
  scale = "none",
  cluster_rows = T,
  cluster_cols = F,
  show_rownames = F,
  show_colnames = F,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # cellwidth = 10, cellheight = 15,
  fontsize = 10
)
graphics.off()

#### d. limma lmfit -----
library(limma)


# differential expression
group <- df_sil1$cluster
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(ssgsea)
fit <- limma::lmFit(ssgsea, design)

lev <- levels(group)
cts <- sapply(lev, function(g) paste0(g, " - (", paste(setdiff(lev, g), collapse = " + "), ")/", length(lev) - 1), USE.NAMES = FALSE)
fr <- limma::eBayes(limma::contrasts.fit(fit, do.call(limma::makeContrasts, c(as.list(cts), list(levels = design)))), robust = TRUE)

# CUTM <- sapply(fr$pv, function(v){
#   t <- apply(fr$p.value, 1, function(x) sum(x < 0.05))
#   rownames(v)[t == ncol(v)]
# })
# names(CUTM) <- levels(group)

q <- apply(fr$p.value, 2, p.adjust, method = "BH") %>% set_colnames(., str_extract(colnames(.), '^\\w+'))
sel <- rowSums(q <= 0.05) > 0

mat <- fr$coefficients[sel, , drop = FALSE] %>% set_colnames(., str_extract(colnames(.), '^\\w+'))
dir.create("results", showWarnings = FALSE, recursive = TRUE)
pheatmap::pheatmap(
  mat,
  annotation_col   = df_sil1 %>% distinct(cluster) %>% set_rownames(.$cluster),
  annotation_colors= ann_colors,    # optional
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = T,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize = 7,
  filename = "heatmap_vs_rest_effect.pdf",
  width = 7, height = 5
)

write.csv(cbind(pathway = rownames(mat), as.data.frame(mat, check.names = FALSE)), "vs_rest_effect_matrix.csv", row.names = FALSE)


# set groups
group <- df_sil1$cluster
Diff <- list()
for(i in seq_along(levels(group))){
  x <- levels(group)[i]
  group_this <- as.character(group)
  group_this[group_this != x] <- 'rest'
  group_this %<>% factor()
  
  design <- model.matrix(~0+group_this)
  colnames(design) <- levels(group_this)
  rownames(design) <- colnames(ssgsea)
  
  compare <- makeContrasts(str_c(colnames(design), collapse = " - "), levels = design)
  fit <- lmFit(ssgsea, design)
  fit2 <- contrasts.fit(fit, compare)
  fit3 <- eBayes(fit2)
  # Diff[[i]] <- topTable(fit3, coef=1, number=Inf) %>% setNames(str_c(colnames(.), '_', str_c(colnames(design), collapse = "_"))) %>% rownames_to_column('pathway')
  Diff[[i]] <- topTable(fit3, coef=1, number=Inf) %>% rownames_to_column('pathway')
}
names(Diff) <- levels(group)
Diff <- plyr::ldply(Diff, .id = 'cluster')

Diff %>% rio::export('pancancer_ssGSEA_hallmark_limmaFit.xlsx')
Diff %>% dplyr::filter(adj.P.Val < 0.05) %>% View()


h1.gs.df <- plyr::ldply(tmp, .id = 'Pathway', function(x){
  data.frame(Pathway.Genes.All = str_c(x, collapse = ';'),
             Genes = str_c(intersect(x, dep_specific2$Genes), collapse = ';')) %>% 
    mutate(Geneset.Size = str_count(Genes, ';') + 1)
}) %>% separate_rows(Genes) %>%
  inner_join(dep_specific2) %>% 
  inner_join(Diff %>% dplyr::rename(Pathway = pathway)) %>% 
  arrange(cluster, adj.P_major)
rio::export(h1.gs.df, 'PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_3873signatures_hallmark1028.xlsx')



ssgsea_median <- aggregate(t(ssgsea), by = list(df_sil1$cluster), median) %>% dplyr::rename(cluster = Group.1)
v_data <- pivot_longer(ssgsea_median, cols = -cluster, names_to = 'pathway', values_to = 'Median') %>% as.data.frame()
p_data <- Diff %>% dplyr::select(pathway, cluster, adj.P.Val)
data <- v_data %>% full_join(p_data, by = c('cluster', 'pathway'))
# data<-merge(v_data,p_data,by.y=c("pathway","cancer"),by.x=c("variable","Group.1"))
# colnames(data)<-c("Hallmark","Subtype","Median","value")
data %<>%
  dplyr::mutate(text = case_when( 
    adj.P.Val < 0.001 ~  "***", # round() 只保留两位小数
    between(adj.P.Val, 0.001, 0.01)  ~ "**",
    between(adj.P.Val, 0.01, 0.05)  ~ "*",
    T ~ ""))

# data %>% rio::export('pancancer_ssGSEA_hallmark_limmaFit_and_median.xlsx')
data %>% rio::export('pancancer_ssGSEA_hallmark_enrich_score_heatmap.xlsx')
data <- rio::import('pancancer_ssGSEA_hallmark_enrich_score_heatmap.xlsx')

mat_heat <- data %>%
  dplyr::select(pathway, cluster, Median) %>%
  pivot_wider(names_from = 'cluster', values_from = 'Median') %>%
  column_to_rownames('pathway') %>% 
  as.matrix()

a <- pheatmap::pheatmap( # to perform a "hc" tree
  mat_heat,
  scale = "none",
  cluster_rows = T,
  cluster_cols = T,
  show_rownames = T,
  show_colnames = T,
  color =colorRampPalette(c("blue", "white","red"))(100),
  cellwidth = 10, cellheight = 10,
  fontsize = 10,
  filename = 'pancancer_ssGSEA_hallmark_enrich_score_heatmap.pdf',
  width = 8, height = 8
)

CUTM <- fr$p.value %>% set_colnames(str_extract(colnames(.), '^\\w+')) %>% 
  as.data.frame() %>% rownames_to_column('Pathway') %>% 
  pivot_longer(cols = -Pathway) %>% filter(value < 0.05) %>% 
  plyr::dlply('name', function(dfsub) dfsub$Pathway) %>% 
  .[lev]

data$refine_text <- apply(data, 1, function(x){
  ifelse(x['pathway'] %in% CUTM[[x['cluster']]],
         x['text'], '')
})

data$new_Var1 <- factor(data$cluster, levels = a$tree_col$labels[a$tree_col$order])
data$new_Var2 <- factor(data$pathway, levels = a$tree_row$labels[a$tree_row$order])



p<-ggplot(data, aes(new_Var1, new_Var2)) + 
  geom_tile(aes(fill = Median), colour = "grey", size = 0.5)+
  scale_fill_gradient2(low = "#5C55AF",mid = "white",high = "#EA2E2D") + 
  geom_text(aes(label=refine_text),col ="black",size = 3) +
  theme_minimal() + # 不要背景
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.ticks=element_blank(),
        axis.title.y=element_blank(), # 去掉 y 轴
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8,vjust =0,color = "black"), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 8,color = "black")) + #调整y轴文字
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","*** p < 0.001","\n\n","Median GSVA scores")) +   # 修改 legend 内容
  scale_x_discrete(position = "top") # 将 X 轴放置在最上面

ggsave('pancancer_ssGSEA_hallmark_enrich_score_heatmap_v2.pdf', p, width = 8, height = 8)


#### e. heatmap --------
# # 两两比较
# pvalue_ls2 <- list()
# design <- model.matrix(~0+group)
# colnames(design) <- levels(group)
# rownames(design) <- colnames(ssgsea)
# 
# #P1
# compare <- makeContrasts(P1-P2,P1-P3,P1-P4,P1-P5,P1-P6,P1-P7,P1-P8, levels = design)
# fit <- lmFit(ssgsea, design)
# fit2 <- contrasts.fit(fit, compare)
# fit3 <- eBayes(fit2)
# pvalue_ls2[[1]] <- fit3$p.value
# 
# #P2
# compare <- makeContrasts(P2-P1,P2-P3,P2-P4,P2-P5,P2-P6,P2-P7,P2-P8, levels = design)
# fit <- lmFit(ssgsea, design)
# fit2 <- contrasts.fit(fit, compare)
# fit3 <- eBayes(fit2)
# pvalue_ls2[[2]] <- fit3$p.value
# 
# #P3
# compare <- makeContrasts(P3-P1,P3-P2,P3-P4,P3-P5,P3-P6,P3-P7,P3-P8, levels = design)
# fit <- lmFit(ssgsea, design)
# fit2 <- contrasts.fit(fit, compare)
# fit3 <- eBayes(fit2)
# pvalue_ls2[[3]] <- fit3$p.value
# 
# #P4
# compare <- makeContrasts(P4-P1,P4-P2,P4-P3,P4-P5,P4-P6,P4-P7,P4-P8, levels = design)
# fit <- lmFit(ssgsea, design)
# fit2 <- contrasts.fit(fit, compare)
# fit3 <- eBayes(fit2)
# pvalue_ls2[[4]] <- fit3$p.value
# 
# #P5
# compare <- makeContrasts(P5-P1,P5-P2,P5-P3,P5-P4,P5-P6,P5-P7,P5-P8, levels = design)
# fit <- lmFit(ssgsea, design)
# fit2 <- contrasts.fit(fit, compare)
# fit3 <- eBayes(fit2)
# pvalue_ls2[[5]] <- fit3$p.value
# 
# #P6
# compare <- makeContrasts(P6-P1,P6-P2,P6-P3,P6-P4,P6-P5,P6-P7,P6-P8, levels = design)
# fit <- lmFit(ssgsea, design)
# fit2 <- contrasts.fit(fit, compare)
# fit3 <- eBayes(fit2)
# pvalue_ls2[[6]] <- fit3$p.value
# 
# #P7
# compare <- makeContrasts(P7-P1,P7-P2,P7-P3,P7-P4,P7-P5,P7-P6,P7-P8, levels = design)
# fit <- lmFit(ssgsea, design)
# fit2 <- contrasts.fit(fit, compare)
# fit3 <- eBayes(fit2)
# pvalue_ls2[[7]] <- fit3$p.value
# 
# #P8
# compare <- makeContrasts(P8-P1,P8-P2,P8-P3,P8-P4,P8-P5,P8-P6,P8-P7, levels = design)
# fit <- lmFit(ssgsea, design)
# fit2 <- contrasts.fit(fit, compare)
# fit3 <- eBayes(fit2)
# pvalue_ls2[[8]] <- fit3$p.value
# 



# combine reslts of `data` and `CUTM`
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)



df_sil1$cluster %<>% factor(levels = levels(data$new_Var1))
df_sil1 %<>% arrange(cluster)

mat_heat <- ssgsea[levels(data$new_Var2), df_sil1$FileName]
pathway_sig <- data %>% filter(refine_text != '') %>% pull(pathway) # only significant pathways
mat_heat <- mat_heat[rownames(mat_heat) %in% pathway_sig, ]

#breaks
bk1 <- quantile(scale(mat_heat), 0.1)
bk2 <- quantile(scale(mat_heat), 0.9)
bk1 <- ifelse(floor(bk1) + 0.5 <= bk1, floor(bk1) + 0.5, floor(bk1))
bk2 <- ifelse(ceiling(bk2) - 0.5 >= bk2, ceiling(bk2) - 0.5, ceiling(bk2))

# my_breaks <- unique(c(seq(-1.5, 0, length.out = 251), 0, seq(0, 3, length.out = 251)))
# my_breaks <- seq(bk1, bk2, by = 0.01)
my_breaks <- seq(-2, 2, by = 0.01)

#colors
my_colors <- rev(c(colorRampPalette(colors = c('red4', 'white'))(length(my_breaks)/2),
                   colorRampPalette(colors = c('white', 'blue4'))(length(my_breaks)/2)))
# my_colors <- rev(c(colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "PiYG")[1:5], 'white'))(length(my_breaks)/2),
#                    colorRampPalette(colors = c('white', RColorBrewer::brewer.pal(11, "PiYG")[7:11]))(length(my_breaks)/2)))



b <- pheatmap::pheatmap( # to perform a "hc" tree
  mat_heat,
  annotation_col = df_sil1 %>% select(cancer_abbr, cluster),
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  # color =colorRampPalette(c("blue", "white","red"))(100),
  color = my_colors, breaks = my_breaks,
  # cellwidth = 10, cellheight = 10,
  fontsize = 10,
  filename = 'pancancer_ssGSEA_hallmark_enrich_score_heatmap_sampleLevel.pdf',
  width = 8, height = 8
)


b <- pheatmap::pheatmap( # to perform a "hc" tree
  mat_heat,
  annotation_col = df_sil1 %>% select(cluster),
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = T,
  cluster_cols = T,
  show_rownames = T,
  show_colnames = F,
  # color =colorRampPalette(c("blue", "white","red"))(100),
  color = my_colors, breaks = my_breaks,
  # cellwidth = 10, cellheight = 10,
  fontsize = 10,
  filename = 'pancancer_ssGSEA_hallmark_enrich_score_heatmap_sampleLevel_v3.pdf',
  width = 8, height = 8
)

b <- pheatmap::pheatmap( # to perform a "hc" tree
  mat_heat,
  annotation_col = df_sil1 %>% select(cluster),
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  # color =colorRampPalette(c("blue", "white","red"))(100),
  color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(df_sil1$cluster)[levels(df_sil1$cluster)])[-length(cumsum(table(df_sil1$cluster)))],
  # cellwidth = 10, cellheight = 10,
  fontsize = 10,
  filename = 'pancancer_ssGSEA_hallmark_enrich_score_heatmap_sampleLevel_v2.pdf',
  width = 8, height = 8
)

b <- pheatmap::pheatmap( # to perform a "hc" tree
  mat_heat,
  annotation_col = df_sil1 %>% select(cluster),
  annotation_colors = ann_colors,
  scale = "column",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  # color =colorRampPalette(c("blue", "white","red"))(100),
  color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(df_sil1$cluster)[levels(df_sil1$cluster)])[-length(cumsum(table(df_sil1$cluster)))],
  # cellwidth = 10, cellheight = 10,
  fontsize = 10,
  filename = 'pancancer_ssGSEA_hallmark_enrich_score_heatmap_sampleLevel_v4.pdf',
  width = 8, height = 8
)

#breaks
quantile(mat_heat)
bk1 <- quantile(mat_heat, 0.1)
bk2 <- quantile(mat_heat, 0.9)

my_breaks <- seq(-0.2, 0.4, by = 0.01)

#colors
# my_colors <- rev(c(colorRampPalette(colors = c('red4', 'white'))(length(my_breaks)/2),
#                    colorRampPalette(colors = c('white', 'blue4'))(length(my_breaks)/2)))
my_colors <- rev(colorRampPalette(colors = c('red4', 'white', 'blue4'))(length(my_breaks)))

b <- pheatmap::pheatmap( # to perform a "hc" tree
  mat_heat,
  annotation_col = df_sil1 %>% select(cluster),
  annotation_colors = ann_colors,
  scale = "column",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  # color =colorRampPalette(c("blue", "white","red"))(100),
  color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(df_sil1$cluster)[levels(df_sil1$cluster)])[-length(cumsum(table(df_sil1$cluster)))],
  # cellwidth = 10, cellheight = 10,
  fontsize = 10,
  filename = 'pancancer_ssGSEA_hallmark_enrich_score_heatmap_sampleLevel_v4.pdf',
  width = 8, height = 8
)

b <- pheatmap::pheatmap( # to perform a "hc" tree
  mat_heat,
  annotation_col = df_sil1 %>% select(cancer_abbr, cluster),
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  # color =colorRampPalette(c("blue", "white","red"))(100),
  color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(df_sil1$cluster)[levels(df_sil1$cluster)])[-length(cumsum(table(df_sil1$cluster)))],
  # cellwidth = 10, cellheight = 10,
  fontsize = 10,
  filename = 'pancancer_ssGSEA_hallmark_enrich_score_heatmap_sampleLevel_v5.pdf',
  width = 11, height = 8
)

my_breaks <- seq(-3, 3, by = 0.01)
my_colors <- rev(colorRampPalette(colors = c('red', 'white', 'blue'))(length(my_breaks)))
pheatmap::pheatmap( # to perform a "hc" tree
  mat_heat[nrow(mat_heat):1, ],
  annotation_col = df_sil1 %>% select(cancer_abbr, cluster),
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = T,
  cluster_cols = F,
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  # color =colorRampPalette(c("blue", "white","red"))(100),
  color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(df_sil1$cluster)[levels(df_sil1$cluster)])[-length(cumsum(table(df_sil1$cluster)))],
  # cellwidth = 10, cellheight = 10,
  fontsize = 7,
  main = 'Pan-cancer ssGSEA hallmark (Z-score)',
  filename = 'pancancer_ssGSEA_hallmark_enrich_score_heatmap_sampleLevel_v6.pdf',
  width = 11, height = 5
)

cancer._color <- as.data.frame(cancer_color) %>%
  rownames_to_column('cancer_abbr') %>% 
  inner_join(df_sil1 %>% distinct(cancer_abbr, cancer.)) %>% 
  pull(cancer_color, cancer.)

ann_colors <- list(
  cancer. = cancer._color,
  # cancer_abbr = cancer_color,
  # cluster = RColorBrewer::brewer.pal(nlevels(df_sil1$cluster), 'Spectral') %>% setNames(levels(df_sil1$cluster)),
  cluster = c('#CA1D1F', '#F4AA63', '#006835', '#7E2F8E', '#2E7CAF') %>% setNames(levels(df_sil1$cluster))#,
  # Type = c(up = 'red3', down = 'blue3')
)

pheatmap::pheatmap( # to perform a "hc" tree
  mat_heat[nrow(mat_heat):1, ] %>% set_rownames(str_remove(rownames(.), '^HALLMARK_')),
  annotation_col = df_sil1 %>% select(cancer., cluster),
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = T,
  cluster_cols = F,
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  # color =colorRampPalette(c("blue", "white","red"))(100),
  color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(df_sil1$cluster)[levels(df_sil1$cluster)])[-length(cumsum(table(df_sil1$cluster)))],
  # cellwidth = 10, cellheight = 10,
  fontsize = 7,
  main = 'Pan-cancer ssGSEA hallmark (Z-score)',
  filename = 'pancancer_ssGSEA_hallmark_enrich_score_heatmap_sampleLevel_v7.pdf',
  width = 9, height = 5
)

### 5.1.4 immune types --------
#### a.LM22 -----
##### ssGSEA --------------
library(magrittr);library(tidyverse);library(GSVA);library(GSEABase)
# with CIBERSORT'S LM22 gene set（547 genes）
# saveRDS(df_sil1, 'df_sil1_5a.rds')
# df_sil1 <- readRDS('df_sil1_5a.rds')
dim(mat_specg) # 3873  1014
dim(df_sil1) # 1014  21

# BiocManager::install("mastR")
# library(mastR)
LM22 <- readxl::read_excel('source/41592_2015_BFnmeth3337_MOESM207_ESM.xls', sheet = 'SuppTable1_DEGs') %>%  as.data.frame()
colnames(LM22) <- LM22[2, ]
LM22 <- LM22[-(1:2), ]
colnames(LM22)[1] <- 'Gene'
rownames(LM22) <- NULL
df_LM22_pansig <- LM22

LM22_genes <- LM22$Gene
length(intersect(LM22_genes, rownames(mat_specg))) # 139

LM22 <- LM22 %>% column_to_rownames('Gene') %>% t() %>% as.data.frame() %>% rownames_to_column('LM22')

gs_LM22 <- plyr::dlply(LM22, 'LM22', function(dfsub){
  dfsub[-1] %>% unlist() %>% .[. == 1] %>% names()
})
# save(gs_LM22, file = 'source/gene_set_LM22.Rdata')
# load('source/gene_set_LM22.Rdata') # gs_LM22


# overlapped genes between LM22 and PUH pan-cancer signatures
rownames(df_LM22_pansig) <- df_LM22_pansig$Gene
df_LM22_pansig <- df_LM22_pansig[intersect(LM22_genes, rownames(mat_specg)), ]
df_LM22_pansig <- dep_specific1 %>% 
  dplyr::rename(Gene = Genes) %>% 
  inner_join(df_LM22_pansig, by = 'Gene')
rownames(df_LM22_pansig) <- df_LM22_pansig$Gene
rio::export(df_LM22_pansig, 'Pancancer_signatures_LM22_overlap.xlsx')

pheatmap::pheatmap(
  df_LM22_pansig %>% select(`B cells naive`:Neutrophils) %>% mutate_all(as.numeric),
  annotation_row = df_LM22_pansig %>% select(Type, cluster),
  annotation_colors = list(cluster = ann_colors$cluster, Type = c(up = 'red3', down = 'blue3')),
  scale = "none",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T, angle_col = 45,
  color =colorRampPalette(c("white","black"))(100),
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  main = 'Pan-cancer signatures & LM22 overlap',
  filename = 'Pancancer_signatures_LM22_overlap.pdf', width = 8, height = 12
)

pheatmap::pheatmap(
  df_LM22_pansig %>% select(`B cells naive`:Neutrophils) %>% mutate_all(as.numeric),
  annotation_row = df_LM22_pansig %>% select(Type, cluster),
  annotation_colors = list(cluster = ann_colors$cluster, Type = c(up = 'red3', down = 'blue3')),
  scale = "none",
  cluster_rows = T,
  cluster_cols = T,
  show_rownames = T,
  show_colnames = T, angle_col = 45,
  color =colorRampPalette(c("white","black"))(100),
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  main = 'Pan-cancer signatures & LM22 overlap',
  filename = 'Pancancer_signatures_LM22_overlap_cluster.pdf', width = 8, height = 12
)



# ssGSEA
ssgsea <- GSVA::gsva(
  param = ssgseaParam(#exprData = mat_specg,
    exprData = scale(mat_specg),
    geneSets = gs_LM22),#, minSize = 15, maxSize = 500
  # expr = mat_specg,
  # gset.idx.list = h1gmt,# min.sz=15, max.sz=500,
  # method = 'ssgsea',
  # kcdf = 'Gaussian', abs.ranking = T
)

# ssgsea <- GSVA::gsva(
#   expr = scale(mat_specg),# a matrix of expression values where rows correspond to genes and columns correspond to samples.
#   gset.idx.list = gs_LM22,
#   method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = T,
# )
# save(ssgsea, file = 'ssgsea_pancancer_LM22.RData')
# save(ssgsea, file = 'ssgsea_pancancer_LM22_sampleScale.RData')
load('ssgsea_pancancer_LM22_sampleScale.RData')

# min-max normalization
ssgsea.1 <- ssgsea
for (i in colnames(ssgsea)) {
  #i <- colnames(ssgsea)[1]
  ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
  
}

pheatmap::pheatmap(
  ssgsea.1,
  annotation_col = df_sil1 %>% dplyr::select(cluster, cancer_abbr),
  annotation_colors = ann_colors,
  scale = "none",
  cluster_rows = T,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # cellwidth = 10, cellheight = 15,
  fontsize = 10
)

a <- pheatmap::pheatmap(
  ssgsea,
  annotation_col = df_sil1 %>% dplyr::select(cluster, cancer_abbr),
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = T,
  cluster_cols = T,
  show_rownames = T,
  show_colnames = F,
  cutree_cols = 4, clustering_method = 'ward.D2', clustering_distance_cols = 'euclidean',
  color =colorRampPalette(c("blue", "white","red"))(100),
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  filename = 'ssgsea_pancancer_LM22_heatmap.pdf', width = 11, height = 8
)


df_sil2 <- as.data.frame(cancer_color) %>% set_names('color') %>% rownames_to_column('cancer_abbr') %>% 
  left_join(df_sil1, .) %>% 
  column_to_rownames('FileName') %>% 
  dplyr::select(cluster, cancer_abbr, class_abbr, color)

# with pheatmap hclust results
fname_ordered <- a$tree_col$labels[a$tree_col$order]
gname_ordered <- a$tree_row$labels[a$tree_row$order]
groups <- cutree(a$tree_col, k = 4)

ann_col <- df_sil2 %>% dplyr::select(cluster, cancer_abbr)
ann_col <- groups %>% as.data.frame() %>% setNames('Immune_type') %>% merge(ann_col, by = 'row.names') %>% column_to_rownames('Row.names')
ann_col$Immune_type %<>% factor(levels = unique(groups[fname_ordered]))
ann_col <- ann_col[names(groups), ]
ann_col %<>% dplyr::select(cancer_abbr, cluster, Immune_type)

ann_colors$Immune_type <- c("red4", "blue4", "green4", "purple4") %>% setNames(unique(groups[fname_ordered])) %>% .[as.character(1:4)]
# ann_colors$cancer_abbr <- cancer_color

pheatmap::pheatmap(
  ssgsea[gname_ordered, fname_ordered],
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(ann_col$Immune_type))[-length(cumsum(table(ann_col$Immune_type)))],
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  main = 'ssGSEA Pan-cancer LM22 scores (Z-score)',
  filename = 'ssgsea_pancancer_LM22_heatmap_v2.pdf', width = 12, height = 9
)








# median of ssGSEA scores
mat_imm_median <- aggregate(t(ssgsea[, rownames(ann_col)]), by = list(ann_col$Immune_type), median) %>% column_to_rownames('Group.1') %>% t()
pheatmap::pheatmap(
  mat_imm_median,
  scale = "row",
  cluster_rows = T,
  cluster_cols = T,
  color =colorRampPalette(c("blue", "white","red"))(100),
  fontsize = 10
)
mat_imm_median %>% as.data.frame() %>% summarise_all(median)
# 4          3         1         2
# 1 0.04414977 0.02659249 0.2222103 0.1184114
# abundance order: 3>4>2>1
abundance_order <- mat_imm_median %>% as.data.frame() %>% summarise_all(median) %>% t() %>% .[, 1] %>% sort(decreasing = T) %>% names()

pheatmap::pheatmap(
  mat_imm_median[, abundance_order],
  scale = "row",
  cluster_rows = T,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F, angle_col = 0,
  annotation_col = data.frame(Immune_type = factor(1:4), row.names = 1:4),
  annotation_colors = list(Immune_type = ann_colors$Immune_type),
  color = colorRampPalette(c("blue", "white","red"))(100),
  # cellwidth = 10, cellheight = 15,
  main = 'ssGSEA Pan-cancer LM22 median scores (Z-score)',
  fontsize = 10,
  filename = 'ssgsea_pancancer_LM22_hclust_median_heatmap.pdf', width = 7.5, height = 6
)

##### hclust ----------
# perform hierarchical clustering using euclidean distance and Ward-linkage

##### hclust -- scale before dist()
# qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# set.seed(10)
# label_colors <- sample(col_vector, 4)
label_colors <- c("red4", "blue4", "green4", "purple4")

# Compute the distance matrix using Euclidean distance:
d <- dist(t(scale(mat_specg)), method = "euclidean") # with rows corresponding to observations and columns to variables
# as.matrix(d) %>% View()
clust <- hclust(d, method = "ward.D2")
groups <- cutree(clust, k = 3)
plot(clust)

dend <- as.dendrogram(clust) # create dendrogram object

pdf('ssgsea_pancancer_LM22_hclust.pdf', width = 3, height = 20)
dend %>%
  dendextend::set("labels_cex", 0.25) %>%
  dendextend::set("branches_k_color", value = label_colors, k = 3) %>%
  # set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3) %>%
  plot(horiz = T, axes = T, ylab = "Euclidean distance & Ward-linkage")
# abline(v = 350, lty = 2)
# colored_bars(cbind(k10[,ncol(k10):1], colorspace::rainbow_hcl(51, c = 70, l  = 50)), pm_heat_dend, rowLabels = c(paste0("k = ", rev(c(5, 10, 20))), "Sample type"), horiz = T)
graphics.off()

pdf('ssgsea_pancancer_LM22_hclust_noLeaf.pdf', width = 20, height = 3)
dend %>%
  dendextend::set("labels_cex", 0.25) %>%
  dendextend::set("branches_k_color", value = label_colors, k=3) %>%
  # set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3) %>%
  plot(horiz = F, axes = T, leaflab = "none", main = "Euclidean distance & Ward-linkage")
# abline(v = 350, lty = 2)
# colored_bars(cbind(k10[,ncol(k10):1], colorspace::rainbow_hcl(51, c = 70, l  = 50)), pm_heat_dend, rowLabels = c(paste0("k = ", rev(c(5, 10, 20))), "Sample type"), horiz = T)
graphics.off()


fname_ordered <- clust$labels[clust$order]
ann_col <- df_sil1 %>% #column_to_rownames('FileName') %>%
  dplyr::select(cluster, cancer_abbr)
ann_col <- groups %>% as.data.frame() %>% setNames('Immune_type') %>% merge(ann_col, by = 'row.names') %>% column_to_rownames('Row.names')
ann_col$Immune_type %<>% factor(levels = unique(groups[clust$labels[clust$order]]))
ann_col <- ann_col[names(groups), ]
ann_col %<>% dplyr::select(cancer_abbr, cluster, Immune_type)

ann_colors$Immune_type <- label_colors %>% setNames(unique(groups[fname_ordered])) %>% .[as.character(1:3)]

# #breaks
# quantile(scale(t(ssgsea)))
# # bk1 <- quantile(scale(t(ssgsea)), 0.1)
# # bk2 <- quantile(scale(t(ssgsea)), 0.9)
# 
# my_breaks <- seq(-2, 2, by = 0.01)
# 
# #colors
# # my_colors <- rev(c(colorRampPalette(colors = c('red4', 'white'))(length(my_breaks)/2),
# #                    colorRampPalette(colors = c('white', 'blue4'))(length(my_breaks)/2)))
# my_colors <- rev(colorRampPalette(colors = c('red3', 'gray', 'blue3'))(length(my_breaks)))


pheatmap::pheatmap(
  ssgsea[, fname_ordered],
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = T,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(ann_col$Immune_type))[-length(cumsum(table(ann_col$Immune_type)))],
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  filename = 'ssgsea_pancancer_LM22_hclust_heatmap.pdf', width = 12, height = 9
)





# median of ssGSEA scores
mat_imm_median <- aggregate(t(ssgsea[, fname_ordered]), by = list(ann_col$Immune_type), median) %>% column_to_rownames('Group.1') %>% t()
pheatmap::pheatmap(
  mat_imm_median,
  scale = "row",
  cluster_rows = T,
  cluster_cols = T,
  show_rownames = T,
  show_colnames = T,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # cellwidth = 10, cellheight = 15,
  fontsize = 10
)
mat_imm_median %>% as.data.frame() %>% summarise_all(median)
# 2          1         3
# 1 0.08664359 0.07792871 0.1040334
abundance_order <- mat_imm_median %>% as.data.frame() %>% summarise_all(median) %>% t() %>% .[, 1] %>% sort(decreasing = T) %>% names()
# abundance order: 1 > 2 > 3



pheatmap::pheatmap(
  mat_imm_median[, abundance_order],
  scale = "none",
  cluster_rows = T,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T, angle_col = 0,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # cellwidth = 10, cellheight = 15,
  main = 'ssGSEA Pan-cancer LM22 median scores',
  fontsize = 10,
  filename = 'ssgsea_pancancer_LM22_hclust_median_heatmap.pdf', width = 7, height = 6
)

pheatmap::pheatmap(
  mat_imm_median[, abundance_order],
  scale = "row",
  cluster_rows = T,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T, angle_col = 0,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # cellwidth = 10, cellheight = 15,
  main = 'ssGSEA Pan-cancer LM22 median scores (Z-score)',
  fontsize = 10,
  filename = 'ssgsea_pancancer_LM22_hclust_median_rowScale_heatmap.pdf', width = 7, height = 6
)


#### sample-protein heatmap
pheatmap::pheatmap(
  mat_specg[df_LM22_pansig$Gene, fname_ordered],
  annotation_row = df_LM22_pansig %>% dplyr::select(Type, cluster),
  annotation_col = ann_col,
  annotation_colors = list(cluster = ann_colors$cluster, Type = c(up = 'red3', down = 'blue3'),
                           cancer_abbr = ann_colors$cancer_abbr, Immune_type = ann_colors$Immune_type),
  scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(ann_col$Immune_type))[-length(cumsum(table(ann_col$Immune_type)))],
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  main = 'Pan-cancer signatures & LM22 overlap expression (Z-score)',
  filename = 'Pancancer_signatures_LM22_overlap_expression.pdf', width = 12, height = 12
)

pheatmap::pheatmap(
  mat_specg[df_LM22_pansig$Gene, fname_ordered],
  annotation_row = df_LM22_pansig %>% dplyr::select(Type, cluster),
  annotation_col = ann_col,
  annotation_colors = list(cluster = ann_colors$cluster, Type = c(up = 'red3', down = 'blue3'),
                           cancer_abbr = ann_colors$cancer_abbr, Immune_type = ann_colors$Immune_type),
  scale = "row",
  cluster_rows = T,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(ann_col$Immune_type))[-length(cumsum(table(ann_col$Immune_type)))],
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  main = 'Pan-cancer signatures & LM22 overlap expression (Z-score)',
  filename = 'Pancancer_signatures_LM22_overlap_expression_cluster.pdf', width = 12, height = 12
)

## save image before change R env
# save.image(file = 'Pancancer_signatures_LM22_generated.RData')

#### b.ImmuneSubtypeClassifier ------------
##### config. 换一个临时的 libpath ------
# # **Right Now, the newest version of xgboost is incompatible. **
# # **Please use a prior version, 1.0.0.1 or 1.0.0.2 should work.**
# 
# .libPaths()
# # [1] "C:/Users/A/AppData/Local/R/win-library/4.2"
# # [2] "C:/Program Files/R/R-4.2.0/library"  
# 
# packageVersion('xgboost') # 1.6.0.1
# installed.packages() %>% as.data.frame() %>% filter(Package == 'xgboost') # installed at [2]
# 
# # withr::with_libpaths(
# #   new = "C:/Users/A/AppData/Local/R/win-library/4.2/custome",
# #   devtools::install_version("xgboost", version = "1.0.0.1", )
# # )
# 
# # devtools::install_version('xgboost', version = '1.0.0.1', lib = 'C:/Users/dell/AppData/Local/R_custome')
# 
# .libPaths(c('C:/Users/dell/AppData/Local/R_custome', .libPaths()))
# 
# install.packages('https://cran.rstudio.com//src/contrib/Archive/xgboost/xgboost_1.0.0.1.tar.gz', version = '1.0.0.1', lib = 'C:/Users/dell/AppData/Local/R_custome')
# packageVersion('xgboost') # 1.0.0.1
# installed.packages() %>% as.data.frame() %>% filter(Package == 'xgboost') # right
# 
# devtools::install_github('CRI-iAtlas/ImmuneSubtypeClassifier', lib = 'C:/Users/dell/AppData/Local/R_custome')
# installed.packages() %>% as.data.frame() %>% filter(Package == 'ImmuneSubtypeClassifier') # right


##### classifier --------
# load('Pancancer_signatures_LM22_generated.RData')

.libPaths(c('C:/Users/dell/AppData/Local/R_custome', .libPaths()))
library(ImmuneSubtypeClassifier)
imm_subtypes <- c('Wound Healing', 'IFN-gamma', 'Inflammatory', 'Lymphocyte Depleted', 'Immunologically Quiet', 'TGF-Beta') %>% set_names(1:6)

# To get a list of the genes needed:
data(ebpp_gene)
dim(ebpp_genes_sig)  ### 485 genes are needed

# To make calls on new data,
Xtest <- as.matrix(mat_specg) # has gene IDs in rownames and sample IDs in column names



# use only specific proteins
res0 <- callEnsemble(X=Xtest, geneids='symbol')

res0 %>% count(BestCall)
#   BestCall   n
# 1        2   2
# 2        4 964

# label signature_name
# s1	  LIexpression_score	
# s2	  CSF1_response	
# s3	  Module3_IFN_score	
# s4	  TGFB_score_21050467	
# s5	  CHANG_CORE_SERUM_RESPONSE_UP	


# # to see where and if gene ID matches have failed:
# match_error <- geneMatchErrorReport(X=Xtest, geneid='symbol')
# match_error$matchError # 0.8969072
# match_error$missingGenes # dim 435 3



# use full matrix
X <- mat.rm.pool[, fname_ordered]
# dfprot %>% filter(Protein %in% rownames(mat.rm.pool)) %>% count(Gene) %>% filter(n > 1) %>% semi_join(dfprot, .) %>% arrange(Gene, Protein) # involved 9 genes
# 
# dfprot1 <- dfprot %>% filter(Protein %in% rownames(mat.rm.pool)) %>% count(Gene) %>% filter(n > 1) %>% anti_join(dfprot, .) %>% arrange(Gene, Protein)
# dfprot2 <- dfprot %>% filter(Protein %in% rownames(mat.rm.pool)) %>% count(Gene) %>% filter(n > 1) %>% semi_join(dfprot, .) %>% arrange(Gene, Protein)
# dfprot2 %<>% filter(Gene != '') %>% filter(str_detect(`Protein Description`, 'isoform', negate = T)) # remove null string and isoforms
# dfprot2 %<>% count(Gene) %>% filter(n == 1) %>% semi_join(dfprot2, .) # remove other genes
# 
# dfprot_icls <- rbind(dfprot1, dfprot2)
# dfprot_icls1 <- dfprot_icls %>% count(Gene) %>% filter(n == 1) %>% semi_join(dfprot_icls, .)
# dfprot_icls2 <- dfprot_icls %>% count(Gene) %>% filter(n > 1) %>% semi_join(dfprot_icls, .) # 不好筛选，干脆都不考虑了，涉及 9 genes
# 
# dfprot_new <- dfprot_icls1 %>% filter(Protein %in% rownames(X))
# length(dfprot_new$Protein) # 11357
# length(dfprot_new$Gene) # 11357
# dim(X) # 11386   436
# 
dfprot_new <- dfprot %>%
  dplyr::rename(Protein = Protein.Group, Gene = Genes) %>% 
  filter(Protein %in% rownames(X))
dfprot_new1 <- dfprot_new %>%
  count(Gene) %>%
  filter(n == 1) %>%
  semi_join(dfprot_new, .)

Xp <- Xg <- X[dfprot_new1$Protein, ]
dim(Xg) # 12721  1014

prot2gene <- dfprot_new$Gene
names(prot2gene) <- dfprot_new$Protein
rownames(Xg) <- prot2gene[rownames(Xg)]
table(table(rownames(Xg)))
# 1 
# 12721

res <- callEnsemble(X = Xg, geneids='symbol')
res$BestType <- imm_subtypes[res$BestCall] %>% factor(levels = imm_subtypes)
res %>% count(BestType)
# BestType   n
# 1         Wound Healing   2
# 2             IFN-gamma 138
# 3          Inflammatory  85
# 4   Lymphocyte Depleted 786
# 5 Immunologically Quiet   2
# 6              TGF-Beta   1
colnames(res)[3:8] <- imm_subtypes



##### heatmap --------
res0$BestType <- imm_subtypes[res0$BestCall] %>% factor(levels = imm_subtypes)
colnames(res0)[3:8] <- imm_subtypes
rownames(res0) <- res0$SampleIDs
res0 <- res0[rownames(ann_col), ]

pheatmap::pheatmap(
  ssgsea[gname_ordered, fname_ordered],
  annotation_col = cbind(ann_col, res0[, 8:3]),
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(ann_col$Immune_type))[-length(cumsum(table(ann_col$Immune_type)))],
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  main = 'ssGSEA Pan-cancer LM22 scores (Z-score)',
  filename = 'ssgsea_pancancer_LM22_heatmap_v2_addImmTypeCls_onlySigP.pdf', width = 12, height = 9
)


rownames(res) <- res$SampleIDs
res <- res[rownames(ann_col), ]

pheatmap::pheatmap(
  ssgsea[gname_ordered, fname_ordered],
  annotation_col = cbind(ann_col, res[, 8:3]),
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(ann_col$Immune_type))[-length(cumsum(table(ann_col$Immune_type)))],
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  main = 'ssGSEA Pan-cancer LM22 scores (Z-score)',
  filename = 'ssgsea_pancancer_LM22_heatmap_v2_addImmTypeCls.pdf', width = 12, height = 9
)

#### c.ESTIMATE -----------
##### ESTIMATE ------
# install.packages('C:/Users/dell/AppData/Local/R_custome/estimate_1.0.11.tar.gz', lib = 'C:/Users/dell/AppData/Local/R_custome')
.libPaths(c('C:/Users/dell/AppData/Local/R_custome', .libPaths()))
library(estimate)
help(package = 'estimate')
# 
# OvarianCancerExpr <- system.file("extdata", "sample_input.txt", package="estimate")
# filterCommonGenes(input.f=OvarianCancerExpr, output.f="OV_10412genes.gct", id="GeneSymbol")
# estimateScore("OV_10412genes.gct", "OV_estimate_score.gct", platform="affymetrix")
# plotPurity(scores="OV_estimate_score.gct", samples="s516", platform="affymetrix")
# 

write.table(Xg, 'pancancer_full_matrix_for_ESTIMATE.txt', row.names = T, quote = F, sep = '\t')
filterCommonGenes(input.f = 'pancancer_full_matrix_for_ESTIMATE.txt', output.f = 'pancancer_full_matrix_for_ESTIMATE.gct', id = 'GeneSymbol')
# [1] "This dataset includes 7836genes."
# [1] "2576genes were mismatched."
estimateScore("pancancer_full_matrix_for_ESTIMATE.gct", "pancancer_full_matrix_for_ESTIMATE_estimate_score.gct", platform="affymetrix")
# [1] "1 gene set: StromalSignature  overlap= 111"
# [1] "2 gene set: ImmuneSignature  overlap= 116"
plotPurity(scores = "pancancer_full_matrix_for_ESTIMATE_estimate_score.gct", samples = 'pancancer_full_matrix_for_ESTIMATE', platform = "affymetrix")

df_estimate <- rio::import('pancancer_full_matrix_for_ESTIMATE_estimate_score.gct', format = '\t')
colnames(df_estimate) <- df_estimate[3, ]
df_estimate <- df_estimate[-(1:3), ]
df_estimate[, -(1:2)] <- df_estimate[, -(1:2)] %>% mutate_all(as.numeric)
df_estimate %<>% dplyr::select(-Description)
rownames(df_estimate) <- NULL
df_estimate %<>% column_to_rownames('NAME') %>%
  t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  mutate(FileName = str_replace_all(FileName, '\\.', '-'))

rio::export(df_estimate, 'pancancer_full_matrix_for_ESTIMATE_estimate_score.xlsx')



##### heatmap update ---------
rownames(res) <- res$SampleIDs
res <- res[rownames(ann_col), ]
rownames(df_estimate) <- df_estimate$FileName
df_estimate <- df_estimate[rownames(ann_col), ]

pheatmap::pheatmap(
  ssgsea[gname_ordered, fname_ordered],
  annotation_col = cbind(ann_col, res[, 8:3], df_estimate[, -1]),
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(ann_col$Immune_type))[-length(cumsum(table(ann_col$Immune_type)))],
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  main = 'ssGSEA Pan-cancer LM22 scores (Z-score)',
  filename = 'ssgsea_pancancer_LM22_heatmap_v2_addImmTypeCls_addESTIMATE.pdf', width = 12, height = 8
)

list(LM22 = ssgsea[gname_ordered, fname_ordered] %>% as.data.frame() %>% rownames_to_column('LM22'),
     ann_col = cbind(ann_col, res[, 8:3], df_estimate[, -1]) %>% rownames_to_column('FileName')) %>% 
  rio::export('ssgsea_pancancer_LM22_heatmap_v2_addImmTypeCls_addESTIMATE.xlsx')

pheatmap::pheatmap(
  ssgsea[gname_ordered, fname_ordered],
  annotation_col = cbind(ann_col, res[, 8:3], df_estimate[, -1]),
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = F,
  color =colorRampPalette(c("blue", "white","red"))(100),
  # color = my_colors, breaks = my_breaks,
  gaps_col = cumsum(table(ann_col$Immune_type))[-length(cumsum(table(ann_col$Immune_type)))],
  cellwidth = 0.1, cellheight = 7,
  fontsize = 7,
  main = 'ssGSEA Pan-cancer LM22 scores (Z-score)',
  filename = 'ssgsea_pancancer_LM22_heatmap_v2_addImmTypeCls_addESTIMATE_v2.pdf', width = 12, height = 12
)

# save.image(file = 'Pancancer_signatures_LM22_imm_est_generated.RData')


#### ComplexHeatmap ------
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
})

mat.pancancer <- rbind(ssgsea[gname_ordered, fname_ordered],
                       mat_heat[nrow(mat_heat):1, ],
                       t(res[, 8:3]),
                       t(df_estimate[, -1])) %>% .[, colnames(mat_heat)]
mat.pancancer.scale <- t(scale(t(mat.pancancer)))

my_breaks <- seq(-5, 5, by = 0.01)
my_colors <- rev(colorRampPalette(colors = c('red', 'white', 'blue'))(length(my_breaks)))
# my_colors <- viridisLite::viridis(length(my_breaks), option = 'G')

# pheatmap::pheatmap(
#   mat.pancancer.scale,
#   # annotation_col = cbind(ann_col, res[, 8:3], df_estimate[, -1]),
#   annotation_col = ann_col %>% select(cancer_abbr, Immune_type, cluster),
#   annotation_colors = ann_colors,
#   scale = "none",
#   cluster_rows = T,
#   cluster_cols = F,
#   clustering_method = 'ward.D2',
#   show_rownames = T,
#   show_colnames = F,
#   # color =colorRampPalette(c("blue", "white","red"))(100),
#   color = my_colors, breaks = my_breaks,
#   # gaps_col = cumsum(table(ann_col$Immune_type))[-length(cumsum(table(ann_col$Immune_type)))],
#   gaps_col = cumsum(table(df_sil1$cluster)[levels(df_sil1$cluster)])[-length(cumsum(table(df_sil1$cluster)))],
#   # cellwidth = 10, cellheight = 15,
#   fontsize = 7,
#   main = 'Pan-cancer ssGSEA hallmark - LM22 scores (Z-score)',
#   # filename = 'pancancer_ssGSEA_scores_all.pdf', width = 12, height = 7
# )

# ComplexHeatmap
ann_use <- ann_col[colnames(mat.pancancer.scale), c("cluster","Immune_type","cancer_abbr")] %>%
  mutate(cluster = factor(cluster, levels = levels(ann_col$cluster)))

col_fun <- circlize::colorRamp2(c(-5, 0, 5), c("red", "white", "blue"))

ha_top <- HeatmapAnnotation(
  df = ann_use,
  col = ann_colors,
  show_annotation_name = TRUE,
  annotation_name_side = "left"
)

ht <- Heatmap(
  mat.pancancer.scale,
  name = "Z-score",
  col = col_fun,
  cluster_rows = TRUE,
  clustering_method_rows = "ward.D2",
  column_split = ann_use$cluster,
  cluster_columns = TRUE,
  clustering_method_columns = "ward.D2",
  cluster_column_slices = FALSE,
  column_gap = unit(1, "mm"),
  show_row_names = TRUE,
  show_column_names = FALSE,
  top_annotation = ha_top,
  heatmap_legend_param = list(title = "Z-score")
)

pdf('pancancer_ssGSEA_scores_all.pdf', width = 20, height = 14)
draw(ht, merge_legend = TRUE)
graphics.off()


#### d.EMT score ------
# load(file = 'Pancancer_signatures_LM22_imm_est_generated.RData')



### 5.1.5 Survival --------
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
})

#### helper functions -----
plot_km <- function(fit_km, data, endpoint_label, time_lab = "Days") {
  ggsurvplot(
    fit_km, data = data,
    pval = TRUE,                 # log-rank p-value on the plot
    risk.table = TRUE,
    conf.int = FALSE,
    surv.median.line = "hv",
    legend.title = "Cluster",
    legend.labs  = paste0("C", levels(data$cluster)),
    xlab = time_lab,
    ylab = paste(endpoint_label, "survival probability"),
    title = paste0(endpoint_label, ": Kaplan–Meier curves"),
    risk.table.title = "Number at risk"
  )
}

save_km <- function(plt, tag) {
  ggsave(paste0("surv_KM_", tag, ".pdf"), plt$plot, width = 7, height = 6)
  ggsave(paste0("surv_RiskTable_", tag, ".pdf"), plt$table, width = 7, height = 3)
}

save_cox_forest <- function(fit_cox, data, tag, main_title) {
  # Forest plot of hazard ratios for model terms
  fp <- ggforest(
    fit_cox, data = data,
    main = main_title,
    fontsize = 1.0
  )
  ggsave(paste0("793patients_5clusters_forest_Cox_", tag, ".pdf"), fp, width = 7, height = 6)
}

save_schoenfeld <- function(cz_obj, tag, main_title) {
  pdf(paste0("793patients_5clusters_schoenfeld_", tag, ".pdf"), width = 7, height = 6)
  plot(cz_obj, main = main_title)  # One panel per coefficient
  dev.off()
}



info_prog0 <- rio::import('input/Prognosis_info_20250922.xlsx')
info_prog_edit <- rio::import('input/修改后 重复病例信息.xlsx')

info_prog0a <- info_prog0 %>% filter(!患者ID %in% info_prog_edit$患者ID)
info_prog0b <- info_prog_edit %>% select(-保留)
info_prog <- rbind(info_prog0a, info_prog0b)

status.converter <- list(
  OS  = c("存活"=0L, "失访"=0L, "死亡"=1L, "不详"=NA_integer_),
  PFS = c(
    "无进展"=0L, "失访"=0L, "死亡"=1L, "疾病进展"=1L,
    "复发后直接手术"=1L, "发现转移后直接化疗"=1L,
    "复发后进行一周期临床试验"=1L, "不详"=NA_integer_
  )
)

info_prog1 <- info_prog %>% mutate(
  `Admission ID` = 患者ID,
  OS.status = OS状态,
  OS.Days = `OS时间(天)`,
  PFS.status = PFS状态,
  PFS.Days = `PFS时间(天)`
) %>%
  select(`Admission ID`:PFS.Days) %>%
  mutate(OS.status = status.converter$OS[OS.status],
         PFS.status = status.converter$PFS[PFS.status]) %>% 
  mutate_at(vars(OS.status, PFS.status), factor) %>% 
  mutate_at(vars(OS.Days, PFS.Days), as.integer)

# df_prog <- info_prog1 %>% count(`Admission ID`) %>% filter(n == 1) %>% 
#   semi_join(info_prog1, .) %>% 
#   inner_join(info_pan %>% distinct(`Admission ID`, patient_ID)) %>% 
#   inner_join(df_sil1)

df_prog <- df_sil1 %>%
  inner_join(info_pan %>% distinct(`Admission ID`, patient_ID)) %>%
  inner_join(info_prog1) %>%
  mutate(
    cluster = str_c('P', str_extract(as.character(cluster), '\\d+$')),
    cluster = factor(cluster, levels = str_c('P', sort(unique(as.numeric(str_remove(cluster, '^P')))))),
    OS_days   = suppressWarnings(as.numeric(OS.Days)),
    OS_event  = suppressWarnings(as.integer(OS.status)),
    PFS_days  = suppressWarnings(as.numeric(PFS.Days)),
    PFS_event = suppressWarnings(as.integer(PFS.status))
  ) %>%
  select(cluster, OS_days, OS_event, PFS_days, PFS_event)

rio::export(df_prog, '793patients_5clusters_prognosis.xlsx')

#### a. OS ------------------
##### KM ----
x_os <- df_prog %>% filter(!is.na(OS_days), OS_days > 0, !is.na(OS_event))
stopifnot(nrow(x_os) > 0L, nlevels(x_os$cluster) >= 2L)

fit_km_os <- survfit(Surv(OS_days, OS_event) ~ cluster, data = x_os)
fit_lr_os <- survdiff(Surv(OS_days, OS_event) ~ cluster, data = x_os)
p_lr_os   <- pchisq(fit_lr_os$chisq, df = length(fit_lr_os$n) - 1L, lower.tail = FALSE)

print(fit_lr_os); cat("OS log-rank p =", signif(p_lr_os, 4), "
")

plt_os <- ggsurvplot(
  fit_km_os, data = x_os, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
  legend.title = "cluster", legend.labs = paste0("Cluster-", levels(x_os$cluster)),
  xlab = "Days", ylab = "OS survival probability"
)
pdf('793patients_5clusters_surv_KM_OS_V1009.pdf', width = 6, height = 6)
print(plt_os)
graphics.off()

fit_cox_os <- coxph(Surv(OS_days, OS_event) ~ cluster, data = x_os)
sm_os <- summary(fit_cox_os)
p_cox_os <- unname(sm_os$logtest["pvalue"])
print(sm_os); cat("OS Cox overall LR p =", signif(p_cox_os, 4), "
")

cz_os <- cox.zph(fit_cox_os)
print(cz_os)

# Forest plot (HRs vs reference level of 'cluster')
save_cox_forest(fit_cox_os, x_os, "OS", "OS: Cox model hazard ratios (HR, 95% CI)")

# PH assumption diagnostic (Schoenfeld)
save_schoenfeld(cz_os, "OS", "OS: Schoenfeld residuals (tests of PH)")
graphics.off()


#### b. PFS --------------
x_pfs <- df_prog %>% filter(!is.na(PFS_days), PFS_days > 0, !is.na(PFS_event))
stopifnot(nrow(x_pfs) > 0L, nlevels(x_pfs$cluster) >= 2L)

fit_km_pfs <- survfit(Surv(PFS_days, PFS_event) ~ cluster, data = x_pfs)
fit_lr_pfs <- survdiff(Surv(PFS_days, PFS_event) ~ cluster, data = x_pfs)
p_lr_pfs   <- pchisq(fit_lr_pfs$chisq, df = length(fit_lr_pfs$n) - 1L, lower.tail = FALSE)

print(fit_lr_pfs); cat("PFS log-rank p =", signif(p_lr_pfs, 4), "
")

plt_pfs <- ggsurvplot(
  fit_km_pfs, data = x_pfs, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
  legend.title = "cluster", legend.labs = paste0("Cluster-", levels(x_pfs$cluster)),
  xlab = "Days", ylab = "PFS survival probability"
)
pdf('793patients_5clusters_surv_KM_PFS_V1009.pdf', width = 6, height = 6)
print(plt_pfs)
graphics.off()

fit_cox_pfs <- coxph(Surv(PFS_days, PFS_event) ~ cluster, data = x_pfs)
sm_pfs <- summary(fit_cox_pfs)
p_cox_pfs <- unname(sm_pfs$logtest["pvalue"])
print(sm_pfs); cat("PFS Cox overall LR p =", signif(p_cox_pfs, 4), "
")

cz_pfs <- cox.zph(fit_cox_pfs)
print(cz_pfs)

# Forest plot (HRs vs reference level of 'cluster')
save_cox_forest(fit_cox_pfs, x_pfs, "PFS", "PFS: Cox model hazard ratipfs (HR, 95% CI)")

# PH assumption diagnpfstic (Schoenfeld)
save_schoenfeld(cz_pfs, "PFS", "PFS: Schoenfeld residuals (tests of PH)")
graphics.off()


#### c. Concise summary table -----
summary_df <- tibble(
  Endpoint   = c("OS", "PFS"),
  N          = c(nrow(x_os), nrow(x_pfs)),
  Events     = c(sum(x_os$OS_event == 1L, na.rm = TRUE), sum(x_pfs$PFS_event == 1L, na.rm = TRUE)),
  Logrank_p  = c(p_lr_os, p_lr_pfs),
  Cox_p      = c(p_cox_os, p_cox_pfs)
)
print(summary_df)
write.csv(summary_df, "793patients_5clusters_survival_summary_V1009.csv", row.names = FALSE)

cat("Done. Saved PDFs and survival_summary.csv in:", normalizePath(getwd(), winslash = "/", mustWork = FALSE), "
")

###################
#add analysis
save_cox_forest <- function(fit_cox, data, tag, main_title) {
  # Forest plot of hazard ratios for model terms
  fp <- ggforest(
    fit_cox, data = data,
    main = main_title,
    fontsize = 1.0
  )
  ggsave(paste0(unique(PTDB$cancer_abbr)[i],"_cluster_",paste(indx$cluster,collapse = "_"),"_",nrow(ind_prog),"patients_forest_Cox_", tag, ".pdf"), fp, width = 7, height = 6)
}

save_schoenfeld <- function(cz_obj, tag, main_title) {
  pdf(paste0(unique(PTDB$cancer_abbr)[i],"_cluster_",paste(indx$cluster,collapse = "_"),"_",nrow(ind_prog),"patients_schoenfeld_", tag, ".pdf"), width = 7, height = 6)
  plot(cz_obj, main = main_title)  # One panel per coefficient
  dev.off()
}


PTDB <- rio::import('./input/PTDB_pancancer_6399_proteins_969_coreSamples_5clusters_edited.xlsx')
PTDB<-PTDB[!is.na(PTDB$compare_cancer_typer_prognosis),]

unique(df_sil1$cancer_abbr)
unique(df_sil1$cluster)



#五组分开做survival

for (i in 1:length(unique(PTDB$cancer_abbr))) {
  indx<-PTDB[PTDB$cancer_abbr==unique(PTDB$cancer_abbr)[i],]
  indx$cluster<-paste0("P",indx$cluster)
  ind_sil1<-df_sil1[df_sil1$cluster%in%indx$cluster&df_sil1$cancer_abbr==unique(PTDB$cancer_abbr)[i],]
  ind_prog <- ind_sil1 %>%
    inner_join(info_pan %>% distinct(`Admission ID`, patient_ID)) %>%
    inner_join(info_prog1) %>%
    mutate(
      cluster = str_c('P', str_extract(as.character(cluster), '\\d+$')),
      cluster = factor(cluster, levels = str_c('P', sort(unique(as.numeric(str_remove(cluster, '^P')))))),
      OS_days   = suppressWarnings(as.numeric(OS.Days)),
      OS_event  = suppressWarnings(as.integer(OS.status)),
      PFS_days  = suppressWarnings(as.numeric(PFS.Days)),
      PFS_event = suppressWarnings(as.integer(PFS.status))
    ) %>%
    select(cluster, OS_days, OS_event, PFS_days, PFS_event)
  
  ###OS
  if(sum(!is.na(ind_prog$OS_days))>0&sum(!is.na(ind_prog$OS_event))>0){
    
  ind_os <- ind_prog %>% filter(!is.na(OS_days), OS_days > 0, !is.na(OS_event))
  stopifnot(nrow(ind_os) > 0L, nlevels(ind_os$cluster) >= 2L)
  
  fit_km_os <- survfit(Surv(OS_days, OS_event) ~ cluster, data = ind_os)
  fit_lr_os <- survdiff(Surv(OS_days, OS_event) ~ cluster, data = ind_os)
  p_lr_os   <- pchisq(fit_lr_os$chisq, df = length(fit_lr_os$n) - 1L, lower.tail = FALSE)
  
  print(fit_lr_os); cat("OS log-rank p =", signif(p_lr_os, 4), "
")
  
  plt_os <- ggsurvplot(
    fit_km_os, data = ind_os, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
    legend.title = "cluster", legend.labs = paste0("Cluster-", levels(ind_os$cluster)),
    xlab = "Days", ylab = "OS survival probability"
  )
  pdf(paste0(unique(PTDB$cancer_abbr)[i],"_cluster_",paste(indx$cluster,collapse = "_"),"_",nrow(ind_prog),'patients_surv_KM_OS_20251208.pdf'), width = 6, height = 6)
  print(plt_os)
  graphics.off()
  
  fit_cox_os <- coxph(Surv(OS_days, OS_event) ~ cluster, data = ind_os)
  sm_os <- summary(fit_cox_os)
  p_cox_os <- unname(sm_os$logtest["pvalue"])
  print(sm_os); cat("OS Cox overall LR p =", signif(p_cox_os, 4), "
")
  
  cz_os <- cox.zph(fit_cox_os)
  print(cz_os)
  
  # Forest plot (HRs vs reference level of 'cluster')
  save_cox_forest(fit_cox_os, ind_os, "OS", "OS: Cox model hazard ratios (HR, 95% CI)")
  
  # PH assumption diagnostic (Schoenfeld)
  save_schoenfeld(cz_os, "OS", "OS: Schoenfeld residuals (tests of PH)")
  graphics.off()
  
  }
  
  #PFS
  if(sum(!is.na(ind_prog$PFS_days))>0&sum(!is.na(ind_prog$PFS_event))>0){
  ind_pfs <- ind_prog %>% filter(!is.na(PFS_days), PFS_days > 0, !is.na(PFS_event))
  stopifnot(nrow(ind_pfs) > 0L, nlevels(ind_pfs$cluster) >= 2L)
  
  fit_km_pfs <- survfit(Surv(PFS_days, PFS_event) ~ cluster, data = ind_pfs)
  fit_lr_pfs <- survdiff(Surv(PFS_days, PFS_event) ~ cluster, data = ind_pfs)
  p_lr_pfs   <- pchisq(fit_lr_pfs$chisq, df = length(fit_lr_pfs$n) - 1L, lower.tail = FALSE)
  
  print(fit_lr_pfs); cat("PFS log-rank p =", signif(p_lr_pfs, 4), "
")
  
  plt_pfs <- ggsurvplot(
    fit_km_pfs, data = ind_pfs, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
    legend.title = "cluster", legend.labs = paste0("Cluster-", levels(ind_pfs$cluster)),
    xlab = "Days", ylab = "PFS survival probability"
  )
  
  pdf(paste0(unique(PTDB$cancer_abbr)[i],"_cluster_",paste(indx$cluster,collapse = "_"),"_",nrow(ind_prog),'patients_surv_KM_PFS_20251208.pdf'), width = 6, height = 6)
  print(plt_pfs)
  graphics.off()
  
  fit_cox_pfs <- coxph(Surv(PFS_days, PFS_event) ~ cluster, data = ind_pfs)
  sm_pfs <- summary(fit_cox_pfs)
  p_cox_pfs <- unname(sm_pfs$logtest["pvalue"])
  print(sm_pfs); cat("PFS Cox overall LR p =", signif(p_cox_pfs, 4), "
")
  
  cz_pfs <- cox.zph(fit_cox_pfs)
  print(cz_pfs)
  
  # Forest plot (HRs vs reference level of 'cluster')
  save_cox_forest(fit_cox_pfs, ind_pfs, "PFS", "PFS: Cox model hazard ratipfs (HR, 95% CI)")
  
  # PH assumption diagnpfstic (Schoenfeld)
  save_schoenfeld(cz_pfs, "PFS", "PFS: Schoenfeld residuals (tests of PH)")
  graphics.off()
  
  
}}

#cox log_rank 

# 1. Config ------
# User paths (ASCII-only)
info_pan <- rio::import('../6_pancancer_harmonizR_QC/output/PUH_pancancer_sample_information_2datesets_1146files.xlsx')
info.all <- rio::import('input/20251009_PUH_sample_information_2997files.xlsx')
# datall <- rio::import('input/imputed_raw_12754proteins_2997samples.csv') %>% column_to_rownames('V1')
datall <- rio::import('../2_batch_effect_evaluation_correction/output/batch_corrected_data_combat3.csv') %>% column_to_rownames('V1')


meta <- info.all %>% filter(sample_type == 'T')
expr <- datall[, meta$FileName]


dim(meta) # 1105   18
dim(expr) # 12754  1105


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
dim(mat_consensus) # 6377 1105
quantile(mat_consensus)
# 0%       25%       50%       75%      100% 
# 3.368983  9.455099 12.306143 14.488626 25.897229 



mat_df<-data.frame(mat_consensus,check.names = FALSE)

mat_df<-data.frame(t(mat_df))

df_sil1 <- readRDS('df_sil1.rds')
rownames(df_sil1) <- df_sil1$FileName
####根据前面的code ,差异分析等用的都是df_sil1中的1014个样本，因此这个cox也用这个数据

mat_df<-mat_df[row.names(mat_df)%in%row.names(df_sil1),]

info_prog0 <- rio::import('input/Prognosis_info_20250922.xlsx')
info_prog_edit <- rio::import('input/修改后 重复病例信息.xlsx')

info_prog0a <- info_prog0 %>% filter(!患者ID %in% info_prog_edit$患者ID)
info_prog0b <- info_prog_edit %>% select(-保留)
info_prog <- rbind(info_prog0a, info_prog0b)

status.converter <- list(
  OS  = c("存活"=0L, "失访"=0L, "死亡"=1L, "不详"=NA_integer_),
  PFS = c(
    "无进展"=0L, "失访"=0L, "死亡"=1L, "疾病进展"=1L,
    "复发后直接手术"=1L, "发现转移后直接化疗"=1L,
    "复发后进行一周期临床试验"=1L, "不详"=NA_integer_
  )
)

info_prog1 <- info_prog %>% mutate(
  `Admission ID` = 患者ID,
  OS.status = OS状态,
  OS.Days = `OS时间(天)`,
  PFS.status = PFS状态,
  PFS.Days = `PFS时间(天)`
) %>%
  select(`Admission ID`:PFS.Days) %>%
  mutate(OS.status = status.converter$OS[OS.status],
         PFS.status = status.converter$PFS[PFS.status]) %>% 
  mutate_at(vars(OS.status, PFS.status), factor) %>% 
  mutate_at(vars(OS.Days, PFS.Days), as.integer)


mat_df$`Admission ID`<-info_pan$`Admission ID`[match(row.names(mat_df),info_pan$FileName)]
#info.all中有第三批数据，当时没有patientID，直接用Admission ID当成 patientID 了，这里把这部分样本拎出来
mat_df_use<-mat_df[!is.na(mat_df$`Admission ID`),]
mat_df_NA<-mat_df[is.na(mat_df$`Admission ID`),]
mat_df_NA$`Admission ID`<-meta$patient_ID[match(row.names(mat_df_NA),meta$FileName)]
mat_df<-rbind(mat_df_use,mat_df_NA)

mat_df$OS.status<-info_prog1$OS.status[match(mat_df$`Admission ID`,info_prog1$`Admission ID`)]
mat_df$OS.Days<-info_prog1$OS.Days[match(mat_df$`Admission ID`,info_prog1$`Admission ID`)]
mat_df$PFS.status<-info_prog1$PFS.status[match(mat_df$`Admission ID`,info_prog1$`Admission ID`)]
mat_df$PFS.Days<-info_prog1$PFS.Days[match(mat_df$`Admission ID`,info_prog1$`Admission ID`)]
#发现有一些没有匹配到Admission ID， 因此没有预后信息，这些样本删掉，不参与后续的cox分析

# setdiff(row.names(mat_df),info.all$FileName)
# 
# info.all$`Admission ID`<-info_pan$`Admission ID`[match(info.all$FileName,info_pan$FileName)]


mat_df1<-mat_df[!is.na(mat_df$OS.status)|!is.na(mat_df$OS.Days)|!is.na(mat_df$PFS.status)|!is.na(mat_df$PFS.Days),]#766

# setdiff(row.names(mat_df),info_pan$FileName)

#OS
mat_df_os<-mat_df1[,c((ncol(mat_df1)-3),(ncol(mat_df1)-2),1:(ncol(mat_df1)-5))]
mat_df_os$OS.status<-as.numeric(mat_df_os$OS.status)
cox_results_OS <- data.frame(
  protein = character(),
  coef = numeric(),
  hr = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
for (i in 3:ncol(mat_df_os)) {
  # 构建Cox模型公式
  formula <- as.formula(paste("Surv(OS.Days, OS.status) ~", names(mat_df_os)[i]))
  
  # 拟合Cox模型
  cox_model <- coxph(formula, data = mat_df_os)
  
  # 提取结果
  cox_summary <- summary(cox_model)
  
  cox_results_OS <- rbind(cox_results_OS, data.frame(
    protein = names(mat_df_os)[i],
    coef = cox_summary$coefficients[1, 1],
    hr = cox_summary$coefficients[1, 2],
    p_value = cox_summary$coefficients[1, 5]
  ))
  
}
# 多重检验校正
cox_results_OS$padj <- p.adjust(cox_results_OS$p_value, method = "BH")
# mat_df_os$OS.status
# mat_df_os$OS.Days
#根据hr值区分上下调，
cox_results_OS$sig_type<-NA
cox_results_OS$sig_type[cox_results_OS$padj<0.05&cox_results_OS$hr>1]<-"up"
cox_results_OS$sig_type[cox_results_OS$padj<0.05&cox_results_OS$hr<1]<-"down"

write.csv(cox_results_OS,"20251024_TPHP_pancancer_cox_OS_add_sig.csv",row.names = F)


#PFS
mat_df_pfs<-mat_df1[,c((ncol(mat_df1)-1),(ncol(mat_df1)),1:(ncol(mat_df1)-5))]
mat_df_pfs$PFS.status<-as.numeric(mat_df_pfs$PFS.status)
# mat_df_pfs
cox_results_PFS <- data.frame(
  protein = character(),
  coef = numeric(),
  hr = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
for (i in 3:ncol(mat_df_pfs)) {
  # 构建Cox模型公式
  formula <- as.formula(paste("Surv(PFS.Days, PFS.status) ~", names(mat_df_pfs)[i]))
  
  # 拟合Cox模型
  cox_model <- coxph(formula, data = mat_df_pfs)
  
  # 提取结果
  cox_summary <- summary(cox_model)
  
  cox_results_PFS <- rbind(cox_results_PFS, data.frame(
    protein = names(mat_df_pfs)[i],
    coef = cox_summary$coefficients[1, 1],
    hr = cox_summary$coefficients[1, 2],
    p_value = cox_summary$coefficients[1, 5]
  ))
  
}
# 多重检验校正
cox_results_PFS$padj <- p.adjust(cox_results_PFS$p_value, method = "BH")
# mat_df_pfs$PFS.status
# mat_df_pfs$PFS.Days
#根据hr值区分上下调，
cox_results_PFS$sig_type<-NA
cox_results_PFS$sig_type[cox_results_PFS$padj<0.05&cox_results_PFS$hr>1]<-"up"
cox_results_PFS$sig_type[cox_results_PFS$padj<0.05&cox_results_PFS$hr<1]<-"down"
write.csv(cox_results_PFS,"20251024_TPHP_pancancer_cox_PFS_add_sig.csv",row.names = F)


#和PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_3873signatures_up_down.xlsx中的结果比较
compair_mat<-rio::import('PTDB_pan-cancer_6377_proteins_1014_coreSamples_5clusters_3873signatures_up_down.xlsx')
compair_mat$match<-paste(compair_mat$Protein,compair_mat$Type,sep="_")
cox_results_OS$match<-paste(cox_results_OS$protein,cox_results_OS$sig_type,sep="_")
cox_results_PFS$match<-paste(cox_results_PFS$protein,cox_results_PFS$sig_type,sep="_")

OS_overlap<-cox_results_OS[cox_results_OS$match%in%compair_mat$match,]
write.csv(OS_overlap,"20251027_TPHP_pancancer_cox_OS_overlap_with_3873signatures_up_down.csv",row.names = F)
OS_overlap<-read.csv("20251027_TPHP_pancancer_cox_OS_overlap_with_3873signatures_up_down.csv",header = T)
OS_overlap$cluster<-compair_mat$cluster[match(OS_overlap$match,compair_mat$match)]
write.csv(OS_overlap,"20251028_TPHP_pancancer_cox_OS_overlap_with_3873signatures_up_down_add_prot_cluster.csv",row.names = F)

PFS_overlap<-cox_results_PFS[cox_results_PFS$match%in%compair_mat$match,]
write.csv(PFS_overlap,"20251027_TPHP_pancancer_cox_PFS_overlap_with_3873signatures_up_down.csv",row.names = F)

PFS_overlap<-read.csv("20251027_TPHP_pancancer_cox_PFS_overlap_with_3873signatures_up_down.csv",header = T)
PFS_overlap$cluster<-compair_mat$cluster[match(PFS_overlap$match,compair_mat$match)]
write.csv(PFS_overlap,"20251028_TPHP_pancancer_cox_PFS_overlap_with_3873signatures_up_down_add__prot_cluster.csv",row.names = F)






