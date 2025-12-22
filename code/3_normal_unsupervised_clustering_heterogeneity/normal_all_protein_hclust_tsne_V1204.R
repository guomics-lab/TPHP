# Working directory
setwd(r"(\\172.16.13.136\tphp\code\3_normal_unsupervised_clustering_heterogeneity)") 

# Source local helpers
source("../source/source_code.R")


# 1.Load environment ------
## Helper functions -------
# Martin Sill's method
# https://github.com/mwsill/mnp_training/blob/master/R/RSpectra_pca.R
require(RSpectra)
require(Rtsne)

prcomp_svds <-
  function(x, retx = TRUE, center = TRUE, scale. = FALSE, tol = NULL, k=10, ...){
    chkDots(...)
    x <- as.matrix(x)
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if(any(sc == 0))
      stop("cannot rescale a constant/zero column to unit variance")
    s <- svds(x, k)
    s$d <- s$d / sqrt(max(1, nrow(x) - 1))
    if (!is.null(tol)) {
      ## we get rank at least one even for a 0 matrix.
      rank <- sum(s$d > (s$d[1L]*tol))
      if (rank < ncol(x)) {
        s$v <- s$v[, 1L:rank, drop = FALSE]
        s$d <- s$d[1L:rank]
      }
    }
    dimnames(s$v) <-
      list(colnames(x), paste0("PC", seq_len(ncol(s$v))))
    r <- list(sdev = s$d, rotation = s$v,
              center = if(is.null(cen)) FALSE else cen,
              scale = if(is.null(sc)) FALSE else sc)
    if (retx) r$x <- x %*% s$v
    class(r) <- "prcomp"
    r
  }

# 2.Data preparation -----
## 2.1 Read files -----
# source("//172.16.13.136/share/members/jiangwenhao/code/myQC.R")
reported_protinfo <- rio::import("../0_process_DIA-NN_data/output/protein_info_from_reports_V1009.csv")
info <- rio::import("../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V6.xlsx")

# DIA-NN protein-group matrices
# Original (log2 expected in file); if raw, we transform below
file_list <- list.files(path = 'input', pattern = 'batch_corrected|12754proteins_2997samples', full.names = T)
names(file_list) <- c(
  'qimp.cmb.mo', 'qimp.cmb', 'qimp.lma', 'qimp.lma.s2'
)
pm_list <- lapply(file_list, function(X){
  message('Reading file ', X, '...\n')
  rio::import(X) %>% column_to_rownames('V1')
})

sapply(pm_list, dim)
# #     qimp.cmb.mo qimp.cmb qimp.lma qimp.lma.s2
# [1,]       12797    12797    12797       12797
# [2,]        2957     2957     2957        2957
sapply(pm_list, function(X) quantile(as.matrix(X)))
#    qimp.cmb.mo  qimp.cmb  qimp.lma qimp.lma.s2
# 0%      5.921994 -1.310841  3.949413    3.161037
# 25%     6.407567  6.388891  6.386089    6.401242
# 50%     6.490974  6.971537  6.892888    7.064453
# 75%    12.365733 12.357216 12.365974   12.446083
# 100%   22.872299 33.539784 24.994231   24.906838

## 2.2 Filter to "N" samples -----
info1 <- info %>% dplyr::filter(sample_type == 'N' & Tissue_heterogeneity_abnormal == FALSE & Low_protein_IDs==FALSE) 
dim(info1)# 514 34
mat_list <- lapply(pm_list, function(X){ # Normal only
  X[, info1$FileName, drop = FALSE]
})

# metadata filter
info1 %>% count(anatomical_classification) %>% arrange(n) %>% View
info2 <- info1 %>% dplyr::select(FileName, tissue_name, anatomical_location, anatomical_classification,
                                 class_abbr)
dim(info2) # 514  5

# class_abbr to label
label_converter <- info2 %>% distinct(class_abbr, anatomical_classification) %>% 
  mutate(label = str_glue("{class_abbr} ({anatomical_classification})")) %>% 
  pull(label, class_abbr)

## 2.3 Median protein matrices -----
mat_med_list <- lapply(mat_list, function(X){
  Xmed <- 2^X %>% t() %>% as.data.frame() %>% rownames_to_column('FileName') %>%
    inner_join(info2 %>% select(FileName, class_abbr)) %>% 
    pivot_longer(cols = -c(FileName, class_abbr)) %>% 
    pivot_wider(id_cols = class_abbr, names_from = 'name', values_from = 'value',
                values_fn = median) %>% 
    as.data.frame() %>% 
    column_to_rownames('class_abbr') %>% 
    log2()
  return(Xmed)
})
saveRDS(mat_med_list, 'output/mat_med_list.rds')



## 2.4 Median matrices (with NA) ----
mat_mask <- rio::import('../1_NA_filter_imputation/output/data_matrix_filtered_05NA_2957_12797.csv')


## 3.1 Hierarchical clustering  -------
### helper functions --------
spearman_based_hclust <- function(
    Xmed, # median matrix (tissue_type x protein)
    main = NULL,
    label_converter = NULL,  # can be a *named character vector* or a *function*
    rho_min = 0.5,
    linkage = "complete", # hclust linkage method
    plot = TRUE
){
  if (!requireNamespace("dendextend", quietly = TRUE))
    stop("Package 'dendextend' is required.")
  if (!requireNamespace("colorspace", quietly = TRUE))
    stop("Package 'colorspace' is required.")
  require(dendextend)
  
  res.cor  <- suppressWarnings(cor(t(Xmed), method = "spearman", use = "pairwise.complete.obs"))
  res.dist <- as.dist(1 - res.cor)              # 1 - Spearman's rho
  res.hcls <- stats::hclust(res.dist, method = linkage)
  res.dend <- as.dendrogram(res.hcls)
  
  # cutoff & cluster count
  h_cut <- 1 - rho_min
  k     <- length(unique(stats::cutree(res.hcls, h = h_cut)))
  
  # labels convert
  old_lab <- labels(res.dend)
  if (is.null(label_converter)) {
    new_lab <- old_lab
  } else if (is.function(label_converter)) {
    new_lab <- label_converter(old_lab)
  } else {
    # assume named character vector; keep originals for any NA
    new_lab <- unname(label_converter[old_lab])
    new_lab[is.na(new_lab)] <- old_lab
  }
  
  set.seed(1); branch_colors <- sample(colorspace::rainbow_hcl(10, c = 70, l  = 50))
  
  # plot
  if (plot) {
    # change right margin
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    cur <- par("mar") # c(bottom, left, top, right)
    par(mar = c(cur[1], cur[2], cur[3], 6))
    par(xpd = NA)
    
    p_dend <- res.dend %>%
      set("labels", new_lab) %>%
      set("labels_cex", 0.7) %>%
      # set("branches_col", value = rep_len(branch_colors, k), h = h_cut) %>% 
      # set("labels_col", value = branch_colors, k) %>%
      set("branches_k_color", value = branch_colors, k) %>%
      color_branches(h = h_cut, col = rep_len(branch_colors, k))
    # k10 <- cutree(res.dend, k = c(5, 10, 20))
    # colored_bars(cbind(k10[,ncol(k10):1], colorspace::rainbow_hcl(51, c = 70, l  = 50)), res.dend, rowLabels = c(paste0("k = ", rev(c(5, 10, 20))), "Sample type"), horiz = T)
    
    plot(p_dend, horiz = TRUE, axes = TRUE,
         xlab = "1-Spearman's rho", main = main)
    # abline(v = h_cut, lty = 2)
  }
  
  list(
    res = list(
      res.cor   = res.cor,
      res.dist  = res.dist,
      res.hcls  = res.hcls,
      res.dend  = res.dend,
      leaf_labels = new_lab
    ),
    params = c(rho_min = rho_min, h_cut = h_cut, k = k, linkage = linkage)
  )
}

pdf('output/TPHP_normal_tissue_cluster_20251204.pdf', width = 3, height = 7.5)
hclust_list <- lapply(names(mat_med_list) %>% setNames(names(mat_med_list)), function(nm){
  message('Analysis of ', nm, '...\n')
  Xmed <- mat_med_list[[nm]]
  spearman_based_hclust(Xmed, main = nm, label_converter = label_converter)
})
graphics.off()
saveRDS(hclust_list, 'normal_tissue_hclust_list.rds')

rio::export(lapply(hclust_list, function(X) X$res$res.cor %>% as.data.frame() %>% rownames_to_column('\\')),
            'output/TPHP_normal_tissue_cluster_20251204_spearmanRho.xlsx')

## 3.2 Hierarchical clustering (median of sample distance matrix) -----


# 4.All protein analysis ------------
## 4.1 PCA-based t-SNE --------
library(RSpectra)
library(Rtsne)
normal_tsne_list <- lapply(names(mat_list), function(nm){
    # nm=2
  message('PCA-based t-SNE of matrix ', nm, '...\n')
  X <- mat_list[[nm]]
  sd<-apply(X, 1, function(v) sd(v, na.rm = T))
  dim(X) # 12754 560
  X<-X[sd>0, ]
  DF <- X %>% t()
  M <- apply(X, 1, function(v) {
    (v - mean(v, na.rm = T)) / sd(v, na.rm = T)
  })


  dim(M) # 514 12779
  m1 <- prcomp(M)
  # summary(m1)
  # m1$x; # pca
  # m1$rotation # eigen$vectors
  # m1$sdev ^ 2 # eigen$values
  
  # plot(cumsum(m1$sdev ^ 2)/sum(m1$sdev ^ 2)) # 500+
  # (cumsum(m1$sdev ^ 2)/sum(m1$sdev ^ 2))[94] # 0.7687
  # which(cumsum(m1$sdev ^ 2)/sum(m1$sdev ^ 2) <= 0.763) %>% max() # 90
  # which(cumsum(m1$sdev ^ 2)/sum(m1$sdev ^ 2) <= 0.9) %>% max() # 246
  # (cumsum(m1$sdev ^ 2)/sum(m1$sdev ^ 2))[420] # 0.97
  
  k <- which(cumsum(m1$sdev ^ 2)/sum(m1$sdev ^ 2) <= 0.9) %>% max()
  
  y <- as.factor(colnames(DF))
  pca <- prcomp_svds(DF, retx = T, center = F, scale. = T, tol = NULL, k = k)
  
  
  # calculate tSNE
  set.seed(2020)
  res <- Rtsne(pca$x,pca=F,max_iter=2500,theta=0,verbose=T, perplexity = 20
               # perplexity must be <= floor((nrow(pca$x) - 1) / 3)
  )
  
  tsne <- res$Y %>% as.data.frame
  colnames(tsne) <- c('tsne1', 'tsne2')
  hc <- hclust(dist(tsne))
  
  res.hc <- tsne %>% 
    mutate(hclust = factor(cutree(hc, 20)),
           hclust = str_c('Cluster ', hclust),
           FileName = colnames(X)) %>% 
    inner_join(info1) %>% 
    set_rownames(.$FileName)
  
  hc.cent <- res.hc %>%
    group_by(hclust) %>%
    dplyr::select(tsne1, tsne2) %>%
    summarize_all(mean)
  
  return(list(res.pca = pca,
              pca.k = k,
              res.tsne = res,
              hc = hc,
              res.hc = res.hc,
              hc.cent = hc.cent))
})

names(normal_tsne_list) <- names(mat_list)
saveRDS(normal_tsne_list, 'output/normal_tsne_list.rds')


for(nm in names(normal_tsne_list)){
  ntsne <- normal_tsne_list[[nm]]
  res.hc <- ntsne$res.hc
  my_colors <- c('#368FC6', '#55BC6D', '#D8894E', '#7C0823', '#9D7DDB', '#BCBD22', '#CF4E9C', '#000000', '#888888')
  my_shapes <- c(
    1, 2, 6, 3, 4, 15, 17, 18, 19
  )
  
  shape_color_pairs <- expand.grid(my_shapes, my_colors)
  
  # anatomical classification
  p1 <- ggplot(res.hc) + 
    geom_point(aes(x = tsne1, y = tsne2, shape = anatomical_classification, color = anatomical_classification), size = 2)+
    scale_shape_manual(values = shape_color_pairs[, 1])+
    scale_color_manual(values = as.character(shape_color_pairs[, 2]))+
    theme_bw() +
    theme(text = element_text(size = 10))
  p2 <- p1 + theme(legend.position = 'none')
  p <- ggarrange(p1, p2, nrow = 1, ncol = 2)
  ggsave(str_c('output/20251204_normalAll_tissue_pca_based_tsne_class_', nm, '.pdf'), p, width = 12, height = 6)
  
  # patient
  p1 <- ggplot(res.hc) + 
    geom_point(aes(x = tsne1, y = tsne2, shape = patient_ID, color = patient_ID), size = 2)+
    scale_shape_manual(values = shape_color_pairs[, 1])+
    scale_color_manual(values = as.character(shape_color_pairs[, 2]))+
    theme_bw() +
    theme(text = element_text(size = 10))
  p2 <- p1 + theme(legend.position = 'none')
  p <- ggarrange(p1, p2, nrow = 1, ncol = 2)
  ggsave(str_c('output/20251204_normalAll_tissue_pca_based_tsne_patient_', nm, '.pdf'), p, width = 12, height = 6)
  
  # Batchg
  res.hc %<>%
    arrange(date, instrument) %>% 
    mutate(Batch = str_c(instrument, str_sub(date, 1, 6)),
           Batch = factor(Batch, levels = unique(Batch)))
  p1 <- ggplot(res.hc) + 
    geom_point(aes(x = tsne1, y = tsne2, color = Batch), size = 2)+
    scale_color_manual(values = RColorBrewer::brewer.pal(n = length(unique(res.hc$Batch)), name = 'Spectral'))+
    theme_bw() +
    theme(text = element_text(size = 10))
  p2 <- p1 + theme(legend.position = 'none')
  p <- ggarrange(p1, p2, nrow = 1, ncol = 2)
  ggsave(str_c('output/20251204_normalAll_tissue_pca_based_tsne_Batch_', nm, '.pdf'), p, width = 12, height = 6)
}

rio::export(lapply(normal_tsne_list, function(X) X$res.hc),
            'output/20251204_normalAll_tissue_pca_based_tsne.xlsx')


