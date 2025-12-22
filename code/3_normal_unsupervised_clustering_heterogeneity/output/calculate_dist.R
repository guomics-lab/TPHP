calculate_pairwise_tissue_distances <- function(mat, index, sample_info) {
  # (输入验证部分保持不变)
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("mat必须是矩阵或数据框")
  }
  
  if (!is.character(index) || length(index) != 1) {
    stop("index必须是单个字符串")
  }
  
  if (!(index %in% colnames(sample_info))) {
    stop(paste0("样本信息必须包含'", index, "'列"))
  }
  
  if (is.null(colnames(mat)) || is.null(rownames(sample_info))) {
    stop("矩阵必须有列名，样本信息必须有行名")
  }
  
  common_samples <- intersect(colnames(mat), rownames(sample_info))
  
  if (length(common_samples) == 0) {
    stop("矩阵列名与样本信息行名没有匹配项")
  }
  
  if (length(common_samples) < ncol(mat) * 0.5) {
    warning(paste("只有", length(common_samples), "个样本匹配，共", 
                  ncol(mat), "个样本"))
  }
  
  # 过滤匹配的样本
  mat <- mat[, common_samples, drop = FALSE]
  sample_info <- sample_info[common_samples, , drop = FALSE]
  
  # 移除包含NA的基因
  na_genes <- apply(mat, 1, function(x) any(is.na(x)))
  if (sum(na_genes) > 0) {
    warning(paste("移除", sum(na_genes), "个包含NA值的基因"))
    mat <- mat[!na_genes, , drop = FALSE]
  }
  
  if (nrow(mat) == 0) {
    stop("移除NA后没有剩余的基因")
  }
  
  # 数据准备
  mat_t <- t(mat)
  tissue_vector <- sample_info[[index]]
  tissues <- sort(unique(tissue_vector)) 
  n_tissues <- length(tissues)
  
  if (n_tissues < 2) {
    stop("需要至少2种组织类型进行比较")
  }
  
  # 计算所有样本间的距离矩阵
  all_dist <- as.matrix(dist(mat_t, method = "euclidean"))
  
  # 组织样本索引准备
  tissue_samples <- list()
  for (tissue in tissues) {
    samples <- which(tissue_vector == tissue)
    if (length(samples) > 0) {
      tissue_samples[[tissue]] <- samples
    } else {
      warning(paste0("组织 '", tissue, "' 没有样本"))
    }
  }
  
  # 初始化结果列表
  results <- list()
  
  # 距离计算核心循环
  for (i in seq_along(tissues)) {
    tissue1 <- tissues[i]
    if (!(tissue1 %in% names(tissue_samples))) next
    
    idx1 <- tissue_samples[[tissue1]]
    samples1 <- rownames(sample_info)[idx1]
    n1 <- length(idx1)
    
    for (j in i:n_tissues) {
      tissue2 <- tissues[j]
      if (!(tissue2 %in% names(tissue_samples))) next
      
      idx2 <- tissue_samples[[tissue2]]
      samples2 <- rownames(sample_info)[idx2]
      n2 <- length(idx2)
      
      # 初始化距离值和比较次数
      dist_values <- numeric(0)
      median_dist <- NA_real_
      
      if (i == j) {
        # 同一组织内的距离（只取上三角）
        if (n1 > 1) {
          dist_values <- all_dist[samples1, samples2][upper.tri(all_dist[samples1, samples2])]
          n_comp <- n1 * (n1 - 1) / 2
        } else {
          # ***修改点：样本数 n1=1 时，不跳过，但没有比较次数***
          n_comp <- 0
          # dist_values 为空，median_dist 保持 NA
        }
      } else {
        # 不同组织间的距离（取完整子矩阵）
        dist_values <- as.vector(all_dist[samples1, samples2])
        n_comp <- n1 * n2
      }
      
      # 仅在有成对距离值时计算中位数
      if (length(dist_values) > 0) {
        median_dist <- median(dist_values, na.rm = TRUE)
      }
      
      # 存储结果
      results[[length(results) + 1]] <- data.frame(
        tissue1 = tissue1,
        tissue2 = tissue2,
        is_within_tissue = (i == j),
        median_distance = median_dist, # 当 n_comp=0 时，这将是 NA
        n_samples_tissue1 = n1,
        n_samples_tissue2 = n2,
        n_comparisons = n_comp,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # 将结果列表合并为数据框
  return(do.call(rbind, results))
}



###################only use tissue_enriched proteins

# #####绘制major组织的hclust
# # 更完整的版本，确保包含所有组织到自身的距离（0）
# library(tidyr)
# library(dplyr)

#  # 计算所有样本间的距离矩阵
# pm<-read.csv("../2_batch_effect_evaluation_correction/output/batch_corrected_data_step2_limma3.csv",row.names = 1,check.names = F, stringsAsFactors = F)
# pm<-as.data.frame(pm)
# info<-read_xlsx("../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V6.xlsx")
# sample_info<-info[(info$sample_type=="N" & info$Tissue_heterogeneity_abnormal == FALSE & info$Low_protein_IDs == FALSE),]
# data<-apply(t(pm[,sample_info$FileName]),c(1,2),as.numeric)
# dim(data)
# head(colnames(data))
# sample_info[grep("sinus",sample_info$anatomical_classification),]
# rownames(sample_info)<-sample_info$FileName
# group_labels<-data.frame(major_category = sample_info$anatomical_classification,
#                          sub_category = sample_info$tissue_name)
# rownames(group_labels)<-sample_info$FileName
# data_enriched<-data[,colnames(data) %in% unique(pm_data$ts.Enrich_rm_contam$Protein)]
# dim(data_enriched)

# all_dist <- calculate_pairwise_tissue_distances(t(data_enriched), "major_category", group_labels)
# dim(all_dist)
# min(all_dist$median_distance,na.rm=T)
# all_dist[is.na(all_dist$median_distance),]$median_distance<-0
# max(all_dist$median_distance)
# dim(all_dist)
# df<-data.frame(tissue1=all_dist$tissue1,
#                  tissue2=all_dist$tissue2,
#                  distance=all_dist$median_distance)

# # 获取所有唯一的组织名称
# all_tissues <- unique(c(df$tissue1, df$tissue2))
# all_tissues
# length(all_tissues)

# # 步骤1：补全对称
# reverse_df <- df %>%
#   rename(tissue1 = tissue2, tissue2 = tissue1)

# complete_df <- bind_rows(df, reverse_df) %>%
#   distinct(tissue1, tissue2, .keep_all = TRUE)

# # 步骤2：添加对角线（组织到自身的距离为0）
# diagonal_df <- data.frame(
#   tissue1 = all_tissues,
#   tissue2 = all_tissues,
#   distance = 0
# )

# final_df <- bind_rows(complete_df, diagonal_df) %>%
#   distinct(tissue1, tissue2, .keep_all = TRUE)

# # 步骤3：转换为宽格式
# distance_matrix <- final_df %>%
#   pivot_wider(
#     names_from = tissue2,
#     values_from = distance
#   ) %>%
#   rename(tissue = tissue1) %>%
#   # 确保列顺序与行顺序一致
#   select(tissue, all_of(all_tissues))


# # 可选：转换为矩阵格式
# matrix_format <- as.matrix(distance_matrix[, -1])
# rownames(matrix_format) <- distance_matrix$tissue

# #画hclust

# distance_matrix_based_hclust <- function(
#   Xmed, # distance matrix (tissue_type x tissue type)
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
  
#   # 确保输入是距离矩阵格式
#   # 如果输入是数据框，转换为矩阵
#   if (is.data.frame(Xmed)) {
#     # 假设第一列是组织名称，其余列是距离值
#     rownames(Xmed) <- Xmed[, 1]
#     Xmed <- as.matrix(Xmed[, -1])
#   }
  
#   # 确保是对称的距离矩阵
#   if (!isSymmetric(Xmed)) {
#     warning("Distance matrix is not symmetric. Making it symmetric by taking the upper triangle.")
#     Xmed[lower.tri(Xmed)] <- t(Xmed)[lower.tri(Xmed)]
#   }
  
#   # 转换为dist对象
#   res.dist <- as.dist(Xmed)
  
#   # 执行层次聚类
#   res.hcls <- stats::hclust(res.dist, method = linkage)
#   res.dend <- as.dendrogram(res.hcls)
  
#   # 确定切割高度和聚类数量
#   if (is.null(dist_cut)) {
#     # 如果没有提供切割高度，使用平均距离
#     dist_cut <- mean(res.dist, na.rm = TRUE)
#   }
#   k <- length(unique(stats::cutree(res.hcls, h = dist_cut)))
  
#   # 标签转换
#   old_lab <- labels(res.dend)
#   if (is.null(label_converter)) {
#     new_lab <- old_lab
#   } else if (is.function(label_converter)) {
#     new_lab <- label_converter(old_lab)
#   } else {
#     # 假设是命名字符向量；对任何NA保留原始标签
#     new_lab <- unname(label_converter[old_lab])
#     new_lab[is.na(new_lab)] <- old_lab
#   }
  
#   # 生成颜色
#   set.seed(1)
#   branch_colors <- sample(colorspace::rainbow_hcl(max(10, k), c = 70, l = 50))
  
#   # 绘图
#   if (plot) {
#     # 改变右边距以适应标签
#     op <- par(no.readonly = TRUE)
#     on.exit(par(op), add = TRUE)
#     cur <- par("mar") # c(bottom, left, top, right)
#     par(mar = c(cur[1], cur[2], cur[3], 8))
#     par(xpd = NA)
    
#     # 创建带颜色的树状图
#     p_dend <- res.dend %>%
#       set("labels", new_lab) %>%
#       set("labels_cex", 0.7) %>%
#       set("branches_k_color", value = branch_colors, k = k) %>%
#       color_branches(h = dist_cut, col = branch_colors[1:k])
    
#     # 绘制树状图
#     plot(p_dend, horiz = TRUE, axes = TRUE,
#          xlab = "Euclidean Distance", main = main)
    
#     # 添加切割线
#     abline(v = dist_cut, lty = 2, col = "red")
#     text(x = dist_cut, y = -0.5, labels = paste("Cutoff =", round(dist_cut, 3)), 
#          pos = 4, cex = 0.8, col = "red")
    
#     # 可选：添加彩色条显示聚类结果
#     cluster_assignments <- cutree(res.hcls, h = dist_cut)
#     colored_bars(colors = branch_colors[cluster_assignments], 
#                  dend = p_dend, 
#                  rowLabels = "Cluster",
#                  horiz = TRUE,
#                  cex.rowLabels = 0.7)
#   }
  
#   # 返回结果
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



# # pdf('output/TPHP_normal_tissue_cluster_20251026.pdf', width = 3, height = 7.5)
# # hclust_list <- lapply(names(mat_med_list) %>% setNames(names(mat_med_list)), function(nm){
# #   message('Analysis of ', nm, '...\n')
# #   Xmed <- mat_med_list[[nm]]
# #   spearman_based_hclust(Xmed, main = nm, label_converter = label_converter)
# # })
# # graphics.off()

# # pdf('output/TPHP_normal_tissue_cluster_20251026.pdf', width = 3, height = 7.5)
# pdf("output/tissue_clster_V6_step2_limma.pdf")
# res<-distance_matrix_based_hclust(Xmed = matrix_format,
#                                main = "Tissue Clustering based on Euclidean Distance",
#                                dist_cut = 300,
#                                linkage = "complete",plot = T)
# dev.off()



