rm(list=ls())
setwd(r"(\\172.16.13.136\tphp\code\5_1_compare_with_RNA\)")
library(readxl)
library(gridExtra)
library(dplyr)
library(tidyr)
library(reshape2)
#0.read.files
# pm_data<-readRDS("//172.16.13.136/tphp/code.20251201.archieved/4_tissue_specificity_analysis/20251031_tissue_Specificity_naive_05minImputated.rds")
tphp_info<-read_xlsx("./input/20251009_PUH_sample_information_3005files_V6.xlsx")
tphp_info <- tphp_info %>% dplyr::filter(sample_type == 'N' & Tissue_heterogeneity_abnormal == FALSE & Low_protein_IDs==FALSE) 
dim(tphp_info)# 514 34



# pm_data<-readRDS("//172.16.13.136/tphp/code/4_tissue_specificity_analysis/20251031_tissue_Specificity_naive_05minImputated.rds")
# tphp_data0<-pm_data$metadata
tphp_data0<-read.csv("./input/batch_corrected_data_step2_limma3.csv",header = T,check.names = F, stringsAsFactors = F,row.names = 1)


# tphp_data0<-read.csv("input/mapped_pg_matrix_1780_14062.csv",check.names = F, stringsAsFactors = F,row.names = 1)
dim(tphp_data0)
tphp_data1<-tphp_data0[,names(tphp_data0)%in%tphp_info$FileName]
# tphp_data1<-tphp_data1[,names()]
tphp_mapping<-read.csv("./input/tphp_mapped_matrix.csv",check.names = F, stringsAsFactors = F,row.names = 1)
tphp_data<-merge(tphp_mapping[,c(1:3)],tphp_data1, by.x = "UniProtKB.Swiss.Prot.ID",by.y = "row.names")
dim(tphp_data) #13742  517
head(colnames(tphp_data))
HPA_data<-read.delim("input/rna_tissue_hpa.tsv",check.names = F,stringsAsFactors = F)
gtex_data<-read.delim("input/rna_tissue_detail_gtex.tsv",check.names = F,stringsAsFactors = F)
# tissue_mapping<-read_xlsx("input/HPA_GTEx_tissue_mapping_clean.xlsx")
tissue_mapping<-rio::import("./input/20251104_tissue_compare.xlsx")
#####
#pancreatic body and pancreatic duct,  pancreatic, these three removed
# tissue_mapping<-tissue_mapping[!tissue_mapping$TPHP_tissue%in%c("pancreatic body","pancreatic duct","pancreas"),]

colnames(tissue_mapping)
tissue_mapping_all<-tissue_mapping[apply(tissue_mapping[,c("tissue_type","HPA_tissue","gtex_tissue")],1,function(x) sum(is.na(x))==0),]
dim(tissue_mapping_all) #42 7
tissue_mapping_all

common_tissues <- unique(tissue_mapping_all$tissue_type)
common_tissues

common_genes <- intersect(tphp_data$Gene.stable.ID, intersect(HPA_data$Gene,gtex_data$Gene))
length(common_genes) #12478

expr_proprietary <- apply(tphp_data[,-c(1:3)],c(1,2),as.numeric)
rownames(expr_proprietary) <- tphp_data$Gene.stable.ID  
expr_published1 <- pivot_wider(HPA_data[,c("Gene","Tissue","nTPM")], names_from = "Tissue", values_from = "nTPM") %>%
    tibble::column_to_rownames("Gene")%>%as.data.frame()
dim(expr_published1)
expr_published2 <-  pivot_wider(gtex_data[,c("Gene","Source tissue","nTPM")], names_from = "Source tissue", values_from = "nTPM") %>%
    tibble::column_to_rownames("Gene")%>%as.data.frame()
dim(expr_published2)
colnames(expr_published2)
colnames(expr_published2)<-tolower(colnames(expr_published2))
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


# Define overlapping proteins and tissues
overlapping_proteins <- common_genes
overlapping_tissues <- common_tissues

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
tphp_info$anatomical_classification[tphp_info$tissue_name=="sigmoid colon"]<-"sigmoid colon"
tphp_info$anatomical_classification[tphp_info$tissue_name=="transverse colon"]<-"transverse colon"

sample_info_proprietary<-tphp_info[,c("FileName","anatomical_classification")]
colnames(sample_info_proprietary) <- c("sample_id", "tissue_type")
unique(sample_info_proprietary$tissue_type)
data_proprietary <- filter_intensity_data(expr_proprietary, sample_info_proprietary, 
                                           overlapping_proteins, unique(as.character(tissue_mapping_all$tissue_type)))

head(data_proprietary$sample_info)
mmm<-data_proprietary$sample_info
unique(mmm$tissue_type)
data_published1 <- filter_published2_data(expr_published1,  overlapping_proteins,unique(as.character(tissue_mapping_all$HPA_tissue)))
data_published2 <- filter_published2_data(expr_published2,  overlapping_proteins,unique(as.character(tissue_mapping_all$gtex_tissue)))
colnames(data_published2) 
dim(data_proprietary$expression)
dim(data_published1)
dim(data_published2)
write.csv(data_published1, "output/data_published1_20251204.csv",row.names = T)  
write.csv(data_published2, "output/data_published2_20251204.csv",row.names = T)
write.csv(data_proprietary$expression, "output/data_proprietary_20251204.csv",row.names = T)
write.csv(data_proprietary$sample_info, "output/data_proprietary_sample_info_20251204.csv",row.names = F)
# =============================================================================
# 2. DATA NORMALIZATION AND PREPROCESSING
# =============================================================================

cat("Data ranges before normalization:\n")
cat(paste("Proprietary (intensity):", round(min(data_proprietary$expression, na.rm = TRUE), 3), 
          "to", round(max(data_proprietary$expression, na.rm = TRUE), 3), "\n"))
cat(paste("Published1 (intensity means):", round(min(data_published1, na.rm = TRUE), 3), 
          "to", round(max(data_published1, na.rm = TRUE), 3), "\n"))
cat(paste("Published2 (intensity means):", round(min(data_published2, na.rm = TRUE), 3), 
          "to", round(max(data_published2, na.rm = TRUE), 3), "\n\n"))

# Normalization strategy:
#  all Log2 transform + Z-score normalization

# Function for Z-score normalization
zscore_normalize <- function(expr_matrix) {
  expr_matrix <- as.matrix(expr_matrix)
  
  # Z-score normalization: (x - mean) / sd for each protein across samples
  normalized <- t(scale(t(expr_matrix), center = TRUE, scale = TRUE))
  
  # Handle proteins with zero variance
  normalized[is.na(normalized)] <- 0
  
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
# median_center <- function(expr_matrix) {
#   expr_matrix <- as.matrix(expr_matrix)
  
#   # Subtract median from each sample (column-wise centering)
#   medians <- apply(expr_matrix, 2, median, na.rm = TRUE)
#   centered <- sweep(expr_matrix, 2, medians, "-")
  
#   return(centered)
# }

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
data_published1_log <- log2(data_published1 + 1)
# Z-score normalization across proteins (row-wise)
data_published1_normalized <- zscore_normalize(data_published1_log)

# Process published2 data (intensity means matrix)
cat("Processing published2 data (intensity means matrix)...\n")
# Log2 transform intensity means
data_published2_log <- log2(data_published2 + 1)
# Z-score normalization across proteins (row-wise)
data_published2_normalized <- zscore_normalize(data_published2_log)

cat("Data ranges after normalization (Z-scores):\n")
cat(paste("Proprietary:", round(min(data_proprietary$expression_normalized, na.rm = TRUE), 3), 
          "to", round(max(data_proprietary$expression_normalized, na.rm = TRUE), 3), "\n"))
cat(paste("Published1:", round(min(data_published1_normalized, na.rm = TRUE), 3), 
          "to", round(max(data_published1_normalized, na.rm = TRUE), 3), "\n"))
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
proprietary_tissue_means <- calc_tissue_means(data_proprietary$expression_normalized, as.data.frame(data_proprietary$sample_info))
dim(proprietary_tissue_means) #12478  20
# class(data_proprietary$sample_info)
names(proprietary_tissue_means)
# =============================================================================
# 4. DATASET SIMILARITY ANALYSIS BY TISSUE TYPE
# =============================================================================
tissue_mapping_all$Prop_Pub1<-"NA"
tissue_mapping_all$Prop_Pub2<-"NA"
tissue_mapping_all$Pub1_Pub2<-"NA"
tissue_mapping_all<-as.data.frame(tissue_mapping_all)
# rownames(tissue_mapping_all)<-tissue_mapping_all$TPHP_anatomical_classification
# tissue_mapping_all
for (tissue_name in overlapping_tissues){

#   tissue_name <-"adrenal gland"
  cat(paste("Analyzing tissue:", tissue_name, "\n"))
  
  # Get tissue-specific means (using normalized data)
  prop_mean <- proprietary_tissue_means[,tissue_name]
#   pub1_mean <- published1_tissue_means[,tissue_name]
  names(prop_mean) <-rownames(proprietary_tissue_means)
#   names(pub1_mean) <-rownames(published1_tissue_means)
  mapped_published1_tissue_name<-unique(tissue_mapping_all[tissue_mapping_all$tissue_type==tissue_name,"HPA_tissue"])%>%as.character()
  mapped_pubblished2_tissue_name<-unique(tissue_mapping_all[tissue_mapping_all$tissue_type==tissue_name,"gtex_tissue"])%>%as.character()
  mapped_pubblished2_tissue_name
  # Check if tissue exists in published1
  if(sum(!is.na(mapped_published1_tissue_name))==0) {
    cat(paste("WARNING: Tissue", tissue_name, "not found in Published1 dataset\n"))
    pub1_mean <- NULL
  } else {
    pub1_mean <- data_published1_normalized[, mapped_published1_tissue_name]
    names(pub1_mean) <-rownames(data_published1_normalized)
  }

  
  # Check if tissue exists in published2
  if(sum(!is.na( mapped_pubblished2_tissue_name))==0) {
    cat(paste("WARNING: Tissue", tissue_name, "not found in Published2 dataset\n"))
    pub2_mean <- NULL
  } else {
    pub2_mean <- data_published2_normalized[,  mapped_pubblished2_tissue_name]
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
    corr_prop_pub1<- cor(prop_mean[common_proteins_prop_pub1], 
                         pub1_mean[common_proteins_prop_pub1], use = "complete.obs", method = "pearson")
    #  tissue_mapping_all[tissue_mapping_all$tissue_name=="tissue_name","Prop_Pub1"] <-  corr_prop_pub1
  }
#   corr_prop_pub1 
  if(length(common_proteins_prop_pub2) >= 10) {
    corr_prop_pub2 <- cor(prop_mean[common_proteins_prop_pub2], 
                         pub2_mean[common_proteins_prop_pub2], use = "complete.obs", method = "pearson")
    # tissue_mapping_all[tissue_mapping_all$tissue_name=="tissue_name","Prop_Pub2"] <-corr_prop_pub2
  }
#    corr_prop_pub2 
  if(length(common_proteins_pub1_pub2) >= 10) {
     corr_pub1_pub2  <- cor(pub1_mean[common_proteins_pub1_pub2], 
                         pub2_mean[common_proteins_pub1_pub2], use = "complete.obs", method = "pearson")   
    #  tissue_mapping_all[tissue_mapping_all$tissue_name=="tissue_name","Pub1_Pub2"]<-corr_pub1_pub2  
  }
#   corr_pub1_pub2  #0.73
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
    library(ggplot2)
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
    
    ggsave(paste0("output/", tissue_name, "_dataset_similarity.png"), 
           combined_plot, width = 15, height = 5, dpi = 300)
  }
  
    tissue_mapping_all[tissue_mapping_all$tissue_type==tissue_name,"Prop_Pub1"] = corr_prop_pub1
    tissue_mapping_all[tissue_mapping_all$tissue_type==tissue_name,"Prop_Pub2"] =   corr_prop_pub2
    tissue_mapping_all[tissue_mapping_all$tissue_type==tissue_name,"Pub1_Pub2"] =  corr_pub1_pub2

}
write.csv(tissue_mapping_all, "output/dataset_correlations_imputated_pearson_20251204.csv")





###############################
# =============================================================================
# 5. DATASET SIMILARITY ANALYSIS BY protein
# =============================================================================

#20251027 changed to overlap gene in overlap tissues, prot vs prot correlation

# # Generate pairwise paired data.frame
# gene_pairs <- as.data.frame(t(combn(common_genes, 2)))
# colnames(gene_pairs) <- c("Gene1", "Gene2")
# 
# 
# gene_pairs$HPA_vs_pub1_cor<-NA
# gene_pairs$HPA_vs_pub2_cor<-NA
data_published1_normalized<-data.frame(data_published1_normalized,check.rows = F,check.names = F)
data_published2_normalized<-data.frame(data_published2_normalized,check.rows = F,check.names = F)

gene_mat<-data.frame(common_genes)
gene_mat$TPHP_vs_pub1_cor_pearson<-NA
gene_mat$TPHP_vs_pub1_P_pearson<-NA
gene_mat$TPHP_vs_pub2_cor_pearson<-NA
gene_mat$TPHP_vs_pub2_P_pearson<-NA
gene_mat$TPHP_vs_pub1_cor_spearman<-NA
gene_mat$TPHP_vs_pub1_P_spearman<-NA
gene_mat$TPHP_vs_pub2_cor_spearman<-NA
gene_mat$TPHP_vs_pub2_P_spearman<-NA

for (i in 1:nrow(gene_mat)){
  
  #   tissue_name <-"adrenal gland"
  cat(paste("Analyzing tissue:", gene_mat$common_genes[i], "\n"))
  
  # Get tissue-specific means (using normalized data)
  prop_mean <- data.frame(t(proprietary_tissue_means[gene_mat$common_genes[i],]))
  names(prop_mean)<-"prot"
  pub1_mean<-data.frame(t(data_published1_normalized[gene_mat$common_genes[i],]))
  pub2_mean <- data.frame(t(data_published2_normalized[gene_mat$common_genes[i],]))
  pub1_mean$tissue_name<-tissue_mapping_all$tissue_type[match(row.names(pub1_mean),tissue_mapping_all$HPA_tissue)]
  pub2_mean$tissue_name<-tissue_mapping_all$tissue_type[match(row.names(pub2_mean),tissue_mapping_all$gtex_tissue)]
  names(pub1_mean)[1]<-"pub1_value"
  names(pub2_mean)[1]<-"pub2_value"
  
  pub1_mean$prop_mean_value<-prop_mean$prot[match(pub1_mean$tissue_name,row.names(prop_mean))]
  # pub1_mean$prop_mean_value<-as.numeric(pub1_mean$prop_mean_value)
  pub2_mean$prop_mean_value<-prop_mean$prot[match(pub2_mean$tissue_name,row.names(prop_mean))]
  # pub2_mean$prop_mean_value<-as.numeric(pub2_mean$prop_mean_value)
  if(sd(pub1_mean$pub1_value)!=0&sd(pub1_mean$prop_mean_value)!=0){
    gene_mat$TPHP_vs_pub1_P_pearson[i]<-cor.test(pub1_mean$prop_mean_value, 
                                             pub1_mean$pub1_value, use = "complete.obs", method = "pearson")$p.value
    gene_mat$TPHP_vs_pub1_cor_pearson[i]<-cor.test(pub1_mean$prop_mean_value, 
                                                   pub1_mean$pub1_value, use = "complete.obs", method = "pearson")$estimate
    
    
   
    gene_mat$TPHP_vs_pub1_P_spearman[i]<-cor.test(pub1_mean$prop_mean_value, 
                                              pub1_mean$pub1_value, use = "complete.obs", method = "spearman")$p.value
    gene_mat$TPHP_vs_pub1_cor_spearman[i]<-cor.test(pub1_mean$prop_mean_value, 
                                                    pub1_mean$pub1_value, use = "complete.obs", method = "spearman")$estimate
    
  }
  
  
  if(sd(pub2_mean$pub2_value)!=0&sd(pub2_mean$prop_mean_value)!=0){
    
  gene_mat$TPHP_vs_pub2_P_pearson[i]<-cor.test(pub2_mean$prop_mean_value, 
                                           pub2_mean$pub2_value, use = "complete.obs", method = "pearson")$p.value
  gene_mat$TPHP_vs_pub2_cor_pearson[i]<-cor.test(pub2_mean$prop_mean_value, 
                                                 pub2_mean$pub2_value, use = "complete.obs", method = "pearson")$estimate
  
  
  
  gene_mat$TPHP_vs_pub2_P_spearman[i]<-cor.test(pub2_mean$prop_mean_value, 
                                           pub2_mean$pub2_value, use = "complete.obs", method = "spearman")$p.value
  gene_mat$TPHP_vs_pub2_cor_spearman[i]<-cor.test(pub2_mean$prop_mean_value, 
                                                  pub2_mean$pub2_value, use = "complete.obs", method = "spearman")$estimate
  }

  
}

gene_mat$TPHP_vs_pub1_Padj_pearson<-p.adjust(gene_mat$TPHP_vs_pub1_P_pearson, method = 'BH')
gene_mat$TPHP_vs_pub2_Padj_pearson<-p.adjust(gene_mat$TPHP_vs_pub2_P_pearson, method = 'BH')
gene_mat$TPHP_vs_pub1_Padj_spearman<-p.adjust(gene_mat$TPHP_vs_pub1_P_spearman, method = 'BH')
gene_mat$TPHP_vs_pub2_Padj_spearman<-p.adjust(gene_mat$TPHP_vs_pub2_P_spearman, method = 'BH')


write.csv(gene_mat, "output/20251205_protenins_correlations_TPHP_vs_HPA_and_gtex.csv",row.names = F)

# 

plot(density(gene_mat$TPHP_vs_pub1_cor_pearson, na.rm = T), col = "red4", lwd = 2, 
     main = "Density Plot of Correlation Coefficients", xlab = "Correlation Coefficient", ylab = "Density",ylim=c(0,1.5))
lines(density(gene_mat$TPHP_vs_pub2_cor_pearson, na.rm = T), col = "blue4", lwd = 2)
lines(density(gene_mat$TPHP_vs_pub1_cor_spearman, na.rm = T), col = "green4", lwd = 2)
lines(density(gene_mat$TPHP_vs_pub2_cor_spearman, na.rm = T), col = "purple4", lwd = 2)
legend("topright", 
       legend = c("TPHP vs pub1 (Pearson)", "TPHP vs pub2 (Pearson)", 
                  "TPHP vs pub1 (Spearman)", "TPHP vs pub2 (Spearman)"),
       col = c("red4", "blue4", "green4", "purple4"), lwd = 2, cex = 0.8)


###################
# Calculate the cor value range and median for four groups under padj<0.05

cor_summary<-data.frame(matrix(NA,nrow = 4,ncol = 4))
names(cor_summary)<-c("cor_type","cor_median","cor_quantile_range_min","cor_quantile_range_max")

cor_summary$cor_type<-c("TPHP_vs_HPA_pearson","TPHP_vs_HPA_spearman","TPHP_vs_gtex_pearson","TPHP_vs_gtex_spearman")

cor_summary$cor_median[1]<-median(gene_mat$TPHP_vs_pub1_cor_pearson[gene_mat$TPHP_vs_pub1_Padj_pearson<0.05&!is.na(gene_mat$TPHP_vs_pub1_Padj_pearson)],na.rm=T)
cor_summary$cor_quantile_range_min[1]<-quantile(gene_mat$TPHP_vs_pub1_cor_pearson[gene_mat$TPHP_vs_pub1_Padj_pearson<0.05&!is.na(gene_mat$TPHP_vs_pub1_Padj_pearson)],na.rm=T)[2]
cor_summary$cor_quantile_range_max[1]<-quantile(gene_mat$TPHP_vs_pub1_cor_pearson[gene_mat$TPHP_vs_pub1_Padj_pearson<0.05&!is.na(gene_mat$TPHP_vs_pub1_Padj_pearson)],na.rm=T)[4]

cor_summary$cor_median[2]<-median(gene_mat$TPHP_vs_pub1_cor_spearman[gene_mat$TPHP_vs_pub1_Padj_spearman<0.05&!is.na(gene_mat$TPHP_vs_pub1_Padj_spearman)],na.rm=T)
cor_summary$cor_quantile_range_min[2]<-quantile(gene_mat$TPHP_vs_pub1_cor_spearman[gene_mat$TPHP_vs_pub1_Padj_spearman<0.05&!is.na(gene_mat$TPHP_vs_pub1_Padj_spearman)],na.rm=T)[2]
cor_summary$cor_quantile_range_max[2]<-quantile(gene_mat$TPHP_vs_pub1_cor_spearman[gene_mat$TPHP_vs_pub1_Padj_spearman<0.05&!is.na(gene_mat$TPHP_vs_pub1_Padj_spearman)],na.rm=T)[4]

cor_summary$cor_median[3]<-median(gene_mat$TPHP_vs_pub2_cor_pearson[gene_mat$TPHP_vs_pub2_Padj_pearson<0.05&!is.na(gene_mat$TPHP_vs_pub2_Padj_pearson)],na.rm=T)
cor_summary$cor_quantile_range_min[3]<-quantile(gene_mat$TPHP_vs_pub2_cor_pearson[gene_mat$TPHP_vs_pub2_Padj_pearson<0.05&!is.na(gene_mat$TPHP_vs_pub2_Padj_pearson)],na.rm=T)[2]
cor_summary$cor_quantile_range_max[3]<-quantile(gene_mat$TPHP_vs_pub2_cor_pearson[gene_mat$TPHP_vs_pub2_Padj_pearson<0.05&!is.na(gene_mat$TPHP_vs_pub2_Padj_pearson)],na.rm=T)[4]

cor_summary$cor_median[4]<-median(gene_mat$TPHP_vs_pub2_cor_spearman[gene_mat$TPHP_vs_pub2_Padj_spearman<0.05&!is.na(gene_mat$TPHP_vs_pub2_Padj_spearman)],na.rm=T)
cor_summary$cor_quantile_range_min[4]<-quantile(gene_mat$TPHP_vs_pub2_cor_spearman[gene_mat$TPHP_vs_pub2_Padj_spearman<0.05&!is.na(gene_mat$TPHP_vs_pub2_Padj_spearman)],na.rm=T)[2]
cor_summary$cor_quantile_range_max[4]<-quantile(gene_mat$TPHP_vs_pub2_cor_spearman[gene_mat$TPHP_vs_pub2_Padj_spearman<0.05&!is.na(gene_mat$TPHP_vs_pub2_Padj_spearman)],na.rm=T)[4]

write.csv(cor_summary,"output/20251204_protenins_correlations_TPHP_vs_HPA_and_gtex_quantile_range_summary.csv",row.names = F)

summary(gene_mat$TPHP_vs_pub1_cor_pearson)
# bbb<-gene_mat[!is.na(gene_mat$TPHP_vs_pub1_Padj_pearson),]
# bbb<-gene_mat[!is.na(gene_mat$TPHP_vs_pub1_cor_pearson),]
quantile(gene_mat$TPHP_vs_pub1_cor_pearson,na.rm=T)
IQR(gene_mat$TPHP_vs_pub1_cor_pearson,na.rm=T)#0.5114363
mmm<-quantile(gene_mat$TPHP_vs_pub1_cor_pearson[gene_mat$TPHP_vs_pub1_Padj_pearson<0.05&!is.na(gene_mat$TPHP_vs_pub1_Padj_pearson)],na.rm=T)
mmm[2]#0.6402637
##############
#go bp 
# Calculate cor p-adjust, p-adjust<0.05 and cor>0 is positive correlation, otherwise p-adjust<0.05 and cor<0 is negative correlation, run enrichment for these two parts separately

#go
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ReactomePA)
library(clusterProfiler)
go_mat<-data.frame(matrix(NA,nrow = 8,ncol = 3))
names(go_mat)<-c("group","sign","gene_list")
go_mat$group<-rep(c("TPHP_vs_HPA_pearson","TPHP_vs_gtex_pearson","TPHP_vs_HPA_spearman","TPHP_vs_gtex_spearman"),each =2)
go_mat$sign<-rep(c("up","down"),4)

#20251204 changed to padj<0.1 to test
go_mat$gene_list[1]<-paste(gene_mat$common_genes[gene_mat$TPHP_vs_pub1_cor_pearson>0&gene_mat$TPHP_vs_pub1_Padj_pearson<0.1&!is.na(gene_mat$TPHP_vs_pub1_cor_pearson)],collapse = ";")
go_mat$gene_list[2]<-paste(gene_mat$common_genes[gene_mat$TPHP_vs_pub1_cor_pearson< 0&gene_mat$TPHP_vs_pub1_Padj_pearson<0.1&!is.na(gene_mat$TPHP_vs_pub1_cor_pearson)],collapse = ";")
go_mat$gene_list[3]<-paste(gene_mat$common_genes[gene_mat$TPHP_vs_pub2_cor_pearson>0&gene_mat$TPHP_vs_pub2_Padj_pearson<0.1&!is.na(gene_mat$TPHP_vs_pub2_cor_pearson)],collapse = ";")
go_mat$gene_list[4]<-paste(gene_mat$common_genes[gene_mat$TPHP_vs_pub2_cor_pearson< 0&gene_mat$TPHP_vs_pub2_Padj_pearson<0.1&!is.na(gene_mat$TPHP_vs_pub2_cor_pearson)],collapse = ";")
go_mat$gene_list[5]<-paste(gene_mat$common_genes[gene_mat$TPHP_vs_pub1_cor_spearman>0&gene_mat$TPHP_vs_pub1_Padj_spearman<0.1&!is.na(gene_mat$TPHP_vs_pub1_cor_spearman)],collapse = ";")
go_mat$gene_list[6]<-paste(gene_mat$common_genes[gene_mat$TPHP_vs_pub1_cor_spearman< 0&gene_mat$TPHP_vs_pub1_Padj_spearman<0.1&!is.na(gene_mat$TPHP_vs_pub1_cor_spearman)],collapse = ";")
go_mat$gene_list[7]<-paste(gene_mat$common_genes[gene_mat$TPHP_vs_pub2_cor_spearman>0&gene_mat$TPHP_vs_pub2_Padj_spearman<0.1&!is.na(gene_mat$TPHP_vs_pub2_cor_spearman)],collapse = ";")
go_mat$gene_list[8]<-paste(gene_mat$common_genes[gene_mat$TPHP_vs_pub2_cor_spearman< 0&gene_mat$TPHP_vs_pub2_Padj_spearman<0.1&!is.na(gene_mat$TPHP_vs_pub2_cor_spearman)],collapse = ";")

for (i in 1:nrow(go_mat)) {
  GO_enrich = enrichGO(gene = unlist(strsplit(go_mat$gene_list[i],";")), # indicates foreground genes, i.e., the gene list to be enriched; [,1] indicates processing the first column of the entrezid_all dataset
                       keyType ="ENSEMBL",#ENSEMBL #UNIPROT
                       OrgDb = "org.Hs.eg.db", ont = "BP", # can input CC/MF/BP/ALL
                       #universe =bg_prot,  # indicates background genes, for de novo species use all assembled unigenes as background; for reference background genes are not needed
                       pvalueCutoff = 1, qvalueCutoff = 1, readable = T)
  
  
  GO_enrich = clusterProfiler::simplify(GO_enrich, cutoff = 0.7)#duplicate
  GO_enrich1  = data.frame(GO_enrich) # convert GO_enrich to data.frame format
  
  GO_enrich1 = subset(GO_enrich1,  pvalue<0.05)
  
  
  
  write.csv(GO_enrich1,paste0("output/20251204_proteins_correlations_",go_mat$group[i],"_",go_mat$sign[i],"_go.csv"),row.names = F)
  
}



  
#############
#20251031 plot the distribution of four groups cor
#pub1:HPA 
#pub2:gtex
#cor>0padj<0.05 is positive, cor<0padj<0.05 is negative, padj>0.05 is notsignificant
library(ggplot2)
df<-gene_mat
df<-read.csv("output/20251204_protenins_correlations_TPHP_vs_HPA_and_gtex.csv",header = T)
tphp_hpa_cor_pearson<-df[,c(1,2,3,10)]
# tphp_hpa_cor_pearson$cor_type<-"pearson"

tphp_hpa_cor_spearman<-df[,c(1,6,7,12)]
# tphp_hpa_cor_spearman$cor_type<-"spearman"

tphp_gtex_cor_pearson<-df[,c(1,4,5,11)]
# tphp_gtex_cor_pearson$cor_type<-"pearson"

tphp_gtex_cor_spearman<-df[,c(1,8,9,13)]
# tphp_gtex_cor_spearman$cor_type<-"spearman"

names(tphp_hpa_cor_pearson)<-names(tphp_hpa_cor_spearman)<-names(tphp_gtex_cor_pearson)<-names(tphp_gtex_cor_spearman)<-c("common_genes","cor","p","padj")
### Only plot spearman, classify cor

tphp_hpa_cor_spearman<-tphp_hpa_cor_spearman[!is.na(tphp_hpa_cor_spearman$cor),]
tphp_hpa_cor_spearman$cor_type[tphp_hpa_cor_spearman$cor>0&tphp_hpa_cor_spearman$padj<0.1]<-"Positive correlation"
tphp_hpa_cor_spearman$cor_type[tphp_hpa_cor_spearman$cor<0&tphp_hpa_cor_spearman$padj<0.1]<-"Negative correlation"
tphp_hpa_cor_spearman$cor_type[tphp_hpa_cor_spearman$padj>=0.1]<-"Notsignificant"
unique(tphp_hpa_cor_spearman$cor_type)
sum(tphp_hpa_cor_spearman$cor>0&tphp_hpa_cor_spearman$padj<0.1)
sum(tphp_hpa_cor_spearman$cor<0&tphp_hpa_cor_spearman$padj<0.1)


min(tphp_hpa_cor_spearman$cor)
max(tphp_hpa_cor_spearman$cor)

a <- ggplot(tphp_hpa_cor_spearman, aes(x = cor, fill = cor_type)) +
  geom_histogram(position = "stack",  # change to stacking
                 alpha = 0.6, bins = 30, color = "white", linewidth = 0.1) +
  labs(title = "Correlation Distribution", x = "Correlation Value", y = "Count", fill = "Correlation Type") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(),
    text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "top"
  ) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.4))+
  coord_cartesian(xlim = c(-1, 1))

a



ggsave("output/20251204_THPA_HPA_spearman_cor_density.pdf",a,width = 8,height = 8)

####

tphp_gtex_cor_spearman<-tphp_gtex_cor_spearman[!is.na(tphp_gtex_cor_spearman$cor),]
tphp_gtex_cor_spearman$cor_type[tphp_gtex_cor_spearman$cor>0&tphp_gtex_cor_spearman$padj<0.1]<-"Positive correlation"
tphp_gtex_cor_spearman$cor_type[tphp_gtex_cor_spearman$cor<0&tphp_gtex_cor_spearman$padj<0.1]<-"Negative correlation"
tphp_gtex_cor_spearman$cor_type[tphp_gtex_cor_spearman$padj>=0.1]<-"Notsignificant"
sum(tphp_gtex_cor_spearman$cor>0&tphp_gtex_cor_spearman$padj<0.1)
sum(tphp_gtex_cor_spearman$cor<0&tphp_gtex_cor_spearman$padj<0.1)
unique(tphp_gtex_cor_spearman$cor_type)

a <- ggplot(tphp_gtex_cor_spearman, aes(x = cor, fill = cor_type)) +
  geom_histogram(position = "stack",  # change to stacking
                 alpha = 0.6, bins = 30, color = "white", linewidth = 0.1) +
  labs(title = "Correlation Distribution", x = "Correlation Value", y = "Count", fill = "Correlation Type") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(),
    text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "top"
  ) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.4))+coord_cartesian(xlim = c(-1, 1))

a

ggsave("output/20251204_THPA_gtex_spearman_cor_density.pdf",a,width = 8,height = 8)


########################
# change to pearson
tphp_hpa_cor_pearson<-tphp_hpa_cor_pearson[!is.na(tphp_hpa_cor_pearson$cor),]
tphp_hpa_cor_pearson$cor_type[tphp_hpa_cor_pearson$cor>0&tphp_hpa_cor_pearson$padj<0.1]<-"Positive correlation"
tphp_hpa_cor_pearson$cor_type[tphp_hpa_cor_pearson$cor<0&tphp_hpa_cor_pearson$padj<0.1]<-"Negative correlation"
tphp_hpa_cor_pearson$cor_type[tphp_hpa_cor_pearson$padj>=0.1]<-"Notsignificant"
unique(tphp_hpa_cor_pearson$cor_type)
sum(tphp_hpa_cor_pearson$cor>0&tphp_hpa_cor_pearson$padj<0.1)
sum(tphp_hpa_cor_pearson$cor<0&tphp_hpa_cor_pearson$padj<0.1)


min(tphp_hpa_cor_pearson$cor)
max(tphp_hpa_cor_pearson$cor)

a <- ggplot(tphp_hpa_cor_pearson, aes(x = cor, fill = cor_type)) +
  geom_histogram(position = "stack",  # change to stacking
                 alpha = 0.6, bins = 30, color = "white", linewidth = 0.1) +
  labs(title = "Correlation Distribution", x = "Correlation Value", y = "Count", fill = "Correlation Type") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(),
    text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "top"
  ) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.4))+
  coord_cartesian(xlim = c(-1, 1))

a



ggsave("output/20251212_THPA_HPA_pearson_cor_density.pdf",a,width = 8,height = 8)

####

tphp_gtex_cor_pearson<-tphp_gtex_cor_pearson[!is.na(tphp_gtex_cor_pearson$cor),]
tphp_gtex_cor_pearson$cor_type[tphp_gtex_cor_pearson$cor>0&tphp_gtex_cor_pearson$padj<0.1]<-"Positive correlation"
tphp_gtex_cor_pearson$cor_type[tphp_gtex_cor_pearson$cor<0&tphp_gtex_cor_pearson$padj<0.1]<-"Negative correlation"
tphp_gtex_cor_pearson$cor_type[tphp_gtex_cor_pearson$padj>=0.1]<-"Notsignificant"
sum(tphp_gtex_cor_pearson$cor>0&tphp_gtex_cor_pearson$padj<0.1)
sum(tphp_gtex_cor_pearson$cor<0&tphp_gtex_cor_pearson$padj<0.1)
unique(tphp_gtex_cor_pearson$cor_type)

a <- ggplot(tphp_gtex_cor_pearson, aes(x = cor, fill = cor_type)) +
  geom_histogram(position = "stack",  # change to stacking
                 alpha = 0.6, bins = 30, color = "white", linewidth = 0.1) +
  labs(title = "Correlation Distribution", x = "Correlation Value", y = "Count", fill = "Correlation Type") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(),
    text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "top"
  ) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.4))+coord_cartesian(xlim = c(-1, 1))

a

ggsave("output/20251212_THPA_gtex_pearson_cor_density.pdf",a,width = 8,height = 8)


#################
# In the summary table, compute the 25%, 50%, and 75% cor without filtering by p-value
gene_mat<-read.csv("output/20251205_protenins_correlations_TPHP_vs_HPA_and_gtex.csv",header=T)
cor_summary<-data.frame(matrix(NA,nrow = 4,ncol = 4))
names(cor_summary)<-c("cor_type","cor_median","cor_quantile_range_min","cor_quantile_range_max")

cor_summary$cor_type<-c("TPHP_vs_HPA_pearson","TPHP_vs_HPA_spearman","TPHP_vs_gtex_pearson","TPHP_vs_gtex_spearman")

cor_summary$cor_median[1]<-median(gene_mat$TPHP_vs_pub1_cor_pearson[!is.na(gene_mat$TPHP_vs_pub1_Padj_pearson)],na.rm=T)
cor_summary$cor_quantile_range_min[1]<-quantile(gene_mat$TPHP_vs_pub1_cor_pearson[!is.na(gene_mat$TPHP_vs_pub1_Padj_pearson)],na.rm=T)[2]
cor_summary$cor_quantile_range_max[1]<-quantile(gene_mat$TPHP_vs_pub1_cor_pearson[!is.na(gene_mat$TPHP_vs_pub1_Padj_pearson)],na.rm=T)[4]

cor_summary$cor_median[2]<-median(gene_mat$TPHP_vs_pub1_cor_spearman[!is.na(gene_mat$TPHP_vs_pub1_Padj_spearman)],na.rm=T)
cor_summary$cor_quantile_range_min[2]<-quantile(gene_mat$TPHP_vs_pub1_cor_spearman[!is.na(gene_mat$TPHP_vs_pub1_Padj_spearman)],na.rm=T)[2]
cor_summary$cor_quantile_range_max[2]<-quantile(gene_mat$TPHP_vs_pub1_cor_spearman[!is.na(gene_mat$TPHP_vs_pub1_Padj_spearman)],na.rm=T)[4]

cor_summary$cor_median[3]<-median(gene_mat$TPHP_vs_pub2_cor_pearson[!is.na(gene_mat$TPHP_vs_pub2_Padj_pearson)],na.rm=T)
cor_summary$cor_quantile_range_min[3]<-quantile(gene_mat$TPHP_vs_pub2_cor_pearson[!is.na(gene_mat$TPHP_vs_pub2_Padj_pearson)],na.rm=T)[2]
cor_summary$cor_quantile_range_max[3]<-quantile(gene_mat$TPHP_vs_pub2_cor_pearson[!is.na(gene_mat$TPHP_vs_pub2_Padj_pearson)],na.rm=T)[4]

cor_summary$cor_median[4]<-median(gene_mat$TPHP_vs_pub2_cor_spearman[!is.na(gene_mat$TPHP_vs_pub2_Padj_spearman)],na.rm=T)
cor_summary$cor_quantile_range_min[4]<-quantile(gene_mat$TPHP_vs_pub2_cor_spearman[!is.na(gene_mat$TPHP_vs_pub2_Padj_spearman)],na.rm=T)[2]
cor_summary$cor_quantile_range_max[4]<-quantile(gene_mat$TPHP_vs_pub2_cor_spearman[!is.na(gene_mat$TPHP_vs_pub2_Padj_spearman)],na.rm=T)[4]

write.csv(cor_summary,"output/20251204_protenins_correlations_TPHP_vs_HPA_and_gtex_quantile_range_summary_ignore_padj.csv",row.names = F)

