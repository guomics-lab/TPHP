# set order
# pm_data<-readRDS("mat_med_na_list.rds")
library(dplyr)
library(magrittr)
#20251201 Only clamp to 3×; do not clamp the highest-expressing tissue’s tissue-enriched proteins

pm_data<-readRDS("./output/20251201_tissue_Specificity_naive_05minImputated_3.rds")
pm_info<-pm_data$metadata

######
pm2 <- pm_data$med.matrix
pm2 <- data.frame(pm2, row.names = pm2$class_abbr)[,-which(names(pm2) == "class_abbr")]
pm2 %<>% log2()


pm_heat_cor <- cor(t(pm2), method = 'spearman')
pm_heat_dist <- as.dist(1 - pm_heat_cor) # 1-Spearman's rho
pm_heat_tree <- hclust(pm_heat_dist, method="complete")
organs_ordered <- pm_heat_tree$labels[pm_heat_tree$order]
organs_ordered %<>% rev() # for mapping to barplot


######

pm2 <- pm_data$med.matrix
pm2 <- data.frame(pm2, row.names = pm2$class_abbr)[,-which(names(pm2) == "class_abbr")]

pm2$anatomical_classification<-row.names(pm2)
pm2<-pm2[,c(ncol(pm2),1:(ncol(pm2)-1))]

df1<-pm_data$ts.df



library(tidyr)
library(dplyr)
library(stringr)
library(magrittr)
library(tibble)

# start
df_drug <- rio::import('./input/drugbank_results_all_20230207.csv')
df_drug %<>% separate_rows(`Target uniprot`) %>% filter(`Target uniprot` != '')


protinfo<-pm_data$ts.tisEnrich_rm_contam
names(protinfo)<-c("tissue_type","UniprotID","Classification","symbol" )
# write.csv(protinfo,"20251029_TPHP_tissue_group_enriched_protein.csv",row.names = F)
# tissue enriched druggable proteins
prot_druggable <- df_drug$`Target uniprot` %>% str_split('; ') %>% unlist() %>% unique() %>% sort() %>% .[. != '']
# length(unique(prot_druggable))

protinfo_tisen <- protinfo %>%
  filter(Classification == 'Tissue enriched') %>% # tissue enriched
  filter(UniprotID %in% all_of(prot_druggable)) %>% # druggable 
  arrange(tissue_type)



##############
#First remove proteins with drugtarget numb<=3


df_tisen_organ <- pm2 %>%
  select(anatomical_classification, all_of(protinfo_tisen$UniprotID)) %>%
  arrange(anatomical_classification)

#get matrix
df_tisen_organ1<-pm2 %>%
  select(anatomical_classification, all_of(protinfo$UniprotID)) %>%
  arrange(anatomical_classification)

df_tisen_organ1$anatomical_classification<-pm_info$anatomical_classification[match(df_tisen_organ1$anatomical_classification,pm_info$class_abbr)]
row.names(df_tisen_organ1)<-df_tisen_organ1$anatomical_classification


pm_1 <- df_tisen_organ1%>% select(-1)
pm_1 %<>% log2()
pm_1 <- apply(pm_1, c(1, 2), as.numeric) %>% t() %>% data.frame()
names(pm_1)<-df_tisen_organ1$anatomical_classification
names(pm_1) <- sub("^([a-z])", "\\U\\1", names(pm_1), perl = TRUE, ignore.case = TRUE)

pm.scale_1 <- apply(pm_1, 1, scale) %>% t() %>% data.frame()
# pm.scale[na_pos] <- NA
names(pm.scale_1) <- names(pm_1)


row.names(pm.scale_1)<-paste(row.names(pm.scale_1),protinfo$symbol[match(row.names(pm.scale_1),protinfo$UniprotID)],sep = "_")


write.csv(pm.scale_1,"./output/20251209_TPHP_All_tissue_enriched_proteins_mat.csv",row.names = T)

# top5 proteins
enrich_top5_ls <- plyr::dlply(protinfo_tisen, 'tissue_type', function(dfsub){
  # dfsub <- protinfo_tisen %>% filter(tissue_type == 'ADG')
  pm <- df_tisen_organ %>% select(all_of(dfsub$UniprotID))
  med1_over_med2 <- apply(pm, 2, function(acol) { # calculate median_1st / median_2nd
    acol_sorted <- sort(acol, decreasing = T)
    return(acol_sorted[1] / acol_sorted[2])
  })
  enrich_top5 <- sort(med1_over_med2, decreasing = T) %>% head(5) %>% names() # top5
  return(enrich_top5)
})

names(df_tisen_organ)

# top X tissue enriched druggable proteins
prot_druggable <- unlist(enrich_top5_ls) # top5
length(prot_druggable) # 206 top5 enriched proteins
# prot_druggable[prot_druggable=="O15554"]
df_top5_drug <- df_drug %>% filter(str_detect(`Target uniprot`, str_c(prot_druggable, collapse = '|')))
dim(df_top5_drug) # 1860   26

df_top5_drug_long <- df_top5_drug %>%
  select(`DrugBank ID`:`Withdrawn`, `Target uniprot`) %>%
  separate_rows(`Target uniprot`, sep = '; ') %>%
  filter(`Target uniprot` != '')

df_top5_drug_target <- plyr::ddply(df_top5_drug_long, '`Target uniprot`', function(dfsub){
  data.frame(`Target uniprot` = dfsub$`Target uniprot`[1],
             Drugs = str_c(dfsub$Name, collapse = '; '),
             DrugNumber = length(unique(dfsub$Name)),
             check.names = F)
})
df_top5_prot <- df_top5_drug_target %>% filter(`Target uniprot` %in% prot_druggable)
# df_top5_prot %>% filter(`Target uniprot` == 'O15554') # 9 drugs for KCNN4

df_top5_prot_over3drug <- df_top5_prot %>% filter(DrugNumber > 3)
df_top5_prot_over3drug %<>% rev()

# df_top5_prot_over3drug %>% select(DrugNumber)
# top5en_prot_over3drug <- df_top5_prot_over3drug$`Target uniprot`

prot_druggable <- df_top5_prot_over3drug$`Target uniprot` # top5
protinfo_tisen <- protinfo %>%
  filter(Classification == 'Tissue enriched') %>% # tissue enriched
  filter(UniprotID %in% all_of(prot_druggable)) %>% # druggable 
  dplyr::mutate(tissue_type = factor(tissue_type, levels = organs_ordered, ordered = T)) %>%
  arrange(tissue_type)

df_tisen_organ <- pm2 %>%
  select(anatomical_classification, all_of(protinfo_tisen$UniprotID)) %>%
  dplyr::mutate(anatomical_classification = factor(anatomical_classification, levels = organs_ordered, ordered = T)) %>%
  arrange(anatomical_classification)
# identical(colnames(df_tisen_organ)[-1], protinfo_tisen$UniprotID) # TRUE


# keepna
pm <- df_tisen_organ %>% select(-1)
pm %<>% log2()
pm <- apply(pm, c(1, 2), as.numeric) %>% t() %>% data.frame()
names(pm) <- df_tisen_organ$anatomical_classification


pm.scale <- apply(pm, 1, scale) %>% t() %>% data.frame()
# pm.scale[na_pos] <- NA
names(pm.scale) <- names(pm)

prot <- colnames(df_tisen_organ)[-1]
row.tissue <- protinfo_tisen$tissue_type

# pm.scale <- apply(pm, 1, scale) %>% t() %>% data.frame()
# names(pm.scale) <- df_tisen$FileName
pm_heat <- pm.scale[prot,]

# concatenate uniprotid and gene names
uniprot_gene <- data.frame(prot = paste0(protinfo_tisen$UniprotID, '_', protinfo_tisen$symbol), row.names = protinfo_tisen$UniprotID)
srt <- row.names(pm_heat)
pm_heat <- merge(pm_heat, uniprot_gene, by = 'row.names')
row.names(pm_heat) <- pm_heat$prot
pm_heat <- pm_heat[, -grep('^Row\\.names$|^prot$', colnames(pm_heat))]
pm_heat <- pm_heat[srt, ]
#####
#Change tissue to full name
tissue_list<-pm_info[,c(7,8)]
tissue_list<-tissue_list%>%distinct()
tissue_list$anatomical_classification<-sub("^([a-z])", "\\U\\1", tissue_list$anatomical_classification, perl = TRUE, ignore.case = TRUE)

names(pm_heat)<-tissue_list$anatomical_classification[match(names(pm_heat),tissue_list$class_abbr)]

#############
# set labels
ann_row <- data.frame(row.tissue=row.tissue,
                      row.names =  row.names(pm_heat))
ann_col <- data.frame(tissue = factor(df_tisen_organ$anatomical_classification),
                      row.names = names(pm_heat))

ann_row$row.tissue<-tissue_list$anatomical_classification[match(ann_row$row.tissue,tissue_list$class_abbr)]
ann_row$row.tissue<-factor(ann_row$row.tissue,levels = unique(ann_row$row.tissue))
ann_col$tissue<-tissue_list$anatomical_classification[match(ann_col$tissue,tissue_list$class_abbr)]



####
# set colors
require(RColorBrewer)
library(ggplot2)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

set.seed(10)
tissue_color <- sample(col_vector, length(unique(ann_col$tissue)))
names(tissue_color) <- unique(ann_col$tissue)
# cancer_color <- c(brewer.pal(10,"Paired")[c(6,8,3)])[1:2]
# names(cancer_color) <- c('1_carcinoma', '2_adjacent')
ann_colors <- list(#cancer=cancer_color,
  tissue=tissue_color,row.tissue=tissue_color)


#####
# pm_heatdf<-cbind(pm_heat,ann_row)
# write.csv(pm_heatdf,"./output/20251209_TPHP_tissue_specific_heatmap_mat.csv",row.names = T)
# plot
dim(pm_heat)
x <- cumsum(table(ann_row$row.tissue))[-length(cumsum(table(ann_row$row.tissue)))]
a <- pheatmap::pheatmap(pm_heat, color = c(brewer.pal(11,"RdYlBu")[9:7],brewer.pal(11,"RdYlBu")[4:2]), fontsize_col = 8,
                        border_color = F,
                        annotation_col = ann_col, annotation_row = ann_row,
                        gaps_col = cumsum(table(ann_col$tissue))[-length(cumsum(table(ann_col$tissue)))],
                        gaps_row = x[!duplicated(x)],
                        fontsize_row=6,
                        annotation_colors = ann_colors,na_col = "grey60",#scale = "row",
                        angle_col = 45,
                        cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = T, 
                        filename = "./output/20251203_TPHP_normal_Ulhen_tissue_top5_enriched_over3_drugs_druggable_heatmap_scale_58protein_64organs_med.matrix.pdf",width=8,height=8)


# add heatmap row legend
row_annotation_matrix <- matrix(1:nrow(ann_row), ncol = 1, 
                                dimnames = list(rownames(ann_row), "Tissue"))

ann_colors1<-tissue_color[names(tissue_color)%in%ann_row$row.tissue]
ann_colors1_list <- list(row.tissue = ann_colors1)
# Now run pheatmap
row_heat <- pheatmap::pheatmap(
  row_annotation_matrix,
  color = rep("white", 1),  # Set a background color
  legend = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = ann_row,
  annotation_colors = ann_colors1_list,  # Use the correct list format
  annotation_names_row = TRUE,
  annotation_legend = TRUE,  # Show legend
  border_color = "white",
  cellwidth = 0.5,  # Very small column width
  cellheight = 6,   # Row height; adjust as needed
  fontsize_row = 6,
  gaps_row = x[!duplicated(x)],  # Use the same gaps
  silent = FALSE
)

# save
ggsave("./output/20251203_row_annotation_tissue.pdf", 
       plot = row_heat$gtable, 
       width = 2,  
       height = 8)




# drug numbers for each protein
df_drugnum <- df_top5_prot_over3drug %>% select(`Target uniprot`, DrugNumber) %>%
  mutate(`Target uniprot` = factor(`Target uniprot`, levels = rev(prot), ordered = T)) %>%
  arrange(`Target uniprot`)

gap_row <- x[!duplicated(x)]
tmp <- ann_row %>% set_rownames(str_remove(rownames(.), '_.+$')) %>% 
  rownames_to_column('Target uniprot')
for(i in seq_along(gap_row)){
  tmp %<>% add_row(`Target uniprot` = str_c('NA', i), .after = gap_row[i])
  gap_row <- gap_row + 1
}
df_drugnum <- tmp %>% slice(nrow(tmp):1) %>% full_join(df_drugnum)
df_drugnum$`Target uniprot` %<>% factor(., levels = .)

fivenum(log10(df_drugnum$DrugNumber)) # 0.602060 0.698970 0.903090 1.301030 2.480007
library(ggplot2)
p <- ggplot(df_drugnum) + 
  geom_tile(aes(1, `Target uniprot`, fill= log10(DrugNumber)), color = 'black') +
  geom_text(aes(1, `Target uniprot`, label = ifelse(DrugNumber > 0, DrugNumber, NA)), size = 5, color = 'black') +
  scale_fill_gradient(limits = c(floor(min(log10(df_drugnum$DrugNumber))), ceiling(max(log10(df_drugnum$DrugNumber)))),
                      # breaks,
                      low="white", high="#FAA449") +
  labs(
    x = '',
    y = 'Protein'
  )+
  scale_x_discrete(expand = c(0, 0))+
  # coord_cartesian(xlim = c(0, length(unique(df_n_long$Organ))) * 1.05,
  #                 ylim = c(0, length(unique(df_n_long$`Target uniprot`)) * 1.3))+
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
p

ggsave('./output/20251203_TPHP_normal_Ulhen_tissue_top5_enriched_over3_drugs_drggable_number_med.matrix.pdf', p, width = 3.5, height = 10)



# ------ number and ratio -----
prot_query<-pm_data$ts.Enrich

df_compare <- rio::import('./input/CELL_2020_TS2_20251031.xlsx')
# df_compare %>% select(ensembl_id:`GTEX_enriched tissue`)

# uniprot to symbol
######

# tphp_ense_map1<-rio::import('./input/tphp_ensemble.txt')
tphp_ense_map<-rio::import('./input/tphp_normal_ensemble_10169.txt')

tphp_ense_map<-tphp_ense_map %>% distinct()

#######
#The table needs to be tidy
prot_query$ensembl_id<-tphp_ense_map$`Gene stable ID`[match(prot_query$Protein,tphp_ense_map$`UniProtKB/Swiss-Prot ID`)]
#####
df_compare_new<-merge(prot_query,df_compare,by="ensembl_id",all.x =T)

######
#20251029 Add HPA RNA-specific data
HPA_RNA_spec<-rio::import("//172.16.13.136/tphp/code.20251201.archieved/4_tissue_specificity_analysis/input/proteinatlas_11cd7f94.tsv")
# HPA_RNA_spec<-HPA_RNA_spec[!is.na(HPA_RNA_spec$`RNA tissue specificity score`),]
df_compare_new$HPA_RNA_tissue_specificity<-HPA_RNA_spec$`RNA tissue specificity`[match(df_compare_new$ensembl_id,HPA_RNA_spec$Ensembl)]
df_compare_new$HPA_RNA_tissue_specific_nTPM<-HPA_RNA_spec$`RNA tissue specific nTPM`[match(df_compare_new$ensembl_id,HPA_RNA_spec$Ensembl)]
########
#20251031 Determine whether tissue specificity is consistent
tissue_compair<-rio::import("./input/20251030_tissue_compare.xlsx")
tphp_inf<-pm_data$metadata
tissue_compair$class_abbr<-tphp_inf$class_abbr[match(tissue_compair$TPHP_anatomical_classification,tphp_inf$anatomical_classification)]
tissue_compair$Wang_tissue<-tolower(tissue_compair$Wang_tissue)
# tissue_compair$gtex_tissue<-tolower(tissue_compair$gtex_tissue)
tissue_compair$HPA_tissue<-tolower(tissue_compair$HPA_tissue)
# tissue_compair$gtex_tissue<-gsub("\\- ","",tissue_compair$gtex_tissue)
tissue_compair$Jiang_tissue_type<-tolower(tissue_compair$Jiang_tissue_type)
df_compare_new$D_Wang_map<-NA
df_compare_new$Jiang_map<-NA
df_compare_new$HPA_map<-NA
df_compare_new$HPA_RNA_map<-NA
for (i in 1:nrow(df_compare_new)) {
  ind_tphp<-df_compare_new$class_abbr[i]
  if(!is.na(df_compare_new$`D.Wang el. al._Enriched_Tissue`[i])){
    ind_w<-tolower(unlist(strsplit(df_compare_new$`D.Wang el. al._Enriched_Tissue`[i],"; ")))
    ind_w_type<-unique(tissue_compair$class_abbr[tissue_compair$Wang_tissue%in%ind_w])
    if(sum(ind_tphp%in%ind_w_type)>0){
      df_compare_new$D_Wang_map[i]<-"TRUE"
      df_compare_new$HPA_map[i]<-"TRUE"
    }else{
      df_compare_new$D_Wang_map[i]<-"FALSE"
      df_compare_new$HPA_map[i]<-"FALSE"
    }
  }
  if(!is.na(df_compare_new$Jiang_category[i])){
    ind_j<-tolower(unlist(strsplit(df_compare_new$`Jiang_enriched tissue`[i],"; ")))
    ind_j_type<-unique(tissue_compair$class_abbr[tissue_compair$Jiang_tissue_type%in%ind_j])
    if(sum(ind_tphp%in%ind_j_type)>0){
      df_compare_new$Jiang_map[i]<-"TRUE"
    }else{
      df_compare_new$Jiang_map[i]<-"FALSE"
    }
  }
  
  if(!is.na(df_compare_new$HPA_RNA_tissue_specificity[i])){
    ind_pr<-unlist(strsplit(df_compare_new$HPA_RNA_tissue_specific_nTPM[i],";"))
    ind_pr<-tolower(unlist(sapply(strsplit(ind_pr,": "),function(e){e[1]})))
    ind_pr_type<-unique(tissue_compair$class_abbr[tissue_compair$HPA_tissue%in%ind_pr])
    if(sum(ind_tphp%in%ind_pr_type)>0){
      df_compare_new$HPA_RNA_map[i]<-"TRUE"
    }else{
      df_compare_new$HPA_RNA_map[i]<-"FALSE"
    }
  }
  
}

# df_compare_new$class_abbr<-tphp_inf$anatomical_classification[match(df_compare_new$class_abbr,tphp_inf$class_abbr)]
# df_compare_new$class_abbr<-sub("^([a-z])", "\\U\\1", df_compare_new$class_abbr, perl = TRUE, ignore.case = TRUE)
rio::export(df_compare_new, './output/20251203_DATABASE_enriched_proteins_compare_update.xlsx')
# length(unique(df_compare_new$Protein))

df_compare_new<-df_compare_new[df_compare_new$`Tissue specificity`=="Tissue enriched",]

df_compare_new %<>% filter(Protein %in% colnames(df_tisen_organ))
# %>% select(-contains('TPHP_')) %>% distinct() # Deduplicate; no duplicates
names(df_compare_new)[3]<-"UniprotID"
df_compare_new %<>%
  dplyr::mutate(UniprotID = factor(UniprotID, levels = colnames(df_tisen_organ)[-1], ordered = T)) %>%
  arrange(UniprotID) %>%
  dplyr::mutate(UniprotID = as.character(UniprotID))
setequal(colnames(df_tisen_organ)[-1], df_compare_new$UniprotID) # TRUE
identical(colnames(df_tisen_organ)[-1], df_compare_new$UniprotID) # TRUE
rio::export(df_compare_new, './output/20251210_DATABASE_enriched_proteins_compare_heatmap_prot.xlsx')
######
#Determine which label in HPA IHC corresponds to which one
df_compare_new$HPA_label<-NA
for(i in 1:nrow(df_compare_new)){
  if(!is.na(df_compare_new$HPA_map[i])&df_compare_new$HPA_map[i]=="TRUE"){
    ind_tiss<-tolower(unlist(strsplit(df_compare_new$`D.Wang el. al._Enriched_Tissue`[i],"; ")))
    ind_label<-unlist(strsplit(df_compare_new$HPA_iH_Level[i],"; "))
    ind_mat<-data.frame(cbind(ind_tiss,ind_label))
    ind_mat$class_abbr<-tissue_compair$class_abbr[match(ind_mat$ind_tiss,tissue_compair$Wang_tissue)]
    ind_mat<-ind_mat[!is.na(ind_mat$class_abbr),]
    df_compare_new$HPA_label[i]<-ind_mat$ind_label[ind_mat$class_abbr==df_compare_new$class_abbr[i]]
  }
}

df_compare_new$HPA_label<-gsub("\\-",NA,df_compare_new$HPA_label)
df_compare_new$HPA_label[is.na(df_compare_new$HPA_label)]<-"Not detected"
df_compare_new$HPA_label<-gsub("Not detected","Not mapped",df_compare_new$HPA_label)


df_lbl <- data.frame(Wang = df_compare_new$`D.Wang el. al_enrichment_category`,
                     Jiang = df_compare_new$Jiang_category,
                     HPA = df_compare_new$HPA_label,
                     HPA_RNA = df_compare_new$HPA_RNA_tissue_specificity) %>%
  add_column(`Target uniprot` = df_compare_new$UniprotID,
             x = factor(df_compare_new$UniprotID, levels = rev(unique(df_compare_new$UniprotID)), ordered = T),
             y = 1,
             .before = 1)


df_lbl$Wang[df_compare_new$D_Wang_map=="FALSE"&!is.na(df_compare_new$D_Wang_map)]<-NA
df_lbl$Jiang[df_compare_new$Jiang_map=="FALSE"&!is.na(df_compare_new$Jiang_map)]<-NA
df_lbl$HPA[df_compare_new$HPA_map=="FALSE"&!is.na(df_compare_new$HPA_map)]<-NA
df_lbl$HPA_RNA[df_compare_new$HPA_RNA_map=="FALSE"&!is.na(df_compare_new$HPA_RNA_map)]<-NA



unique(HPA_RNA_spec$`RNA tissue specificity`)

df_lbl <- df_drugnum %>% slice(nrow(df_drugnum):1) %>% 
  full_join(df_lbl)
df_lbl$x <- df_lbl$`Target uniprot`
df_lbl$x %<>% factor(., levels = rev(.))
df_lbl$Wang[is.na(df_lbl$Wang)]<-"Not mapped"
df_lbl$Jiang[is.na(df_lbl$Jiang)]<-"Not mapped"
df_lbl$HPA[is.na(df_lbl$HPA)]<-"Not mapped"
df_lbl$HPA_RNA[is.na(df_lbl$HPA_RNA)]<-"Not mapped"
df_lbl$Wang<-gsub("Mixed","Not mapped",df_lbl$Wang)
df_lbl$Wang<-gsub("Expressed in all","Not mapped",df_lbl$Wang)

unique(df_lbl$Wang)
unique(df_lbl$Jiang)
unique(df_lbl$HPA)
unique(df_lbl$HPA_RNA)

p1 <- ggplot(df_lbl) +
  aes(x = x, y = y, fill = Wang) +
  geom_point(shape = 21, size = 5) +
  scale_fill_manual(values = c('Tissue enriched' = '#984EA3', 'Group enriched' = '#E5C494', 'Not mapped' = 'white')) +
  coord_flip() +
  theme_void() +
  theme(axis.title.y = element_text(size = 12, hjust = 0.5, color = 'black', angle = 90),
        axis.text.y = element_text(size = 10, color = 'black', angle = 0),
        axis.ticks.y = element_blank()
  )

p2 <- ggplot(df_lbl) +
  aes(x = x, y = y, fill = Jiang) +
  geom_point(shape = 21, size = 5) +
  scale_fill_manual(values = c('prt_specific' = '#C17768', 'prt_enriched_not_spec' = '#D6C990', 'Not mapped' = 'white', 'prt_hk' = '#416F88')) +
  coord_flip() +
  theme_void() +
  theme(axis.title.y = element_text(size = 12, hjust = 0.5, color = 'black', angle = 90),
        axis.text.y = element_text(size = 10, color = 'black', angle = 0),
        axis.ticks.y = element_blank()
  )

p3 <- ggplot(df_lbl) +
  aes(x = x, y = y, fill = HPA) +
  geom_point(shape = 21, size = 5) +
  scale_fill_manual(values = c('High' = '#C17768', 'Medium' = '#D6C990', 'Not mapped' = 'white')) +
  coord_flip() +
  theme_void() +
  theme(axis.title.y = element_text(size = 12, hjust = 0.5, color = 'black', angle = 90),
        axis.text.y = element_text(size = 10, color = 'black', angle = 0),
        axis.ticks.y = element_blank()
  )

p4 <- ggplot(df_lbl) +
  aes(x = x, y = y, fill = HPA_RNA) +
  geom_point(shape = 21, size = 5) +
  scale_fill_manual(values = c('Tissue enriched' = '#C17768', 'Tissue enhanced' = '#D6C990', 'Not mapped' = 'white','Group enriched' = '#E5C494')) +
  coord_flip() +
  theme_void() +
  theme(axis.title.y = element_text(size = 12, hjust = 0.5, color = 'black', angle = 90),
        axis.text.y = element_text(size = 10, color = 'black', angle = 0),
        axis.ticks.y = element_blank()
  )


p <- ggpubr::ggarrange(p1, p2, p3,p4, nrow = 1)
ggsave('./output/20251203_TPHP_normal_Ulhen_tissue_top5_enriched_over3_drugs_druggable_compare.pdf', p, width = 12, height = 10)#fig 3B need

# dev.off()

# df_top5_prot_over3drug_stat
df_top5_prot_over3drug_stat <- plyr::ddply(df_top5_prot_over3drug, '`Target uniprot`', function(dfsub){
  # dfsub=df_top5_prot_over3drug[1,]
  drugs <- unique(str_split(dfsub$Drugs, '; ')[[1]])
  df_tmp <- df_top5_drug %>%
    filter(Name %in% drugs, `Target uniprot` == dfsub$`Target uniprot`[1]) %>% 
    distinct()
  
  status <- df_tmp %>%
    select(Approved:Investigational) %>%
    apply(1, function(arow){
      names(arow[arow]) %>% sort() %>% str_c(collapse = ';')
    })
  
  actions <- df_tmp %>% pull(`Target actions`)
  data.frame(DrugName = drugs,
             DrugStatus = status,
             DrugActions = actions)
})


dfbar_status <- plyr::ddply(df_top5_prot_over3drug_stat, '`Target uniprot`', function(dfsub){
  dfret <- dfsub %>%
    select(DrugStatus) %>%
    separate_rows(DrugStatus) %>%
    filter(DrugStatus != '') %>%
    count(DrugStatus)
  return(dfret)
})
dfbar_status$`Target uniprot` %<>% factor(., levels = rev(protinfo_tisen$UniprotID), ordered = T)
dfbar_status$DrugStatus %<>% factor(., levels = rev(c('Approved', 'Investigational', 'Experimental', 'Illicit')), ordered = T)
length(unique(dfbar_status$DrugStatus)) # 4 kinds of drug actions

dfbar_actions <- plyr::ddply(df_top5_prot_over3drug_stat, '`Target uniprot`', function(dfsub){
  dfret <- dfsub %>%
    select(DrugActions) %>%
    separate_rows(DrugActions, sep = '; ') %>%
    filter(DrugActions != '') %>%
    count(DrugActions)
  if(nrow(dfret) == 0){
    dfret <- data.frame(DrugActions = 'NA',
                        n = 1)
  }
  return(dfret)
})
dfbar_actions$`Target uniprot` %<>% factor(., levels = rev(protinfo_tisen$UniprotID), ordered = T)
dfbar_actions$DrugActions[dfbar_actions$DrugActions != 'NA'] %<>% str_to_sentence()
dfbar_actions$DrugActions %<>% factor(., levels = rev(c(setdiff(sort(dfbar_actions$DrugActions), 'NA'), 'NA')), ordered = T)
length(unique(dfbar_actions$DrugActions)) # 45 kinds of drug actions

length(unique(dfbar_status$`Target uniprot`))
length(unique(dfbar_actions$`Target uniprot`))


dfbar_status$`Target uniprot` %<>% as.character()
dfbar_status <- df_drugnum %>% slice(nrow(df_drugnum):1) %>% 
  full_join(dfbar_status)
dfbar_status$`Target uniprot` %<>% factor(., levels = rev(unique(.)))

p1 <- ggplot() +
  geom_bar(data = dfbar_status, aes(x = `Target uniprot`, y = n, fill = DrugStatus, group = DrugStatus), position="fill", stat="identity") +
  geom_text(data = df_top5_prot_over3drug_stat %>% count(`Target uniprot`, name = 'DrugNumber'),
            aes(x = `Target uniprot`, y = 1.01, label = DrugNumber),
            color = "black", size = 3, angle = 0, hjust = 0) + 
  labs(y = 'Ratio') +
  scale_fill_brewer(palette = "PuBu", direction = 1) +
  coord_flip() +
  theme_minimal() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())

#########
dfbar_status_1<-dfbar_status[!is.na(dfbar_status$row.tissue),]

library(dplyr)

dfbar_status_1 <- dfbar_status_1 %>%
  group_by(`Target uniprot`) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()
write.csv(dfbar_status_1,"./output/20251210_TPHP_normal_Ulhen_tissue_top5_enriched_over3_drugs_tables4_10.csv",row.names = F)








#######################
#Among actions that are >10% within 60% of targets — after filtering by percent > 10%, the count n > 0.6*length(prot_druggable)
length(prot_druggable) # 58
dfbar_actions_percent <- plyr::ddply(df_top5_prot_over3drug_stat, '`Target uniprot`', function(dfsub){
  dfret <- dfsub %>%
    select(DrugActions) %>%
    separate_rows(DrugActions, sep = '; ') %>%
    filter(DrugActions != '') %>%
    count(DrugActions)
  dfret$n <- dfret$n / sum(dfret$n)
  if(nrow(dfret) == 0){
    dfret <- data.frame(DrugActions = 'NA',
                        n = 1)
  }
  return(dfret)
})


dfbar_actions_percent %>%
  filter(n > 0.1) %>%
  count(DrugActions) %>%
  filter(n > 0.6*length(prot_druggable) )
#   DrugActions  n
# 1   inhibitor 45
# 2   substrate 45

drug_action_selected <- dfbar_actions_percent %>%
  filter(n > 0.1) %>%
  count(DrugActions) %>%
  arrange(desc(n)) %>%
  slice(1:12) %>%
  pull(DrugActions)

length(drug_action_selected) # 12
drug_action_selected %<>% str_to_sentence()


dfbar_actions_select1 <- dfbar_actions %>% filter(DrugActions %in% drug_action_selected)
dfbar_actions_select2 <- dfbar_actions %>% filter(!(DrugActions %in% drug_action_selected))
dfbar_actions_select2$DrugActions <- 'Others'
dfbar_actions_select <- rbind(dfbar_actions_select1, dfbar_actions_select2)
dfbar_actions_select$DrugActions %<>% factor(., levels = c(drug_action_selected, 'Others'), ordered = T)

df_color <- rio::import('//172.16.13.136/share/members/jiangwenhao/TPHP/input/PUH_tissue_colorset_20230210.xlsx')
tissue_color <- str_c('#', df_color$color[1:51])

action_colors <- str_c('#', df_color$color)[1:length(levels(dfbar_actions_select$DrugActions))]
names(action_colors) <- levels(dfbar_actions_select$DrugActions)
# action_colors['NA'] <- '#DDDDDD'


dfbar_actions_select$`Target uniprot` %<>% as.character()
dfbar_actions_select <- df_drugnum %>% slice(nrow(df_drugnum):1) %>% 
  full_join(dfbar_actions_select)
dfbar_actions_select$`Target uniprot` %<>% factor(., levels = rev(unique(.)))


######################
dfbar_actions_select_1<-dfbar_actions_select[!is.na(dfbar_actions_select$row.tissue),]
dfbar_actions_select_1 <- dfbar_actions_select_1 %>%
  group_by(`Target uniprot`) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()
write.csv(dfbar_actions_select_1,"./output/20251210_TPHP_normal_Ulhen_tissue_top5_enriched_over3_drugs_tables4_11.csv",row.names = F)

##################

p3 <- ggplot() +
  geom_bar(data = dfbar_actions_select, aes(x = `Target uniprot`, y = n, fill = DrugActions, group = DrugActions), position="fill", stat="identity") +
  labs(y = 'Ratio') +
  scale_fill_manual(values = action_colors) +
  coord_flip() +
  theme_minimal() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())
ggsave('./output/TPHP_normal_Ulhen_tissue_top5_enriched_over3_drugs_druggable_actions12_20251203.pdf', p3, width = 4, height = 10)

p <- ggpubr::ggarrange(p1, p3, nrow = 1)
ggsave('./output/TPHP_normal_Ulhen_tissue_top5_enriched_over3_drugs_druggable_status_actions12_20251203.pdf', p, width = 8, height = 10)#fig 3B

##################
#Obtain table data

df_top5_prot_over3drug_stat_1 <- plyr::ddply(df_top5_prot_over3drug, '`Target uniprot`', function(dfsub){
  # dfsub=df_top5_prot_over3drug[1,]
  drugs <- unique(str_split(dfsub$Drugs, '; ')[[1]])
  df_tmp <- df_top5_drug %>%
    filter(Name %in% drugs, `Target uniprot` == dfsub$`Target uniprot`[1]) %>% 
    distinct()
  
  status <- df_tmp %>%
    select(Approved:Investigational) %>%
    apply(1, function(arow){
      names(arow[arow]) %>% sort() %>% str_c(collapse = ';')
    })
  
  actions <- df_tmp %>% pull(`Target actions`)
  DrugBank_ID<-df_tmp %>%pull(`DrugBank ID`)
  data.frame(DrugName = drugs,
             DrugStatus = status,
             DrugActions = actions,
             DrugID=DrugBank_ID)
})

df_top5_prot_over3drug_stat_1<-df_top5_prot_over3drug_stat_1[,-2]


df_merged <- df_top5_prot_over3drug_stat_1 %>%
  group_by(`Target uniprot`) %>%
  summarise(across(everything(), ~ paste(na.omit(.), collapse = "-")))

df_merged$DrugNumber<-df_top5_prot_over3drug$DrugNumber[match(df_merged$`Target uniprot`,df_top5_prot_over3drug$`Target uniprot`)]
write.csv(df_merged,"./output/20251210_TPHP_normal_Ulhen_tissue_top5_enriched_over3_drugs_tables4_12.csv",row.names = F)

################
#all drug

# top5 proteins
enrich_top5_ls <- plyr::dlply(protinfo_tisen, 'tissue_type', function(dfsub){
  # dfsub <- protinfo_tisen %>% filter(tissue_type == 'ADG')
  pm <- df_tisen_organ %>% select(all_of(dfsub$UniprotID))
  med1_over_med2 <- apply(pm, 2, function(acol) { # calculate median_1st / median_2nd
    acol_sorted <- sort(acol, decreasing = T)
    return(acol_sorted[1] / acol_sorted[2])
  })
  enrich_top5 <- sort(med1_over_med2, decreasing = T) %>% head(5) %>% names() # top5
  return(enrich_top5)
})

names(df_tisen_organ)

# top X tissue enriched druggable proteins
prot_druggable <- unique(protinfo$UniprotID) # top5
length(prot_druggable) # 206 top5 enriched proteins
# prot_druggable[prot_druggable=="O15554"]
df_all_drug <- df_drug %>% filter(str_detect(`Target uniprot`, str_c(prot_druggable, collapse = '|')))
dim(df_all_drug) # 1860   26

df_all_drug_long <- df_all_drug %>%
  select(`DrugBank ID`:`Withdrawn`, `Target uniprot`) %>%
  separate_rows(`Target uniprot`, sep = '; ') %>%
  filter(`Target uniprot` != '')

df_all_drug_target <- plyr::ddply(df_all_drug_long, '`Target uniprot`', function(dfsub){
  data.frame(`Target uniprot` = dfsub$`Target uniprot`[1],
             Drugs = str_c(dfsub$Name, collapse = '; '),
             DrugNumber = length(unique(dfsub$Name)),
             check.names = F)
})

df_all_prot <- df_all_drug_target %>% filter(`Target uniprot` %in% prot_druggable)
# df_top5_prot %>% filter(`Target uniprot` == 'O15554') # 9 drugs for KCNN4


df_all_prot_over_stat <- plyr::ddply(df_all_prot, '`Target uniprot`', function(dfsub){
  # dfsub=df_top5_prot_over3drug[1,]
  drugs <- unique(str_split(dfsub$Drugs, '; ')[[1]])
  df_tmp <- df_all_drug %>%
    filter(Name %in% drugs, `Target uniprot` == dfsub$`Target uniprot`[1]) %>% 
    distinct()
  
  status <- df_tmp %>%
    select(Approved:Investigational) %>%
    apply(1, function(arow){
      names(arow[arow]) %>% sort() %>% str_c(collapse = ';')
    })
  
  actions <- df_tmp %>% pull(`Target actions`)
  DrugBank_ID<-df_tmp %>%pull(`DrugBank ID`)
  data.frame(DrugName = drugs,
             DrugStatus = status,
             DrugActions = actions,
             DrugID=DrugBank_ID)
})
write.csv(df_all_prot_over_stat,"./output/20251210_TPHP_normal_Ulhen_tissue_top5_enriched_over3_drugs_tables4_13.csv",row.names = F)


df_all_prot_over_stat<-df_all_prot_over_stat[,-2]


df_merged_all <- df_all_prot_over_stat %>%
  group_by(`Target uniprot`) %>%
  summarise(across(everything(), ~ paste(na.omit(.), collapse = "-")))

df_merged_all$DrugNumber<-df_all_prot$DrugNumber[match(df_merged_all$`Target uniprot`,df_all_prot$`Target uniprot`)]
write.csv(df_merged_all,"./output/20251210_TPHP_normal_Ulhen_tissue_top5_enriched_over3_drugs_tables4_14.csv",row.names = F)

all_dfbar_status <- plyr::ddply(df_all_prot_over_stat, '`Target uniprot`', function(dfsub){
  dfret <- dfsub %>%
    select(DrugStatus) %>%
    separate_rows(DrugStatus) %>%
    filter(DrugStatus != '') %>%
    count(DrugStatus)
  return(dfret)
})

all_dfbar_status$`Target uniprot` %<>% factor(., levels = rev(protinfo$UniprotID), ordered = T)
all_dfbar_status$DrugStatus %<>% factor(., levels = rev(c('Approved', 'Investigational', 'Experimental', 'Illicit')), ordered = T)
length(unique(all_dfbar_status$DrugStatus)) # 4 kinds of drug actions

write.csv(all_dfbar_status,"./output/20251210_TPHP_normal_Ulhen_tissue_top5_enriched_over3_drugs_tables4_15.csv",row.names = F)

all_dfbar_actions <- plyr::ddply(df_all_prot_over_stat, '`Target uniprot`', function(dfsub){
  dfret <- dfsub %>%
    select(DrugActions) %>%
    separate_rows(DrugActions, sep = '; ') %>%
    filter(DrugActions != '') %>%
    count(DrugActions)
  if(nrow(dfret) == 0){
    dfret <- data.frame(DrugActions = 'NA',
                        n = 1)
  }
  return(dfret)
})

all_dfbar_actions$`Target uniprot` %<>% factor(., levels = rev(protinfo$UniprotID), ordered = T)
all_dfbar_actions$DrugActions[all_dfbar_actions$DrugActions != 'NA'] %<>% str_to_sentence()
all_dfbar_actions$DrugActions %<>% factor(., levels = rev(c(setdiff(sort(all_dfbar_actions$DrugActions), 'NA'), 'NA')), ordered = T)



write.csv(all_dfbar_actions,"./output/20251210_TPHP_normal_Ulhen_tissue_top5_enriched_over3_drugs_tables4_16.csv",row.names = F)





