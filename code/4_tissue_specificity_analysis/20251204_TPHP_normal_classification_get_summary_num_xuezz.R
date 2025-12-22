library(tidyr)
library(dplyr)
library(stringr)
library(magrittr)
library(tibble)
pm_data<-readRDS("./output/20251201_tissue_Specificity_naive_05minImputated_3.rds")
######
pm2 <- pm_data$med.matrix
pm2 <- data.frame(pm2, row.names = pm2$class_abbr)[,-which(names(pm2) == "class_abbr")]

pm_info<-pm_data$metadata
# pm2 <- pm_data$imp
# pm2[is.na(pm2)]<-min(pm2,na.rm=T)*0.5
pm2$anatomical_classification<-row.names(pm2)
pm2<-pm2[,c(ncol(pm2),1:(ncol(pm2)-1))]



sect_mat<-rio::import('./input/20230710_tissue_comparison_edited.xlsx')
pm_info$Inner_class<-sect_mat$Inner_class[match(pm_info$tissue_name,sect_mat$Detailed_tissue_type)]
pm_info$Outer_class<-sect_mat$Outer_class[match(pm_info$tissue_name,sect_mat$Detailed_tissue_type)]
pm_info$Rough_tissue_type<-sect_mat$Rough_tissue_type[match(pm_info$tissue_name,sect_mat$Detailed_tissue_type)]
pm_info$Detailed_tissue_type<-sect_mat$Detailed_tissue_type[match(pm_info$tissue_name,sect_mat$Detailed_tissue_type)]
pm_info<-pm_info[,-c(35:10203)]


pm_info<-pm_info[!is.na(pm_info$Inner_class),]
pm_info$class_finall<-pm_info$Inner_class
pm_info$class_finall[pm_info$anatomical_classification==pm_info$Detailed_tissue_type]<-pm_info$Outer_class[pm_info$anatomical_classification==pm_info$Detailed_tissue_type]
# pm_info_use<-pm_info[pm_info$class_finall=="Only in PUH" ,]



tp_info<-pm_info[,c(8,39)]
tp_info<-tp_info%>%distinct()


tmp<-pm_data$ts.Enrich_rm_contam
names(tmp)<-c("tissue_type","UniprotID","Classification","symbol")
df_compare <- rio::import('./input/CELL_2020_TS2_20251031.xlsx')
# df_compare %>% select(ensembl_id:`GTEX_enriched tissue`)

# uniprot to symbol
######

tphp_ense_map<-rio::import('./input/tphp_normal_ensemble_10169.txt')
tphp_ense_map<-tphp_ense_map %>% distinct()
#######
#The table needs to be tidy

tmp$ensembl_id<-tphp_ense_map$`Gene stable ID`[match(tmp$UniprotID,tphp_ense_map$`UniProtKB/Swiss-Prot ID`)]
nnn<-tmp[is.na(tmp$ensembl_id),]#307
#####

df_compare_new<-merge(tmp,df_compare,by="ensembl_id",all.x =T)



setdiff(tp_info$class_abbr,df_compare_new$tissue_type)

HPA_RNA_spec<-rio::import("./input/proteinatlas_11cd7f94.tsv")
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
tissue_compair$Jiang_tissue_type<-gsub("\\- ","",tissue_compair$Jiang_tissue_type)
tissue_compair$Jiang_tissue_type<-tolower(tissue_compair$Jiang_tissue_type)
df_compare_new$D_Wang_map<-NA
df_compare_new$Jiang_map<-NA
df_compare_new$HPA_map<-NA
df_compare_new$HPA_RNA_map<-NA
df_compare_new$class_abbr<-df_compare_new$tissue_type
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
    ind_j<-tolower(unlist(strsplit(df_compare_new$`Jiang enriched tissue`[i],"; ")))
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

df_lbl <- data.frame(Wang = df_compare_new$`D.Wang el. al_enrichment_category`,
                     Jiang = df_compare_new$Jiang_category,
                     HPA_RNA = df_compare_new$HPA_RNA_tissue_specificity) %>%
  add_column(`Target uniprot` = df_compare_new$UniprotID,
             tissue_type=df_compare_new$tissue_type,
             x = factor(df_compare_new$UniprotID, levels = rev(unique(df_compare_new$UniprotID)), ordered = T),
             y = 1,
             .before = 1)


df_lbl$Wang[df_compare_new$D_Wang_map=="FALSE"&!is.na(df_compare_new$D_Wang_map)]<-NA
df_lbl$Jiang[df_compare_new$Jiang_map=="FALSE"&!is.na(df_compare_new$Jiang_map)]<-NA
df_lbl$HPA_RNA[df_compare_new$HPA_RNA_map=="FALSE"&!is.na(df_compare_new$HPA_RNA_map)]<-NA


df_lbl$Wang[!df_lbl$Wang %in% c("Group enriched", "Tissue enriched")] <- NA

df_lbl$Jiang[!df_lbl$Jiang %in% c("prt_specific", "prt_enriched_not_spec")] <- NA
df_lbl$HPA_RNA[!df_lbl$HPA_RNA %in% c("Group enriched", "Tissue enriched")] <- NA


tab_tis<-data.frame(table(tmp$tissue_type,tmp$Classification))
tab_tis_mat<-tab_tis[tab_tis$Var2=="Tissue enriched",]
tab_gro_mat<-tab_tis[tab_tis$Var2=="Group enriched",]
tp_info$tissue_enriched_prot<-tab_tis_mat$Freq[match(tp_info$class_abbr,tab_tis_mat$Var1)]
tp_info$group_enriched_prot<-tab_gro_mat$Freq[match(tp_info$class_abbr,tab_gro_mat$Var1)]
df_lbl$tiss_group<-paste(df_lbl$tissue_type,df_lbl$`Target uniprot`,sep="_")
tmp$tiss_group<-paste(tmp$tissue_type,tmp$UniprotID,sep = "_")
df_lbl$enrich_type<-tmp$Classification[match(df_lbl$tiss_group,tmp$tiss_group)]

for (i in 1:nrow(tp_info)) {
  tp_info$te_mapped_with_Jiang[i]<-length(unique(df_lbl$`Target uniprot`[df_lbl$enrich_type=="Tissue enriched"&df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$Jiang)]))
  tp_info$ge_mapped_with_Jiang[i]<-length(unique(df_lbl$`Target uniprot`[df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$Jiang)]))
  tp_info$te_mapped_with_Wang[i]<-length(unique(df_lbl$`Target uniprot`[df_lbl$enrich_type=="Tissue enriched"&df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$Wang)]))
  tp_info$ge_mapped_with_Wang[i]<-length(unique(df_lbl$`Target uniprot`[df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$Wang)]))
  tp_info$te_mapped_with_Jiang_Wang[i]<-length(unique(df_lbl$`Target uniprot`[df_lbl$enrich_type=="Tissue enriched"&df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$Wang)|
                                                                              df_lbl$enrich_type=="Tissue enriched"&df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$Jiang)]))
  tp_info$ge_mapped_with_Jiang_Wang[i]<-length(unique(df_lbl$`Target uniprot`[df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$Wang)|
                                                                              df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$Jiang)]))
  
  tp_info$te_mapped_with_HPA_RNA[i]<-length(unique(df_lbl$`Target uniprot`[df_lbl$enrich_type=="Tissue enriched"&df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$HPA_RNA)]))
  tp_info$ge_mapped_with_HPA_RNA[i]<-length(unique(df_lbl$`Target uniprot`[df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$HPA_RNA)]))
  
  tp_info$te_mapped_with_all3[i]<-length(unique(df_lbl$`Target uniprot`[df_lbl$enrich_type=="Tissue enriched"&df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$Wang)|
                                                                          df_lbl$enrich_type=="Tissue enriched"&df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$Jiang)|
                                                                          df_lbl$enrich_type=="Tissue enriched"&df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$HPA_RNA)]))
  tp_info$ge_mapped_with_all3[i]<-length(unique(df_lbl$`Target uniprot`[df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$Wang)|
                                                                          df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$Jiang)|
                                                                          df_lbl$tissue_type==tp_info$class_abbr[i]&!is.na(df_lbl$HPA_RNA)]))
  
}
####20251127 Add, for each tissue, the count that is tissue-enriched only at the HPA RNA level
ddd<-df_lbl[grepl("Tissue enriched",df_lbl$HPA_RNA),]

ddd<-ddd[!grepl("Tissue enriched",ddd$Wang),]
ddd<-ddd[!grepl("prt_specific",ddd$Jiang),]
ddd<-ddd[!grepl("Tissue enriched",ddd$enrich_type),]
mmm<-data.frame(table(ddd$tissue_type))
tp_info$te_only_in_PHA_RNA<-mmm$Freq[match(tp_info$class_abbr,mmm$Var1)]

tp_info$te_only_in_PHA_RNA[is.na(tp_info$te_only_in_PHA_RNA)]<-0

# write.csv(tp_info,"./output/20251203_TPHP_tissue_specific_type_mapping_with_Jiang_Wang_HPA_RNA_V2_add_only_HPA_RNA_tissue_enrich.csv",row.names = F)
# tp_info<-read.csv("./output/20251203_TPHP_tissue_specific_type_mapping_with_Jiang_Wang_HPA_RNA_V2_add_only_HPA_RNA_tissue_enrich.csv",header = T)
#Change tissue to full name
tissue_list<-pm_info[,c(7,8)]
tissue_list<-tissue_list%>%distinct()
tissue_list$anatomical_classification<-sub("^([a-z])", "\\U\\1", tissue_list$anatomical_classification, perl = TRUE, ignore.case = TRUE)
tp_info$class_abbr<-tissue_list$anatomical_classification[match(tp_info$class_abbr,tissue_list$class_abbr)]
write.csv(tp_info,"./output/20251203_TPHP_tissue_specific_type_mapping_with_Jiang_Wang_HPA_RNA_V2_add_only_HPA_RNA_tissue_enrich_change_class_abbr.csv",row.names = F)
# 


















       