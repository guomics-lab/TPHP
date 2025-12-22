rm(list = ls())
#cell 2020
library(stringr)
library(tibble)
cell_bd<-readxl::read_xlsx("//172.16.13.136/tphp/code/5_2_compare_with_previous/input/mmc2.xlsx",sheet = 2,skip = 2)
gtex<-read.delim("//172.16.13.136/tphp/code/5_2_compare_with_previous/input/biobank.search.1755523013528.tsv",sep = "\t")
gtex_info<-data.frame(`Individual ID`= str_match(gtex$sampleId, "GTEX-(.*?)-")[,2],
age=gtex$ageBracket,
sex=gtex$sex,check.names = F)%>%unique()
intersect(cell_bd$`Individual ID`,gtex_info$`Individual ID`)
cell_info<-merge(cell_bd,gtex_info,by="Individual ID")
write.csv(cell_info,"//172.16.13.136/tphp/code/5_2_compare_with_previous/output/2020cell_individual_age_sex.csv")


#2025cell age and sex
#wide to long
liugh_bd<-readxl::read_xlsx("//172.16.13.136/tphp/code/5_2_compare_with_previous/input/2025c_liugh.xlsx",sheet = 1,skip = 2)
# View(liugh_bd)
colnames(liugh_bd)
liugh_info<-reshape2::melt(liugh_bd[-nrow(liugh_bd),],id.vars=c("NO.","Age (y)","Gender","Tissue number", "Group"))
# View(liugh_info)
liugh_info<-liugh_info[!is.na(liugh_info$value),]
View(liugh_info)

write.csv(liugh_info,"//172.16.13.136/tphp/code/5_2_compare_with_previous/output/2025cell_individual_age_sex.csv")

#tphp age and sex
tphp_bd<-readxl::read_xlsx("//172.16.13.136/tphp/code/0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V6.xlsx")
tphp_bd<-tphp_bd[,-grep("Age|Gender",colnames(tphp_bd))]
tphp_patient<-readxl::read_xlsx('//172.16.13.136/tphp/code/Supplementary_tables/20251209_nature_revision_version/Supplementary Table 1. Sample origin information.xlsx')
# # View(tphp_patient)
colnames(tphp_bd)
tphp_info<-merge(tphp_bd, tphp_patient,by.x="patient_ID",by.y="Patient_ID")
# View(tphp_info)
write.csv(tphp_info,"//172.16.13.136/tphp/code/5_2_compare_with_previous/output/tphp_individual_age_sex.csv")


#cell 2020 mapping tissue type
tissue_2020<-readxl::read_xlsx("//172.16.13.136/tphp/code/5_2_compare_with_previous/input/20230710_cell_tissue_list_edit.xlsx")
rownames(tissue_2020)<-tissue_2020$tissue_type
cell_info$tissue_type<-tissue_2020[cell_info$`Tissue Name`,]$organ
# View(cell_info)
liugh_info$tissue_type<-stringr::str_to_lower(liugh_info$variable)
# View(liugh_info)
intersect(unique(liugh_info$tissue_type),unique(tphp_info$tissue_name))
intersect(unique(cell_info$tissue_type),unique(tphp_info$tissue_name))
intersect(unique(liugh_info$tissue_type),unique(tphp_info$anatomical_classification))
intersect(unique(cell_info$tissue_type),unique(tphp_info$anatomical_location))

a<-intersect(unique(liugh_info$tissue_type),unique(tphp_info$tissue_name))
b<-intersect(unique(liugh_info$tissue_type),unique(tphp_info$anatomical_classification))
union(a,b) #10
tphp_map_2025<-tphp_info[(tphp_info$anatomical_classification %in% union(a,b) ) |(tphp_info$tissue_name %in% union(a,b) ), ]
c<-intersect(unique(cell_info$tissue_type),unique(tphp_info$tissue_name))
d<-intersect(unique(cell_info$tissue_type),unique(tphp_info$anatomical_classification))#22
d#22
tphp_map_2020<-tphp_info[(tphp_info$anatomical_classification %in% union(c,d) ) |(tphp_info$tissue_name %in% union(c,d) ), ]
tphp_map<-tphp_info[which(tphp_info$FileName %in% intersect(tphp_map_2020$FileName, tphp_map_2025$FileName)),]
# View(tphp_map)
dim(tphp_map)  #[1] 63 34
table(tphp_map$anatomical_classification)
# adrenal gland         heart         liver          lung      pancreas          skin        spleen 
# 8            18            12            14             2             5             4 


#map with anatomical classification
write.csv(tphp_map,"//172.16.13.136/tphp/code/5_2_compare_with_previous/output/20251210_tphp_mapping_previous_study_info.xlsx")

###map proteins
res<-readRDS("../4_tissue_specificity_analysis/output/20251201_tissue_Specificity_naive_05minImputated_3.rds")
tphp_pm <- apply(res$metadata[,-c(1:34)],c(1,2),as.numeric)%>%t() %>% as.data.frame()%>%rownames_to_column("UniProtKB/Swiss-Prot ID")
tphp_ensemble<-read.delim("input/tphp_normal_ensemble_10169.txt",sep = "\t",check.names = F,stringsAsFactors = F)  
tphp_ensemble2<-unique(tphp_ensemble[,c("UniProtKB/Swiss-Prot ID","Gene stable ID")])
tphp_merge<-merge(tphp_ensemble2,tphp_pm, by="UniProtKB/Swiss-Prot ID",all.x=T)  #10971
# tphp_merge<-tphp_merge[!is.na(tphp_merge$`UniProtKB/Swiss-Prot ID`) & !is.na(tphp_merge$`Gene stable ID`),]


cell2020_pm<-readxl::read_xlsx("input/CELL_2020_TS2.xlsx",sheet=6,skip = 3) #12627
# sum(duplicated(cell2020_pm$gene.id))
# [1] 0


cell2025pm<-readxl::read_xlsx("input/2025c_liugh.xlsx",sheet = 4,skip = 2)
liugh_ensemble<-read.delim("input/liugh_ensemble.txt",sep = "\t",check.names = F,stringsAsFactors = F)
# sum(liugh_ensemble$`HGNC symbol`=="")
# [1] 199
# > sum(liugh_ensemble$`Gene name`=="")
# [1] 0
liugh_ensemble2<-unique(liugh_ensemble[,c("Gene name","Gene stable ID")])
ligh_merge<-merge(liugh_ensemble2,cell2025pm, by="Gene name",by.y="Protein",all.x=T) #13320


write.csv(tphp_merge,"output/tphp_protein_normalized_pm_10169_ensemblID.csv",row.names = F)
write.csv(cell2020_pm,"output/cell2020_protein_normalized_pm_ensemblID.csv",row.names = F)
write.csv(ligh_merge,"output/liugh2025_protein_normalized_pm_ensemblID.csv",row.names = F)
