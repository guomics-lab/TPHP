setwd("//172.16.13.136/tphp/code/0_3_sample_info_mapping_data_matrix/")

###0.read data###################
pm0<- read.csv("input/mapped_pg_matrix_2856_13609.csv", row.names = 1)
dim(pm0) #2856 13609
head(pm0[,1:10] )
labels<- readxl::read_xlsx("../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V5_2.xlsx")
pool0<- read.csv( "input/pool_pg_matrix_149_13609.csv", row.names = 1 )
dim(pool0) #149 13609
head(pool0[,1:10])
str(labels)
#manual check timsTOF data QC-remove sample with abnoammly low protein identication number 
library(dplyr)
check_files<- readxl::read_xlsx("input/QC2856_identity_20251201_source_check.xlsx")
str(check_files)
lower_files<- check_files %>% filter(Is.Lower.Ingroup.major   == "TRUE" & Is.Lower.Ingroup.detailed == "TRUE" ) 
dim(lower_files) #48 files
# colnames(check_files)
# dim(lower_files)    #[1] 88
info<- labels %>%
    rename(tissue_type_major = anatomical_classification,
                    tissue_type_detailed = tissue_name) %>%
                    mutate(Low_protein_IDs  = ifelse( FileName %in% lower_files$FileName ,  TRUE, FALSE)) %>%
    mutate(imputation_group = paste(sample_type, tissue_type_major, sep="_"))                       

all_info<-labels %>%mutate(Low_protein_IDs  = ifelse( FileName %in% lower_files$FileName ,  TRUE, FALSE))
writexl::write_xlsx(all_info, "../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V6.xlsx" )
#filter out low group   
info1<-info[info$Low_protein_IDs!=TRUE,]   
dim(info1) #[1] 2957 34
View(info1)
write.csv(info1, "output/batch_info_tab_batchserver.csv", row.names = F  )

###1. prepare paired data matrix for QC######################
identical(colnames(pm0), colnames(pool0)) #T
s_pm0 <-as.data.frame(rbind(pm0, pool0))  
s_pm<-t(s_pm0[labels$FileName,]  ) 
colnames(s_pm) <- labels$FileName 
rownames(s_pm) <- colnames(s_pm0)
dim(s_pm) 
pm<-s_pm[,info1$FileName] 
dim(pm) #[1] 13609  2957
pm_batchserver<-pm %>%as.data.frame %>% tibble::rownames_to_column('prob')
write.csv(pm_batchserver, "output/all_mapped_pg_matrix_13609_2957_batchserver.csv", row.names = F)
