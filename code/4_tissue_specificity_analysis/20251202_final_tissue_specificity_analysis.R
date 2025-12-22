rm(list = ls())
library(tidyverse)
library(magrittr)
library(readxl)
library(writexl)
library(tidyr)
library(dplyr)
library(stringr)
# Working directory
setwd("//172.16.13.136/tphp/code/4_tissue_specificity_analysis/")

# Source local helpers
source("../source/source_code.R")

# 2023-02-07
## Ulhen's classification method
# df<-read.csv("//172.16.13.136/tphp/code/4_tissue_specificity_analysis/input/data_matrix_filtered_05NA_2997samples_12754proteins.csv",row.names = 1,check.names = F)
df<-read.csv("//172.16.13.136/tphp/code/4_tissue_specificity_analysis/input/data_matrix_filtered_05NA_2957_12797.csv",row.names = 1,check.names = F)
# df1<-df
head(df[,1:10])
df<-apply(df,c(1,2),as.numeric) 
ai<-read_xlsx("//172.16.13.136/tphp/code/4_tissue_specificity_analysis/input/20251009_PUH_sample_information_3005files_V6.xlsx")
ai<-as.data.frame(ai)
str(ai)
df1<-df[rownames(df)%in%as.character(ai$FileName),]

# di<-inner_join(ai,df[,c(1,16:ncol(df))],by="FileName")
# table(di$patient_ID)
ai0<-ai%>% dplyr::filter((sample_type == 'N') & (Tissue_heterogeneity_abnormal == FALSE) & (Low_protein_IDs ==FALSE)  & (patient_ID %in% c("DL_1","DL_2","DL_3","DL_4","DL_7","JD_1","JD_2","JD_3","JD_4")))
df0<-df[ai0$FileName,]
dim(df0) #12754 pros,  466 samples
#remove high missing values 
source("//172.16.13.136/tphp/code/2_batch_effect_evaluation_correction/src_claude_v1_vsc.R",encoding = "utf-8")
df_filter<-removeHighNAProteins(df0, ai0,
                                group_columns = c("anatomical_classification"),
                                na_threshold = 0.5,group_threshold = 1.0)  
# Group by anatomical_classification: 64 groups
# Total: 12754 proteins; removed: 2585 proteins
# (NA rate >= 0.50 in >= 100.0% of groups)
# Removed proteins (first 10): A4FU69, A5PLK6, A6NHZ5, A8K0R7, A8MU46, B4DJY2, B7ZAP0, F5H4A9, O60259, O60830...

df0<-df_filter$filtered_data
dim(df0) #[1]   466 10169
# df1<-apply(df0[,-c(1:18)],c(1,2),as.numeric)
min(df0,na.rm = T) # 2.091316
max(df0,na.rm=T) # 28.62249
df1<-cbind(ai0,df0)

rio::export(df1, '20251031_TPHP_normalData_for_classification_naive_466_10169.xlsx')

####calculate median 
####impute NA with 0.5*min
pm1<-df0
pm1[is.na(pm1)]<-0.5*min(pm1,na.rm=T)
dim(pm1)
pm2<-aggregate(pm1,by=list(ai0$class_abbr),median)

########ulhen tissue specificity analysis
X <-pm2
colnames(X)[1]<-"class_abbr"
clas <- ulhen_class(X, fct = 4)

#summary ulhen class result
message('Specificity visulization of matrix ', nm)
res.clas1 <- clas[, apply(clas, 2, function(y) !all(y == 'not detected'))]
ts.sub <- res.clas1 %>% # tissue specificity
    rename(class_abbr = tissue_type) %>% 
    pivot_longer(cols = -class_abbr, names_to = 'Protein', values_to = 'Tissue specificity') %>% 
    mutate(`Tissue specificity` = str_to_sentence(`Tissue specificity`)) %>% 
    inner_join(dfprot %>% rename(Protein = Protein.Group))
ts.sub.Enrich <- ts.sub %>% filter(str_detect(`Tissue specificity`, 'enriched'))
ts.sub.tisEnrich <- ts.sub.Enrich %>% filter(`Tissue specificity` == 'Tissue enriched')
ret <- list(metadata = df1,
            med.matrix=X,
            ts.df=clas, 
            ts.all = ts.sub, 
            ts.Enrich = ts.sub.Enrich, 
            ts.tisEnrich = ts.sub.tisEnrich)
saveRDS(ret, '20251126_tissue_Specificity_naive_05minImputated_3.rds')

#####20251202 remove contaminant proteins
contam<-read_xlsx("../0_QC/input/599blk_protein_0.1FreqPeptide.xlsx") #266
ret<-readRDS("INPUT/20251126_tissue_Specificity_naive_05minImputated_3.rds")
ret$ts.Enrich_rm_contam<-ret$ts.Enrich %>%
    filter(!Protein %in% contam$Protein.ID)
ret$ts.tisEnrich_rm_contam<-ret$ts.tisEnrich %>%
    filter(!Protein %in% contam$Protein.ID)
setdiff(ret$ts.tisEnrich$Protein,ret$ts.tisEnrich_rm_contam$Protein)
saveRDS(ret, 'output/20251201_tissue_Specificity_naive_05minImputated_3.rds')
dim(ret$ts.df)
