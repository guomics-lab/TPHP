setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../source/source_code.R")



# 1. TCGA, CPTAC ----------------------------------------------------------
res <- readxl::read_xlsx("PUH_TCGA_CPTAC_DEA_nofilterP_chr_20251202.xlsx")
res %<>% mutate_at(vars(PUH_gAdj, TCGA_log2FC, CPTAC_log2FC), as.numeric)
df_TCGA <- res %>%
  filter(!is.na(organ)) %>% 
  filter((abs(PUH_gAdj) > 0.5) & (abs(TCGA_log2FC) > log2(1.6))) %>% 
  # filter(abs(TCGA_consistency) == 1) %>% 
  select(protein, organ, PUH_gAdj, TCGA_log2FC, TCGA_consistency)

df_CPTAC <- res %>% 
  filter(!is.na(organ)) %>% 
  filter((abs(PUH_gAdj) > 0.5) & (abs(CPTAC_log2FC) > log2(1.2))) %>%
  # filter(abs(CPTAC_consistency) == 1) %>% 
  select(protein, organ, PUH_gAdj, CPTAC_log2FC, CPTAC_consistency)

# organ name to cancer name
ht_organ2cancer <- c('Glioblastoma', 'Cervical carcinoma', 'Colon carcinoma', 'Esophageal carcinoma', 'Fallopian tube carcinoma', 'Gallbladder carcinoma', 'Renal carcinoma', 'Hepatocellular carcinoma', 'Lung carcinoma', 'Diffused large B-cell carcinoma', 'Breast carcinoma', 'Muscle tumor', 'Ovarian carcinoma', 'Pancreas carcinoma', 'Male genitalia carcinoma', 'Prostate carcinoma', 'Rectum carcinoma', 'Gastrointestinal stromal tumors', 'Gastric carcinoma', 'Testis carcinoma', 'Laryngocarcinoma', 'Thymoma and thymic carcinoma', 'Thyroid carcinoma', 'Laryngocarcinoma', 'Endometrial carcinoma')
names(ht_organ2cancer) <- c('brain', 'cervix_uteri', 'colon', 'esophagus', 'fallopian tube', 'gall bladder', 'kidney', 'liver', 'lung', 'lymph node', 'mammary gland', 'muscle', 'ovary', 'pancreas', 'penis', 'prostate', 'rectum', 'small intestine', 'stomach', 'testis', 'throat', 'thymus', 'thyroid', 'tongue', 'uterus')
df_TCGA$organ <- ht_organ2cancer[df_TCGA$organ]
df_CPTAC$organ <- ht_organ2cancer[df_CPTAC$organ]


df_TCGA %<>% add_column(TCGA.Mapped = TRUE, .before = 1)
df_CPTAC %<>% add_column(CPTAC.Mapped = TRUE, .before = 1)

df_comb <- full_join(df_TCGA, df_CPTAC)
df_comb %>% rio::export('comb_TCGA_CPTAC_filtered_v1202.xlsx')




df_comb_chr <- rio::import('PUH_TCGA_CPTAC_DEA_nofilterP_chr_20251202.xlsx')



# 2. HPA ------------------------------------------------------------------
meta <- rio::import('output/T_NT_info_check.xlsx')
df_dep <- rio::import('../6_3paired_tumor_nontumor/output/tumor_nontumor_analysis_reports_v1202.xlsx', sheet = 3)
df_dep %<>% rename(Uniprot = protein)
df_dep$dysregulation <- df_dep$direction
df_dep %<>% inner_join(meta %>% distinct(cancer, cancer_abbr, anatomical_classification))

# read HPA supplementary table 6
hpa_tbls <- read_excel_allsheets('input/Table S6.xlsx')

# map PUH organ names to HPA cancer names
hpa2tphp_df <- tibble::tribble(
  ~tphp,                               ~hpa,
  "Breast carcinoma",                  "Breast cancer",
  "Cervical carcinoma",                "Cervical cancer",
  "Colon carcinoma",                   "Colorectal cancer",
  "Diffused large B-cell carcinoma",   NA,
  "Endometrial carcinoma",             "Endometrial cancer",
  "Esophageal carcinoma",              NA,
  "Fallopian tube carcinoma",          NA,
  "Gallbladder carcinoma",             NA,
  "Gastric carcinoma",                 "Stomach cancer",
  "Gastrointestinal stromal tumors",   NA,
  "Glioblastoma",                      "Glioma",
  "Hepatocellular carcinoma",          "Liver cancer",
  "Laryngocarcinoma",                  "Head and neck cancer",
  "Lung carcinoma",                    "Lung cancer",
  "Male genitalia carcinoma",          NA,
  "Muscle tumor",                      NA,
  "Ovarian cancer",                    "Ovarian cancer",
  "Pancreas carcinoma",                "Pancreatic cancer",
  "Prostate carcinoma",                "Prostate cancer",
  "Rectum carcinoma",                  "Colorectal cancer",
  "Renal carcinoma",                   "Renal cancer",
  "Testis carcinoma",                  "Testis cancer",
  "Thymoma and thymic carcinoma",      NA,
  "Thyroid carcinoma",                 "Thyroid cancer",
  "Tongue carcinoma",                  NA,
  NA,                                  "Melanoma",
  NA,                                  "Urothelial cancer"
)
hpa2tphp <- hpa2tphp_df %>% drop_na() %>% pull(tphp, hpa)


df_hpa <- plyr::ldply(hpa_tbls, .id = 'cancer')
df_hpa %<>%
  mutate(cancer = hpa2tphp[cancer]) %>% 
  filter(cancer %in% unique(df_dep$cancer))
length(unique(df_hpa$cancer)) # 16


# # read uniprot data
df_uniprot <- rio::import('../6_3paired_tumor_nontumor/input/uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cgene_names_-2023.06.16-09.12.08.30.tsv')

# combine
dfa <- df_uniprot %>% select(From, Entry) %>% rename(EnsemblIDs = From, Uniprot = Entry) %>% 
  left_join(df_hpa, .) %>% 
  inner_join(df_dep, .)
dfb <- dfprot %>% rename(Symbols = Genes, Uniprot = Protein.Group) %>% 
  left_join(df_hpa, .) %>% 
  inner_join(df_dep, .)
df <- rbind(dfa, dfb %>% anti_join(dfa))
dim(df) # 52568    26

# 标记
df$MATCHED <- (df$`Prognostic Types` == 'Unfavorable' & df$dysregulation == 'Up') | (df$`Prognostic Types` == 'Favorable' & df$dysregulation == 'Down')

# 列重排&筛选
df %<>% select(cancer, cancer_abbr, Uniprot, dysregulation, `Prognostic Types`, MATCHED, everything())

# 筛选log-rank p value
df %<>% filter(`log-rank P Values` < 0.05)


df %>% count(cancer)
# cancer    n
# 1          Breast carcinoma 1828
# 2        Cervical carcinoma 1882
# 3           Colon carcinoma 2384
# 4     Endometrial carcinoma 1479
# 5         Gastric carcinoma 1741
# 6              Glioblastoma  188
# 7  Hepatocellular carcinoma 2280
# 8          Laryngocarcinoma 1359
# 9            Lung carcinoma  894
# 10           Ovarian cancer 1764
# 11       Pancreas carcinoma 1754
# 12       Prostate carcinoma  546
# 13         Rectum carcinoma 3858
# 14          Renal carcinoma 1199
# 15         Testis carcinoma  897
# 16        Thyroid carcinoma  257

df %>% count(dysregulation)
#   dysregulation     n
# 1          Down  8698
# 2            Up 15612

df %>% count(`Prognostic Types`)
# Prognostic Types     n
# 1        Favorable 11722
# 2      Unfavorable 12588

length(unique(df$Uniprot)) # 7336



# 3. combine --------------------------------------------------------------
# dfprot <- rio::import('//172.16.13.136/TPHP/TPL/libs/20220616_fraglib/protein.tsv')
# df_comb %<>% left_join(dfprot, by = c('protein' = 'Protein ID'))

# df_comb %<>% left_join(df_dep %>% distinct(Uniprot, Genes) %>% rename(protein = Uniprot))
df_comb <- df_comb_chr %>% left_join(df_dep %>% distinct(Uniprot, Genes) %>% rename(protein = Uniprot))

df %<>% rename(HPA_cancer = cancer)
df %<>% add_column(protein = df$Uniprot, .before = 1)

tmpcol <- intersect(colnames(df), colnames(df_comb)) # "protein"  "cancer_abbr" "Genes" "anatomical_classification"
df_comb[1,tmpcol]; df[1,tmpcol]
unique(df$cancer_abbr) # 16
unique(df_comb$cancer_abbr) # 17
# > df$cancer_abbr %>% setdiff(df_comb$cancer_abbr, .)
# [1] "ESCA" "GBCA" "TOCA"
# > df$cancer_abbr %>% setdiff(df_comb$cancer_abbr)
# [1] "ENCA" "TGCT"

df_final_outer <- df_comb %>% full_join(df) %>% select(cancer_abbr, anatomical_classification, protein, Genes, everything())
df_final <- df_comb %>% inner_join(df) %>% select(cancer_abbr, anatomical_classification, protein, Genes, everything())
df_final1 <- df_final %>% filter(MATCHED)

list(
  `PUH-TCGA&CPTAC-HPA` = df_final1 %>% filter(TCGA_significant == 'TRUE', CPTAC_significant == 'TRUE') %>%
    arrange(cancer_abbr, protein),
  `PUH-TCGA|CPTAC-HPA` = df_final1 %>% arrange(cancer_abbr, protein),
  `PUH-TCGA-CPTAC-HPA_concat` = df_final_outer %>% arrange(cancer_abbr, protein)
) %>% rio::export('PUH_TCGA_CPTAC_HPA_mapping_v1202.xlsx')




