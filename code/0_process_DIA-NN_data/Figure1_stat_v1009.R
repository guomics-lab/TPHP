setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('../source/source_code.R')


pm0_human <- rio::import(r"(..\..\results\rlt_append_RCA_G3\report.pro_matrix_full_0.001Q_PGQ0.01_wide_human.tsv)")
df <- pm0_human %>% column_to_rownames('Protein.Group')
df[df == 0]<-NA 
df %>% min(na.rm = T) # 12.77158
df1<-df<- df %>%t()%>% as.data.frame()%>%
  tibble::rownames_to_column(var = "FileName")
dim(df1) # 3262 13610
labels <- readRDS('20251009_PUH_sample_information_3005files_qc.RDS')

# change cancer labels
labels %<>% mutate(
  cancer = ifelse(tissue_name == 'reactive hyperplastic lymph nodes', 'Diffused large B-cell carcinoma', cancer),
  cancer_abbr = ifelse(tissue_name == 'reactive hyperplastic lymph nodes', 'DLBCL', cancer_abbr)
) %>% 
  mutate(cancer_subtype = cancer_abbr, .after = class_abbr) %>% 
  mutate(cancer_abbr = ifelse(str_detect(cancer_abbr, '^BRCA'), 'BRCA', cancer_abbr),
         cancer_abbr = ifelse(str_detect(cancer_abbr, 'CCOC|HGSOC'), 'OC', cancer_abbr)) %>% 
  mutate(cancer = ifelse(str_detect(cancer, '^Breast carcinoma'), 'Breast carcinoma', cancer),
         cancer = ifelse(str_detect(cancer, 'ovarian cancer$'), 'Ovarian cancer', cancer))

df_info <- inner_join(labels, df1, by = 'FileName')

# stat of sub-matrices
FN_sub_stat <- df_info %>% 
  filter(sample_type %in% c('F', 'N')) %>% 
  group_by(anatomical_classification, class_abbr) %>% 
  group_modify(~ {
    subm <- .x %>%
      select(-(1:end)) %>%
      removeColsNa()
    
    tibble(
      n_files = nrow(subm),
      n_proteins = ncol(subm)
    )
  }) %>%
  ungroup()

TNT_sub_stat <- df_info %>% 
  filter(sample_type %in% c('T', 'NT')) %>% 
  group_by(cancer, cancer_abbr) %>% 
  group_modify(~ {
    subm <- .x %>%
      select(-(1:end)) %>%
      removeColsNa()
    
    tibble(
      n_files = nrow(subm),
      n_proteins = ncol(subm)
    )
  }) %>%
  ungroup()

# special detailed classes
cls_special <- c('blood', 'brain', 'ear', 'eye', 'nose')
N_sub_stat_special <- df_info %>% 
  filter(sample_type %in% c('N')) %>% 
  filter(anatomical_classification %in% cls_special) %>% 
  group_by(anatomical_classification, class_abbr, tissue_name) %>% 
  group_modify(~ {
    subm <- .x %>%
      select(-(1:end)) %>%
      removeColsNa()
    
    tibble(
      n_files = nrow(subm),
      n_proteins = ncol(subm)
    )
  }) %>%
  ungroup()

# special aggregated classes
# FO, MO, EMBRYO
dat1 <- df_info %>% 
  filter(sample_type %in% c('N')) %>% 
  filter(class_abbr %in% c('FT', 'UTE', 'VA')) %>% 
  mutate(class_agg = 'FO', .before = class_abbr)
dat2 <- df_info %>% 
  filter(sample_type %in% c('N')) %>% 
  filter(class_abbr %in% c('CG', 'EPIDI', 'SEV', 'SED')) %>% 
  mutate(class_agg = 'MO', .before = class_abbr)
dat3 <- df_info %>% 
  filter(sample_type %in% c('F')) %>% 
  mutate(class_agg = 'EMBRYO', .before = class_abbr)

FN_sub_stat_agg <- rbind(dat1, dat2, dat3) %>% 
  group_by(class_agg) %>% 
  group_modify(~ {
    subm <- .x %>%
      select(-(1:end)) %>%
      removeColsNa()
    
    tibble(
      n_files = nrow(subm),
      n_proteins = ncol(subm)
    )
  }) %>%
  ungroup()

list(tumor.nontumor = TNT_sub_stat,
     fetal.normal = FN_sub_stat,
     normal.detailed = N_sub_stat_special,
     fetal.normal.agg = FN_sub_stat_agg) %>% 
  .write_excel('Figure1_stat_v1009.xlsx')





# Check data ------
df_info %>%
  filter(sample_type %in% c('F', 'N')) %>% 
  filter(str_detect(anatomical_classification, 'brain')) %>% 
  View()

df_info %>%
  filter(sample_type %in% c('F', 'N')) %>% 
  filter(str_detect(tissue_name, 'brain')) %>% 
  View()

df_info %>%
  filter(sample_type %in% c('F', 'N')) %>% 
  filter(str_detect(tissue_name, 'penis')) %>% 
  View()

df_info %>% 
  filter(sample_type %in% c('F')) %>% View()

