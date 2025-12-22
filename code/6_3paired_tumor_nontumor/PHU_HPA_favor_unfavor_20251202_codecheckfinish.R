source("../source/source_code.R")

meta <- rio::import('output/T_NT_info_check.xlsx')
df_dep <- rio::import('../6_3paired_tumor_nontumor/output/tumor_nontumor_analysis_reports_v1202.xlsx', sheet = 3)
df_dep %<>% rename(Uniprot = protein)
df_dep$dysregulation <- df_dep$direction
df_dep %<>% inner_join(meta %>% distinct(cancer, cancer_abbr, anatomical_classification))

# read HPA supplementary table 6
hpa_tbls <- read_excel_allsheets('input/Table S6.xlsx')

# map PUH organ names to HPA cancer names
# tphp2hpa <- c('', 'cervical', 'colorectal', '', 'liver', 'lung', '', 'breast', '', 'ovarian', 'pancreatic', '', 'colorectal', '', 'stomach', 'testis', 'head_and_neck', 'thyroid', '', 'endometrial')
# names(tphp2hpa) <- c('brain', 'cervix uteri', 'colon', 'gall bladder', 'liver', 'lung', 'lymph node', 'mammary gland', 'muscle', 'ovary', 'pancreas', 'penis', 'rectum', 'small intestine', 'stomach', 'testis', 'throat', 'thyroid', 'tongue', 'uterus')
# df_dep$cancer <- tphp2hpa[df_dep$organ]
# df_dep$cancer %<>% str_c(., ' cancer') %>% str_replace_all('_', '') %>% str_to_sentence()
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
length(unique(df_hpa$cancer))


# # read uniprot data
# # df_uniprot <- rio::import('Y:/members/jiangwenhao/TPHP/20220908/drug/uniprot/uniprot-download_true_fields_accession_2Cid_2Cprotein_name_2Cgene_na-2022.09.21-14.18.23.22.xlsx')
# df_uniprot <- rio::import('uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cgene_names_-2023.06.16-09.12.08.30.tsv')
df_uniprot <- rio::import('input/uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cgene_names_-2023.06.16-09.12.08.30.tsv')

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
# df %<>% select(-p_t, -p_wilcox, -pAdj_wilcox) # 去掉不涉及的列
df %<>% select(cancer, cancer_abbr, Uniprot, dysregulation, `Prognostic Types`, MATCHED, everything())

# 筛选log-rank p value
df %<>% filter(`log-rank P Values` < 0.05)
rio::export(df, 'PHU_dysregulation_HPA_favor_unfavor_overlap_20251202.xlsx')





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

X <- df_dep %>% count(cancer, name = 'x')
Y <- df %>% count(cancer, name = 'y')
Z <- X %>% inner_join(Y, by = 'cancer') %>% mutate(ratio = round(100 * y / x, 2))
rio::export(Z, 'PUH_dep_HPA_FavorUnfavor_ratio_20251202.xlsx')




# # sankey plot
# colnames(df)
# # [1] "cancer"             "organ"              "Uniprot"            "dysregulation"     
# # [5] "Prognostic Types"   "MATCHED"            "log2FC"             "pAdj_t"            
# # [9] "EnsemblIDs"         "Symbols"            "Mean Expression"    "Sample Numbers"    
# # [13] "Expression Cutoffs" "log-rank P Values" 
library(ggalluvial)
# Sheet1 <- df %>% count(cancer)
Sheet1 <- Z %>% select(cancer, y, ratio) %>% mutate(n = str_glue("{y}, {ratio}%")) %>% select(cancer, n)
Sheet2 <- df %>% count(dysregulation)
Sheet3 <- df %>% count(dysregulation, `Prognostic Types`)
Sheet4 <- df %>% count(`Prognostic Types`)
df_n <- list(Sheet1, Sheet2, Sheet4) %>% plyr::ldply(function(dfsub){
  colnames(dfsub) <- c('from', 'to')
  return(dfsub)
})
df_n %<>% mutate(to = str_glue("{from} (n = {to})"))
ht_fromto <- df_n$to
names(ht_fromto) <- df_n$from

df_alluvium <- df %>% count(cancer, dysregulation, `Prognostic Types`)
df_alluvium %<>% mutate_at(vars(-n), function(x) ht_fromto[x])

# set.seed(2023)
# my_colors <- sample(colorRampPalette(c("#2266BB", "#66BB22", "#BB2266"),bias=1)(n=16))
my_colors <- c('#2B2B6B', '#9D7DDB', '#F6D151', '#D8894E', '#E58989', '#BA6CA4', '#D34A8F', '#103D47',
               '#8080BC', '#55BC6D', '#7C0823', '#2578A0', "#2266BB", "#66BB22", "#BB2266", "#8D7341")
p <- ggplot(
  df_alluvium,
  aes(
    y     = n,
    axis1 = cancer,
    axis2 = dysregulation,
    axis3 = `Prognostic Types`
  )
) +
  # flows
  geom_alluvium(
    aes(fill = cancer),
    alpha    = 0.9,
    width    = 0.5,
    knot.pos = 0.4,
    reverse  = FALSE
  ) +
  # nodes
  geom_stratum(
    width   = 0.5,
    alpha   = 0.55,
    reverse = FALSE,
    color   = "grey20",
    size    = 0.25
  ) +
  # labels on nodes: name (count)
  geom_text(
    stat     = "stratum",
    aes(label = sprintf("%s\n(n = %s)", after_stat(stratum), scales::comma(after_stat(count)))),
    reverse  = FALSE,
    size     = 3,
    lineheight = 0.9,
    family   = "sans"
  ) +
  scale_fill_manual(values = my_colors) +
  guides(fill = "none") +
  scale_x_continuous(
    breaks = 1:3,
    labels = c(
      "Cancer type",
      "Dysregulation (PUH dataset)",
      "Prognostic types (HPA dataset)"
    )
  ) +
  ggtitle("PUH dysregulated – HPA favorable/unfavorable proteins") +
  theme_bw(base_size = 11, base_family = "sans") +
  theme(
    panel.grid      = element_blank(),
    panel.background= element_blank(),
    axis.line       = element_blank(),
    axis.title      = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0, size = 12),
    plot.margin     = margin(6, 6, 6, 6)
  )


ggsave('PUH_dep_HPA_FavorUnfavor_sankey_20251202.pdf', p, height = 10, width = 7)




