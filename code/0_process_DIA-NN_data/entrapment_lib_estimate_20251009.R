# Reference: https://doi.org/10.1038/s41592-025-02719-x
##  FDP.combined <- N.ε * (1 + 1/r) / (N.τ + N.ε)
# rstudioapi::getActiveDocumentContext()$path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('../source/source_code.R')
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(cowplot)
library(pheatmap)
library(broom)
library(arrow)

TPHP_HOME <- '//192.168.99.100/TPHP/'


fdp_est <- function(N.τ, N.ε, r) {
  FDP.combined <- N.ε * (1 + 1/r) / (N.τ + N.ε)
  return(FDP.combined)
}


# # 0.Entrapment library prepare ----------
lib_tphp <- rio::import('//192.168.99.100/TPHP/library/TPHPlib_frag1025_swissprot_final.tsv')
lib_entrap <- rio::import('//192.168.99.100/TPHP/library/tphpSwiss1025_mpwIGC_IGC88_library.tsv') %>%
  filter(!(PeptideSequence %in% lib_tphp$PeptideSequence)) %>% 
  mutate(AverageExperimentalRetentionTime = NA)
# length(unique(lib_tphp$PeptideSequence)) # 484,391
# length(unique(lib_entrap$PeptideSequence)) # 558,085   # 1,042,476
# setdiff(names(lib_tphp), names(lib_entrap)) # AverageExperimentalRetentionTime

pro_tphp <- rio::import('../../TPL/libs/20220616_fraglib/protein.tsv')
real.proteotipic <- pro_tphp %>% filter(`Indistinguishable Proteins` == '') %>% pull(`Protein ID`)
false.proteotipic <- pro_tphp %>% filter(`Indistinguishable Proteins` != '') %>% pull(`Protein ID`)

dim(lib_tphp) # 16722511       16
dim(lib_entrap) # 16041520       18


X1 <- lib_tphp %>% distinct(ModifiedPeptideSequence, PrecursorCharge) %>% 
  unite(col = 'Precursor.Id', ModifiedPeptideSequence, PrecursorCharge, sep = '', remove = F)
nrow(X1) # 689,568
X2 <- lib_entrap %>% distinct(ModifiedPeptideSequence, PrecursorCharge) %>% 
  unite(col = 'Precursor.Id', ModifiedPeptideSequence, PrecursorCharge, sep = '', remove = F)
nrow(X2) # 739,318    # 1,428,886

length(unique(lib_tphp$ProteinId)) # 15,332
length(unique(lib_entrap$ProteinId)) # 256,259    # 271,591


r_pr <- length(unique(lib_entrap$PeptideSequence)) / length(unique(lib_tphp$PeptideSequence)) # 1.152137    # 2.152137
r_pep <- nrow(X2) / nrow(X1) # 1.072147       # 2.072147
r_pg <- length(unique(lib_entrap$ProteinId)) / length(unique(lib_tphp$ProteinId)) # 16.714    # 17.714

tbl1 <- data.frame(name = c('r_pr', 'r_pep', 'r_pg'),
                   value = c(r_pr, r_pep, r_pg))

# 1.Data readin ------
tphp <- arrow::read_parquet(file.path(TPHP_HOME, 'results/rlt_append_RCA_G3/report.pro_matrix_full_0005LibQval_PGMaxLFQ.parquet'))


# TO BE CONTINUE
# info1 <- readxl::read_excel('20250811_PUH_sample_information_1900files_v11.xlsx') %>% as.data.frame()
# info1 %<>% mutate(ID = sapply(FileName, function(x){
#   str_extract(x, '^(\\w+?)\\d+.+_(\\d+)$', group = c(1, 2)) %>%
#     str_c(collapse = '_')
# }), .before = 1)
df_abbr <- read.delim('//172.16.13.136/share/members/jiangwenhao/TPHP/input/sample_types_abbr_20230113.txt', stringsAsFactors = F, check.names = F, na.strings = '')
h <- c(df_abbr$Abbr, 'Pool') %>% setNames(c(df_abbr$Entire, 'pool'))
# h <- hash::hash(keys = df_abbr$'Entire', values = df_abbr$'Abbr')


info1.qc <- rio::import('../6_pancancer_harmonizR_QC/input/20250811_PUH_sample_information_1900files_v11_forQC.xlsx') %>%
  mutate(DateTime_parsed = as.POSIXct(DateTime, format = "%Y-%m-%d_%H-%M-%S"),
         Date = as.Date(DateTime_parsed),
         # YearMonth = ifelse(is.na(Date), NA_character_, format(Date, "%Y-%m")),
         Date_char = ifelse(is.na(Date), NA_character_, format(DateTime_parsed, "%Y%m%d%H%M%S")),
         acq_order = dplyr::dense_rank(
           dplyr::if_else(is.na(Date_char), paste0("Z", row_number()), paste0("A", Date_char))
         ))

info2 <- rio::import('20250812_PUH_RCA_sample_and_pool_1179files_info.xlsx')
info2.sample <- rio::import('TPHP_RCA_sample_info_20250807.xlsx') %>% 
  mutate(patient_ID = `住院号`, ana_cn = `癌症类型`, sample_type = Cancer_type,
         sample_id = str_replace(BatchID, '\\-', '_')) %>% 
  distinct(sample_id, patient_ID, ana_cn) %>% 
  full_join(info2 %>% distinct(tissue_name, anatomical_classification, patient_ID)) %>% 
  mutate(tissue_name = ifelse(ana_cn == '胶质母细胞瘤', 'glioblastoma', tissue_name),
         tissue_name = ifelse(ana_cn == '反应性增生淋巴结', 'reactive hyperplastic lymph nodes', tissue_name),
         anatomical_classification = ifelse(ana_cn == '胶质母细胞瘤', 'brain', anatomical_classification),
         anatomical_classification = ifelse(ana_cn == '反应性增生淋巴结', 'lymph node', anatomical_classification)) %>% 
  drop_na()

info.pc <- rio::import('../6_pancancer_harmonizR_QC/output/PUH_pancancer_sample_information_2datesets_1146files.xlsx') %>%
  select(-batch, -Date) # including bladder cancer
# dim(info1.qc %>% inner_join(info.pc))
labelT <- rio::import('../0_process_DIA-NN_data/input/TPHP_cancer_abbr_labels.xlsx')
labelTs <- rio::import('../0_process_DIA-NN_data/input/TPHP_cancer_patients1006.xlsx')

#new
dnn.pg <- rio::import('../../results/rlt_append_RCA_G3/report.pg_matrix.tsv')
info.new <- data.frame(
  FileName = colnames(dnn.pg %>% select(-(1:N.Proteotypic.Sequences))) %>% str_remove_all('.+\\\\|.+/|\\.d$')#unique(tphp$Run)
)
info.new1 <- info.new %>% filter(str_detect(FileName, 'RCA', negate = T)) %>% 
  inner_join(info1.qc) %>% 
  select(FileName, DIA_ID:anatomical_classification, patient_ID, DateTime:acq_order)
info.new2 <- info.new %>% filter(str_detect(FileName, 'RCA')) %>% 
  mutate(sample_id = str_extract(FileName, 'b\\d+_(\\d+|pool)')) %>% 
  left_join(info2.sample)
info.new2a <- info.new2 %>% filter(!str_detect(sample_id, 'pool')) %>% 
  filter(!is.na(patient_ID)) %>%
  full_join(info2.sample %>% select(-ana_cn)) %>% 
  inner_join(info2 %>% distinct(sample_id, patient_ID, sample_type))
info.new2b <- info.new2 %>% filter(str_detect(sample_id, 'pool')) %>%
  mutate(sample_type = 'p') %>% 
  mutate_at(vars(tissue_name:anatomical_classification), function(x) 'pool')
info.new2 <- rbind(info.new2a, info.new2b) %>%
  mutate(new = 'data3',
         instrument = str_extract(FileName, '[A-Z]+'),
         date.in.filename = str_extract(FileName, '\\d+'),
         year = str_sub(date.in.filename, 1, 4),
         month = str_sub(date.in.filename, 5, 6)) %>% 
  mutate(end = str_extract(FileName, '\\d+$')) %>% 
  arrange(end)# %>% 
  # arrange(date.in.filename) %>%
  # mutate(acq_order = 1:nrow(.),
  #        acq_order = ifelse(is.na(FileName), NA, acq_order))

df_info <- plyr::rbind.fill(info.new1, info.new2) %>% 
  mutate(tissue_name = str_replace_all(tissue_name, intToUtf8('0x00A0'), ' '),
         anatomical_location = str_replace_all(anatomical_location, intToUtf8('0x00A0'), ' ') %>% 
           str_remove_all('[\\r\\n]+') %>% str_trim())
# df_info %>% filter(tissue_name == 'reactive hyperplastic lymph nodes')

info.all1 <- df_info %>% filter(!sample_type %in% c('T', 'NT'))
info.all2a <- df_info %>% filter(sample_type %in% c('T', 'NT')) %>%
  filter(patient_ID %in% (df_info %>% filter(tissue_name == 'reactive hyperplastic lymph nodes') %>% pull(patient_ID)))
info.all2b <- df_info %>% filter(sample_type %in% c('T', 'NT')) %>%
  filter(patient_ID %in% labelTs$patient_ID) %>%
  # semi_join(labelTs %>% select(patient_ID)) %>% 
  select(-tissue_name) %>% 
  inner_join(labelTs %>% distinct(patient_ID, tissue_name))
info.all2 <- rbind(info.all2a, info.all2b)

rm.pool <- 'CAD20250801yuel_TPHP_RCA_90minDIA_b8_pool2_Slot1-3_1_13194' # it's an insufficient injection, and has been re-injected
info.all <- rbind(info.all1, info.all2) %>% 
  left_join(labelT) %>%
  arrange(date.in.filename, Date_char, end) %>%
  mutate(acq_order = 1:nrow(.)) %>% 
  relocate(cancer, cancer_abbr, .after = anatomical_classification) %>% 
  mutate(class_abbr = h[anatomical_classification], .after = anatomical_classification) %>% 
  filter(!FileName %in% rm.pool)
# info.all %>% filter(tissue_name == 'reactive hyperplastic lymph nodes')


info.all.clean <- info.all %>% filter(!tissue_name %in% c('bladder carcinoma', 'skin carcinoma')) #3005
rio::export(info.all.clean, '20251009_PUH_sample_information_3005files_qc.xlsx')
saveRDS(info.all.clean, '20251009_PUH_sample_information_3005files_qc.RDS')
# info.all %>% filter(is.na(cancer_abbr)) %>% dplyr::filter(sample_type %in% c('T', 'NT')) %>% 
#   rio::export('T_NT_info_check.xlsx')
info.all.mini <- info.all.clean %>%
  mutate(date = ifelse(!is.na(date), date, date.in.filename)) %>%
  select(FileName, sample_type:patient_ID,
         DIA_ID, sample_id, instrument:trans, new)
rio::export(info.all.mini, '20251009_PUH_sample_information_3005files.xlsx')



tphp %<>%
  mutate(Is.pr.Target = Precursor.Id %in% X1$Precursor.Id,
         Is.pg.Target = Protein.Group %in% unique(lib_tphp$ProteinId))

# 2.FDP estimate -------
## global pr FDP
tmp <- tphp %>% distinct(Precursor.Id, Is.pr.Target) %>% dplyr::count(Is.pr.Target) %>% pull(n)
global.pr.FDP.combined <- fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pr)

## global pg FDP
tmp <- tphp %>% distinct(Protein.Group, Is.pg.Target) %>% dplyr::count(Is.pg.Target) %>% pull(n)
global.pg.FDP.combined <- fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pg)

## specific pr FDP
runSpecific.pr.FDP.combined <- tphp %>% plyr::ddply('Run', function(dfsub){
  tmp <- dfsub %>% distinct(Protein.Group, Is.pr.Target) %>% dplyr::count(Is.pr.Target) %>% pull(n)
  if(length(tmp) == 1) tmp <- c(0, tmp)
  data.frame(
    Run = dfsub$Run[1],
    runSpecific.pr.FDP.combined = fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pr)
  )
}) %>% arrange(desc(runSpecific.pr.FDP.combined))

## specific pg FDP
runSpecific.pg.FDP.combined <- tphp %>% plyr::ddply('Run', function(dfsub){
  tmp <- dfsub %>% distinct(Protein.Group, Is.pg.Target) %>% dplyr::count(Is.pg.Target) %>% pull(n)
  if(length(tmp) == 1) tmp <- c(0, tmp)
  data.frame(
    Run = dfsub$Run[1],
    runSpecific.pg.FDP.combined = fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pg)
  )
}) %>% arrange(desc(runSpecific.pg.FDP.combined))

# identity
df_identity <- tphp %>% plyr::ddply('Run', function(dfsub){
  npg <- dfsub %>%
    distinct(Protein.Group, Is.pg.Target) %>%
    nrow()
  npr <- dfsub %>%
    distinct(Precursor.Id, Is.pr.Target) %>%
    nrow()
  data.frame(`# precursors` = npr, `# protein groups` = npg, check.names = F)
})

tbl2 <- data.frame(
  global.pr.FDP.combined = global.pr.FDP.combined,
  global.pg.FDP.combined = global.pg.FDP.combined
)

dfbox <- rbind(runSpecific.pg.FDP.combined %>% setNames(c('Run', 'FDP')) %>% mutate(Type = 'pg'),
               runSpecific.pr.FDP.combined %>% setNames(c('Run', 'FDP')) %>% mutate(Type = 'pr'))
p_runFDP <- ggplot(dfbox) +
  aes(x = Type, y = FDP, fill = Type) +
  geom_boxplot(width = 0.6, position = position_dodge()) +
  # geom_text(data = dfbox %>% filter(FDP > 0.01),
  #           aes(label = round(FDP, 4)), hjust = -0.5) +
  labs(x = str_glue("DIA runs\n(n={length(unique(dfbox$Run))})"), y = 'Run-specific FDP') +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0, 0.05, 0)) +
  # scale_fill_manual(values = target_decoy_color) +
  ggsci::scale_fill_npg() +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  theme(text = element_text(size = 12))
ggsave('entrapment_estimate_combined_runspecific.pdf', p_runFDP, width = 5, height = 4)


bad_runs <- runSpecific.pr.FDP.combined %>%
  filter(runSpecific.pr.FDP.combined > 0.01) %>% 
  pull(Run)
bad_runs_info <- info.all.clean %>%
  filter(str_detect(FileName, str_c(bad_runs, collapse = '|')))

bad_runs_info_annotate <- bad_runs_info %>%
  distinct(sample_type, tissue_name) %>%
  semi_join(info.all.clean, .) %>%
  inner_join(df_identity %>% rename(FileName = Run)) %>%
  inner_join(runSpecific.pr.FDP.combined %>% rename(FileName = Run)) %>%
  inner_join(runSpecific.pg.FDP.combined %>% rename(FileName = Run)) %>%
  arrange(sample_type, tissue_name, runSpecific.pr.FDP.combined, runSpecific.pg.FDP.combined)


tbl3 <- df_identity %>%
  inner_join(runSpecific.pr.FDP.combined) %>% 
  inner_join(runSpecific.pg.FDP.combined)

# OUTPUT -------
list(r.values = tbl1,
     global.FDP = tbl2,
     # run.specific.pr.FDP = runSpecific.pr.FDP.combined,
     # run.specific.pg.FDP = runSpecific.pg.FDP.combined,
     identity = tbl3,
     bad.runs = bad_runs_info,
     bad.run.groups = bad_runs_info_annotate
     ) %>% 
  rio::export('entrapment_lib_estimate_source_data_20251009.xlsx')


# save.image('entrapment_lib_estimate_20251009.RData')


# try different thresholds --------
### benchmark -----
# rep<-read_parquet(r"(..\..\results\report-CAF20250725_b12_13\report.parquet)", 
#                   col_select = c("Precursor.Id", "Protein.Group", "Run","Proteotypic","Protein.Group","Q.Value","Global.Q.Value",
#                                  "PG.Q.Value","Global.PG.Q.Value","Precursor.Quantity","Lib.Q.Value","Lib.PG.Q.Value","PG.MaxLFQ"
#                   ))
# dim(rep)#[1]  74883157        12
# pep<-rep %>% 
#   filter(Proteotypic == 1 , 
#          # Lib.Q.Value < 0.005,
#          Global.Q.Value <5e-4 , 
#          Q.Value< 5e-4 ,  
#          Precursor.Quantity != 0
#          ) 
# 
# pep_sum<-length(unique(pep$Precursor.Id))
# pep_sum  # 333400
# 
# pro<-pep %>% 
#   filter(#Lib.PG.Q.Value < 0.005,
#          Global.PG.Q.Value<5e-4,
#          PG.Q.Value < 5e-4
#          # Global.Q.Value <0.01 , 
#          # Q.Value< 0.01 ,  
#          # Proteotypic == 1 , 
#   ) 
# pro_sum<-length(unique(pro$Protein.Group))
# pro_sum  # 13608
# 
# thresholds <- c(1e-2, 5e-3, 1e-3, 5e-4, 4e-4, 3e-4, 2e-4)
# # reports_test_thresholds1
# 
# reports_test_thresholds <- lapply(thresholds, function(thre){
#   cat(thre, '...\n')
#   tphp <- rep %>%
#     filter(Proteotypic == 1, Precursor.Quantity != 0,
#            PG.Q.Value < thre, Global.PG.Q.Value < thre,
#            Q.Value < thre, Global.Q.Value < thre) %>% 
#     mutate(Is.pr.Target = Precursor.Id %in% X1$Precursor.Id,
#            Is.pg.Target = Protein.Group %in% unique(lib_tphp$ProteinId))
#   ## global pr FDP
#   tmp <- tphp %>% distinct(Precursor.Id, Is.pr.Target) %>% dplyr::count(Is.pr.Target) %>% pull(n)
#   global.pr.FDP.combined <- fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pr)
#   pr.N.τ <- tmp[2]
#   pr.N.ε <- tmp[1]
#   
#   ## global pg FDP
#   tmp <- tphp %>% distinct(Protein.Group, Is.pg.Target) %>% dplyr::count(Is.pg.Target) %>% pull(n)
#   global.pg.FDP.combined <- fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pg)
#   pg.N.τ <- tmp[2]
#   pg.N.ε <- tmp[1]
#   
#   ## specific pr FDP
#   runSpecific.pr.FDP.combined <- tphp %>% plyr::ddply('Run', function(dfsub){
#     tmp <- dfsub %>% distinct(Protein.Group, Is.pr.Target) %>% dplyr::count(Is.pr.Target) %>% pull(n)
#     if(length(tmp) == 1) tmp <- c(0, tmp)
#     data.frame(
#       Run = dfsub$Run[1],
#       runSpecific.pr.FDP.combined = fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pr)
#     )
#   }) %>% arrange(desc(runSpecific.pr.FDP.combined))
#   
#   ## specific pg FDP
#   runSpecific.pg.FDP.combined <- tphp %>% plyr::ddply('Run', function(dfsub){
#     tmp <- dfsub %>% distinct(Protein.Group, Is.pg.Target) %>% dplyr::count(Is.pg.Target) %>% pull(n)
#     if(length(tmp) == 1) tmp <- c(0, tmp)
#     data.frame(
#       Run = dfsub$Run[1],
#       runSpecific.pg.FDP.combined = fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pg)
#     )
#   }) %>% arrange(desc(runSpecific.pg.FDP.combined))
#   
#   # identity
#   df_identity <- tphp %>% plyr::ddply('Run', function(dfsub){
#     npg <- dfsub %>%
#       distinct(Protein.Group, Is.pg.Target) %>%
#       nrow()
#     npr <- dfsub %>%
#       distinct(Precursor.Id, Is.pr.Target) %>%
#       nrow()
#     data.frame(`# precursors` = npr, `# protein groups` = npg, check.names = F)
#   })
#   
#   tbl2 <- data.frame(
#     global.pr.FDP.combined = global.pr.FDP.combined,
#     global.pg.FDP.combined = global.pg.FDP.combined
#   )
#   
#   dfbox <- rbind(runSpecific.pg.FDP.combined %>% setNames(c('Run', 'FDP')) %>% mutate(Type = 'pg'),
#                  runSpecific.pr.FDP.combined %>% setNames(c('Run', 'FDP')) %>% mutate(Type = 'pr'))
#   
#   tbl3 <- df_identity %>%
#     inner_join(runSpecific.pr.FDP.combined) %>% 
#     inner_join(runSpecific.pg.FDP.combined)
#   
#   global.identity <- list(
#     pg.N.ε = pg.N.ε,
#     pr.N.τ = pr.N.τ,
#     pr.N.ε = pr.N.ε,
#     pg.N.τ = pg.N.τ
#   )
#   
#   ret <- list(global.FDP = tbl2,
#               run.specific.FDP = dfbox,
#               global.identity = global.identity,
#               identity = tbl3)
#   return(ret)
# })
# names(reports_test_thresholds) <- thresholds
# 
# sum.tbl2 <- plyr::ldply(reports_test_thresholds, function(ret){
#   ret$global.FDP
# }, .id = 'Q.cutoff')
# 
# sum.box <- plyr::ldply(reports_test_thresholds, function(ret){
#   ret$run.specific.FDP
# }, .id = 'Q.cutoff')
# 
# sum.identity.glb <- plyr::ldply(reports_test_thresholds, function(ret){
#   as.data.frame(ret$global.identity)
# }, .id = 'Q.cutoff')
# 
# sum.tbl3 <- plyr::ldply(reports_test_thresholds, function(ret){
#   ret$identity
# }, .id = 'Q.cutoff')
# 
# # ggplot(sum.tbl3) +
# #   aes(x = `# protein groups`, y = Q.cutoff, color = Q.cutoff) +
# #   geom_boxplot(width = 0.6) +
# #   # labs(x = str_glue("DIA runs\n(n={length(unique(dfbox$Run))})"), y = 'Run-specific FDP') +
# #   # scale_fill_manual(values = target_decoy_color) +
# #   ggsci::scale_color_npg() +
# #   theme_classic() +
# #   theme(axis.text.x = element_blank()) +
# #   theme(text = element_text(size = 12))
# # 
# # ggplot(sum.box) +
# #   aes(x = FDP, y = Q.cutoff, color = Q.cutoff) +
# #   geom_boxplot()
# 
# Qcutoff.selection <- sum.identity.glb %>% inner_join(sum.tbl2) %>% 
#   relocate(pg.N.ε, .after = pg.N.τ) %>% 
#   rename(global.pr.FDP = global.pr.FDP.combined,
#          global.pg.FDP = global.pg.FDP.combined)
# Qcutoff.run.specific <- sum.tbl3
# 
# # pheatmap(apply(Qcutoff.selection, 1:2, as.numeric),
# #          display_numbers = T, number_format = '%.2e', color = NA, number_color = 'black',
# #          cluster_rows = F, cluster_cols = F)
# 
# plot_Qcut.global <- Qcutoff.selection %>%
#   rename(`pr.FDP\n(global)` = global.pr.FDP,
#          `pg.FDP\n(global)` = global.pg.FDP) %>% 
#   pivot_longer(cols = -Q.cutoff) %>% 
#   ggplot() +
#   aes(x = name, y = Q.cutoff, color = Q.cutoff) +
#   geom_text(aes(label = round(value, 4))) +
#   labs(x = '') +
#   ggsci::scale_color_npg() +
#   theme_bw() +
#   theme(text = element_text(size = 12),
#         # axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
#         panel.grid = element_blank(), legend.position = 'none')
# 
# plot_Qcut <- sum.tbl3 %>%
#   rename(pr.FDP = runSpecific.pr.FDP.combined,
#          pg.FDP = runSpecific.pg.FDP.combined) %>% 
#   pivot_longer(cols = c(pr.FDP, pg.FDP,
#                         `# precursors`, `# protein groups`)) %>% 
#   ggplot() +
#   facet_wrap(~name, nrow = 1, scales = 'free_x') +
#   aes(x = value, y = Q.cutoff, color = Q.cutoff) +
#   geom_boxplot(width = 0.6) +
#   # labs(x = str_glue("DIA runs\n(n={length(unique(dfbox$Run))})"), y = 'Run-specific FDP') +
#   # scale_fill_manual(values = target_decoy_color) +
#   ggsci::scale_color_npg() +
#   theme_classic() +
#   theme(text = element_text(size = 12))
# 
# ggsave('entrapment_analysis_Q.cutoff_selection.pdf',
#        ggarrange(plot_Qcut.global, plot_Qcut, nrow = 1, widths = c(1, 2)),
#        width = 12, height = 5)
# 
# 
# 
# list(Qcutoff.selection = Qcutoff.selection,
#      Qcutoff.run.specific = Qcutoff.run.specific) %>% 
#   rio::export('entrapment_analysis_Q.cutoff.xlsx') # Also named as 'sF2B_entrapment_analysis_Qselection.xlsx'
# 
# 
# 
### removed data -----
# data_info1 <- rio::import('input/20250811_PUH_sample_information_1900files_v11.xlsx')
# data_ea <- read_excel_allsheets('input/sF2B_entrapment_analysis_Qselection.xlsx')
# data_ea$Qcutoff.selection %<>% mutate(Q.cutoff = factor(Q.cutoff, sort(as.numeric(Q.cutoff))))
# data_ea$Qcutoff.run.specific %<>% mutate(Q.cutoff = factor(Q.cutoff, sort(unique(as.numeric(Q.cutoff)))))
# 
# data_removed <- plyr::ddply(data_ea$Qcutoff.run.specific, 'Q.cutoff', function(tb){
#   data.frame(FileName = data_info1$FileName %>% setdiff(tb$Run))
# }) %>% inner_join(data_info1)
# rio::export(data_info1, 'QValues_filter_data_removed.xlsx')
# 
# data_remain <- plyr::ddply(data_ea$Qcutoff.run.specific, 'Q.cutoff', function(tb){
#   data.frame(FileName = data_info1$FileName %>% intersect(tb$Run))
# }) %>% inner_join(data_info1)
# 
# data_removed %>% group_by(Q.cutoff) %>% dplyr::count()
# data_remain %>% group_by(Q.cutoff) %>% distinct(sample_type, tissue_name) %>% dplyr::count()
# data_remain %>% group_by(Q.cutoff) %>% distinct(anatomical_classification) %>% dplyr::count()
# View(rep %>% filter(Run %in% data_removed$FileName[data_removed$Q.cutoff == 4e-4])) # In the 120 files, min(PG.Q) > 4e-4
# 
# 
# tmp.rep <- rep %>% filter(Run %in% data_removed$FileName[data_removed$Q.cutoff == 4e-4])
# length(unique(tmp.rep$Run)) # 120
# 
# # tmp.rep.stat <- tmp.rep %>%
# #   select(Run, Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value) %>% 
# #   ggpubr::get_summary_stats(digits = 10)
# 
# tmp.rep.stat <- tmp.rep %>% 
#   select(Run, Q.Value, Global.Q.Value, PG.Q.Value, Global.PG.Q.Value) %>%
#   pivot_longer(cols = -Run) %>% 
#   group_by(Run, name) %>% 
#   summarise(
#     across(everything(), list(
#       mean = ~mean(.x, na.rm = TRUE),
#       sd = ~sd(.x, na.rm = TRUE),
#       median = ~median(.x, na.rm = TRUE),
#       min = ~min(.x, na.rm = TRUE),
#       max = ~max(.x, na.rm = TRUE),
#       n = ~sum(!is.na(.x))
#     ))
#   ) %>% ungroup()
# 
# tmp.rep.stat %>% group_by(Run, value_min)
# 
# tmp.rep %>%
#   filter(Proteotypic == 1, Precursor.Quantity != 0,
#          PG.Q.Value < 5e-3, Global.PG.Q.Value < 4e-4,
#          Q.Value < 4e-4, Global.Q.Value < 4e-4) %>% 
#   distinct(Run) %>% pull() %>% length() # 48
# 
# a <- rep %>% filter(Proteotypic != 1)
# aa <- unique(a$Protein.Group)
# table(aa %in% lib_tphp$ProteinId) # ALL not in TPHP library
# table(str_detect(aa, ';')) # ALL containing ;
# 
# 
# b <- rep %>% filter(Proteotypic == 1)
# bb <- unique(b$Protein.Group)
# # table(bb %in% lib_tphp$ProteinId)
# table(str_detect(bb, ';')) # ALL not containing ;
# 
# bb %>% intersect(real.proteotipic) %>% length() # 14938
# 
# 
# tmp.rep %>%
#   filter(Proteotypic == 1, Precursor.Quantity != 0,
#          !(Protein.Group %in% false.proteotipic),
#          PG.Q.Value < 0.01, Global.PG.Q.Value < 4e-4,
#          Q.Value < 4e-4, Global.Q.Value < 4e-4) %>% 
#   distinct(Run) %>% pull() %>% length() # 120
# 

## benchmark2 -----
rep<-read_parquet(r"(..\..\results\rlt_append_RCA_G3\report.parquet)", 
                  col_select = c("Precursor.Id", "Protein.Group", "Genes", "Run","Proteotypic","Protein.Group","Q.Value","Global.Q.Value",
                                 "PG.Q.Value","Global.PG.Q.Value","Precursor.Quantity","Lib.Q.Value","Lib.PG.Q.Value","PG.MaxLFQ"
                  )
)
dim(rep)#[1]  260848055        13

thresholds <- c(1e-2, 5e-3, 1e-3, 5e-4, 4e-4, 3e-4, 2e-4, 1e-4)
# reports_test_thresholds1

reports_test_thresholds <- lapply(thresholds, function(thre){
  cat(thre, '...\n')
  tphp <- rep %>%
    filter(Proteotypic == 1, Precursor.Quantity != 0,
           !(Protein.Group %in% false.proteotipic),
           PG.Q.Value < 0.01, Global.PG.Q.Value < thre,
           Q.Value < thre, Global.Q.Value < thre) %>% 
    mutate(Is.pr.Target = Precursor.Id %in% X1$Precursor.Id,
           Is.pg.Target = Protein.Group %in% unique(lib_tphp$ProteinId))
  ## global pr FDP
  tmp <- tphp %>% distinct(Precursor.Id, Is.pr.Target) %>% dplyr::count(Is.pr.Target) %>% pull(n)
  global.pr.FDP.combined <- fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pr)
  pr.N.τ <- tmp[2]
  pr.N.ε <- tmp[1]
  
  ## global pg FDP
  tmp <- tphp %>% distinct(Protein.Group, Is.pg.Target) %>% dplyr::count(Is.pg.Target) %>% pull(n)
  global.pg.FDP.combined <- fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pg)
  pg.N.τ <- tmp[2]
  pg.N.ε <- tmp[1]
  
  ## specific pr FDP
  runSpecific.pr.FDP.combined <- tphp %>% plyr::ddply('Run', function(dfsub){
    tmp <- dfsub %>% distinct(Protein.Group, Is.pr.Target) %>% dplyr::count(Is.pr.Target) %>% pull(n)
    if(length(tmp) == 1) tmp <- c(0, tmp)
    data.frame(
      Run = dfsub$Run[1],
      runSpecific.pr.FDP.combined = fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pr)
    )
  }) %>% arrange(desc(runSpecific.pr.FDP.combined))
  
  ## specific pg FDP
  runSpecific.pg.FDP.combined <- tphp %>% plyr::ddply('Run', function(dfsub){
    tmp <- dfsub %>% distinct(Protein.Group, Is.pg.Target) %>% dplyr::count(Is.pg.Target) %>% pull(n)
    if(length(tmp) == 1) tmp <- c(0, tmp)
    data.frame(
      Run = dfsub$Run[1],
      runSpecific.pg.FDP.combined = fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pg)
    )
  }) %>% arrange(desc(runSpecific.pg.FDP.combined))
  
  # identity
  df_identity <- tphp %>% plyr::ddply('Run', function(dfsub){
    npg <- dfsub %>%
      distinct(Protein.Group, Is.pg.Target) %>%
      nrow()
    npr <- dfsub %>%
      distinct(Precursor.Id, Is.pr.Target) %>%
      nrow()
    data.frame(`# precursors` = npr, `# protein groups` = npg, check.names = F)
  })
  
  tbl2 <- data.frame(
    global.pr.FDP.combined = global.pr.FDP.combined,
    global.pg.FDP.combined = global.pg.FDP.combined
  )
  
  dfbox <- rbind(runSpecific.pg.FDP.combined %>% setNames(c('Run', 'FDP')) %>% mutate(Type = 'pg'),
                 runSpecific.pr.FDP.combined %>% setNames(c('Run', 'FDP')) %>% mutate(Type = 'pr'))
  
  tbl3 <- df_identity %>%
    inner_join(runSpecific.pr.FDP.combined) %>% 
    inner_join(runSpecific.pg.FDP.combined)
  
  global.identity <- list(
    pg.N.ε = pg.N.ε,
    pr.N.τ = pr.N.τ,
    pr.N.ε = pr.N.ε,
    pg.N.τ = pg.N.τ
  )
  
  ret <- list(global.FDP = tbl2,
              run.specific.FDP = dfbox,
              global.identity = global.identity,
              identity = tbl3)
  return(ret)
})
names(reports_test_thresholds) <- thresholds

sum.tbl2 <- plyr::ldply(reports_test_thresholds, function(ret){
  ret$global.FDP
}, .id = 'Q.cutoff')

sum.box <- plyr::ldply(reports_test_thresholds, function(ret){
  ret$run.specific.FDP
}, .id = 'Q.cutoff')

sum.identity.glb <- plyr::ldply(reports_test_thresholds, function(ret){
  as.data.frame(ret$global.identity)
}, .id = 'Q.cutoff')

sum.tbl3 <- plyr::ldply(reports_test_thresholds, function(ret){
  ret$identity
}, .id = 'Q.cutoff')

# ggplot(sum.tbl3) +
#   aes(x = `# protein groups`, y = Q.cutoff, color = Q.cutoff) +
#   geom_boxplot(width = 0.6) +
#   # labs(x = str_glue("DIA runs\n(n={length(unique(dfbox$Run))})"), y = 'Run-specific FDP') +
#   # scale_fill_manual(values = target_decoy_color) +
#   ggsci::scale_color_npg() +
#   theme_classic() +
#   theme(axis.text.x = element_blank()) +
#   theme(text = element_text(size = 12))
# 
# ggplot(sum.box) +
#   aes(x = FDP, y = Q.cutoff, color = Q.cutoff) +
#   geom_boxplot()

Qcutoff.selection <- sum.identity.glb %>% inner_join(sum.tbl2) %>% 
  relocate(pg.N.ε, .after = pg.N.τ) %>% 
  rename(global.pr.FDP = global.pr.FDP.combined,
         global.pg.FDP = global.pg.FDP.combined)
Qcutoff.run.specific <- sum.tbl3

# pheatmap(apply(Qcutoff.selection, 1:2, as.numeric),
#          display_numbers = T, number_format = '%.2e', color = NA, number_color = 'black',
#          cluster_rows = F, cluster_cols = F)

### XXX old figure XXX -----
# plot_Qcut.global <- Qcutoff.selection %>%
#   rename(`pr.FDP\n(global)` = global.pr.FDP,
#          `pg.FDP\n(global)` = global.pg.FDP) %>% 
#   pivot_longer(cols = -Q.cutoff) %>% 
#   ggplot() +
#   aes(x = name, y = Q.cutoff, color = Q.cutoff) +
#   geom_text(aes(label = round(value, 4))) +
#   labs(x = '') +
#   ggsci::scale_color_npg() +
#   theme_bw() +
#   theme(text = element_text(size = 12),
#         # axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
#         panel.grid = element_blank(), legend.position = 'none')
# 
# plot_Qcut <- sum.tbl3 %>%
#   rename(pr.FDP = runSpecific.pr.FDP.combined,
#          pg.FDP = runSpecific.pg.FDP.combined) %>% 
#   pivot_longer(cols = c(pr.FDP, pg.FDP,
#                         `# precursors`, `# protein groups`)) %>% 
#   ggplot() +
#   facet_wrap(~name, nrow = 1, scales = 'free_x') +
#   aes(x = value, y = Q.cutoff, color = Q.cutoff) +
#   geom_boxplot(width = 0.6) +
#   # labs(x = str_glue("DIA runs\n(n={length(unique(dfbox$Run))})"), y = 'Run-specific FDP') +
#   # scale_fill_manual(values = target_decoy_color) +
#   ggsci::scale_color_npg() +
#   theme_classic() +
#   theme(text = element_text(size = 12))
# 
# ggsave('entrapment_analysis_Q.cutoff_selection.pdf',
#        ggarrange(plot_Qcut.global, plot_Qcut, nrow = 1, widths = c(1, 2)),
#        width = 12, height = 5)

### new figure ----
## Left panel: plot_Qcut.global, cowplot style & English annotations
plot_Qcut.global <-
  Qcutoff.selection %>%
  rename(`pr.FDP` = global.pr.FDP,
         `pg.FDP` = global.pg.FDP) %>%
  pivot_longer(cols = -Q.cutoff) %>%
  mutate(name = factor(name, levels = c('pr.FDP', 'pg.FDP', 'pr.N.τ', 'pr.N.ε', 'pg.N.τ', 'pg.N.ε'))) %>%
  ggplot(aes(x = name, y = Q.cutoff, color = Q.cutoff)) +
  geom_text(aes(label = round(value, 4))) +
  labs(title = NULL, x = NULL, y = 'Q cutoff') +
  ggsci::scale_color_npg(name = 'Q cutoff') +
  cowplot::theme_cowplot(font_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = 'none'
  )

# ---- Right panel: plot_Qcut
Qcutoff.run.specific.long <- Qcutoff.run.specific %>%
  rename(pr.FDP = runSpecific.pr.FDP.combined,
         pg.FDP = runSpecific.pg.FDP.combined) %>%
  pivot_longer(cols = c(pr.FDP, pg.FDP, `# precursors`, `# protein groups`),
               names_to = 'metric', values_to = 'value')

metrics_order <- c('pr.FDP', 'pg.FDP', '# precursors', '# protein groups')
make_plot <- function(d, useBoxplot = TRUE, show_legend = FALSE, show_ytitle = FALSE) {
  # Precompute scalar parameters to avoid using ifelse() inside theme()/labs()
  y_title <- if (show_ytitle) "Q cutoff" else NULL
  y_title_theme <- if (show_ytitle) element_text() else element_text(color = NA)
  y_text_theme  <- if (show_ytitle) element_text() else element_blank()
  leg_pos <- if (show_legend) "right" else "none"
  the_title <- unique(d$metric)[1]  # Ensure scalar string for plot title
  geom <- if (useBoxplot) geom_boxplot(width = 0.6, outlier.stroke = 0.4) else geom_violin(aes(fill = Q.cutoff), width = 0.8)
  
  ggplot(d, aes(x = value, y = Q.cutoff, color = Q.cutoff)) +
    geom +
    ggsci::scale_color_npg(name = "Q cutoff") +
    ggsci::scale_fill_npg(name = "Q cutoff") +
    labs(title = NULL, x = the_title, y = y_title) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = leg_pos,
      axis.title.y    = y_title_theme,
      axis.text.y     = y_text_theme,
      plot.margin  = margin(t = 5.5, r = 2, b = 5.5, l = 2) # Reduce left and right plot margins (units in pt)
    )
}

plots_all <- list(
  make_plot(filter(Qcutoff.run.specific.long, metric == 'pr.FDP'), show_ytitle = T, useBoxplot = F),
  make_plot(filter(Qcutoff.run.specific.long, metric == 'pg.FDP'), useBoxplot = F),
  make_plot(filter(Qcutoff.run.specific.long, metric == '# precursors')),
  make_plot(filter(Qcutoff.run.specific.long, metric == '# protein groups'))
)

plot_with_leg <- make_plot(filter(Qcutoff.run.specific.long, metric == metrics_order[1]), show_legend = TRUE)
legend_g <- cowplot::get_legend(plot_with_leg)

row_plots <- cowplot::plot_grid(
  plotlist = plots_all, nrow = 1, # align = "hv",
  align = "h", axis = "l",
  rel_widths = c(5, 5, 5, 5)
)
plot_Qcut <- cowplot::plot_grid(row_plots, legend_g, nrow = 1, rel_widths = c(1, 0.1))
plot_Qcut


## ---- Compose both with cowplot
# Define a small helper function: generate a centered, bold top title
make_title <- function(txt, size = 12, vpad = 0.08) {
  # vpad is only used to adjust the relative height proportion in rel_heights (optional)
  cowplot::ggdraw() + 
    cowplot::draw_label(
      txt, x = 0.5, y = 1,
      hjust = 0.5, vjust = 1,
      fontface = "bold", size = size
    )
}

left_block  <- cowplot::plot_grid(
  make_title("Global"),         # Left-side title (custom text)
  plot_Qcut.global,             # Your left-side plot
  ncol = 1, rel_heights = c(0.05, 1)  # 0.10 can be tuned to adjust title height proportion
)

right_block <- cowplot::plot_grid(
  make_title("Run-specific"),        # Right-side title (custom text)
  plot_Qcut,                        # Your right-side plot
  ncol = 1, rel_heights = c(0.05, 1)
)

combined_cow <-
  cowplot::plot_grid(
    left_block,                     # Left module with title
    right_block,                    # Right module with title
    nrow = 1,
    align = 'v',           # vertical alignment of panels
    axis = 'lr',           # keep left/right y-axes aligned
    rel_widths = c(1, 2)
  )

ggsave('output/entrapment_analysis_Q.cutoff_selection.pdf',
       combined_cow,
       width = 12, height = 5)


### table -------
list(
  Qcutoff.selection = Qcutoff.selection %>% inner_join(
    data.frame(plyr::daply(Qcutoff.run.specific, 'Q.cutoff', nrow)) %>%
      setNames('n.files') %>% rownames_to_column('Q.cutoff'),
    .
  ),
  Qcutoff.run.specific = Qcutoff.run.specific %>% rename(FileName = Run)# %>% full_join(data_info1)
) %>% 
  rio::export('output/entrapment_analysis_Q.cutoff.xlsx') # Also named as 'sF2B_entrapment_analysis_Qselection.xlsx'



# EOF -------------
cat(thre, '...\n')
tphp <- rep %>%
  filter(Proteotypic == 1, Precursor.Quantity != 0,
         !(Protein.Group %in% false.proteotipic),
         PG.Q.Value < 0.01, Global.PG.Q.Value < 4e-4,
         Q.Value < 5e-4, Global.Q.Value < 4e-4) %>% 
  mutate(Is.pr.Target = Precursor.Id %in% X1$Precursor.Id,
         Is.pg.Target = Protein.Group %in% unique(lib_tphp$ProteinId))
## global pr FDP
tmp <- tphp %>% distinct(Precursor.Id, Is.pr.Target) %>% dplyr::count(Is.pr.Target) %>% pull(n)
global.pr.FDP.combined <- fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pr)
pr.N.τ <- tmp[2]
pr.N.ε <- tmp[1]

## global pg FDP
tmp <- tphp %>% distinct(Protein.Group, Is.pg.Target) %>% dplyr::count(Is.pg.Target) %>% pull(n)
global.pg.FDP.combined <- fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pg)
pg.N.τ <- tmp[2]
pg.N.ε <- tmp[1]

## specifi pr FDP
runSpecific.pr.FDP.combined <- tphp %>% plyr::ddply('Run', function(dfsub){
  tmp <- dfsub %>% distinct(Protein.Group, Is.pr.Target) %>% dplyr::count(Is.pr.Target) %>% pull(n)
  if(length(tmp) == 1) tmp <- c(0, tmp)
  data.frame(
    Run = dfsub$Run[1],
    runSpecific.pr.FDP.combined = fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pr)
  )
}) %>% arrange(desc(runSpecific.pr.FDP.combined))

## specific pg FDP
runSpecific.pg.FDP.combined <- tphp %>% plyr::ddply('Run', function(dfsub){
  tmp <- dfsub %>% distinct(Protein.Group, Is.pg.Target) %>% dplyr::count(Is.pg.Target) %>% pull(n)
  if(length(tmp) == 1) tmp <- c(0, tmp)
  data.frame(
    Run = dfsub$Run[1],
    runSpecific.pg.FDP.combined = fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pg)
  )
}) %>% arrange(desc(runSpecific.pg.FDP.combined))

# identity
df_identity <- tphp %>% plyr::ddply('Run', function(dfsub){
  npg <- dfsub %>%
    distinct(Protein.Group, Is.pg.Target) %>%
    nrow()
  npr <- dfsub %>%
    distinct(Precursor.Id, Is.pr.Target) %>%
    nrow()
  data.frame(`# precursors` = npr, `# protein groups` = npg, check.names = F)
})

tbl2 <- data.frame(
  global.pr.FDP.combined = global.pr.FDP.combined,
  global.pg.FDP.combined = global.pg.FDP.combined
)

dfbox <- rbind(runSpecific.pg.FDP.combined %>% setNames(c('Run', 'FDP')) %>% mutate(Type = 'pg'),
               runSpecific.pr.FDP.combined %>% setNames(c('Run', 'FDP')) %>% mutate(Type = 'pr'))

tbl3 <- df_identity %>%
  inner_join(runSpecific.pr.FDP.combined) %>% 
  inner_join(runSpecific.pg.FDP.combined)

global.identity <- list(
  pr.N.τ = pr.N.τ,
  pr.N.ε = pr.N.ε,
  pg.N.τ = pg.N.τ,
  pg.N.ε = pg.N.ε
)

ret <- list(global.FDP = tbl2,
            run.specific.FDP = dfbox,
            global.identity = global.identity,
            identity = tbl3)



# # select 4e-4 as Cutoff of all Q.Values except PG.Q < 0.01 -----
# ## Q.Values filter ---------
# thre <- 4e-4
# tphp <- rep %>%
#   filter(Proteotypic == 1, Precursor.Quantity != 0,
#          !(Protein.Group %in% false.proteotipic),
#          PG.Q.Value < 0.01, Global.PG.Q.Value < thre,
#          Q.Value < thre, Global.Q.Value < thre) %>% 
#   mutate(Is.pr.Target = Precursor.Id %in% X1$Precursor.Id,
#          Is.pg.Target = Protein.Group %in% unique(lib_tphp$ProteinId))
# 
# write_parquet(tphp, r"(..\..\results\rlt_append_RCA_G3\report.pro_matrix_full_0.0004Q_PGQ0.01.parquet)")
# 
# ## save original matrix ------
# tphp_final<-tphp%>% 
#   select(Protein.Group,Run,PG.MaxLFQ) %>% as.data.frame() %>% 
#   # pivot_wider(names_from = Run, values_from = PG.MaxLFQ) 
#   group_by(Run,Protein.Group) %>%
#   summarise(mean = mean(PG.MaxLFQ), 
#             count = n(),
#             min=min(PG.MaxLFQ), 
#             max=max(PG.MaxLFQ),   
#             .groups = 'drop')
# head(tphp_final)
# dim(tphp_final)#[1] 11097294       6
# identical(tphp_final$max, tphp_final$min)#True
# sum(tphp_final$count==1)== nrow(tphp_final)#False
# print(hist(tphp_final$count, breaks = 100) )
# min(tphp_final$count) #1
# max(tphp_final$count) #3299
# # mapping<-
# pm0<-tphp_final %>% 
#   select(Protein.Group,Run,min) %>% 
#   pivot_wider(names_from = Run, values_from = min,values_fill=NA) %>% 
#   as.data.frame()
# dim(pm0)#[1] 13391 1958
# write.table(pm0, file = r"(..\..\results\rlt_append_RCA_G3\report.pro_matrix_full_0.0004Q_PGQ0.01_wide.tsv)", sep = "\t",row.names = F,quote = F)
# 
# pm0_human<-pm0 %>% filter(Protein.Group %in% lib_tphp$ProteinId) %>% filter(!str_detect(Protein.Group, 'iRT'))
# dim(pm0_human) # 13245 1985
# write.table(pm0_human, file = r"(..\..\results\rlt_append_RCA_G3\report.pro_matrix_full_0.0004Q_PGQ0.01_wide_human.tsv)", sep = "\t",row.names = F,quote = F)
# 
# 
# ## mapping sample info ------
# # read protein matrix after FDP-control from entrap analysis 
# df <- pm0_human %>% column_to_rownames('Protein.Group')
# head(colnames(df))  #file_name
# head(rownames(df))    #protein uniprot ID
# # "DIA-NN: Sometimes DIA-NN will report a zero as the best estimate for a precursor or protein quantity. Such zero quantities are omitted from protein/gene matrices. "
# df[df == 0]<-NA 
# df %>% min(na.rm = T) # 18.4645 # 2.840796  #entrapVersion: 12.88134    #2022 version: 37.4647
# df1<-df<- df %>%t()%>% as.data.frame()%>%
#   tibble::rownames_to_column(var = "FileName")
# dim(df1)  # 1957 13246   # mild filter version: 1957 14063
# 
# labels <- readxl::read_excel('20250811_PUH_sample_information_1900files_v11.xlsx')
# df_info <- inner_join(labels, df1, by = 'FileName')
# 
# pm0_human_mapped<-df_info %>% column_to_rownames('FileName') %>% 
#   select(-(1:new)) %>% as.matrix()
# is.numeric(pm0_human_mapped) # TRUE
# dim(pm0_human_mapped) # 1900 13245
# 
# pm_tphp1 <- pm0_human_mapped[df_info$FileName[df_info$sample_type != 'p'], ]
# dim(pm_tphp1)
# write.csv(pm_tphp1, 'mapped_pg_matrix_1780_13245.csv', row.names = T)                         
# 
# pm_tphp2 <- pm0_human_mapped[df_info$FileName[df_info$sample_type == 'p'], ] # pooling
# dim(pm_tphp2 ) #120 13245
# head(pm_tphp2[,1:10])
# write.csv(pm_tphp2, 'pool_pg_matrix_120_13245.csv', row.names = T)

# select 1e-3 as Cutoff of all Q.Values except PG.Q < 0.01 -----
## Q.Values filter ---------
thre <- 1e-3
tphp <- rep %>%
  filter(Proteotypic == 1, Precursor.Quantity != 0,
         !(Protein.Group %in% false.proteotipic),
         PG.Q.Value < 0.01, Global.PG.Q.Value < thre,
         Q.Value < thre, Global.Q.Value < thre) %>% 
  mutate(Is.pr.Target = Precursor.Id %in% X1$Precursor.Id,
         Is.pg.Target = Protein.Group %in% unique(lib_tphp$ProteinId))
dim(tphp) # 212616368        15
write_parquet(tphp, r"(..\..\results\rlt_append_RCA_G3\report.pro_matrix_full_0.001Q_PGQ0.01.parquet)")

## save original matrix ------
tphp_final<-tphp%>% 
  select(Protein.Group,Run,PG.MaxLFQ) %>% as.data.frame() %>% 
  # pivot_wider(names_from = Run, values_from = PG.MaxLFQ) 
  group_by(Run,Protein.Group) %>%
  summarise(mean = mean(PG.MaxLFQ), 
            count = n(),
            min=min(PG.MaxLFQ), 
            max=max(PG.MaxLFQ),   
            .groups = 'drop')
head(tphp_final)
dim(tphp_final)#[1] 20306085        6
identical(tphp_final$max, tphp_final$min)#True
sum(tphp_final$count==1)== nrow(tphp_final)#False
print(hist(tphp_final$count, breaks = 100) )
min(tphp_final$count) #1
max(tphp_final$count) #3509
# mapping<-
pm0<-tphp_final %>% 
  select(Protein.Group,Run,min) %>% 
  pivot_wider(names_from = Run, values_from = min,values_fill=NA) %>% 
  as.data.frame()
dim(pm0)#[1] 13840  3263
write.table(pm0, file = r"(..\..\results\rlt_append_RCA_G3\report.pro_matrix_full_0.001Q_PGQ0.01_wide.tsv)", sep = "\t",row.names = F,quote = F)

pm0_human<-pm0 %>% filter(Protein.Group %in% lib_tphp$ProteinId) %>% filter(!str_detect(Protein.Group, 'iRT'))
dim(pm0_human) # 13609  3263
write.table(pm0_human, file = r"(..\..\results\rlt_append_RCA_G3\report.pro_matrix_full_0.001Q_PGQ0.01_wide_human.tsv)", sep = "\t",row.names = F,quote = F)


## mapping sample info ------
# read protein matrix after FDP-control from entrap analysis 
df <- pm0_human %>% column_to_rownames('Protein.Group')
head(colnames(df))  #file_name
head(rownames(df))    #protein uniprot ID
# "DIA-NN: Sometimes DIA-NN will report a zero as the best estimate for a precursor or protein quantity. Such zero quantities are omitted from protein/gene matrices. "
df[df == 0]<-NA 
df %>% min(na.rm = T) # 12.77158 # 16.40772 # 18.4645 # 2.840796  #entrapVersion: 12.88134    #2022 version: 37.4647
df1<-df<- df %>%t()%>% as.data.frame()%>%
  tibble::rownames_to_column(var = "FileName")
dim(df1) # 3262 13610  # 1957 13545   # mild filter version: 1957 14063

# labels <- readxl::read_excel('20250811_PUH_sample_information_1900files_v11.xlsx')
# labels <- rio::import('20251009_PUH_sample_information_3112files.xlsx')
labels <- readRDS('20251009_PUH_sample_information_3005files_qc.RDS')
df_info <- inner_join(labels, df1, by = 'FileName')

pm0_human_mapped<-df_info %>% column_to_rownames('FileName') %>% 
  select(-(1:end)) %>% as.matrix()
is.numeric(pm0_human_mapped) # TRUE
# pm0_human_mapped <- pm0_human_mapped[, colSums(!is.na(pm0_human_mapped)) > 0]
dim(pm0_human_mapped) # 3005 13609 # 3997 13609 # 3141 13613  # 1900 13544


pm_tphp1 <- pm0_human_mapped[df_info$FileName[df_info$sample_type != 'p'], ]
dim(pm_tphp1) # 2856 13609
write.csv(pm_tphp1, 'mapped_pg_matrix_2856_13609.csv', row.names = T)

pm_tphp2 <- pm0_human_mapped[df_info$FileName[df_info$sample_type == 'p'], ] # pooling
dim(pm_tphp2) #149 13609
head(pm_tphp2[,1:10])
write.csv(pm_tphp2, 'pool_pg_matrix_149_13609.csv', row.names = T)

# 
# pm_tphp <- rbind(pm_tphp1, pm_tphp2)
# dim(pm_tphp) # 3109 13609
# length(intersect(rownames(pm_tphp), df_info$FileName)) # 3109

# save protien info table -----
# tphp1 <- arrow::read_parquet(r"(..\..\results\rlt_combine_entraplib2_merge_G3\report.parquet)")
# tphp2 <- arrow::read_parquet(r"(..\..\results\report-CAF20250725_b12_13\report.parquet)")
# pt1 <- tphp1 %>% distinct(Protein.Group, Genes)
# pt2 <- tphp2 %>% distinct(Protein.Group, Genes)
# prot_info <- pt1 %>% rbind(pt2) %>% distinct()
prot_info <- rep %>% distinct(Protein.Group, Genes)
rio::export(prot_info, 'output/protein_info_from_reports_V1009.csv')



# Post-filter FDP estimate -------
tphp_post <- tphp %>% filter(Run %in% labels$FileName)
dim(tphp_post) # 198558513        15

## global pr FDP
tmp <- tphp_post %>% distinct(Precursor.Id, Is.pr.Target) %>% dplyr::count(Is.pr.Target) %>% pull(n)
global.pr.FDP.combined <- fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pr)

## global pg FDP
tmp <- tphp_post %>% distinct(Protein.Group, Is.pg.Target) %>% dplyr::count(Is.pg.Target) %>% pull(n)
global.pg.FDP.combined <- fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pg)

## specific pr FDP
runSpecific.pr.FDP.combined <- tphp_post %>% plyr::ddply('Run', function(dfsub){
  tmp <- dfsub %>% distinct(Protein.Group, Is.pr.Target) %>% dplyr::count(Is.pr.Target) %>% pull(n)
  if(length(tmp) == 1) tmp <- c(0, tmp)
  data.frame(
    Run = dfsub$Run[1],
    runSpecific.pr.FDP.combined = fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pr)
  )
}) %>% arrange(desc(runSpecific.pr.FDP.combined))

## specific pg FDP
runSpecific.pg.FDP.combined <- tphp_post %>% plyr::ddply('Run', function(dfsub){
  tmp <- dfsub %>% distinct(Protein.Group, Is.pg.Target) %>% dplyr::count(Is.pg.Target) %>% pull(n)
  if(length(tmp) == 1) tmp <- c(0, tmp)
  data.frame(
    Run = dfsub$Run[1],
    runSpecific.pg.FDP.combined = fdp_est(N.τ = tmp[2], N.ε = tmp[1], r = r_pg)
  )
}) %>% arrange(desc(runSpecific.pg.FDP.combined))

# identity
df_identity <- tphp_post %>% plyr::ddply('Run', function(dfsub){
  npg <- dfsub %>%
    distinct(Protein.Group, Is.pg.Target) %>%
    nrow()
  npr <- dfsub %>%
    distinct(Precursor.Id, Is.pr.Target) %>%
    nrow()
  data.frame(`# precursors` = npr, `# protein groups` = npg, check.names = F)
})

tbl2 <- data.frame(
  global.pr.FDP.combined = global.pr.FDP.combined,
  global.pg.FDP.combined = global.pg.FDP.combined
)

dfbox <- rbind(runSpecific.pg.FDP.combined %>% setNames(c('Run', 'FDP')) %>% mutate(Type = 'pg'),
               runSpecific.pr.FDP.combined %>% setNames(c('Run', 'FDP')) %>% mutate(Type = 'pr'))
p_runFDP <- ggplot(dfbox) +
  aes(x = Type, y = FDP, fill = Type) +
  geom_boxplot(width = 0.6, position = position_dodge()) +
  # geom_text(data = dfbox %>% filter(FDP > 0.01),
  #           aes(label = round(FDP, 4)), hjust = -0.5) +
  labs(x = str_glue("DIA runs\n(n={length(unique(dfbox$Run))})"), y = 'Run-specific FDP') +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0, 0.05, 0)) +
  # scale_fill_manual(values = target_decoy_color) +
  ggsci::scale_fill_npg() +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  theme(text = element_text(size = 12))
ggsave('entrapment_estimate_combined_runspecific_001Q_01PGQ_3005Run.pdf', p_runFDP, width = 5, height = 4)


bad_runs <- runSpecific.pr.FDP.combined %>%
  filter(runSpecific.pr.FDP.combined > 0.01) %>% 
  pull(Run)
bad_runs_info <- df_info %>%
  filter(str_detect(FileName, str_c(bad_runs, collapse = '|')))

bad_runs_info_annotate <- bad_runs_info %>%
  distinct(sample_type, tissue_name) %>%
  semi_join(df_info, .) %>%
  inner_join(df_identity %>% rename(FileName = Run)) %>%
  inner_join(runSpecific.pr.FDP.combined %>% rename(FileName = Run)) %>%
  inner_join(runSpecific.pg.FDP.combined %>% rename(FileName = Run)) %>%
  arrange(sample_type, tissue_name, runSpecific.pr.FDP.combined, runSpecific.pg.FDP.combined)


tbl3 <- df_identity %>%
  inner_join(runSpecific.pr.FDP.combined) %>% 
  inner_join(runSpecific.pg.FDP.combined)

# counter
tphp_post_human <- tphp_post %>%
  filter(Protein.Group %in% lib_tphp$ProteinId) %>%
  filter(!str_detect(Protein.Group, 'iRT'))
pro.num <- tphp_post_human %>% distinct(Run, Protein.Group) %>% 
  count(Run, name = '# proteins')
pr.num <- tphp_post_human %>% distinct(Run, Precursor.Id) %>% 
  count(Run, name = '# precursors')
pep.num <- tphp_post_human %>% distinct(Run, Precursor.Id) %>% 
  mutate(Stripped.Sequence = str_remove_all(Precursor.Id, '\\d+$|\\(.+?\\)')) %>% 
  distinct(Run, Stripped.Sequence) %>% 
  count(Run, name = '# peptides')
df_identity <- pro.num %>% 
  inner_join(pep.num) %>% 
  inner_join(pr.num)



rio::export(df_identity, 'identity_3005Run_V1009.xlsx')


# OUTPUT -------
list(r.values = tbl1,
     global.FDP = tbl2,
     # run.specific.pr.FDP = runSpecific.pr.FDP.combined,
     # run.specific.pg.FDP = runSpecific.pg.FDP.combined,
     identity = tbl3,
     bad.runs = bad_runs_info,
     bad.run.groups = bad_runs_info_annotate
) %>% 
  rio::export('entrapment_lib_estimate_001Q_01PGQ_3005Run_source_data_20251009.xlsx')
