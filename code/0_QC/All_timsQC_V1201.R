rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = T)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggsci)
source('../source/source_code.R')

# 1.read data ----------
## 1.0 data reading -------
info5 <- rio::import('../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V5.xlsx')
meta <- info5 %>% filter(sample_type != 'p')
smp <- rio::import('../0_process_DIA-NN_data/input/mapped_pg_matrix_2856_13609.csv') %>%
  column_to_rownames('V1')

df_smp <- smp %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(meta, .)

## 1.1 remove low-identity (< Q1-1.5*IQR) -------
# est.identity <- df_smp %>%
#   plyr::ddply(c('sample_type', 'tissue_name'), function(dfsub){
#     cat(dfsub$sample_type[1], dfsub$tissue_name[1], '...\r')
#     X <- dfsub %>%
#       column_to_rownames('FileName') %>%
#       select(-(1:ncol(meta))) %>% t()
#     # X <- dat_raw[, dfsub$FileName, drop = FALSE]
#     ret <- data.frame(FileName = dfsub$FileName,
#                       `# proteins` = colSums(!is.na(X)),
#                       n_detailed_sample = nrow(dfsub),
#                       check.names = F)
#     ret$outlier.lower.ingroup <- get_outliers(ret$`# proteins`)[1]
#     return(ret)
#   })
# est.identity %<>%
#   group_by(sample_type, tissue_name) %>%
#   summarise(identity.mean = mean(`# proteins`),
#             identity.min = min(`# proteins`),
#             identity.max = max(`# proteins`),
#             identity.range.length = identity.max - identity.min,
#             .groups = 'drop') %>%
#   inner_join(est.identity, .) %>%
#   mutate(Group = str_c(sample_type, ' - ', tissue_name)) %>% 
#   arrange(identity.mean) %>%
#   mutate(Group = factor(Group, levels = unique(Group)),
#          Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup) %>% 
#   inner_join(meta %>% select(FileName, sample_id, Rep_type, tissue_name, sample_type, anatomical_classification))
# group_order <- est.identity %>% arrange(`# proteins`) %>% distinct(Group) %>% pull()
# est.identity %<>% mutate(Group = factor(Group, group_order))
# rio::export(est.identity, 'output/QC2856_identity_20251201_source.xlsx')

# plot_ident <- ggplot(est.identity) +
#   aes(x = `# proteins`, y = Group) +
#   geom_boxplot(data = est.identity %>% filter(identity.range.length <= 1000),
#                color = 'black', outlier.color = '#c23190', outlier.size = 3) +
#   geom_boxplot(data = est.identity %>% filter(identity.range.length > 1000),
#                color = 'red4', outlier.color = '#c23190', outlier.size = 3) +
#   geom_point(aes(color = tissue_name), size = 1.2) +
#   labs(x = "# proteins", y = "Sample ID", subtitle = "Protein identification") +
#   # ggsci::scale_color_aaas(name = 'Sample type') +
#   theme_bw() +
#   theme(text = element_text(size = 10, color = 'black'))
# ggsave('output/QC2856_identity_20251201.pdf', plot_ident, width = 40, height = 18)
est.identity0 <- plyr::ddply(df_smp, c("sample_type","class_abbr"), function(d){
  cat(d$sample_type[1], d$class_abbr[1], "...\n")
  
  X <- d %>%
    tibble::column_to_rownames("FileName") %>%
    dplyr::select(-(1:ncol(meta))) %>%
    t()
  
  nprot <- unname(colSums(!is.na(X))[d$FileName])
  
  ret <- data.frame(
    FileName=d$FileName, sample_type=d$sample_type, class_abbr=d$class_abbr,
    tissue_name=d$tissue_name, `# proteins`=nprot,
    n_major_sample=nrow(d),
    stringsAsFactors=FALSE, check.names=FALSE
  )
  
  # major-level stats (replicated to each file row in this class)
  ret <- ret %>% dplyr::mutate(
    outlier.lower.class = get_outliers(`# proteins`)[1],
    identity.mean.class = mean(`# proteins`),
    identity.min.class  = min(`# proteins`),
    identity.max.class  = max(`# proteins`),
    identity.range.length.class = identity.max.class - identity.min.class
  )
  
  # detailed-level stats (tissue within this class), then join by shared key without "by="
  ts <- plyr::ddply(ret, "tissue_name", function(t){
    data.frame(
      tissue_name=t$tissue_name[1],
      n_detailed_sample=nrow(t),
      outlier.lower.tissue=get_outliers(t$`# proteins`)[1],
      identity.mean.tissue=mean(t$`# proteins`),
      identity.min.tissue=min(t$`# proteins`),
      identity.max.tissue=max(t$`# proteins`),
      identity.range.length.tissue=max(t$`# proteins`) - min(t$`# proteins`),
      stringsAsFactors=FALSE, check.names=FALSE
    )
  })
  
  ret <- ret %>%
    dplyr::left_join(ts) %>%
    dplyr::mutate(
      Is.Lower.Ingroup.tissue = `# proteins` < outlier.lower.tissue,
      Is.Lower.Ingroup.class  = `# proteins` < outlier.lower.class
    )
  
  ret
})

est.identity <- est.identity0 %>%
  dplyr::mutate(
    Group_major    = stringr::str_c(sample_type," - ",class_abbr),
    Group_detailed = stringr::str_c(sample_type," - ",class_abbr," - ",tissue_name)
  ) %>%
  dplyr::arrange(identity.mean.tissue, identity.mean.class, `# proteins`) %>%
  dplyr::mutate(
    Group_major    = factor(Group_major,    levels=unique(Group_major)),
    Group_detailed = factor(Group_detailed, levels=unique(Group_detailed))
  ) %>%
  dplyr::left_join(meta %>% dplyr::select(FileName, sample_id, Rep_type, tissue_name, sample_type, class_abbr))

est.identity %<>% select(FileName, sample_id, Rep_type, Group_major, Is.Lower.Ingroup.class, Group_detailed, Is.Lower.Ingroup.tissue,
                         everything())
# one table only (each FileName row carries both tissue- and class-level estimates)
.write_excel(est.identity, "output/QC2856_identity_20251201_source.xlsx")


plot_ident <- ggplot(est.identity) +
  aes(x=`# proteins`, y=Group_major) +
  geom_boxplot(data=est.identity %>% dplyr::filter(identity.range.length.tissue <= 1000),
               color="black", outlier.shape=NA) +
  geom_boxplot(data=est.identity %>% dplyr::filter(identity.range.length.tissue > 1000),
               color="red4", outlier.shape=NA) +
  geom_point(data=est.identity %>% dplyr::filter(!Is.Lower.Ingroup.tissue),
             aes(color=tissue_name), size=1.2, alpha=0.9) +
  geom_point(data=est.identity %>% dplyr::filter(Is.Lower.Ingroup.tissue),
             color="#c23190", size=2.5) +
  geom_point(data=est.identity %>% dplyr::filter(Is.Lower.Ingroup.class),
             color="#c23190", size=5.0) +
  labs(x="# proteins", y="Sample", subtitle="Protein identification") +
  theme_bw() + theme(text=element_text(size=10, color="black"))

ggsave("output/QC2856_identity_20251201.pdf", plot_ident, width=40, height=18)


## 1.2 remove high-missing proteins (>50%), impute_group = sample_type + detailed_tissue ----
# keep proteins with present% >= 50% in at least one impute_group
mat.smp <- t(smp)
impute_group <- interaction(meta$sample_type, meta$tissue_name, drop = TRUE)
pres <- !is.na(mat.smp.keep)

idx_list <- split(seq_len(ncol(mat.smp.keep)), impute_group)
present_by_group <- lapply(idx_list, function(ii) rowMeans(pres[, ii, drop = FALSE]))
present_max <- Reduce(pmax, present_by_group)

keep_prot <- present_max >= 0.50
mat.smp.keep <- mat.smp[keep_prot, , drop = FALSE]

## 1.3 impute NA (global_min-based) ----
mat.smp.keep_log2 <- log2(mat.smp.keep)

min_log2 <- min(mat.smp.keep_log2, na.rm = TRUE)
n_na <- sum(is.na(mat.smp.keep_log2))

set.seed(2022)
noise_sd <- 0.001
jit <- rnorm(ceiling(n_na * 1.1), mean = 1, sd = noise_sd)
jit <- jit[jit > 0]
imp <- min_log2 + log2(0.5) + log2(jit[seq_len(n_na)])

mat.smp.keep_log2[is.na(mat.smp.keep_log2)] <- imp


## 1.4 batch effect correct (ComBat) ----
# protect tissue_name; mean.only = FALSE
mod <- model.matrix(~ tissue_name, data = meta)

mat.smp.keep_log2_combat <- sva::ComBat(
  dat = mat.smp.keep_log2,
  batch = meta$new,
  mod = mod,
  par.prior = TRUE,
  prior.plots = TRUE,
  mean.only = FALSE
)
dim(mat.smp.keep_log2_combat) # 13064 2856
saveRDS(mat.smp.keep_log2_combat, 'output/mapped_pg_matrix_2856_13064_combat.csv')


# 2.Pooling ------
## 2.0 data reading -------
info5 <- rio::import('../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V5.xlsx')
meta.pool <- info5 %>% filter(sample_type == 'p')
pool <- rio::import('../0_process_DIA-NN_data/input/pool_pg_matrix_149_13609.csv') %>%
  column_to_rownames('V1')

df_pool <- pool %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(meta.pool, .)

# mat.smp.correct <- rio::import('../2_batch_effect_evaluation_correction/output/batch_corrected_data_step2_limma3.csv') %>% 
#   column_to_rownames('V1')

## 2.1 data preprocess -----
pool1 <- rio::import('../1_NA_filter_imputation/output/data_matrix_filtered_05NA_2957_12797.csv') %>% 
  column_to_rownames('V1') %>% .[meta.pool$FileName, ]
mat.pool <- t(pool1) %>% removeRowsAllNa() %>% removeColsAllNa() %>% removeRowsNa(threshold = 0.5)
mat.pool.q <- preprocessCore::normalize.quantiles(mat.pool, copy = T, keep.names = T)
mat.pool.q.log <- log2(mat.pool.q + 1)  

sample_na<-apply(mat.pool.q.log, 2, function(x) sum(is.na(x))) 
which(sample_na==nrow(mat.pool.q.log))  #0
pro_na<-apply(mat.pool.q.log, 1, function(x) sum(is.na(x))) 
which(pro_na==ncol(mat.pool.q.log)) #0
min.pool <- min(mat.pool.q.log, na.rm = TRUE)
imputation <- log2(0.5) + min.pool

mat.pool.imp <- mat.pool.q.log; mat.pool.imp[is.na(mat.pool.imp)] <- imputation

# two-steps limma removeBatchEffect
# step1
step1_corrected <- limma::removeBatchEffect(
  x = mat.pool.imp,
  batch = meta.pool$new,
  design = model.matrix(~ instrument + anatomical_classification, 
                        data = meta.pool)
)
# step2
step2_corrected <- limma::removeBatchEffect(
  x = step1_corrected,
  batch = meta.pool$instrument,
  design = model.matrix(~ anatomical_classification, 
                        data = meta.pool)
)
mat.pool.correct <- step2_corrected

df_est_pool <- mat.pool.correct %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>%
  inner_join(meta.pool, .)

df_est_pool_filter <- df_est_pool %>% filter(anatomical_classification != 'nail')


# 3.sF2-3 QC -------
## 3.1 rep data ------
rep_check <- meta %>% filter(Is.Rep)
# rep_check %>% count(sample_id) %>% count(n)

mat.smp.correct <- rio::import('../2_batch_effect_evaluation_correction/output/batch_corrected_data_step2_limma3.csv') %>% 
  column_to_rownames('V1')
df_est_rep <- mat.smp.correct %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(rep_check, .)
df_est_rep %<>% count(sample_id) %>% filter(n > 1) %>% semi_join(df_est_rep, .)


rm.fid <- c()

rep.identity <- #df_est_rep %>% 
  mat.smp %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(rep_check, .) %>% 
  plyr::ddply('sample_id', function(dfsub){
  cat(dfsub$sample_id[1], '...\r')
  X <- dfsub %>%
    column_to_rownames('file_id') %>%
    select(-(1:tail(colnames(rep_check), 1))) %>% t()
  ret <- data.frame(file_id = dfsub$file_id,
                    `# proteins` = colSums(!is.na(X)),
                    check.names = F)
  ret$outlier.lower.ingroup <- get_outliers(ret$`# proteins`)[1]
  return(ret)
})
rep.identity %<>%
  group_by(sample_id) %>%
  summarise(identity.mean = mean(`# proteins`),
            identity.min = min(`# proteins`),
            identity.max = max(`# proteins`),
            identity.range.length = identity.max - identity.min) %>%
  inner_join(rep.identity, .) %>%
  arrange(identity.mean) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id)),
         Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup) %>% 
  inner_join(meta %>% select(file_id, Rep_type, class_abbr, sample_type))
rep.sid_order <- rep.identity %>% arrange(`# proteins`) %>% distinct(sample_id) %>% pull()
rep.identity %<>% 
  mutate(sample_id = factor(sample_id, rep.sid_order),
         Rep_type = ifelse(sample_type == 'p', 'p', Rep_type),
         Rep_type = ifelse(!is.na(Rep_type), Rep_type, 'Reference'),
         Rep_type = factor(Rep_type, c('Reference', 'b', 't', 'p')))

low.rep.ident <- rep.identity %>% filter(identity.range.length > 1000) %>%
  group_by(sample_id) %>% arrange(`# proteins`) %>% slice(1) %>% 
  pull(file_id)
rm.fid %<>% append(low.rep.ident)
rm.sid <- df_est_rep %>% filter(file_id %in% rm.fid) %>% distinct(sample_id) %>% pull() %>% setdiff(c('p'))


rep.pearson <- plyr::ddply(#df_est_rep %>%
  # filter(!(sample_id %in% rm.sid),
  #        !(file_id %in% rm.fid)),
  df_est_rep,
  'sample_id', function(dfsub){
    cat(dfsub$sample_id[1], '...\r')
    X <- dfsub %>%
      column_to_rownames('file_id') %>%
      select(-(1:tail(colnames(rep_check), 1))) %>% t()
    corX <- cor(X, use = 'pairwise.complete.obs', method = 'pearson')
    corX_pair <- corX %>% as.data.frame() %>%
      rownames_to_column('ID1') %>%
      pivot_longer(cols = -ID1, names_to = 'ID2') %>%
      unite('pair', ID1, ID2, sep = '-') %>%
      pull(pair) %>%
      matrix(nrow = nrow(corX), ncol = ncol(corX),
             byrow = T, dimnames = dimnames(corX))
    ret <- data.frame(ID.pair = corX_pair %>% .[upper.tri(., diag = F)],
                      pearson.r = corX %>% .[upper.tri(., diag = F)],
                      check.names = F)
    ret$outlier.lower.ingroup <- get_outliers(ret$pearson.r)[1]
    return(ret)
  })
rep.pearson %<>%
  group_by(sample_id) %>%
  summarise(pearson.r.mean = mean(pearson.r)) %>%
  inner_join(rep.pearson, .) %>%
  arrange(pearson.r.mean) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id)))

boxplot(rep.pearson$pearson.r)

# After checking, the Rb1_76_t1 should be removed
df_est_rep_filter <- df_est_rep %>% filter(!file_id %in% c('Rb1_76_t1'))
df_est_rep$FileName[df_est_rep$file_id == 'Rb1_76_t1']
# "N20250708yuel_TPHP_RCA_90mDIA_b1_76_Slot2-27_1_34640"

## 3.2 combine  ------
df_est <- plyr::rbind.fill(df_est_pool_filter, df_est_rep_filter)
# saveRDS(df_est, 'df_est.rds')
# df_est <- readRDS('df_est.rds')

# rep data sorting
rep1 <- df_est %>% count(sample_id) %>% filter(n == 1) %>% anti_join(df_est, .) # 658
rep1 %<>% mutate(Rep_type = ifelse(sample_type == 'p', 'p', Rep_type)) %>%
  mutate(Rep_type = ifelse(is.na(Rep_type), 'o', Rep_type))
dim(rep1) #  658 12831

# Check wrong labels: 'o' annotated as 't'
rep1_type_check <- rep1 %>%
  mutate(Rep_type_check = str_extract(file_id, '([bt])\\d$', group = 1),
         Rep_type_check = ifelse(sample_type == 'p', 'p', Rep_type_check),
         Rep_type_check = ifelse(is.na(Rep_type_check), 'o', Rep_type_check)) %>% 
  select(FileName, sample_id, file_id, Rep_type, Rep_type_check)
o.fid <- rep1_type_check %>% filter(Rep_type != Rep_type_check) %>% pull(file_id)

rep1 %<>%
  mutate(Rep_type = ifelse(file_id %in% o.fid, 'o', Rep_type)) %>%
  mutate(Type = c(o = 'Original', t = 'Technical', b = 'Biological', p = 'Pooling')[Rep_type], .before = new)

# Separate 't' and 'b'
rep1p <- rep1 %>% filter(Rep_type == 'p')
rep1t <- rep1 %>% filter(Rep_type %in% c('o', 't'))
rep1t <- rep1t %>% count(sample_id) %>% filter(n > 1) %>% semi_join(rep1t, .)
rep1b <- rep1 %>% filter(Rep_type %in% c('o', 'b'))
rep1b <- rep1b %>% count(sample_id) %>% filter(n > 1) %>% semi_join(rep1b, .)

# some sample_ids without 'o' (only t and b) -- transform t to o
rbind(rep1b, rep1t, rep1p) %>% distinct(file_id) %>% nrow() # 659
rbind(rep1b, rep1t, rep1p) %>% distinct(file_id) %>% anti_join(rep1, .) %>% select(FileName, sample_id:cancer_abbr, DateTime:new) %>% View() # NULL
rbind(rep1b, rep1t, rep1p) %>% distinct(file_id) %>% nrow() %>% identical(nrow(rep1)) # TRUE (659)

# envelope data_list (actually a data.frame)
rep1_list <- list(Pooling = rep1p, Technical = rep1t, Biological = rep1b) %>% plyr::ldply()

# annotations (n/pairs of samples)
n_pool <- rep1_list %>% filter(Rep_type == 'p') %>% distinct(file_id) %>% nrow()
n_trep_pairs <- rep1_list %>% filter(Rep_type == 't') %>% distinct(sample_id) %>% nrow()
n_brep_pairs <- rep1_list %>% filter(Rep_type == 'b') %>% distinct(sample_id) %>% nrow()

# saveRDS(rep1_list, 'rep1_list.rds')
# rep1_list <- load('rep1_list.rds')


## 3.3 sF2, all data ------
### c.CV: pool and rep -------
dfcv <- plyr::ddply(rep1_list, c('.id', 'sample_id'), function(data_tmp){
  cat(data_tmp$.id[1], data_tmp$sample_id[1], '...\r')
  expr <- data_tmp %>% select(-(1:Low_protein_IDs)) %>% removeColsAllNa()
  ret <- data.frame(
    Protein = colnames(expr),
    mean = apply(2^expr, 2, mean, na.rm = T),
    sd = apply(2^expr, 2, sd, na.rm = T),
    CV = apply(2^expr, 2, cv, na.rm = T)
  )
  return(ret)
}) %>% rename(Type = .id)
dfcv %<>% mutate(Type = factor(Type, levels = c('Pooling', 'Technical', 'Biological')))

dfcv_wide <- dfcv %>% pivot_wider(id_cols = Type:sample_id, names_from = Protein, values_from = CV) # to be saved
# dfcv <- dfcv_wide %>% pivot_longer(cols = -(Type:sample_id), names_to = 'Protein', values_to = 'CV')



#### output -----------
p_cv <- ggplot(dfcv) +
  aes(x = Type, y = CV, fill = Type) +
  # geom_jitter(color = '#CCCCCC', size=0.1, alpha=0.3,
  #             position = position_jitterdodge(jitter.width = 1, seed = 0)) +
  geom_violin(aes(color = Type), width = 0.9) +
  # geom_boxplot(width=0.1, color = '#000000', outlier.size = 0.5) +
  stat_summary(fun = median, geom = "point", size = 3, color = "black") +
  labs(x = '', y = 'Coefficient of Variation') +
  scale_color_manual(values = c("#8DA0CBFF", "#FC8D62FF", "#66C2A5FF")) +
  scale_fill_manual(values = c("#8DA0CB99", "#FC8D6299", "#66C2A599")) +
  theme_classic() +
  theme(text = element_text(size = 15, color = "black"),
        legend.position='none')

p_cv <- p_cv +
  annotate(
    'text',
    x = 'Pooling',
    y = max(na.omit(dfcv$CV)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcv$CV[dfcv$Type == 'Pooling'], na.rm = T))}\n",
      "(n = {n_pool})"
    ),
    color="black", size = 5, angle = 0
  ) +
  annotate(
    'text',
    x = 'Technical',
    y = max(na.omit(dfcv$CV)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcv$CV[dfcv$Type == 'Technical'], na.rm = T))}\n",
      "(pairs = {n_trep_pairs})"
    ),
    color="black", size = 5, angle = 0
  ) +
  annotate(
    'text',
    x = 'Biological',
    y = max(na.omit(dfcv$CV)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcv$CV[dfcv$Type == 'Biological'], na.rm = T))}\n",
      "(pairs = {n_brep_pairs})"
    ),
    color="black", size = 5, angle = 0
  )

ggsave('output/20251201TPHP_sF2c_CV_violin.pdf', p_cv, width = 6, height = 6)



### d.correlation: pool and rep -------
# dfcor <- plyr::ddply(rep1, 'sample_id', function(data_tmp){
#   cor_mat <- data_tmp %>%
#     column_to_rownames('file_id') %>% select(-(1:new)) %>%
#     t() %>% log2() %>% cor(use = 'pairwise.complete.obs', method = 'pearson')
#   ret <- data.frame(pearson.r = sort(cor_mat[upper.tri(cor_mat)], decreasing = T)) %>% 
#     mutate(rank = 1:nrow(.))
#   return(ret)
# })
# dfcor %<>% mutate(sample_id = ifelse(sample_id == 'p', 'Pooling', sample_id),
#                   Type = ifelse(sample_id == 'Pooling', sample_id, 'Replicates'))
dfcor <- plyr::ddply(rep1_list, c('.id', 'sample_id'), function(data_tmp){
  cat(data_tmp$.id[1], data_tmp$sample_id[1], '...\r')
  cor_mat <- data_tmp %>%
    column_to_rownames('file_id') %>% select(-(1:Low_protein_IDs)) %>%
    t() %>% log2() %>% cor(use = 'pairwise.complete.obs', method = 'pearson')
  
  vec <- cor_mat[upper.tri(cor_mat)]
  inds <- which(upper.tri(cor_mat), arr.ind = T) # indices of upper-triangle elements
  
  ret <- data.frame(pearson.r = vec,
                    file_id1 = rownames(cor_mat)[inds[, 'row']],
                    file_id2 = colnames(cor_mat)[inds[, 'col']]) %>% 
    arrange(desc(pearson.r)) %>% 
    mutate(rank = 1:nrow(.))
  return(ret)
}) %>% rename(Type = .id)
dfcor %<>% mutate(Type = factor(Type, levels = c('Pooling', 'Technical', 'Biological')))


#### output -----------
p_cor <- ggplot(dfcor)+
  aes(x = Type, y = pearson.r, fill = Type) +
  # geom_jitter(color = '#CCCCCC', size=0.1, alpha=0.3,
  #             position = position_jitterdodge(jitter.width = 1, seed = 0)) +
  geom_violin(aes(color = Type), width = 0.9) +
  geom_boxplot(width=0.1, color = '#000000', outlier.size = 0.5) +
  # stat_summary(fun=mean, geom="point", size=2, color="#000000") +
  labs(x = '', y = 'Pearson.r') +
  scale_color_manual(values = c("#8DA0CBFF", "#FC8D62FF", "#66C2A5FF")) +
  scale_fill_manual(values = c("#8DA0CB99", "#FC8D6299", "#66C2A599")) +
  theme_classic() +
  theme(text = element_text(size = 15,color = "black"), legend.position='none')

p_cor <- p_cor +
  annotate(
    'text',
    x = 'Pooling',
    y = max(na.omit(dfcor$pearson.r)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcor$pearson.r[dfcor$Type == 'Pooling'], na.rm = T))}\n",
      "(n = {n_pool})"
    ),
    color="black", size = 5, angle = 0
  ) +
  annotate(
    'text',
    x = 'Technical',
    y = max(na.omit(dfcor$pearson.r)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcor$pearson.r[dfcor$Type == 'Technical'], na.rm = T))}\n",
      "(pairs = {n_trep_pairs})"
    ),
    color="black", size = 5, angle = 0
  ) +
  annotate(
    'text',
    x = 'Biological',
    y = max(na.omit(dfcor$pearson.r)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcor$pearson.r[dfcor$Type == 'Biological'], na.rm = T))}\n",
      "(pairs = {n_brep_pairs})"
    ),
    color="black", size = 5, angle = 0
  )
ggsave('output/20251201TPHP_sF2d_pearson_violin.pdf', p_cor, width = 6, height = 6)


### save rds -----
qc.all <- list(cv = dfcv, cor = dfcor)
# saveRDS(qc.all, 'qc.all.rds')

ggsave('output/20251201TPHP_sF2cd.pdf',
       ggpubr::ggarrange(p_cv, p_cor), width = 12, height = 6)

list(all.cv = dfcv_wide, all.pearson = dfcor) %>% rio::export('output/20251201TPHP_sF2cd.xlsx')

## 3.4 sF3, T/NT data ------
rep2_list <- rep1_list %>% filter(!sample_type %in% c('F', 'N'), !FileName %in% pool_N)

# annotations (n/pairs of samples)
n_pool <- rep2_list %>% filter(Rep_type == 'p') %>% distinct(file_id) %>% nrow() # 92
n_trep_pairs <- rep2_list %>% filter(Rep_type == 't') %>% distinct(sample_id) %>% nrow() # 163
n_brep_pairs <- rep2_list %>% filter(Rep_type == 'b') %>% distinct(sample_id) %>% nrow() # 13

### c.CV: pool and rep -------
dfcv <- plyr::ddply(rep2_list, c('.id', 'sample_id'), function(data_tmp){
  cat(data_tmp$.id[1], data_tmp$sample_id[1], '...\r')
  expr <- data_tmp %>% select(-(1:Low_protein_IDs))
  ret <- data.frame(
    Protein = colnames(expr),
    mean = apply(2^expr, 2, mean, na.rm = T),
    sd = apply(2^expr, 2, sd, na.rm = T),
    CV = apply(2^expr, 2, cv, na.rm = T)
  )
  return(ret)
}) %>% rename(Type = .id)
dfcv %<>% mutate(Type = factor(Type, levels = c('Pooling', 'Technical', 'Biological')))

dfcv_wide <- dfcv %>% pivot_wider(id_cols = Type:sample_id, names_from = Protein, values_from = CV) # to be saved


#### output -----------
p_cv <- ggplot(dfcv) +
  aes(x = Type, y = CV, fill = Type) +
  # geom_jitter(color = '#CCCCCC', size=0.1, alpha=0.3,
  #             position = position_jitterdodge(jitter.width = 1, seed = 0)) +
  geom_violin(aes(color = Type), width = 0.6) +
  # geom_boxplot(width=0.1, color = '#000000', outlier.size = 0.5) +
  stat_summary(fun = median, geom = "point", size = 3, color = "black") +
  labs(x = '', y = 'Coefficient of Variation') +
  scale_color_manual(values = c("#8DA0CBFF", "#FC8D62FF", "#66C2A5FF")) +
  scale_fill_manual(values = c("#8DA0CB99", "#FC8D6299", "#66C2A599")) +
  theme_classic() +
  theme(text = element_text(size = 15, color = "black"),
        legend.position='none')

p_cv <- p_cv +
  annotate(
    'text',
    x = 'Pooling',
    y = max(na.omit(dfcv$CV)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcv$CV[dfcv$Type == 'Pooling'], na.rm = T))}\n",
      "(n = {n_pool})"
    ),
    color="black", size = 5, angle = 0
  ) +
  annotate(
    'text',
    x = 'Technical',
    y = max(na.omit(dfcv$CV)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcv$CV[dfcv$Type == 'Technical'], na.rm = T))}\n",
      "(pairs = {n_trep_pairs})"
    ),
    color="black", size = 5, angle = 0
  ) +
  annotate(
    'text',
    x = 'Biological',
    y = max(na.omit(dfcv$CV)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcv$CV[dfcv$Type == 'Biological'], na.rm = T))}\n",
      "(pairs = {n_brep_pairs})"
    ),
    color="black", size = 5, angle = 0
  )

ggsave('output/20251201TPHP_sF3c_CV_violin.pdf', p_cv, width = 6, height = 6)



### d.correlation: pool and rep -------
dfcor <- plyr::ddply(rep2_list, c('.id', 'sample_id'), function(data_tmp){
  cat(data_tmp$.id[1], data_tmp$sample_id[1], '...\r')
  cor_mat <- data_tmp %>%
    column_to_rownames('file_id') %>% select(-(1:Low_protein_IDs)) %>%
    t() %>% log2() %>% cor(use = 'pairwise.complete.obs', method = 'pearson')
  
  vec <- cor_mat[upper.tri(cor_mat)]
  inds <- which(upper.tri(cor_mat), arr.ind = T) # indices of upper-triangle elements
  
  ret <- data.frame(pearson.r = vec,
                    file_id1 = rownames(cor_mat)[inds[, 'row']],
                    file_id2 = colnames(cor_mat)[inds[, 'col']]) %>% 
    arrange(desc(pearson.r)) %>% 
    mutate(rank = 1:nrow(.))
  return(ret)
}) %>% rename(Type = .id)
dfcor %<>% mutate(Type = factor(Type, levels = c('Pooling', 'Technical', 'Biological')))


#### output -----------
p_cor <- ggplot(dfcor)+
  aes(x = Type, y = pearson.r, fill = Type) +
  # geom_jitter(color = '#CCCCCC', size=0.1, alpha=0.3,
  #             position = position_jitterdodge(jitter.width = 1, seed = 0)) +
  geom_violin(aes(color = Type), width = 0.9) +
  geom_boxplot(width=0.1, color = '#000000', outlier.size = 0.5) +
  # stat_summary(fun=mean, geom="point", size=2, color="#000000") +
  labs(x = '', y = 'Pearson.r') +
  scale_color_manual(values = c("#8DA0CBFF", "#FC8D62FF", "#66C2A5FF")) +
  scale_fill_manual(values = c("#8DA0CB99", "#FC8D6299", "#66C2A599")) +
  theme_classic() +
  theme(text = element_text(size = 15,color = "black"), legend.position='none')

p_cor <- p_cor +
  annotate(
    'text',
    x = 'Pooling',
    y = max(na.omit(dfcor$pearson.r)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcor$pearson.r[dfcor$Type == 'Pooling'], na.rm = T))}\n",
      "(n = {n_pool})"
    ),
    color="black", size = 5, angle = 0
  ) +
  annotate(
    'text',
    x = 'Technical',
    y = max(na.omit(dfcor$pearson.r)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcor$pearson.r[dfcor$Type == 'Technical'], na.rm = T))}\n",
      "(pairs = {n_trep_pairs})"
    ),
    color="black", size = 5, angle = 0
  ) +
  annotate(
    'text',
    x = 'Biological',
    y = max(na.omit(dfcor$pearson.r)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcor$pearson.r[dfcor$Type == 'Biological'], na.rm = T))}\n",
      "(pairs = {n_brep_pairs})"
    ),
    color="black", size = 5, angle = 0
  )
ggsave('output/20251201TPHP_sF3d_pearson_violin.pdf', p_cor, width = 6, height = 6)

### save rds -----
qc.tnt <- list(cv = dfcv, cor = dfcor)
# saveRDS(qc.tnt, 'qc.tnt.rds')

ggsave('output/20251201TPHP_sF3cd.pdf',
       ggpubr::ggarrange(p_cv, p_cor), width = 12, height = 6)

list(T.NT.cv = dfcv_wide, T.NT.pearson = dfcor) %>% rio::export('output/20251201TPHP_sF3cd.xlsx')






# Final: Source data ---------
list(sF2c_CV = dfcv_wide,
     sF2d_pearson = dfcor) %>%
  rio::export('output/source_data_sF2_V1009.xlsx')


# XXXXX --------------
pool <- rio::import('../0_process_DIA-NN_data/input/pool_pg_matrix_149_13609.csv') %>%
  column_to_rownames('V1')
dat_raw <- rbind(smp, pool) %>% t() %>% log2() %>% as.data.frame()

pool.matrices <- readRDS('pool.matrices.rds')
pool_cmb <- pool.matrices$combat
# dat_imp <- rio::import('input/imputed_raw_12754proteins_3005samples.csv') %>% 
#   column_to_rownames('V1') %>% log2()
dat_cmb <- rio::import('../2_batch_effect_evaluation_correction/output/batch_corrected_data_combat3.csv') %>%
  column_to_rownames('V1')
# mat_mask <- rio::import('../1_NA_filter_imputation/log2_imputed_qrilc_mask.csv') %>% 
#   column_to_rownames('V1') %>% 
#   as.matrix()
# mat_mask <- mat_mask[rownames(dat_imp), colnames(dat_imp)]
# dat_raw <- dat_imp
# dat_raw[mat_mask] <- NA


# log2 transform
# dat <- plyr::rbind.fill(smp, pool) %>%
#   t() %>% log2() %>% .[, info$FileName]
# dim(dat) # 13609  3005
# identical(colnames(dat), info$FileName) # TRUE
# 
# dat_imp <- dat_imp[, info$FileName]
# identical(colnames(dat_imp), info$FileName) # TRUE
smp1 <- dat_cmb[, setdiff(colnames(dat_cmb), colnames(pool))] %>% as.data.frame() %>% rownames_to_column('protein')
pool1 <- pool_cmb %>% as.data.frame() %>% rownames_to_column('protein')
dat <- smp1 %>% full_join(pool1) %>% column_to_rownames('protein') %>%
  .[, intersect(info$FileName, colnames(.))]
dim(dat) # 12754  2961
identical(colnames(dat), intersect(info$FileName, colnames(dat))) # TRUE
# 
# dat_imp <- dat_imp[, info$FileName]
# identical(colnames(dat_imp), info$FileName) # TRUE

# sample_id, Rep, file_id, etc.
info1 <- info %>%
  mutate(
    file_id = ifelse(!is.na(file_id), file_id, sample_id)
  )

not_rep <- info1 %>% count(sample_id) %>% filter(n == 1) %>% semi_join(info1, .) %>% 
  mutate(Is.Rep = FALSE)
rep <- info1 %>% count(sample_id) %>% filter(n > 1) %>% semi_join(info1, .) %>% 
  mutate(Is.Rep = TRUE) %>% 
  mutate(Rep_type = if_else(new == "data3" & str_detect(FileName, "brep"),
                            "b", Rep_type))

rep1p <- rep %>% filter(sample_type == 'p') %>% 
  mutate(file_id = ifelse(!is.na(file_id), file_id,
                          str_c('pool', str_extract(FileName, '^[A-Z]+'),
                                str_extract(FileName, '\\d+$'))))

rep1a <- rep %>% filter(sample_type != 'p') %>%
  group_by(sample_id) %>% slice(1) %>%
  ungroup() %>%
  mutate(
    Rep = NA, Rep_type = NA,
    file_id = ifelse(!is.na(file_id), str_remove(file_id, '_[bt]\\d+$'), sample_id)
  )
rep1b <- rep %>% filter(sample_type != 'p') %>%
  group_by(sample_id) %>% slice(-1) %>%
  mutate(
    Rep = ifelse(is.na(Rep), sample_id, Rep),
    Rep_type = ifelse(is.na(Rep_type), 't', Rep_type), # artificially checked
    file_id = str_c(sample_id, '_', Rep_type, row_number())
  ) %>% 
  ungroup()

info2 <- rbind(not_rep, rep1a, rep1b, rep1p) #%>% 
  # arrange(acq_order) %>% 
  # select(FileName, file_id, sample_id, sample_type:patient_ID, Is.Rep, Rep,
  #        Rep_type, instrument, batch, trans, date, year, month, yearmonth, new)

meta <- info2
rio::export(meta, '20251009_PUH_sample_information_3005files_meta.xlsx')


# info5 %<>% mutate(file_id = ifelse(!is.na(file_id), file_id, sample_id))
# rio::export(info5, 'input/20251009_PUH_sample_information_3005files_V6.xlsx')


# save(dat, smp, pool, dat_imp, meta, file = 'All_timsQC_V1116_input.RData')
save(dat, dat_raw, meta, file = 'All_timsQC_V1116_input.RData')


# 2.identity -----
df_ident <- data.frame(FileName = colnames(dat_raw),
                       `# proteins` = colSums(!is.na(dat_raw)),
                       check.names = F) %>% 
  filter(FileName %in% meta$FileName)

# 3.replicates --------
rep_check <- meta %>% filter(Is.Rep)

df_est_rep <- dat %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(rep_check, .)


rm.fid <- c()

rep.identity <- plyr::ddply(df_est_rep, 'sample_id', function(dfsub){
  cat(dfsub$sample_id[1], '...\r')
  # X <- dfsub %>%
  #   column_to_rownames('file_id') %>%
  #   select(-(1:new)) %>% t()
  X <- dat_raw[, dfsub$FileName]
  ret <- data.frame(file_id = dfsub$file_id,
                    `# proteins` = colSums(!is.na(X)),
                    check.names = F)
  ret$outlier.lower.ingroup <- get_outliers(ret$`# proteins`)[1]
  return(ret)
})
rep.identity %<>%
  group_by(sample_id) %>%
  summarise(identity.mean = mean(`# proteins`),
            identity.min = min(`# proteins`),
            identity.max = max(`# proteins`),
            identity.range.length = identity.max - identity.min) %>%
  inner_join(rep.identity, .) %>%
  arrange(identity.mean) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id)),
         Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup) %>% 
  inner_join(meta %>% select(file_id, Rep_type, class_abbr, sample_type))
rep.sid_order <- rep.identity %>% arrange(`# proteins`) %>% distinct(sample_id) %>% pull()
rep.identity %<>% 
  mutate(sample_id = factor(sample_id, rep.sid_order),
         Rep_type = ifelse(sample_type == 'p', 'p', Rep_type),
         Rep_type = ifelse(!is.na(Rep_type), Rep_type, 'Reference'),
         Rep_type = factor(Rep_type, c('Reference', 'b', 't', 'p')))

low.rep.ident <- rep.identity %>% filter(identity.range.length > 1000) %>%
  group_by(sample_id) %>% arrange(`# proteins`) %>% slice(1) %>% 
  pull(file_id)
rm.fid %<>% append(low.rep.ident)
rm.sid <- df_est_rep %>% filter(file_id %in% rm.fid) %>% distinct(sample_id) %>% pull() %>% setdiff(c('p'))


rep.pearson <- plyr::ddply(#df_est_rep %>%
                             # filter(!(sample_id %in% rm.sid),
                             #        !(file_id %in% rm.fid)),
                           df_est_rep,
                           'sample_id', function(dfsub){
                             cat(dfsub$sample_id[1], '...\r')
                             X <- dfsub %>%
                               column_to_rownames('file_id') %>%
                               select(-(1:new)) %>% t()
                             corX <- cor(X, use = 'pairwise.complete.obs', method = 'pearson')
                             corX_pair <- corX %>% as.data.frame() %>%
                               rownames_to_column('ID1') %>%
                               pivot_longer(cols = -ID1, names_to = 'ID2') %>%
                               unite('pair', ID1, ID2, sep = '-') %>%
                               pull(pair) %>%
                               matrix(nrow = nrow(corX), ncol = ncol(corX),
                                      byrow = T, dimnames = dimnames(corX))
                             ret <- data.frame(ID.pair = corX_pair %>% .[upper.tri(., diag = F)],
                                               pearson.r = corX %>% .[upper.tri(., diag = F)],
                                               check.names = F)
                             ret$outlier.lower.ingroup <- get_outliers(ret$pearson.r)[1]
                             return(ret)
                           })
rep.pearson %<>%
  group_by(sample_id) %>%
  summarise(pearson.r.mean = mean(pearson.r)) %>%
  inner_join(rep.pearson, .) %>%
  arrange(pearson.r.mean) %>%
  mutate(sample_id = factor(sample_id, levels = unique(sample_id)))

## ---- Correlation heatmaps by .id (simplified, robust) --------------------
# Input expected:
#   - rep1_list: either a named list of data.frames (e.g., list(Pooling=..., Technical=..., Biological=...))
#                 or a single data.frame that already contains a column `.id`
#   - new: integer; the first `new` columns are metadata; numeric assay columns start at (new+1)
# Output:
#   - One PDF and one PNG heatmap per `.id` under pheatmap_cor_by_id/

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})


### ### ---------- 0) Normalize input to a single data.frame with `.id` -----------
### normalize_rep1 <- function(x) {
###   if (is.list(x) && !inherits(x, "data.frame")) {
###     # Named list -> bind; `.id` from list names
###     dplyr::bind_rows(x, .id = ".id")
###   } else if (inherits(x, "data.frame")) {
###     stopifnot(".id" %in% names(x))
###     x
###   } else {
###     stop("`input` must be a named list of data.frames or a data.frame containing `.id`.")
###   }
### }
### 
### 
### # rep1_list_raw <- dat_raw %>% t() %>% as.data.frame() %>% 
### #   rownames_to_column('FileName') %>% 
### #   inner_join(rep1_list %>% select(1:new), .)
### # df_all <- normalize_rep1(rep1_list_raw)
### 
### dat_cmb_mask <- dat_cmb
### dat_cmb_mask[mat_mask] <- NA
### rep1_list_cmb<- dat_cmb_mask %>% t() %>% as.data.frame() %>% 
###   rownames_to_column('FileName') %>% 
###   inner_join(rep1_list %>% select(1:new), .)
### df_all <- normalize_rep1(rep1_list_cmb)
### 
### 
### 
### ### ---------- 1) Correlation matrix builder ---------------------------------
### # - Ensures unique file_id per group
### # - Selects numeric assay block after the first `new` columns
### # - log2 transform with a small pseudocount to avoid -Inf
### build_cor_mat <- function(data_sub, new, pseudocount = 1) {
###   ds <- data_sub %>%
###     dplyr::distinct(file_id, .keep_all = TRUE)
###   
###   # extract numeric matrix
###   mat <- ds %>%
###     dplyr::select(-(1:new)) %>%
###     dplyr::mutate(across(everything(), as.numeric)) %>%
###     as.matrix()
###   
###   # rows = file_id, cols = assay features -> transpose to have cols = file_id for cor()
###   rownames(mat) <- ds$file_id
###   mat <- t(mat)
###   
###   # # transform
###   # mat <- log2(mat + pseudocount)
###   # mat[!is.finite(mat)] <- NA
###   
###   # if fewer than 2 file_id, return NA matrix to be skipped later
###   if (ncol(mat) < 2) return(matrix(NA_real_, nrow = 1, ncol = 1,
###                                    dimnames = list("insufficient", "insufficient")))
###   
###   cor(mat, use = "pairwise.complete.obs", method = "pearson")
### }
### 
### ### ---------- 2) Split by .id and compute correlations ----------------------
### df_list  <- split(df_all, df_all$.id)        # names: .id values
### cor_list <- purrr::imap(df_list, ~ build_cor_mat(.x, new = new))
### 
### ### ---------- 3) Build annotations per .id (de-duplicate file_id) -----------
### .type_lvls <- c("Pooling", "Technical", "Biological")
### 
### ann_list <- purrr::imap(cor_list, function(cor_mat, id_name) {
###   ids <- colnames(cor_mat)
###   if (is.null(ids)) ids <- character(0)
###   
###   ann <- df_all %>%
###     dplyr::filter(.id == id_name, file_id %in% ids) %>%
###     dplyr::group_by(file_id) %>%
###     dplyr::summarise(
###       # sample_id = paste(unique(sample_id), collapse = ";"),
###       Type = first(.id),
###       instrument = instrument,
###       yearmonth = yearmonth,
###       new = new,
###       .groups = "drop"
###     ) %>%
###     dplyr::mutate(Type = factor(Type, levels = .type_lvls)) %>%
###     as.data.frame()
###   
###   rownames(ann) <- ann$file_id
###   ann$file_id <- NULL
###   
###   list(
###     ann_col = ann[ids, , drop = FALSE],
###     ann_row = ann[rownames(cor_mat), , drop = FALSE]
###   )
### })
### 
### ### ---------- 4) Plot & export ----------------------------------------------
### out_dir <- "pheatmap_cor_by_id_cmb"
### if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
### 
### purrr::iwalk(cor_list, function(cor_mat, id_name) {
###   # skip groups with insufficient samples
###   if (ncol(cor_mat) < 2 || all(!is.finite(cor_mat))) {
###     warning(sprintf("%s: skipped heatmap (insufficient distinct file_id or all-NA)", id_name))
###     return(invisible(NULL))
###   }
###   
###   ann <- ann_list[[id_name]]
###   pdf_file <- file.path(out_dir, sprintf("cor_%s.pdf", id_name))
###   png_file <- file.path(out_dir, sprintf("cor_%s.png", id_name))
###   
###   grDevices::pdf(pdf_file, width = 8, height = 8)
###   pheatmap::pheatmap(
###     cor_mat,
###     annotation_col = ann$ann_col,
###     annotation_row = ann$ann_row,
###     annotation_colors = list(
###       new = ggsci::pal_npg()(length(unique(ann$ann_col$new))) %>% setNames(sort(unique(ann$ann_col$new)))
###     ),
###     main = sprintf("Pearson correlation (file_id) — %s", id_name),
###     clustering_distance_rows = "euclidean",
###     clustering_distance_cols = "euclidean",
###     clustering_method = "complete",
###     show_rownames = FALSE,# display_numbers = T, number_color = '#222222',
###     show_colnames = FALSE,
###     border_color = NA, fontsize = 7,
###     annotation_legend = TRUE,
###     legend = TRUE
###   )
###   grDevices::dev.off()
###   
###   # grDevices::png(png_file, width = 2000, height = 2000, res = 300)
###   # pheatmap::pheatmap(
###   #   cor_mat,
###   #   annotation_col = ann$ann_col,
###   #   annotation_row = ann$ann_row,
###   #   main = sprintf("Pearson correlation (file_id) — %s", id_name),
###   #   clustering_distance_rows = "euclidean",
###   #   clustering_distance_cols = "euclidean",
###   #   clustering_method = "complete",
###   #   show_rownames = FALSE,
###   #   show_colnames = FALSE,
###   #   border_color = NA,
###   #   annotation_legend = TRUE,
###   #   legend = TRUE
###   # )
###   # grDevices::dev.off()
### })
### 
### 
### 

# 4.samples -----
df_est <- dat %>% t() %>% as.data.frame() %>%
  rownames_to_column('FileName') %>% 
  inner_join(meta, .)

est.identity <- df_est %>%
  filter(sample_id != 'p') %>% 
  filter(!(file_id %in% rm.fid)) %>% 
  plyr::ddply(c('sample_type', 'class_abbr'), function(dfsub){
    cat(dfsub$sample_type[1], dfsub$class_abbr[1], '...\r')
    # X <- dfsub %>%
    #   column_to_rownames('FileName') %>%
    #   select(-(1:date)) %>% t()
    X <- dat_raw[, dfsub$FileName, drop = FALSE]
    ret <- data.frame(FileName = dfsub$FileName,
                      `# proteins` = colSums(!is.na(X)),
                      check.names = F)
    ret$outlier.lower.ingroup <- get_outliers(ret$`# proteins`)[1]
    return(ret)
  })
est.identity %<>%
  group_by(sample_type, class_abbr) %>%
  summarise(identity.mean = mean(`# proteins`),
            identity.min = min(`# proteins`),
            identity.max = max(`# proteins`),
            identity.range.length = identity.max - identity.min,
            .groups = 'drop') %>%
  inner_join(est.identity, .) %>%
  mutate(Group = str_c(sample_type, ' - ', class_abbr)) %>% 
  arrange(identity.mean) %>%
  mutate(Group = factor(Group, levels = unique(Group)),
         Is.Lower.Ingroup = `# proteins` < outlier.lower.ingroup) %>% 
  inner_join(meta %>% select(FileName, Rep_type, class_abbr, sample_type, tissue_name))
group_order <- est.identity %>% arrange(`# proteins`) %>% distinct(Group) %>% pull()
est.identity %<>% mutate(Group = factor(Group, group_order))



rm.fid
# [1] "Rb13_40"    "Rb1_76_t1"  "Rb9_53_b1"  "Rb4_81_b1"  "Rb12_60_b1" "poolN7467"  "DIA_31_b1"  "Rb11_9"     "DIA_29_b1"  "Rb9_83"    
# [11] "SCA_500_t1" "CA_DIA_202"
meta %>% pull(FileName, file_id) %>% .[rm.fid] %>% unname()
# [1] "CAD20250829yuel_TPHP_RCA_90minDIA_b13_40_Slot2-1_1_13637"       "N20250708yuel_TPHP_RCA_90mDIA_b1_76_Slot2-27_1_34640"          
# [3] "N20250806yuel_TPHP_RCA_90mDIA_b9_53_brep_Slot2-51_1_35244"      "CAD20250715yuel_TPHP_RCA_90minDIA_b4_81_brep_Slot2-50_1_12979" 
# [5] "CAD20250819yuel_TPHP_RCA_90minDIA_b12_60_brep_Slot2-50_1_13560" "N20210825yuel_nail_pool_Slot1-6_1_7467"                        
# [7] "K20200716yuel_TPHP_DIA_353_Slot2-4_1_1215"                      "N20250814yuel_TPHP_RCA_90mDIA_b11_9_Slot1-24_1_35300"          
# [9] "K20200713yuel_TPHP_DIA_351_Slot2-2_1_1189"                      "N20250806yuel_TPHP_RCA_90mDIA_b9_83_Slot2-44_1_35232"          
# [11] "N20220306yuel_TPHP_DIA_SCA_500_Slot1-24_1_11492"                "N20201111yuel_TPHP_202_Slot1-32_1_2086" 

rm.sid
# [1] "CA_DIA_202" "Rb11_9"     "Rb13_40"    "Rb9_83"     "DIA_29"     "DIA_31"     "Rb12_60"    "Rb1_76"     "Rb4_81"     "Rb9_53"    
# [11] "SCA_500"

# Output of 1-4 -------
plot_ident_rep <- ggplot(rep.identity) +
  aes(x = `# proteins`, y = sample_id) +
  geom_boxplot(data = rep.identity %>% filter(identity.range.length <= 1000),
               color = 'black', outlier.color = '#c23190', outlier.size = 3) +
  geom_boxplot(data = rep.identity %>% filter(identity.range.length > 1000),
               color = 'red4', outlier.color = '#c23190', outlier.size = 3) +
  geom_point(aes(color = class_abbr, shape = sample_type), size = 1.2) +
  labs(x = "# proteins", y = "Sample ID", subtitle = "Protein identification") +
  # ggsci::scale_color_aaas(name = 'Replicates type') +
  scale_color_manual(values = organ_color) +
  theme_bw() +
  theme(text = element_text(size = 10, color = 'black'))
plot_pearson_rep <- ggplot(rep.pearson) +
  aes(x = pearson.r, y = sample_id) +
  # geom_boxplot(color = '#000000') +
  geom_point(color = '#000000', size = 1.2) +
  labs(x = "Pearson's r", y = "Sample ID", subtitle = "Correlation - same.sample.id") +
  theme_bw() +
  theme(text = element_text(size = 10, color = 'black'))
plot_rep <- ggarrange(plot_ident_rep, plot_pearson_rep, nrow = 1, ncol = 2, widths = c(4, 3))
ggsave('All_tims_rep_20251116.pdf', plot_rep, width = 10, height = 25)


plot_ident <- ggplot(est.identity) +
  aes(x = `# proteins`, y = Group) +
  geom_boxplot(data = est.identity %>% filter(identity.range.length <= 1000),
               color = 'black', outlier.color = '#c23190', outlier.size = 3) +
  geom_boxplot(data = est.identity %>% filter(identity.range.length > 1000),
               color = 'red4', outlier.color = '#c23190', outlier.size = 3) +
  geom_point(aes(color = tissue_name), size = 1.2) +
  labs(x = "# proteins", y = "Sample ID", subtitle = "Protein identification") +
  # ggsci::scale_color_aaas(name = 'Sample type') +
  theme_bw() +
  theme(text = element_text(size = 10, color = 'black'))
ggsave('All_tims_identity_20251116.pdf', plot_ident, width = 40, height = 18)


list(rep.identity = rep.identity,
     rep.pearson = rep.pearson,
     sample.identity = est.identity) %>% 
  rio::export('All_tims_QC_source_20251116.xlsx')

df_est_rep %>% filter(!(file_id %in% rm.fid)) %>%
  rio::export('output/All_tims_replicates_20251116.xlsx')


# 5.sF2 ---------
meta <- rio::import('20251009_PUH_sample_information_3005files_meta.xlsx')
# meta_filter <- rio::import('PUH_simple_meta_2990files_10labels.xlsx') %>% filter(new != 'data3')
# 'N20210825yuel_nail_pool_Slot1-6_1_7467' %in% meta_filter$FileName # FALSE
meta <- meta %>%
  mutate(transport.capillary = trans,
         month = str_sub(date, 1, 8),
         tissue_name = ifelse(is.na(tissue_name), 'pool', tissue_name)) %>% 
  select(FileName, file_id, sample_id, sample_type, tissue_name, anatomical_classification,
         instrument, yearmonth, batch, transport.capillary) 

load('All_timsQC_V1116_input.RData')
rep1 <- rio::import('output/All_tims_replicates_20251116.xlsx') # i.e. df_est_rep

# dim(dat) # 12754  2961
# # dim(pool) # 150 13609
# dim(pool1) # 5733  150
# dim(rep1) # 653 12788
# dim(meta) # 2970   10
# 
# dat <- dat[, intersect(colnames(dat), meta$FileName)]
# # pool <- t(pool)[, intersect(colnames(pool), meta$FileName)]
# all1 <- cbind(dat, pool)
# dim(dat) # 13544  1692
# dim(pool) # 13544   119
# dim(all1) # 13544  1811

all1 <- dat_raw
# for quantitative analysis
all1.fill <- dat
# all1.fill <- rio::import('../2_batch_effect_evaluation_correction/input/quantile_log2_imputed_12754proteins_3005samples.csv') %>% 
#   column_to_rownames('V1')
# dim(all1.fill) # 12754  3005
# all1.fill <- all1.fill[intersect(rownames(all1), rownames(all1.fill)), intersect(colnames(all1), colnames(all1.fill))]
# dim(all1.fill) # 12754  3005

# rep1 data sorting
rep1 <- rep1 %>% filter(FileName %in% meta$FileName) # 662
rep1 <- rep1 %>% count(sample_id) %>% filter(n == 1) %>% anti_join(rep1, .) # 651
rep1 %<>% mutate(Rep_type = ifelse(sample_type == 'p', 'p', Rep_type)) %>%
  mutate(Rep_type = ifelse(is.na(Rep_type), 'o', Rep_type))
dim(rep1) #  642 12788

# Check wrong labels: 'o' annotated as 't'
rep1_type_check <- rep1 %>%
  mutate(Rep_type_check = str_extract(file_id, '([bt])\\d$', group = 1),
         Rep_type_check = ifelse(sample_type == 'p', 'p', Rep_type_check),
         Rep_type_check = ifelse(is.na(Rep_type_check), 'o', Rep_type_check)) %>% 
  select(FileName, sample_id, file_id, Rep_type, Rep_type_check)
o.fid <- rep1_type_check %>% filter(Rep_type != Rep_type_check) %>% pull(file_id)

rep1 %<>%
  mutate(Rep_type = ifelse(file_id %in% o.fid, 'o', Rep_type)) %>%
  mutate(Type = c(o = 'Original', t = 'Technical', b = 'Biological', p = 'Pooling')[Rep_type], .before = new)

# Separate 't' and 'b'
rep1p <- rep1 %>% filter(Rep_type == 'p')
rep1t <- rep1 %>% filter(Rep_type %in% c('o', 't'))
rep1t <- rep1t %>% count(sample_id) %>% filter(n > 1) %>% semi_join(rep1t, .)
rep1b <- rep1 %>% filter(Rep_type %in% c('o', 'b'))
rep1b <- rep1b %>% count(sample_id) %>% filter(n > 1) %>% semi_join(rep1b, .)

# some sample_ids without 'o' (only t and b) -- transform t to o
rbind(rep1b, rep1t, rep1p) %>% distinct(file_id) %>% nrow() # 651
rbind(rep1b, rep1t, rep1p) %>% distinct(file_id) %>% anti_join(rep1, .) %>% select(FileName, sample_id:cancer_abbr, DateTime:new) %>% View() # NULL

rbind(rep1b, rep1t, rep1p) %>% distinct(file_id) %>% nrow() %>% identical(nrow(rep1)) # TRUE (651)

# envelope data_list (actually a data.frame)
rep1_list <- list(Pooling = rep1p, Technical = rep1t, Biological = rep1b) %>% plyr::ldply()

# annotations (n/pairs of samples)
n_pool <- rep1_list %>% filter(Rep_type == 'p') %>% distinct(file_id) %>% nrow()
n_trep_pairs <- rep1_list %>% filter(Rep_type == 't') %>% distinct(sample_id) %>% nrow()
n_brep_pairs <- rep1_list %>% filter(Rep_type == 'b') %>% distinct(sample_id) %>% nrow()




pm_list <- list(all.samples = all1.fill, pool.only = all1.fill[, meta$FileName[meta$sample_type == 'p']])
meta_df <- meta %>% select(-FileName)

## save.image -----
# save.image(file = 'All_timsQC_v1116_sF2_data.RData')
load('All_timsQC_v1116_sF2_data.RData')

## c.CV: pool and rep -------
# rep1_list <- list(Pooling = rep1 %>% filter(sample_type == 'p'),
#                   Replicates = rep1 %>% filter(sample_type != 'p'))
# dfcv <- plyr::ldply(rep1_list, function(data_tmp){
#   ret <- data.frame(
#     Protein = data_tmp %>% select(-(1:new)) %>% colnames(),
#     mean = data_tmp %>% 
#       select(-(1:new)) %>% 
#       apply(2, mean, na.rm = T),
#     sd = data_tmp %>% 
#       select(-(1:new)) %>% 
#       apply(2, sd, na.rm = T),
#     CV = data_tmp %>% 
#       select(-(1:new)) %>% 
#       apply(2, cv)
#   )
#   return(ret)
# }, .id = 'Type')
dfcv <- plyr::ddply(rep1_list, c('.id', 'sample_id'), function(data_tmp){
  cat(data_tmp$.id[1], data_tmp$sample_id[1], '...\r')
  expr <- data_tmp %>% select(-(1:Low_protein_IDs))
  ret <- data.frame(
    Protein = colnames(expr),
    mean = apply(2^expr, 2, mean, na.rm = T),
    sd = apply(2^expr, 2, sd, na.rm = T),
    CV = apply(2^expr, 2, cv, na.rm = T)
  )
  return(ret)
}) %>% rename(Type = .id)
dfcv %<>% mutate(Type = factor(Type, levels = c('Pooling', 'Technical', 'Biological')))

dfcv_wide <- dfcv %>% pivot_wider(id_cols = Type:sample_id, names_from = Protein, values_from = CV) # to be saved
# dfcv <- dfcv_wide %>% pivot_longer(cols = -(Type:sample_id), names_to = 'Protein', values_to = 'CV')



### output -----------
p_cv <- ggplot(dfcv) +
  aes(x = Type, y = CV, fill = Type) +
  # geom_jitter(color = '#CCCCCC', size=0.1, alpha=0.3,
  #             position = position_jitterdodge(jitter.width = 1, seed = 0)) +
  geom_violin(aes(color = Type), width = 0.6) +
  # geom_boxplot(width=0.1, color = '#000000', outlier.size = 0.5) +
  stat_summary(fun = median, geom = "point", size = 3, color = "black") +
  labs(x = '', y = 'Coefficient of Variation') +
  scale_color_manual(values = c("#8DA0CBFF", "#FC8D62FF", "#66C2A5FF")) +
  scale_fill_manual(values = c("#8DA0CB99", "#FC8D6299", "#66C2A599")) +
  theme_classic() +
  theme(text = element_text(size = 15, color = "black"),
        legend.position='none')

p_cv <- p_cv +
  annotate(
    'text',
    x = 'Pooling',
    y = max(na.omit(dfcv$CV)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcv$CV[dfcv$Type == 'Pooling'], na.rm = T))}\n",
      "(n = {n_pool})"
    ),
    color="black", size = 5, angle = 0
  ) +
  annotate(
    'text',
    x = 'Technical',
    y = max(na.omit(dfcv$CV)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcv$CV[dfcv$Type == 'Technical'], na.rm = T))}\n",
      "(pairs = {n_trep_pairs})"
    ),
    color="black", size = 5, angle = 0
  ) +
  annotate(
    'text',
    x = 'Biological',
    y = max(na.omit(dfcv$CV)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcv$CV[dfcv$Type == 'Biological'], na.rm = T))}\n",
      "(pairs = {n_brep_pairs})"
    ),
    color="black", size = 5, angle = 0
  )

ggsave('output/20251116TPHP_sF2c_CV_violin.pdf', p_cv, width = 6, height = 6)



## d.correlation: pool and rep -------
# dfcor <- plyr::ddply(rep1, 'sample_id', function(data_tmp){
#   cor_mat <- data_tmp %>%
#     column_to_rownames('file_id') %>% select(-(1:new)) %>%
#     t() %>% log2() %>% cor(use = 'pairwise.complete.obs', method = 'pearson')
#   ret <- data.frame(pearson.r = sort(cor_mat[upper.tri(cor_mat)], decreasing = T)) %>% 
#     mutate(rank = 1:nrow(.))
#   return(ret)
# })
# dfcor %<>% mutate(sample_id = ifelse(sample_id == 'p', 'Pooling', sample_id),
#                   Type = ifelse(sample_id == 'Pooling', sample_id, 'Replicates'))
dfcor <- plyr::ddply(rep1_list, c('.id', 'sample_id'), function(data_tmp){
  cat(data_tmp$.id[1], data_tmp$sample_id[1], '...\r')
  cor_mat <- data_tmp %>%
    column_to_rownames('file_id') %>% select(-(1:Low_protein_IDs)) %>%
    t() %>% log2() %>% cor(use = 'pairwise.complete.obs', method = 'pearson')
  
  vec <- cor_mat[upper.tri(cor_mat)]
  inds <- which(upper.tri(cor_mat), arr.ind = T) # indices of upper-triangle elements
  
  ret <- data.frame(pearson.r = vec,
                    file_id1 = rownames(cor_mat)[inds[, 'row']],
                    file_id2 = colnames(cor_mat)[inds[, 'col']]) %>% 
    arrange(desc(pearson.r)) %>% 
    mutate(rank = 1:nrow(.))
  return(ret)
}) %>% rename(Type = .id)
dfcor %<>% mutate(Type = factor(Type, levels = c('Pooling', 'Technical', 'Biological')))


### output -----------
p_cor <- ggplot(dfcor)+
  aes(x = Type, y = pearson.r, fill = Type) +
  # geom_jitter(color = '#CCCCCC', size=0.1, alpha=0.3,
  #             position = position_jitterdodge(jitter.width = 1, seed = 0)) +
  geom_violin(aes(color = Type), width = 0.9) +
  geom_boxplot(width=0.1, color = '#000000', outlier.size = 0.5) +
  # stat_summary(fun=mean, geom="point", size=2, color="#000000") +
  labs(x = '', y = 'Pearson.r') +
  scale_color_manual(values = c("#8DA0CBFF", "#FC8D62FF", "#66C2A5FF")) +
  scale_fill_manual(values = c("#8DA0CB99", "#FC8D6299", "#66C2A599")) +
  theme_classic() +
  theme(text = element_text(size = 15,color = "black"), legend.position='none')

p_cor <- p_cor +
  annotate(
    'text',
    x = 'Pooling',
    y = max(na.omit(dfcor$pearson.r)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcor$pearson.r[dfcor$Type == 'Pooling'], na.rm = T))}\n",
      "(n = {n_pool})"
    ),
    color="black", size = 5, angle = 0
  ) +
  annotate(
    'text',
    x = 'Technical',
    y = max(na.omit(dfcor$pearson.r)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcor$pearson.r[dfcor$Type == 'Technical'], na.rm = T))}\n",
      "(pairs = {n_trep_pairs})"
    ),
    color="black", size = 5, angle = 0
  ) +
  annotate(
    'text',
    x = 'Biological',
    y = max(na.omit(dfcor$pearson.r)) * 1.01,
    label = str_glue(
      "Median: {sprintf('%.3f', median(dfcor$pearson.r[dfcor$Type == 'Biological'], na.rm = T))}\n",
      "(pairs = {n_brep_pairs})"
    ),
    color="black", size = 5, angle = 0
  )
ggsave('output/20251116TPHP_sF2d_pearson_violin.pdf', p_cor, width = 6, height = 6)



## e.reduce dimension: pool and all.samples -------
id_col <- 'file_id'
date_col <- 'DateTime'
var_col <- colnames(meta_df) %>% setdiff(c(id_col, date_col, 'sample_id'))
seed <- 10

res_dr_ls <- list()
for(i in seq_along(names(pm_list))){
  cat(i, '...\n')
  nm <- names(pm_list)[i]
  X <- pm_list[[nm]]%>% set_colnames(., unname(pull(meta, file_id, FileName)[colnames(.)]))
  
  res_dr_ls[[i]] <- beca.DR(X, meta_df, id_col, var_col, date_col, seed)
}
names(res_dr_ls) <- names(pm_list)

export_plots_separated_legends(res_dr_ls$all.samples$plots_pca, "output/20251009TPHP_pca_plots_all.pdf", width = 6, height = 6)
export_plots_separated_legends(res_dr_ls$all.samples$plots_tsne, "output/20251009TPHP_tsne_plots_all.pdf", width = 6, height = 6)
export_plots_separated_legends(res_dr_ls$all.samples$plots_umap, "output/20251009TPHP_umap_plots_all.pdf", width = 6, height = 6)

export_plots_separated_legends(res_dr_ls$pool.only$plots_pca, "output/20251009TPHP_pca_plots_pool.pdf", width = 6, height = 6)
export_plots_separated_legends(res_dr_ls$pool.only$plots_tsne, "output/20251009TPHP_tsne_plots_pool.pdf", width = 6, height = 6)
export_plots_separated_legends(res_dr_ls$pool.only$plots_umap, "output/20251009TPHP_umap_plots_pool.pdf", width = 6, height = 6)

# save(res_dr_ls$all.samples, file = 'output/res_dr_all1_13255x1811.RData')
save(res_dr_ls, file = 'output/res_dr_ls_all1_12754 x 2847+150.RData')


# Pool Matrix Heatmap --------
## annotations ----
ann_col.new <- info1 %>% column_to_rownames('FileName') %>%
  filter(sample_type == 'p') %>% 
  select(yearmonth, instrument, new)
ann_color.new <- list(
  new = ggsci::pal_d3()(length(unique(ann_col.new$new))) %>% setNames(sort(unique(ann_col.new$new))),
  instrument = ggsci::pal_futurama()(length(unique(ann_col.new$instrument))) %>% setNames(sort(unique(ann_col.new$instrument))),
  yearmonth = colorRampPalette(colors = c(RColorBrewer::brewer.pal(11, "Spectral")))(length(unique(ann_col.new$yearmonth))) %>% setNames(sort(unique(ann_col.new$yearmonth)))
)



## raw --------
pool_raw <- t(pool)

pdf('pool_matrix_raw.pdf', width = 18, height = 8)
pheatmap::pheatmap(
  log2(pool_raw)[rowSums(!is.na(pool_raw)) == ncol(pool_raw), ],
  scale = 'row',
  annotation_col = ann.new,
  # # annotation_row = ann$ann_row,
  annotation_colors = ann_color.new,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = FALSE,# display_numbers = T, number_color = '#222222',
  show_colnames = TRUE,
  border_color = NA, fontsize = 7,
  annotation_legend = TRUE,
  legend = TRUE
)
graphics.off()


## cmb --------
pool_cmb_mask <- dat_cmb_mask[, meta$FileName[meta$sample_type == 'p']]

pdf('pool_matrix_cmb.pdf', width = 18, height = 8)
pheatmap::pheatmap(
  pool_cmb_mask[rowSums(!is.na(pool_cmb_mask)) == ncol(pool_cmb_mask), ],
  scale = 'row',
  annotation_col = ann.new,
  # # annotation_row = ann$ann_row,
  annotation_colors = ann_color.new,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = FALSE,# display_numbers = T, number_color = '#222222',
  show_colnames = TRUE,
  border_color = NA, fontsize = 7,
  annotation_legend = TRUE,
  legend = TRUE
)
graphics.off()

## imp --------
dat_imp_mask <- dat_imp
dat_imp_mask[mat_mask] <- NA
pool_imp_mask <- dat_imp_mask[, meta$FileName[meta$sample_type == 'p']]

pdf('pool_matrix_imp.pdf', width = 18, height = 8)
pheatmap::pheatmap(
  log2(pool_imp_mask )[rowSums(!is.na(pool_imp_mask )) == ncol(pool_imp_mask ), ],
  scale = 'row',
  annotation_col = ann.new,
  # # annotation_row = ann$ann_row,
  annotation_colors = ann_color.new,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = FALSE,# display_numbers = T, number_color = '#222222',
  show_colnames = TRUE,
  border_color = NA, fontsize = 7,
  annotation_legend = TRUE,
  legend = TRUE
)
graphics.off()











