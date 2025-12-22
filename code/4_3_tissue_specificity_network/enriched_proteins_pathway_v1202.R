# 0.config ------------
# rm(list = ls())
# pacman::p_unload(pacman::p_loaded(), character.only = T)

wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
source('../source/source_code.R')
level_color <- c(RNA = 'green4', Protein = 'purple4', Both = 'blue4')

## helper functions ----------
my.enrich.network <- function(enrich.table, size.scales = c(10, 370), seed = 1, plotOnly = FALSE, simMethod = 'jaccard', clustMethod = 'markov', clustNameMethod = 'pagerank', drawEllipses = TRUE, fontSize = 4, nodeSize = 'Pathway.Size', colorBy = 'p.adjust', colorType = 'pval', minClusterSize = 2){
  require(aPEAR)
  
  set.seed(seed)
  enw <- enrichmentNetwork(
    enrich.table, plotOnly = plotOnly,
    clustMethod = clustMethod, # method for detecting pathway clusters. Available methods: 'markov', 'hier' and 'spectral'. Using 'spectral' method requires that you have Spectrum package installed.
    simMethod = simMethod, # method for calculating similarity matrix between the pathways. Available methods: 'jaccard', 'cosine' and 'cor'
    clustNameMethod = clustNameMethod, # method for selecting cluster names. Available methods: 'pagerank', 'hits' and 'none'
    drawEllipses = drawEllipses, fontSize = fontSize, nodeSize = nodeSize,
    colorBy = colorBy, colorType = colorType,
    minClusterSize = minClusterSize
  )
  enw$plot <- enw$plot +
    scale_color_viridis_c(name = 'Log10(p.adjust)', option = 'C', end = 0.8, limits = c(-10, -2)) +
    scale_size_continuous(name = "Pathway size", range = c(2, 6), limits = size.scales)
  return(enw)
}


# 1.Read data --------
info5 <- rio::import('../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V6.xlsx')
info5 %<>% filter(!Tissue_heterogeneity_abnormal, !Low_protein_IDs, !FileName %in% pancancer_low_pearson)

# Inputs
res_ts <- readRDS('../4_tissue_specificity_analysis/output/20251201_tissue_Specificity_naive_05minImputated_3.rds')

# df_enrich <- res_ts$ts.Enrich
df_tis_en <- res_ts$ts.tisEnrich
df_tis_en_ <- df_tis_en %>% inner_join(info5 %>% distinct(class_abbr, anatomical_classification))
# ts.med <- res_ts$med.matrix

# Basic check of tissue specificity classes
count_tbl <- df_tis_en %>% dplyr::count(`Tissue specificity`, name = 'n')

# Tissue-enriched proteins (UniProt IDs)
enriched_uniprot <- df_tis_en$Protein %>% unique()
length(enriched_uniprot) # 1723

# Metascape enrichment
writeClipboard(enriched_uniprot)

# split by tissue classification
enrich_list <- split(df_tis_en$Protein, df_tis_en$class_abbr)


# class colors
class_color <- c(organ_color, cancer_color) %>% unname() %>% 
  .[1:length(sort(unique(res_ts$ts.all$class_abbr)))] %>% 
  setNames(sort(unique(res_ts$ts.all$class_abbr)))




# 2.Enriched proteins, common enrich -----
library(clusterProfiler)
library(org.Hs.eg.db)

## 2.1 compareCluster with Enrich-list ------
### gobp -------
ego_enrich <- compareCluster(
  geneClusters = enrich_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = FALSE
)
# saveRDS(ego_enrich, 'N_Enrich_ego_compare.rds')
# ego_enrich <- readRDS('N_Enrich_ego_compare.rds')
# dotplot(ego_enrich, color = "p.adjust", showCategory = 2, font.size = 10)

### kegg -----
tmp <- clusterProfiler::bitr_kegg(enriched_uniprot, fromType='uniprot', toType='ncbi-proteinid', organism='hsa') %>% 
  setNames(c('Protein', 'proteinid'))
# 1.54% of input gene IDs are fail to map...
df_tis_en1 <- tmp %>% inner_join(df_tis_en)
tmp1 <- df_tis_en1 %>% distinct(class_abbr, proteinid)

enrich_pid_list <- split(tmp1$proteinid, tmp1$class_abbr)
# saveRDS(enrich_pid_list, 'enrich_pid_list.rds')
# enrich_pid_list <- readRDS('enrich_pid_list.rds')
ekegg_enrich <- compareCluster(enrich_pid_list, fun = enrichKEGG,
                               pAdjustMethod = "BH",
                               keyType = 'ncbi-proteinid',
                               minGSSize = 10, maxGSSize = 500,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2)
# saveRDS(ekegg_enrich, 'N_Enrich_ekegg_compare.rds')
# ekegg_enrich <- readRDS('N_Enrich_ekegg_compare.rds')
# dotplot(ekegg_enrich)


pacman::p_unload(clusterProfiler, enrichplot, ReactomePA, AnnotationDbi, BiocGenerics, IRanges, S4Vectors, Biobase)

## 2.2 compute generality per pathway term ----------------------
### gobp ----
ego_enrich_df <- ego_enrich %>%
  filter(pvalue < 0.05, qvalue < 0.05) %>% # in fact redundant
  as.data.frame() %>%
  separate_rows(geneID, sep = '/') %>% 
  rename(class_abbr = Cluster,
         Protein = geneID) %>% 
  mutate(Pathway.Size = as.numeric(str_extract(BgRatio, '^\\d+'))) %>% 
  inner_join(df_tis_en %>% distinct(class_abbr, Protein), .)

go_enrich_summary <- ego_enrich_df %>%
  group_by(ID, Description, Pathway.Size) %>%
  summarise(
    generality = n_distinct(class_abbr),
    Specificity = length(unique(res_ts$ts.all$class_abbr)) - generality,
    # n_up = n_distinct(class_abbr[direction == "Up"]),
    # n_down = n_distinct(class_abbr[direction == "Down"]),
    # geneID in compareCluster is "/"-delimited character of SYMBOLs per cancer
    Protein.list = list(unique(Protein[!is.na(Protein)])),
    RichFactor = stats::median(RichFactor, na.rm = TRUE),
    FoldEnrichment = stats::median(FoldEnrichment, na.rm = TRUE),
    zScore = stats::median(zScore, na.rm = TRUE),
    pvalue = stats::median(pvalue, na.rm = TRUE),
    p.adjust = stats::median(p.adjust, na.rm = TRUE),
    qvalue = stats::median(qvalue, na.rm = TRUE),
    Count = stats::median(Count, na.rm = TRUE),
    geneID = str_c(sort(unique(Protein)), collapse = '/'), # actually proteinid
    .groups = "drop"
  ) %>%
  mutate(
    # directionality = n_up - n_down,
    n_protein = lengths(Protein.list),
    label = if_else(generality == 1 & qvalue < 1e-10,
                    Description, NA_character_)
  )

# each organ select top5
selected_goid <- ego_enrich_df %>% 
  semi_join(go_enrich_summary %>% filter(!is.na(label))) %>% 
  group_by(class_abbr) %>% distinct(ID, Description, p.adjust) %>%
  arrange(p.adjust) %>% slice(1:5) %>% pull(ID)
go_enrich_summary %<>% mutate(label = ifelse(ID %in% selected_goid, label, NA))

go_enrich_summary %<>% left_join(
  ego_enrich_df %>%
    filter(Description %in% go_enrich_summary$Description[go_enrich_summary$generality == 1]) %>%
    distinct(Description, class_abbr)
) %>% as.data.frame()

# Uniform jitter position (applied to both points and labels)
pos_jit <- position_jitter(width = 0.4, height = 0.4, seed = 42)

p_go_enrich <- ggplot(go_enrich_summary) +
  aes(
    x     = -log10(p.adjust),
    y     = Specificity,
    fill  = class_abbr,
    size  = Pathway.Size#,
    # alpha = Pathway.Size      
  ) +
  # 
  geom_point(
    alpha = 0.5,
    shape    = 21,             # 
    colour   = "grey25",       # 
    stroke   = 0.25,
    position = pos_jit
  ) +
  # Vertical center dashed line
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    colour     = "grey60"
  ) +
  # Labels: share the same jitter as the scatter points + mandatory short connecting lines
  ggrepel::geom_text_repel(
    aes(label = label),
    position           = pos_jit,   # Align with points
    color              = "black",
    size               = 2.5,
    min.segment.length = 0,         # 
    segment.size       = 0.25,
    segment.color      = "grey30",
    segment.alpha      = 0.8,
    max.overlaps       = Inf,
    seed               = 42         # ggrepel set seed
  ) +
  # size: 拉开差异
  scale_size_continuous(
    name  = "Proteins in GO",
    range = c(0.2, 5)               #
  ) +
  # alpha: 大size → 更透明
  # scale_alpha_continuous(
  #   range = c(0.45, 0.9),           # 0.35(大点) ~ 0.9(小点更实)
  #   trans = "reverse",
  #   guide = "none"                  # 
  # ) +
  # fill: viridis
  # scale_fill_viridis_c(
  #   name      = "-Log10(p.adjust.median)",
  #   end       = 0.8,
  #   direction = -1
  # ) +
  scale_fill_manual(values = class_color) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  guides(
    fill = guide_legend(override.aes = list(size = 4))
  ) +
  labs(
    x        = "-Log10(p.adjust.median)",
    y        = "Specificity (# tissues - # tissues enriched)",
    title    = "Common GO functions among 64 tissues",
    subtitle = "Labeled in unique tissue"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )
ggsave(filename = "N_tissue_enriched_GO_bubble_plot.pdf", plot = p_go_enrich, width = 9, height = 6)


#### v2 plot ------
# remove Count < 3
go_enrich_summary_filter <- ego_enrich_df %>%
  filter(Count >= 3) %>% 
  group_by(ID, Description, Pathway.Size) %>%
  summarise(
    generality = n_distinct(class_abbr),
    Specificity = length(unique(res_ts$ts.all$class_abbr)) - generality,
    # geneID in compareCluster is "/"-delimited character of SYMBOLs per cancer
    Protein.list = list(unique(Protein[!is.na(Protein)])),
    RichFactor = stats::median(RichFactor, na.rm = TRUE),
    FoldEnrichment = stats::median(FoldEnrichment, na.rm = TRUE),
    zScore = stats::median(zScore, na.rm = TRUE),
    pvalue = stats::median(pvalue, na.rm = TRUE),
    p.adjust = stats::median(p.adjust, na.rm = TRUE),
    qvalue = stats::median(qvalue, na.rm = TRUE),
    Count = stats::median(Count, na.rm = TRUE),
    geneID = str_c(sort(unique(Protein)), collapse = '/'), # actually proteinid
    .groups = "drop"
  ) %>%
  mutate(
    # directionality = n_up - n_down,
    n_protein = lengths(Protein.list),
    label = if_else(generality == 1 & qvalue < 1e-10,
                    Description, NA_character_)
  )

# each organ select top5
selected_goid <- ego_enrich_df %>% 
  semi_join(go_enrich_summary_filter %>% filter(!is.na(label))) %>% 
  group_by(class_abbr) %>% distinct(ID, Description, p.adjust) %>%
  arrange(p.adjust) %>% slice(1:5) %>% pull(ID)
go_enrich_summary_filter %<>% mutate(label = ifelse(ID %in% selected_goid, label, NA))

go_enrich_summary_filter %<>% left_join(
  ego_enrich_df %>%
    filter(Description %in% go_enrich_summary_filter$Description[go_enrich_summary_filter$generality == 1]) %>%
    distinct(Description, class_abbr)
) %>% as.data.frame()

# select top5 organs
selected_organ_go <- go_enrich_summary_filter %>%
  filter(generality == 1) %>%
  arrange(qvalue) %>% distinct(class_abbr) %>% pull() %>% head(5)

qvalue.scales <- go_enrich_summary_filter %>%
  filter(generality == 1, class_abbr %in% selected_organ_go) %>%
  arrange(qvalue) %>% group_by(class_abbr) %>% slice(1:5) %>%
  pull(qvalue) %>% quantile()
qvalue.cutoff <- min(1e-10, qvalue.scales['100%']) # 0.001

go_enrich_summary_v2 <- ego_enrich_df %>%
  group_by(ID, Description, Pathway.Size) %>%
  summarise(
    generality = n_distinct(class_abbr),
    Specificity = length(unique(res_ts$ts.all$class_abbr)) - generality,
    Protein.list = list(unique(Protein[!is.na(Protein)])),
    RichFactor = stats::median(RichFactor, na.rm = TRUE),
    FoldEnrichment = stats::median(FoldEnrichment, na.rm = TRUE),
    zScore = stats::median(zScore, na.rm = TRUE),
    pvalue = stats::median(pvalue, na.rm = TRUE),
    p.adjust = stats::median(p.adjust, na.rm = TRUE),
    qvalue = stats::median(qvalue, na.rm = TRUE),
    Count = stats::median(Count, na.rm = TRUE),
    geneID = str_c(sort(unique(Protein)), collapse = '/'), # actually proteinid
    .groups = "drop"
  ) %>%
  mutate(
    # directionality = n_up - n_down,
    n_protein = lengths(Protein.list),
    label = if_else(generality == 1 & qvalue < qvalue.cutoff,
                    Description, NA_character_)
  )

# each organ select top5
selected_goid_v2 <- ego_enrich_df %>% 
  semi_join(go_enrich_summary_v2 %>% filter(!is.na(label))) %>% 
  group_by(class_abbr) %>% distinct(ID, Description, p.adjust) %>%
  arrange(p.adjust) %>% slice(1:5) %>% pull(ID)
go_enrich_summary_v2 %<>% mutate(label = ifelse(ID %in% selected_goid_v2, label, NA))

go_enrich_summary_v2 %<>% left_join(
  ego_enrich_df %>%
    filter(Description %in% go_enrich_summary_v2$Description[go_enrich_summary_v2$generality == 1]) %>%
    distinct(Description, class_abbr)
) %>%
  mutate(class_abbr = factor(class_abbr, selected_organ_go)) %>% 
  as.data.frame()



data_colored <- go_enrich_summary_v2 %>%
  filter(class_abbr %in% selected_organ_go)
data_grey <- go_enrich_summary_v2 %>% anti_join(data_colored) %>% 
  mutate(class_abbr = 'Others')
data <- data_colored %>% rbind(data_grey)

# class_color_v2 <- c(class_color[61:64], 'grey') %>% setNames(c(unique(data_colored$class_abbr), 'Others'))
# class_color_v2 <- c("#E64B35", "#2B2B6B", "#00A087", "#D8894E", "grey") %>% setNames(c(unique(data_colored$class_abbr), 'Others'))
class_color_v2 <- c(class_color[as.character(unique(data_colored$class_abbr))], '#EEEEEE') %>%
  setNames(c(as.character(unique(data_colored$class_abbr)), 'Others'))

p_go_enrich <- ggplot(data) +
  aes(
    x     = -log10(p.adjust),
    y     = Specificity,
    fill  = class_abbr,
    size  = Pathway.Size#,
    # alpha = Pathway.Size       # 
  ) +
  # 
  geom_point( # for legend sequence
    alpha = 0.8,
    shape    = NA,
    size = 0,
    colour   = 'black',
    stroke   = 0.25,
    data = data
  ) +
  geom_point(
    alpha = 0.5,
    shape    = 21,
    colour   = "#EEEEEE",
    stroke   = 0.25,
    position = pos_jit,
    data = data %>% semi_join(data_grey)
  ) +
  geom_point(
    alpha = 0.8,
    shape    = 21,
    colour   = 'black',
    stroke   = 0.25,
    position = pos_jit,
    data = data %>% semi_join(data_colored)
  ) +
  # # 
  # geom_vline(
  #   xintercept = 0,
  #   linetype   = "dashed",
  #   colour     = "grey60"
  # ) +
  # 
  ggrepel::geom_text_repel(
    aes(label = label), data = data %>% semi_join(data_colored),
    position           = pos_jit,   # 
    color              = "black",
    size               = 2,
    min.segment.length = 0,         # 
    segment.size       = 0.25,
    segment.color      = "grey30",
    segment.alpha      = 0.8,
    max.overlaps       = Inf,
    seed               = 42         # 
  ) +
  #
  scale_size_continuous(
    name  = "Proteins in GOBP",
    range = c(0.5, 6)               # 
  ) +
  
  # scale_alpha_continuous(
  #   range = c(0.45, 0.9),           # 
  #   trans = "reverse",
  #   guide = "none"                  # 
  # ) +
  # fill: viridis
  # scale_fill_viridis_c(
  #   name      = "-Log10(p.adjust.median)",
  #   end       = 0.8,
  #   direction = -1
  # ) +
  # scale_fill_manual(values = c(class_color[unique(data_colored$class_abbr)], 'Others' = 'grey')) +
  scale_fill_manual(values = class_color_v2) +
  guides(fill = guide_legend("Tissue", order = 1, override.aes = list(size = 6))) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    x        = "-Log10(p.adjust.median)",
    y        = "Specificity (# tissues - # tissues enriched)",
    title    = "Common GOBP functions among 64 tissues",
    subtitle = "Labeled in unique tissue"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )
ggsave(filename = "N_tissue_enriched_GOBP_bubble_plot_20251202_v2.pdf", plot = p_go_enrich, width = 9, height = 6)


#### Network GOBP -----
library(aPEAR)

go_enrich_summary_filter %>% arrange(pvalue) %>% pull(class_abbr) %>% unique()
# [1] "LI"    "CE"    NA      "EYE-l" "TE"    "HEA"   "MU-sk" "ADG"   "GB"    "THY"   "LU"    "SI"    "EAR-e"
# [14] "EYE-s" "BM"    "CL"    "TDN"   "OV"    "TOS"   "EAR-t" "EYE-c" "EAR-s" "EYE-i" "KI"    "EAR-c" "EPIG" 
# [27] "SG"    "CG"    "MAG"   "PA"    "PR"    "SPI"   "TON"   "WEA"   "ST"    "EAR-o"
organs <- selected_organ_go
size.scales <- go_enrich_summary_filter %>% filter(generality == 1, class_abbr %in% organs) %>% pull(Pathway.Size) %>% quantile(c(0, 1))

# network analysis
enw.list.jaccard.hier <- lapply(organs %>% setNames(organs), function(organ){
  cat('Analysing ', organ, '...\n')
  go_id <- go_enrich_summary_filter %>%
    filter(generality == 1, class_abbr == organ) %>% pull(ID)
  minClusterSize <- min(10, ceiling(length(go_id) / 3))
  my.enrich.network(enrich.table = go_enrich_summary_filter %>%
                      filter(class_abbr %in% organ, ID %in% go_id,
                             p.adjust < 1e-2),
                    size.scales = size.scales, simMethod = 'jaccard', clustMethod = 'hier', minClusterSize = minClusterSize)
})
# saveRDS(enw.list.jaccard.hier, 'enw.list.jaccard.hier.rds')

organ.labels <- c('LI', 'CE', 'TE', 'EYE-l', 'HEA')
plot.network.jaccard.hier <- ggpubr::ggarrange(
  plotlist = lapply(enw.list.jaccard.hier, function(x) x$plot)[organ.labels],
  labels = organ.labels, common.legend = TRUE,
  ncol = 3, nrow = 2, heights = c(2, 0.5)
)
ggsave('N_tissue_enriched_GOBP_network_plot_v1202.pdf', plot.network.jaccard.hier, width = 5*3, height = 5*2)

enrich.network <- plyr::ldply(enw.list.jaccard.hier,
                              function(x) x$clusters, .id = 'class_abbr') %>%
  select(class_abbr, Cluster, ID) %>% 
  rename(Description = ID) %>% 
  # mutate(Cluster = str_remove(Cluster, '\\.1$')) %>% 
  inner_join(go_enrich_summary_filter) %>% 
  mutate(class_abbr = ifelse(class_abbr == 'CE', 'BRAIN', class_abbr)) %>% 
  arrange(class_abbr, Cluster, qvalue) %>% 
  select(-label) %>%
  group_by(class_abbr, Cluster) %>% # add new.cluster labels
  mutate(p.adjust.median = median(p.adjust),
         Cluster.new = Description[ which.min(qvalue) ],
         .before = Cluster) %>%
  ungroup() %>% 
  arrange(class_abbr, p.adjust.median, p.adjust)
# enrich.network %>% filter(str_detect(Cluster, 'postsynapse assembly'))

list(
  GOBP = ego_enrich[] %>% 
    mutate(Cluster = as.character(Cluster),
           Cluster = ifelse(Cluster == 'CE', 'BRAIN', Cluster)) %>% 
    arrange(Cluster),
  GOBP.summary = go_enrich_summary_filter %>%
    mutate(class_abbr = ifelse(class_abbr == 'CE', 'BRAIN', class_abbr)) %>% 
    group_by(Protein.list) %>%
    mutate(UniprotID = str_c(unlist(Protein.list), collapse = '/'),
           .after = geneID) %>% ungroup() %>%
    arrange(qvalue) %>% 
    select(-Protein.list, -label),
  enrich.network = enrich.network
) %>% 
  .write_excel('N_tissue_enriched_GOBP_network_plot_v1202.xlsx')


### kegg ----
ekegg_enrich_df <- ekegg_enrich %>%
  filter(pvalue < 0.05, qvalue < 0.05,
         # !category %in% c('Human Diseases'#, 'Organismal Systems'),
         !subcategory %in% c('Infectious disease: viral')
  ) %>% # in fact redundant
  as.data.frame() %>%
  separate_rows(geneID, sep = '/') %>% 
  rename(class_abbr = Cluster,
         proteinid = geneID) %>% 
  mutate(Pathway.Size = as.numeric(str_extract(BgRatio, '^\\d+'))) %>% 
  inner_join(df_tis_en1 %>% distinct(class_abbr, Protein, proteinid), .)

kegg_enrich_summary <- ekegg_enrich_df %>%
  group_by(ID, Description, Pathway.Size) %>%
  summarise(
    generality = n_distinct(class_abbr),
    Specificity = length(unique(res_ts$ts.all$class_abbr)) - generality,
    Protein.list = list(unique(Protein[!is.na(Protein)])),
    RichFactor = stats::median(RichFactor, na.rm = TRUE),
    FoldEnrichment = stats::median(FoldEnrichment, na.rm = TRUE),
    zScore = stats::median(zScore, na.rm = TRUE),
    pvalue = stats::median(pvalue, na.rm = TRUE),
    p.adjust = stats::median(p.adjust, na.rm = TRUE),
    qvalue = stats::median(qvalue, na.rm = TRUE),
    Count = stats::median(Count, na.rm = TRUE),
    geneID = str_c(sort(unique(proteinid)), collapse = '/'), # actually proteinid
    .groups = "drop"
  ) %>%
  mutate(
    # directionality = n_up - n_down,
    n_protein = lengths(Protein.list),
    label = if_else(generality == 1 & qvalue < 1e-10,
                    Description, NA_character_)
  ) %>% 
  as.data.frame()
kegg_enrich_summary %<>% left_join(
  ekegg_enrich_df %>%
    filter(Description %in% kegg_enrich_summary$Description[kegg_enrich_summary$generality == 1]) %>%
    distinct(category, subcategory, Description, class_abbr)
)

# 
pos_jit <- position_jitter(width = 0.4, height = 0.4, seed = 42)

p_kegg_enrich <- ggplot(kegg_enrich_summary) +
  aes(
    x     = -log10(p.adjust),
    y     = Specificity,
    fill  = class_abbr,
    size  = Pathway.Size#,
    # alpha = Pathway.Size      
  ) +
  # 
  geom_point(
    alpha = 0.5,
    shape    = 21,            
    colour   = "grey25",       
    stroke   = 0.25,
    position = pos_jit
  ) +
  # 
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    colour     = "grey60"
  ) +
  # 
  ggrepel::geom_text_repel(
    aes(label = label),
    position           = pos_jit,   # 
    color              = "black",
    size               = 2.5,
    min.segment.length = 0,         # 
    segment.size       = 0.25,
    segment.color      = "grey30",
    segment.alpha      = 0.8,
    max.overlaps       = Inf,
    seed               = 42         # 
  ) +
  #
  scale_size_continuous(
    name  = "Proteins in KEGG",
    range = c(0.5, 6)               # 
  ) +
  # 
  # scale_alpha_continuous(
  #   range = c(0.45, 0.9),           
  #   trans = "reverse",
  #   guide = "none"                  
  # ) +
  # fill: viridis
  # scale_fill_viridis_c(
  #   name      = "-Log10(p.adjust.median)",
  #   end       = 0.8,
  #   direction = -1
  # ) +
  scale_fill_manual(values = class_color) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  guides(
    fill = guide_legend(override.aes = list(size = 4))
  ) +
  labs(
    x        = "-Log10(p.adjust.median)",
    y        = "Specificity (# tissues - # tissues enriched)",
    title    = "Common KEGG functions among 64 tissues",
    subtitle = "Labeled in unique tissue"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )
ggsave(filename = "N_tissue_enriched_KEGG_bubble_plot_v1202.pdf", plot = p_kegg_enrich, width = 9, height = 6)
rio::export(kegg_enrich_summary, 'kegg_enrich_summary_v1202.xlsx')
saveRDS(kegg_enrich_summary, 'kegg_enrich_summary_v1202.rds')

#### v2 plot ------
data_colored <- kegg_enrich_summary %>%
  filter(class_abbr %in% c('LI', 'CE', 'ADG', 'MU-sk'))
data_grey <- kegg_enrich_summary %>% anti_join(data_colored) %>% 
  mutate(class_abbr = 'Others')
data <- data_colored %>% rbind(data_grey)

# class_color_v2 <- c(class_color[61:64], 'grey') %>% setNames(c(unique(data_colored$class_abbr), 'Others'))
class_color_v2 <- c("#E64B35", "#2B2B6B", "#00A087", "#D8894E", "grey") %>% setNames(c(unique(data_colored$class_abbr), 'Others'))

p_kegg_enrich <- ggplot(data) +
  aes(
    x     = -log10(p.adjust),
    y     = Specificity,
    fill  = class_abbr,
    size  = Pathway.Size#,
    # alpha = Pathway.Size      
  ) +
  
  geom_point(
    alpha = 0.8,
    shape    = 21,            
    colour   = "grey25",       
    stroke   = 0.25,
    position = pos_jit,
    data = data
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    colour     = "grey60"
  ) +
  
  ggrepel::geom_text_repel(
    aes(label = label),
    position           = pos_jit,   
    color              = "black",
    size               = 2,
    min.segment.length = 0,         
    segment.size       = 0.25,
    segment.color      = "grey30",
    segment.alpha      = 0.8,
    max.overlaps       = Inf,
    seed               = 42        
  ) +
  
  scale_size_continuous(
    name  = "Proteins in KEGG",
    range = c(0.5, 6)               
  ) +
  
  # scale_alpha_continuous(
  #   range = c(0.45, 0.9),           
  #   trans = "reverse",
  #   guide = "none"                 
  # ) +
  # fill: viridis
  # scale_fill_viridis_c(
  #   name      = "-Log10(p.adjust.median)",
  #   end       = 0.8,
  #   direction = -1
  # ) +
  # scale_fill_manual(values = c(class_color[unique(data_colored$class_abbr)], 'Others' = 'grey')) +
  scale_fill_manual(values = class_color_v2) +
  guides(fill = guide_legend("Tissue", order = 1, override.aes = list(size = 6))) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    x        = "-Log10(p.adjust.median)",
    y        = "Specificity (# tissues - # tissues enriched)",
    title    = "Common KEGG functions among 64 tissues",
    subtitle = "Labeled in unique tissue"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )
ggsave(filename = "N_tissue_enriched_KEGG_bubble_plot_20251202_v2.pdf", plot = p_kegg_enrich, width = 9, height = 6)


##### zoom in ---------
# p_kegg_enrich <- ggplot(data) +
#   aes(
#     x     = -log10(med_padj),
#     y     = Specificity,
#     fill  = class_abbr,
#     size  = Pathway.Size#,
#     # alpha = Pathway.Size      
#   ) +
#   
#   geom_point(
#     alpha = 0.8,
#     shape    = 21,             
#     colour   = "grey25",      
#     stroke   = 0.25,
#     position = pos_jit,
#     data = data
#   ) +
#   
#   geom_vline(
#     xintercept = 0,
#     linetype   = "dashed",
#     colour     = "grey60"
#   ) +
#   
#   ggrepel::geom_text_repel(
#     aes(label = label),
#     position           = pos_jit,  
#     color              = "black",
#     size               = 2,
#     min.segment.length = 0,         
#     segment.size       = 0.25,
#     segment.color      = "grey30",
#     segment.alpha      = 0.8,
#     max.overlaps       = Inf,
#     seed               = 42         
#   ) +
#   
#   scale_size_continuous(
#     name  = "Proteins in KEGG",
#     range = c(0.5, 6)               
#   ) +
#   
#   # scale_alpha_continuous(
#   #   range = c(0.45, 0.9),           
#   #   trans = "reverse",
#   #   guide = "none"                  
#   # ) +
#   # fill: viridis
#   # scale_fill_viridis_c(
#   #   name      = "-Log10(p.adjust.median)",
#   #   end       = 0.8,
#   #   direction = -1
#   # ) +
#   # scale_fill_manual(values = c(class_color[unique(data_colored$class_abbr)], 'Others' = 'grey')) +
#   scale_fill_manual(values = c(class_color[61:64], 'grey') %>% setNames(c(unique(data_colored$class_abbr), 'Others'))) +
#   guides(fill = guide_legend("Tissue", order = 1, override.aes = list(size = 6))) +
#   scale_x_continuous(breaks = scales::pretty_breaks(), limits = c(8, NA)) +
#   scale_y_continuous(breaks = c(63), limits = c(62.5, NA)) +
#   labs(
#     x        = "-Log10(p.adjust.median)",
#     y        = "Specificity (# tissues - # tissues enriched)",
#     title    = "Common KEGG functions among 64 tissues",
#     subtitle = "Labeled in unique tissue"
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(face = "bold")
#   )
# ggsave(filename = "N_tissue_enriched_KEGG_bubble_plot_v2_zoomIn.pdf", plot = p_kegg_enrich, width = 6, height = 3)

#### Network KEGG -----

library(aPEAR)

# kegg_enrich_summary %>% arrange(pvalue) %>% pull(class_abbr) %>% unique()
# [1] NA      "LI"    "CE"    "MU-sk" "ADG"   "SG"    "PA"    "THY"   "TDN"   "PTH"   "EYE-l" "BLA"   "SI"    "EYE-c" "TON"  
# [16] "EPIDI" "TOS"   "CL"    "WEA"   "ST"    "LN"    "BU"    "MAG"   "EAR-t" "VA"    "CG"  
organs <- c('LI', 'CE', 'MU-sk', 'ADG')

kegg_id <- kegg_enrich_summary %>% filter(generality == 1, class_abbr == 'LI') %>% pull(ID)
set.seed(1)
enw1 <- enrichmentNetwork(
  kegg_enrich_summary %>% filter(ID %in% kegg_id), plotOnly = FALSE,
  clustMethod = 'hier', # method for detecting pathway clusters. Available methods: 'markov', 'hier' and 'spectral'. Using 'spectral' method requires that you have Spectrum package installed.
  simMethod = 'cosine', # method for calculating similarity matrix between the pathways. Available methods: 'jaccard', 'cosine' and 'cor'
  clustNameMethod = 'pagerank', # method for selecting cluster names. Available methods: 'pagerank', 'hits' and 'none'
  drawEllipses = TRUE, fontSize = 4, nodeSize = 'Pathway.Size',
  colorBy = 'p.adjust', colorType = 'pval',# pCutoff = -5,
  minClusterSize = 2
)
p1 <- enw1$plot +
  scale_color_viridis_c(name = 'Log10(p.adjust)', option = 'C', end = 0.8, limits = c(-10, -2)) +
  scale_size_continuous(name = "Pathway size", range = c(2, 6), limits = c(10, 370))


kegg_id <- kegg_enrich_summary %>% filter(generality == 1, class_abbr == 'CE') %>% pull(ID)
# pheatmap(aPEAR::pathwaySimilarity(
#   kegg_enrich_summary %>% filter(ID %in% kegg_id), 'geneID',
#   method = 'cosine'
# ))
set.seed(1)
enw2 <- enrichmentNetwork(
  kegg_enrich_summary %>% filter(ID %in% kegg_id), plotOnly = FALSE,
  clustMethod = 'hier', # method for detecting pathway clusters. Available methods: 'markov', 'hier' and 'spectral'. Using 'spectral' method requires that you have Spectrum package installed.
  simMethod = 'cosine', # method for calculating similarity matrix between the pathways. Available methods: 'jaccard', 'cosine' and 'cor'
  clustNameMethod = 'pagerank', # method for selecting cluster names. Available methods: 'pagerank', 'hits' and 'none'
  drawEllipses = TRUE, fontSize = 4, nodeSize = 'Pathway.Size',
  colorBy = 'p.adjust', colorType = 'pval',# pCutoff = -5,
  minClusterSize = 2
)
p2 <- enw2$plot +
  scale_color_viridis_c(name = 'Log10(p.adjust)', option = 'C', end = 0.8, limits = c(-10, -2)) +
  scale_size_continuous(name = "Pathway size", range = c(2, 6), limits = c(10, 370))


kegg_id <- kegg_enrich_summary %>% filter(generality == 1, class_abbr == 'MU-sk') %>% pull(ID)
set.seed(1)
enw3 <- enrichmentNetwork(
  kegg_enrich_summary %>% filter(ID %in% kegg_id), plotOnly = FALSE,
  clustMethod = 'hier', # method for detecting pathway clusters. Available methods: 'markov', 'hier' and 'spectral'. Using 'spectral' method requires that you have Spectrum package installed.
  simMethod = 'cosine', # method for calculating similarity matrix between the pathways. Available methods: 'jaccard', 'cosine' and 'cor'
  clustNameMethod = 'pagerank', # method for selecting cluster names. Available methods: 'pagerank', 'hits' and 'none'
  drawEllipses = TRUE, fontSize = 4, nodeSize = 'Pathway.Size',
  colorBy = 'p.adjust', colorType = 'pval',# pCutoff = -5,
  minClusterSize = 2
)
p3 <- enw3$plot +
  scale_color_viridis_c(name = 'Log10(p.adjust)', option = 'C', end = 0.8, limits = c(-10, -2)) +
  scale_size_continuous(name = "Pathway size", range = c(2, 6), limits = c(10, 370))


kegg_id <- kegg_enrich_summary %>% filter(generality == 1, class_abbr == 'ADG') %>% pull(ID)
set.seed(1)
enw4 <- enrichmentNetwork(
  kegg_enrich_summary %>% filter(ID %in% kegg_id), plotOnly = FALSE,
  clustMethod = 'hier', # method for detecting pathway clusters. Available methods: 'markov', 'hier' and 'spectral'. Using 'spectral' method requires that you have Spectrum package installed.
  simMethod = 'cosine', # method for calculating similarity matrix between the pathways. Available methods: 'jaccard', 'cosine' and 'cor'
  clustNameMethod = 'pagerank', # method for selecting cluster names. Available methods: 'pagerank', 'hits' and 'none'
  drawEllipses = TRUE, fontSize = 4, nodeSize = 'Pathway.Size',
  colorBy = 'p.adjust', colorType = 'pval',# pCutoff = -5,
  minClusterSize = 2
)
p4 <- enw4$plot +
  scale_color_viridis_c(name = 'Log10(p.adjust)', option = 'C', end = 0.8, limits = c(-10, -2)) +
  scale_size_continuous(name = "Pathway size", range = c(2, 6), limits = c(10, 370))


plot.network <- ggpubr::ggarrange(p1, p2, p3, p4, labels = organs, common.legend = TRUE)
ggsave('N_tissue_enriched_KEGG_network_plot_20251119.pdf', plot.network, width = 10, height = 10)

enrich.network <- list(enw1$clusters, enw2$clusters, enw3$clusters, enw4$clusters) %>%
  setNames(organs) %>% 
  plyr::ldply(.id = 'class_abbr') %>% 
  select(class_abbr, Cluster, ID) %>% 
  rename(Description = ID) %>% 
  inner_join(kegg_enrich_summary) %>% 
  mutate(class_abbr = ifelse(class_abbr == 'CE', 'BRAIN', class_abbr)) %>% 
  arrange(class_abbr, Cluster, qvalue) %>% 
  select(-label)

list(
  KEGG = ekegg_enrich[] %>% 
    mutate(Cluster = as.character(Cluster),
           Cluster = ifelse(Cluster == 'CE', 'BRAIN', Cluster)) %>% 
    arrange(Cluster),
  KEGG.summary = kegg_enrich_summary %>%
    mutate(class_abbr = ifelse(class_abbr == 'CE', 'BRAIN', class_abbr)) %>% 
    group_by(Protein.list) %>%
    mutate(UniprotID = str_c(unlist(Protein.list), collapse = '/'),
           .after = geneID) %>% ungroup() %>%
    arrange(qvalue) %>% 
    select(-Protein.list, -label),
  enrich.network = enrich.network
) %>% 
  .write_excel('N_tissue_enriched_KEGG_network_plot_20251119.xlsx')



### output ------
# list(GOBP = go_enrich_summary %>% select(-Protein.list),
#      KEGG = kegg_enrich_summary %>% select(-Protein.list)) %>% 
#   rio::export('N_tissue_enriched_bubble_plot.xlsx')
# 
# list(GOBP = go_enrich_summary,
#      KEGG = kegg_enrich_summary) %>% 
#   saveRDS('N_tissue_enriched_bubble_plot.rds')
# # Enrich.common.pathway <- readRDS('common_Enrich_bubble_plot.rds')


# 3. HPA RNA tissue-enriched genes -----
hpa_rna_raw <- rio::import('input/proteinatlas_0ac3e69c.tsv')

hpa_tis_en <- hpa_rna_raw %>%
  filter(`RNA tissue specificity` == "Tissue enriched") %>%
  dplyr::select(Gene, `Gene synonym`, Uniprot, contains('RNA tissue')) %>% 
  mutate(HPA.enriched.tissue = `RNA tissue specific nTPM` %>% 
           str_remove(': [\\d\\.]+$') %>% 
           str_remove(' 1$'))
# writeClipboard(unique(sort(hpa_tis_en$HPA.enriched.tissue)))
# writeClipboard(unique(sort(df_tis_en_$anatomical_classification)))

# enriched_uniprot_hpa <- hpa_tis_en$Uniprot %>% str_split(', ') %>%
#   unlist() %>% unique()
# 
# tmp <- clusterProfiler::bitr_kegg(enriched_uniprot_hpa, fromType='uniprot', toType='ncbi-proteinid', organism='hsa') %>% 
#   setNames(c('Protein', 'proteinid'))
# # 3.86% of input gene IDs are fail to map...
# hpa_tis_en1 <- hpa_tis_en %>%
#   rename(Protein = Uniprot) %>%
#   separate_rows(Protein, sep = ', ') %>%
#   inner_join(tmp, .)
# tmp1 <- hpa_tis_en1 %>% distinct(HPA.enriched.tissue, proteinid)
# 
# hpa_enrich_pid_list <- split(tmp1$proteinid, tmp1$HPA.enriched.tissue)
# # saveRDS(hpa_enrich_pid_list, 'hpa_enrich_pid_list.rds')
# # enrich_pid_list <- readRDS('hpa_enrich_pid_list.rds')

# id mapping
hpa_id_map <- clusterProfiler::bitr(
  unique(hpa_tis_en$Gene),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

df_hpa_tis_en1 <- hpa_tis_en %>%
  inner_join(hpa_id_map, by = c("Gene" = "SYMBOL")) %>%
  distinct(HPA.enriched.tissue, Gene, ENTREZID)

tmp_rna <- df_hpa_tis_en1 %>%
  distinct(HPA.enriched.tissue, ENTREZID)

enrich_gene_list_rna <- split(tmp_rna$ENTREZID, tmp_rna$HPA.enriched.tissue)
# saveRDS(enrich_gene_list_rna, 'enrich_gene_list_rna.rds')
# enrich_gene_list_rna <- readRDS('enrich_gene_list_rna.rds')

## KEGG ------------
### kegg compareCluster -----
library(clusterProfiler)
library(org.Hs.eg.db)

# hpa_ekegg_enrich <- compareCluster(hpa_enrich_pid_list, fun = enrichKEGG,
#                                pAdjustMethod = "BH",
#                                keyType = 'ncbi-proteinid',
#                                minGSSize = 10, maxGSSize = 500,
#                                pvalueCutoff = 0.05,
#                                qvalueCutoff = 0.2)
# # saveRDS(hpa_ekegg_enrich, 'hpa_Enrich_ekegg_compare.rds')
# # hpa_ekegg_enrich <- readRDS('hpa_Enrich_ekegg_compare.rds')
ekegg_enrich_rna <- compareCluster(
  enrich_gene_list_rna,
  fun           = enrichKEGG,
  keyType       = "ncbi-geneid",
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
# saveRDS(ekegg_enrich_rna, 'ekegg_enrich_rna_compare.rds')
# ekegg_enrich_rna <- readRDS('ekegg_enrich_rna_compare.rds')
enrichplot::dotplot(ekegg_enrich_rna, showCategory = 1) +
  theme(axis.text.x = element_text(angle = 45))


pacman::p_unload(clusterProfiler, enrichplot, ReactomePA, AnnotationDbi, BiocGenerics, IRanges, S4Vectors, Biobase)

### kegg bubbles -----
ekegg_enrich_rna_df <- ekegg_enrich_rna %>%
  filter(pvalue < 0.05, qvalue < 0.05,
         # !category %in% c('Human Diseases'#, 'Organismal Systems'),
         !subcategory %in% c('Infectious disease: viral')
  ) %>% # in fact redundant
  as.data.frame() %>%
  separate_rows(geneID, sep = '/') %>% 
  rename(HPA.enriched.tissue = Cluster,
         ENTREZID = geneID) %>% 
  mutate(Pathway.Size = as.numeric(str_extract(BgRatio, '^\\d+'))) %>% 
  inner_join(df_hpa_tis_en1 %>% distinct(HPA.enriched.tissue, Gene, ENTREZID), .)

kegg_enrich_rna_summary <- ekegg_enrich_rna_df %>%
  group_by(ID, Description, Pathway.Size) %>%
  summarise(
    generality = n_distinct(HPA.enriched.tissue),
    Specificity = length(unique(ekegg_enrich_rna_df$HPA.enriched.tissue)) - generality,
    Gene.list = list(unique(Gene[!is.na(Gene)])),
    RichFactor = stats::median(RichFactor, na.rm = TRUE),
    FoldEnrichment = stats::median(FoldEnrichment, na.rm = TRUE),
    zScore = stats::median(zScore, na.rm = TRUE),
    pvalue = stats::median(pvalue, na.rm = TRUE),
    p.adjust = stats::median(p.adjust, na.rm = TRUE),
    qvalue = stats::median(qvalue, na.rm = TRUE),
    Count = stats::median(Count, na.rm = TRUE),
    geneID = str_c(sort(unique(ENTREZID)), collapse = '/'),
    .groups = "drop"
  ) %>%
  mutate(
    # directionality = n_up - n_down,
    n_protein = lengths(Gene.list),
    label = if_else(generality == 1 & p.adjust < 1e-8,
                    Description, NA_character_)
  ) %>% 
  as.data.frame()
kegg_enrich_rna_summary %<>% left_join(
  ekegg_enrich_rna_df %>%
    filter(Description %in% kegg_enrich_rna_summary$Description[kegg_enrich_rna_summary$generality == 1]) %>%
    distinct(category, subcategory, Description, HPA.enriched.tissue), .
)

pos_jit <- position_jitter(width = 0.4, height = 0.4, seed = 42)

rna_class_color <- class_color[1:length(unique(ekegg_enrich_rna_df$HPA.enriched.tissue))] %>%
  setNames(sort(unique(ekegg_enrich_rna_df$HPA.enriched.tissue)))
p_kegg_enrich_rna <- ggplot(kegg_enrich_rna_summary) +
  aes(x = -log10(p.adjust), y = Specificity, fill = HPA.enriched.tissue, size = Pathway.Size) +
  geom_point(alpha = 0.5, shape = 21, colour = "grey25", stroke = 0.25, position = pos_jit) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  ggrepel::geom_text_repel(aes(label = label), position = pos_jit, color = "black", size = 2,
                           min.segment.length = 0, segment.size = 0.25, segment.color = "grey30",
                           segment.alpha = 0.8, max.overlaps = Inf, seed = 42) +
  scale_size_continuous(name = "Proteins in KEGG", range = c(0.5, 6)) +
  scale_fill_manual(values = rna_class_color) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(x = "-Log10(p.adjust.median)",
       y = "Specificity (# tissues - # tissues enriched)",
       title = "Common KEGG functions among 27 tissues",
       subtitle = "Labeled in unique tissue") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "HPA_tissue_enriched_KEGG_bubble_plot_20251119.pdf",
       plot = p_kegg_enrich_rna, width = 9, height = 6)
.write_excel(kegg_enrich_rna_summary, 'kegg_enrich_rna_summary_20251119.xlsx')
saveRDS(kegg_enrich_rna_summary, 'kegg_enrich_rna_summary_20251119.rds')

#### v2 plot ------
rna_data_colored <- kegg_enrich_rna_summary %>%
  filter(HPA.enriched.tissue %in% c('liver', 'retina', 'salivary gland', 'adrenal gland', 'lymphoid tissue', 'skin', 'heart muscle'))
rna_data_grey <- kegg_enrich_rna_summary %>% anti_join(rna_data_colored) %>% 
  mutate(HPA.enriched.tissue = 'Others')
rna_data <- rna_data_colored %>% rbind(rna_data_grey)

rna_class_color_v2 <- rna_class_color[unique(rna_data_colored$HPA.enriched.tissue)] %>%
  append(c(Others = "grey"))

p_kegg_enrich_rna <- ggplot(rna_data) +
  aes(x = -log10(p.adjust), y = Specificity, fill = HPA.enriched.tissue, size = Pathway.Size) +
  geom_point(alpha = 0.8, shape = 21, colour = "grey25", stroke = 0.25, position = pos_jit, data = rna_data) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  ggrepel::geom_text_repel(aes(label = label), position = pos_jit, color = "black", size = 2,
                           min.segment.length = 0, segment.size = 0.25, segment.color = "grey30",
                           segment.alpha = 0.8, max.overlaps = Inf, seed = 42) +
  scale_size_continuous(name = "Proteins in KEGG", range = c(0.5, 6)) +
  scale_fill_manual(values = rna_class_color_v2) +
  guides(fill = guide_legend("Tissue", order = 1, override.aes = list(size = 6))) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(x = "-Log10(p.adjust.median)",
       y = "Specificity (# tissues - # tissues enriched)",
       title = "Common KEGG functions among 27 tissues",
       subtitle = "Labeled in unique tissue") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "HPA_tissue_enriched_KEGG_bubble_plot_20251119_v2.pdf",
       plot = p_kegg_enrich_rna, width = 9, height = 6)



# 4.Aggregate HPA and TPHP ------
## 4.1 data prepare ------
tissue_id_mapping <- rio::import('input/tphp_hpa_tissue_match.xlsx')
tissue_id_mapping1 <- tissue_id_mapping %>%
  rename(anatomical_classification = TPHP.enriched.tissue) %>% 
  left_join(info5 %>% distinct(anatomical_classification, class_abbr)) %>%
  mutate(class_abbr = ifelse(is.na(class_abbr), HPA.enriched.tissue, class_abbr))
kegg_enrich_summary <- readRDS('kegg_enrich_summary_20251119.rds')
kegg_enrich_rna_summary <- readRDS('kegg_enrich_rna_summary_20251119.rds')

# format FEATURE.list
kegg_enrich_summary$Protein.list <- sapply(kegg_enrich_summary$Protein.list, function(prots){
  str_c(prots, collapse = ', ')
})
kegg_enrich_rna_summary$Gene.list <- sapply(kegg_enrich_rna_summary$Gene.list, function(genes){
  str_c(genes, collapse = ', ')
})

# Add Protein.list in HPA-table
hpa_gene2prot <- hpa_tis_en %>% 
  mutate(Uniprot = ifelse(Uniprot == '', '-', Uniprot)) %>%
  pull(Uniprot, Gene)
kegg_enrich_rna_summary$Protein.list <- sapply(kegg_enrich_rna_summary$Gene.list, function(genes){
  vec_genes <- str_split(genes, ', ')[[1]]
  str_c(hpa_gene2prot[vec_genes], collapse = ', ')
})

# Aggregate two tables
summ1 <- kegg_enrich_rna_summary %>%
  filter(generality == 1) %>% 
  right_join(tissue_id_mapping1, .) %>%
  select(class_abbr, category, subcategory, Description, ID,
         Pathway.Size, Protein.list, geneID, RichFactor:Count) %>% 
  mutate(LevelSrc = 'RNA', .before = 1)
summ2 <- kegg_enrich_summary %>%
  filter(generality == 1) %>% 
  right_join(info5 %>% distinct(anatomical_classification, class_abbr), .) %>%
  select(class_abbr, category, subcategory, Description, ID,
         Pathway.Size, Protein.list, geneID, RichFactor:Count) %>% 
  mutate(LevelSrc = 'Protein', .before = 1)
kegg_enrich_agg <- bind_rows(summ1, summ2) %>%
  group_by(class_abbr, category, subcategory, Description, ID, Pathway.Size) %>%
  summarise(
    Level = case_when(
      n_distinct(LevelSrc) == 2 ~ "Both",
      first(LevelSrc) == "Protein" ~ "Protein",
      TRUE ~ "RNA"
    ),
    RichFactor      = median(RichFactor,      na.rm = TRUE),
    FoldEnrichment  = median(FoldEnrichment,  na.rm = TRUE),
    zScore          = median(zScore,          na.rm = TRUE),
    pvalue          = median(pvalue,          na.rm = TRUE),
    p.adjust        = median(p.adjust,        na.rm = TRUE),
    qvalue          = median(qvalue,          na.rm = TRUE),
    Count           = median(Count,           na.rm = TRUE),
    geneID          = str_c(
      str_c(unique(geneID[LevelSrc == "Protein"]), collapse = ","),
      str_c(unique(geneID[LevelSrc == "RNA"]),     collapse = ","),
      sep = ".//."
    ),
    Protein.list1   = first(Protein.list[LevelSrc == "Protein"]),
    Protein.list2   = first(Protein.list[LevelSrc == "RNA"]),
    .groups = "drop"
  ) %>%
  mutate(
    Protein.list.shared = map2_chr(Protein.list1, Protein.list2, ~{
      p1 <- if (!is.na(.x)) str_split(.x, ",\\s*")[[1]] else character()
      p2 <- if (!is.na(.y)) str_split(.y, ",\\s*")[[1]] else character()
      paste(intersect(p1, p2), collapse = ", ")
    })
  ) %>% 
  left_join(tissue_id_mapping1) %>% 
  mutate(class_abbr = ifelse(class_abbr == 'CE', 'BRAIN', class_abbr)) %>% 
  as.data.frame()

# saveRDS(kegg_enrich_agg, 'kegg_enrich_agg.rds')


## 4.2 network ------
quantile(kegg_enrich_agg$Pathway.Size)
# 0%  25%  50%  75% 100% 
# 10   42   76  122  483 

## 4.2.1 choose a seed tissue (class_abbr) and pathway set ----------
vec_class    <- sort(unique(kegg_enrich_agg$class_abbr))
target_class <- vec_class[4]   # e.g. "ADG"; change if needed

kegg_id <- kegg_enrich_agg %>%
  filter(class_abbr == target_class) %>%
  pull(ID) %>%
  unique()

kegg_net_tbl <- kegg_enrich_agg %>%
  filter(class_abbr == target_class, ID %in% kegg_id) %>%
  mutate(
    n_Protein = ifelse(Level == 'Protein', 1, 0),
    n_RNA = ifelse(Level == 'RNA', 1, 0),
    n_Both = ifelse(Level == 'Both', 1, 0),
  ) %>% 
  arrange(p.adjust) %>% 
  as.data.frame()

set.seed(1)
enw1 <- enrichmentNetwork(
  kegg_net_tbl,
  plotOnly        = FALSE,
  clustMethod = 'hier', # method for detecting pathway clusters. Available methods: 'markov', 'hier' and 'spectral'. Using 'spectral' method requires that you have Spectrum package installed.
  simMethod = 'cosine', # method for calculating similarity matrix between the pathways. Available methods: 'jaccard', 'cosine' and 'cor'
  clustNameMethod = 'pagerank', # method for selecting cluster names. Available methods: 'pagerank', 'hits' and 'none'
  drawEllipses = TRUE, fontSize = 4, nodeSize = 'Pathway.Size',
  colorBy = 'p.adjust', colorType = 'pval',# pCutoff = -5,
  minClusterSize = 2
)

cls1 <- enw1$clusters
p1 <- enw1$plot +
  scale_color_viridis_c(name = 'Log10(p.adjust)', option = 'C', end = 0.8, limits = c(-10, -2)) +
  scale_size_continuous(name = "Pathway size", range = c(2, 6), limits = c(10, 490))

q1 <- ggplot_build(p1)
q1$data[[7]]$shape <- 21
q1$data[[7]]$fill <- q1$data[[7]]$colour

hash_table <- kegg_net_tbl %>% pull(Level, Description)
q1$data[[7]]$colour <- unname(level_color[hash_table[q1$data[[7]]$ID]])

q <- ggplot_gtable(q1)
p <- ggplotify::as.ggplot(q)

print(p)


# Test1 ------
suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(scatterpie)
  library(ggrepel)
})

# ----------------------------
# 1) Vertices: keep ALL attributes, enforce 1 row per ID
#    Decision: take the first row per ID
# ----------------------------
vertices_df <- kegg_net_tbl %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup() %>%
  rename(name = ID)

# Safety: igraph requires unique, non-NA names
if (any(is.na(vertices_df$name))) stop("Some vertex IDs are NA.")
if (any(duplicated(vertices_df$name))) stop("Duplicate vertex IDs remain after slice(1).")

# ----------------------------
# 2) Gene sets per ID (union across rows)
# ----------------------------
gene_tbl <- kegg_net_tbl %>%
  group_by(ID) %>%
  summarise(
    geneID = paste(geneID, collapse = ","),
    .groups = "drop"
  ) %>%
  filter(ID %in% vertices_df$name) %>%
  arrange(match(ID, vertices_df$name))

gene_sets <- gene_tbl$geneID |>
  gsub("\\.//\\.", ",", x = _) |>
  strsplit("[,/]+") |>
  lapply(function(x) unique(trimws(x[nzchar(trimws(x))])))

names(gene_sets) <- gene_tbl$ID
n_path <- length(gene_sets)

# ----------------------------
# 3) Similarity matrix (Jaccard)
# ----------------------------
sim_mat <- matrix(0, nrow = n_path, ncol = n_path, dimnames = list(names(gene_sets), names(gene_sets)))

for (i in seq_len(n_path)) {
  gi <- gene_sets[[i]]
  for (j in i:n_path) {
    gj <- gene_sets[[j]]
    inter <- length(intersect(gi, gj))
    uni <- length(unique(c(gi, gj)))
    s <- if (uni == 0) 0 else inter / uni
    sim_mat[i, j] <- s
    sim_mat[j, i] <- s
  }
}

# ----------------------------
# 4) Edge list from similarity matrix
# ----------------------------
sim_cutoff <- 0.2  # adjust density here

edge_idx <- which(sim_mat > sim_cutoff, arr.ind = TRUE)
edge_idx <- edge_idx[edge_idx[, 1] < edge_idx[, 2], , drop = FALSE]

edge_df <- data.frame(
  from   = rownames(sim_mat)[edge_idx[, 1]],
  to     = colnames(sim_mat)[edge_idx[, 2]],
  weight = sim_mat[edge_idx],
  stringsAsFactors = FALSE
)

# ----------------------------
# 5) Graph + layout
# ----------------------------
g <- graph_from_data_frame(edge_df, vertices = vertices_df, directed = FALSE)

set.seed(1)
lay <- ggraph::create_layout(g, layout = "fr")

nodes_plot <- as.data.frame(lay) %>%
  mutate(
    Protein = .data$n_Protein,
    RNA     = .data$n_RNA,
    Both    = .data$n_Both,
    Total   = Protein + RNA + Both,
    # radius in layout units; tuned for ~20–80 nodes. adjust 0.02/0.08 if needed
    r       = ifelse(Total > 0, (sqrt(Total) / max(sqrt(Total))) * 0.08, 0.02)
  )

# ----------------------------
# 6) Plot: edges + pie nodes
# ----------------------------
p <- ggraph(lay) +
  geom_edge_link(aes(width = weight, alpha = weight), colour = "grey65") +
  scale_edge_width(range = c(0.2, 2), guide = "none") +
  scale_edge_alpha(range = c(0.2, 0.8), guide = "none") +
  
  scatterpie::geom_scatterpie(
    data = nodes_plot,
    aes(x = x, y = y, r = r),
    cols = c("Protein", "RNA", "Both"),
    color = "white",
    size = 0.25
  ) +
  scale_fill_manual(values = c(Protein = "#1f77b4", RNA = "#2ca02c", Both = "#ff7f0e")) +
  
  ggrepel::geom_text_repel(
    data = nodes_plot,
    aes(x = x, y = y, label = Description),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  
  theme_void() +
  guides(fill = guide_legend(title = "Pie"))

p

# Test2 ------

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)

## 4.2.1 choose a seed tissue (class_abbr) and pathway set ----------
vec_class    <- sort(unique(kegg_enrich_agg$class_abbr))
target_class <- vec_class[4]   # e.g. "ADG"; change if needed

kegg_id <- kegg_enrich_agg %>%
  filter(class_abbr == target_class) %>%
  pull(ID) %>%
  unique()

## 4.2.2 aggregate across tissues to one row per pathway for aPEAR
## here Level pies are based on counts of tissues in each Level category
kegg_net_tbl <- kegg_enrich_agg %>%
  filter(class_abbr == target_class, ID %in% kegg_id) %>%
  mutate(
    n_Protein = ifelse(Level == 'Protein', 1, 0),
    n_RNA = ifelse(Level == 'RNA', 1, 0),
    n_Both = ifelse(Level == 'Both', 1, 0),
  ) %>% 
  arrange(p.adjust) %>% 
  as.data.frame()

## 4.2.3 run aPEAR to obtain similarity matrix + clusters ------------
set.seed(1)
enw1 <- enrichmentNetwork(
  kegg_net_tbl,
  plotOnly        = FALSE,
  clustMethod     = "hier",     # same as your previous settings
  simMethod       = "cosine",
  clustNameMethod = "pagerank",
  drawEllipses    = FALSE,
  fontSize        = 4,
  nodeSize        = "Pathway.Size",
  colorBy         = "p.adjust",
  colorType       = "pval",
  minClusterSize  = 2
)

sim_mat  <- enw1$sim
clusters <- enw1$clusters       # named by pathway ID

## 4.2.4 edge list from similarity matrix ----------------------------
edges <- sim_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("from") %>%
  tidyr::pivot_longer(-from, names_to = "to", values_to = "weight") %>%
  filter(from != to, weight > 0)   # you can tighten this (e.g. weight > 0.3)

## 4.2.5 node attributes + clusters ----------------------------------
nodes <- kegg_net_tbl %>%
  mutate(
    Cluster = clusters[ID],
    Cluster = factor(Cluster)
  )

## 4.2.6 igraph + layout ---------------------------------------------
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

set.seed(1)
g_layout <- ggraph::create_layout(g, layout = "fr") %>%
  as.data.frame()   # x, y, name, and vertex attributes

## merge node attributes into layout
nodes_plot <- g_layout %>%
  left_join(
    nodes,
    by = c("name" = "ID")
  ) %>%
  mutate(
    Protein = n_Protein,
    RNA     = n_RNA,
    Both    = n_Both
  )

## edge layout for geom_segment
edges_plot <- edges %>%
  left_join(
    g_layout %>% select(name, x, y),
    by = c("from" = "name")
  ) %>%
  rename(x = x, y = y) %>%
  left_join(
    g_layout %>% select(name, xend = x, yend = y),
    by = c("to" = "name")
  )

## 4.2.7 pie-network plot --------------------------------------------
p_pie <- ggplot() +
  ## edges
  geom_segment(
    data = edges_plot,
    aes(x = x, y = y, xend = xend, yend = yend, size = weight),
    alpha  = 0.3,
    colour = "grey60"
  ) +
  scale_size_continuous(range = c(0.1, 1.2), guide = "none") +
  
  ## pie nodes: radius ~ sqrt(Pathway.Size) (or change to sqrt(total_n))
  geom_scatterpie(
    data = nodes_plot,
    aes(
      x = x,
      y = y,
      r = sqrt(Pathway.Size) * 0.03
    ),
    cols  = c("Protein", "RNA", "Both"),
    color = NA
  ) +
  scale_fill_manual(
    name   = "Level",
    values = c(
      Protein = "#E64B35",  # Protein-only enriched pathways (in tissues)
      RNA     = "#4DBBD5",  # RNA-only enriched
      Both    = "#00A087"   # Both protein & RNA enriched
    )
  ) +
  coord_equal() +
  theme_void() +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold")
  ) +
  facet_wrap(~ Cluster) +
  ggtitle(paste0("KEGG network with Level pies (seed tissue: ", target_class, ")"))

print(p_pie)

ggsave(
  filename = paste0("TPHP_HPA_tissue_enriched_KEGG_pie_network_", target_class, ".pdf"),
  plot     = p_pie,
  width    = 9,
  height   = 9
)


# Test3 ------
## 4.2 network with pie nodes (using only aPEAR::enrichmentNetwork) ------

library(aPEAR)
library(igraph)
library(ggraph)
library(scatterpie)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)

id2clusters <- enw1$clusters %>% pull(Cluster, ID)   # named by pathway ID
base_plot <- enw1$plot
kegg_net_tbl$Cluster <- id2clusters[kegg_net_tbl$Description]
kegg_net_tbl$Cluster <- factor(kegg_net_tbl$Cluster)

## 4.2.4 build similarity matrix manually (Jaccard on gene sets) -----
##     -> edges for our own igraph/ggraph layout

# turn ".//." into "," first, then split on "," or "/"
kegg_net_tbl2 <- kegg_net_tbl %>%
  group_by(ID, Description, Cluster) %>%   # keep other grouping vars as appropriate
  summarise(
    geneID = paste(geneID, collapse = ","),  # merge gene lists
    .groups = "drop"
  ) %>%
  distinct(ID, .keep_all = TRUE)

# Rebuild gene_sets and similarity matrix using the deduplicated table
gene_sets <- kegg_net_tbl2$geneID %>%
  gsub("\\.//\\.", ",", .) %>%
  strsplit("[,/]+")

names(gene_sets) <- kegg_net_tbl2$ID
n_path <- length(gene_sets)

sim_mat <- matrix(
  0,
  nrow = n_path,
  ncol = n_path,
  dimnames = list(kegg_net_tbl2$ID, kegg_net_tbl2$ID)
)

for (i in seq_len(n_path)) {
  gi <- unique(gene_sets[[i]])
  for (j in i:n_path) {
    gj <- unique(gene_sets[[j]])
    inter <- length(intersect(gi, gj))
    union <- length(unique(c(gi, gj)))
    s <- if (union == 0) 0 else inter / union
    sim_mat[i, j] <- s
    sim_mat[j, i] <- s
  }
}

sim_cutoff <- 0.2

edge_idx <- which(sim_mat > sim_cutoff, arr.ind = TRUE)
edge_idx <- edge_idx[edge_idx[, "row"] < edge_idx[, "col"], , drop = FALSE]

edge_df <- data.frame(
  from   = rownames(sim_mat)[edge_idx[, "row"]],
  to     = colnames(sim_mat)[edge_idx[, "col"]],
  weight = sim_mat[edge_idx],
  row.names = NULL
)

vertices_df <- kegg_net_tbl2 %>%
  rename(name = ID)

g <- graph_from_data_frame(edge_df, vertices = vertices_df, directed = FALSE)
set.seed(1)
g_layout <- ggraph::create_layout(g, layout = "fr")

nodes_plot <- as.data.frame(g_layout) %>%
  mutate(
    Protein = n_Protein,
    RNA     = n_RNA,
    Both    = n_Both
  )

edges_plot <- edge_df %>%
  left_join(nodes_plot %>% select(name, x, y),  by = c("from" = "name")) %>%
  rename(x = x, y = y) %>%
  left_join(nodes_plot %>% select(name, xend = x, yend = y), by = c("to" = "name"))

## 4.2.7 pie-network plot --------------------------------------------

p_pie <- ggplot() +
  ## edges
  geom_segment(
    data  = edges_plot,
    aes(x = x, y = y, xend = xend, yend = yend, size = weight),
    colour = "grey70",
    alpha  = 0.4
  ) +
  scale_size_continuous(range = c(0.1, 1.5), guide = "none") +
  
  ## pie nodes: radius ~ sqrt(Pathway.Size)
  geom_scatterpie(
    data = nodes_plot,
    aes(
      x = x,
      y = y,
      r = sqrt(Pathway.Size) * 0.03
    ),
    cols  = c("Protein", "RNA", "Both"),
    color = NA
  ) +
  scale_fill_manual(
    name   = "Level (tissue count)",
    values = c(
      Protein = "#E64B35",  # Protein-only enriched tissues
      RNA     = "#4DBBD5",  # RNA-only enriched tissues
      Both    = "#00A087"   # enriched in both Protein & RNA
    )
  ) +
  coord_equal() +
  facet_wrap(~ Cluster) +
  theme_void() +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold", size = 9)
  ) +
  ggtitle(
    paste0(
      "KEGG network with Level pies (seed tissue: ",
      target_class, ")"
    )
  )

print(p_pie)

ggsave(
  filename = paste0("N_tissue_enriched_KEGG_pie_network_", target_class, ".pdf"),
  plot     = p_pie,
  width    = 9,
  height   = 9
)






