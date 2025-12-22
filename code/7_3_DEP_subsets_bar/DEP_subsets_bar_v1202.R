# 0. Packages ----
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Reproducibility & parallel
set.seed(1L)
options(stringsAsFactors = FALSE)
n_cores <- min(64L, parallel::detectCores())
BiocParallel::register(BiocParallel::SnowParam(workers = max(1L, n_cores - 2L)))

# Working dir & helpers
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('../source/source_code.R')
library(data.table)


# 1.Read data --------
# DEP
df_dep <- rio::import('../6_3paired_tumor_nontumor/output/tumor_nontumor_analysis_reports_v1202.xlsx', sheet = 3)

# DEP subsets
dia_target <- read.csv('../6_3paired_tumor_nontumor/TPHP_CA_targetlist_3.csv', check.names = F)
tumor_specific <- dia_target %>%
  pivot_longer(cols = everything(), names_to = 'cancer_abbr', values_to = 'protein') %>% 
  filter(protein != '') %>% 
  left_join(df_dep)

tumor_ledep <- rio::import('../7_1_LEDEP/output/LEDEP_218protein.xlsx')

tumor_enriched <- rio::import('../7_2_T_vs_N/output/tumor_enriched_q1Filter.xlsx', sheet = 'Tumor.enriched')


# 2.Statastics ----
stat_spec <- tumor_specific %>% count(cancer_abbr, direction) %>% mutate(DEP.Type = 'Tumor.specific')
stat_ledep <- tumor_ledep %>% count(cancer_abbr, direction) %>% mutate(DEP.Type = 'LEDEP')
stat_enriched <- tumor_enriched %>% count(cancer_abbr) %>% mutate(direction = 'Up', DEP.Type = 'Tumor.enriched')
df_stat <- rbind(stat_spec, stat_ledep, stat_enriched) %>% 
  mutate(DEP.Type = factor(DEP.Type, c('Tumor.specific', 'LEDEP', 'Tumor.enriched')))

# ggplot(df_stat) +
#   facet_wrap(~ DEP.Type, scales = 'free') +
#   aes(x = n, y = cancer_abbr, fill = direction) +
#   geom_bar(color = NA, stat = 'identity', position = position_stack()) +
#   geom_text(aes(label = n), position = position_stack(vjust = 0.5, reverse = FALSE)) +
#   labs(x = '', y = 'Cancer') +
#   scale_fill_manual(values = c(Up = "#CF3F54", Down = "#03A2B3")) +
#   scale_x_log10(breaks = c(1, 10, 100, 1000),
#                 expand = expansion(mult = c(0, 0))) +
#   theme_classic(base_size = 10) +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size = 10, face = "bold"))

# v1
# plot.bars1 <- ggplot(df_stat) +
#   facet_wrap(~ DEP.Type, scales = "free") +
#   aes(x = n, y = cancer_abbr, fill = direction) +
#   geom_bar(color = NA, alpha = 0.95, stat = "identity", position = position_stack()) +
#   geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
#   labs(x = "", y = "Cancer") +
#   scale_fill_manual(values = sample_color[1:2] %>% setNames(c('Up', 'Down'))) +
#   scale_x_continuous(
#     trans  = scales::pseudo_log_trans(base = 10, sigma = 5),  # allows 0 baseline. sigma=5 to avoid 4+5>=10
#     breaks = c(1, 10, 30, 100, 300, 1000),
#     expand = expansion(mult = c(0, 0))
#   ) +
#   theme_classic(base_size = 10) +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size = 10, face = "bold"))
# ggsave('DEP_subsets_bar_v1.pdf', plot.bars1, width = 8, height = 4)

# # v1 narrow
# plot.bars1 <- ggplot(df_stat) +
#   facet_wrap(~ DEP.Type, scales = "free", space = "free_x") +
#   aes(x = n, y = cancer_abbr, fill = direction) +
#   geom_bar(color = NA, alpha = 0.95, stat = "identity", position = position_stack()) +
#   geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
#   labs(x = "", y = "Cancer") +
#   scale_fill_manual(values = sample_color[1:2] %>% setNames(c('Up', 'Down'))) +
#   scale_x_continuous(
#     trans  = scales::pseudo_log_trans(sigma = 5),  # allows 0 baseline. sigma=5 to avoid 4+5>=10
#     breaks = c(1, 10, 30, 100, 300, 1000),
#     expand = expansion(mult = c(0, 0))
#   ) +
#   theme_classic(base_size = 8) +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size = 10, face = "bold"))
# ggsave('DEP_subsets_bar_v1_narrow.pdf', plot.bars1, width = 5, height = 4)

# # v2 free_x only
# plot.bars2 <- ggplot(df_stat) +
#   facet_wrap(~ DEP.Type, scales = "free_x") +
#   aes(x = n, y = cancer_abbr, fill = direction) +
#   geom_bar(color = NA, alpha = 0.95, stat = "identity", position = position_stack()) +
#   geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
#   labs(x = "", y = "Cancer") +
#   scale_fill_manual(values = sample_color[1:2] %>% setNames(c('Up', 'Down'))) +
#   scale_x_continuous(
#     trans  = scales::pseudo_log_trans(base = 10, sigma = 5),  # allows 0 baseline. sigma=5 to avoid 4+5>=10
#     breaks = c(1, 10, 100, 1000),
#     expand = expansion(mult = c(0, 0))
#   ) +
#   theme_classic(base_size = 10) +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size = 10, face = "bold"))
# ggsave('DEP_subsets_bar_v2.pdf', plot.bars2, width = 8, height = 4)

# # v3
# df_stat_up <- df_stat %>% filter(direction == 'Up')
# df_stat_down <- df_stat %>% filter(direction == 'Down')
# 
# df_stat_all <- df_stat %>% 
#   group_by(cancer_abbr, DEP.Type) %>% 
#   summarise(direction = 'Down', # psudo-down, actually "All"
#             n = sum(n), .groups = 'drop')
# 
# df_stat_psudo <- df_stat %>% filter(direction == 'Up') %>% 
#   rbind(df_stat_all) %>% 
#   group_by(cancer_abbr, DEP.Type) %>% 
#   reframe(n_pos = n[direction == 'Down'] * n[direction == 'Up']) %>% 
#   inner_join(df_stat_down, .)
# 
# ggplot(df_stat) +
#   facet_wrap(~ DEP.Type, scales = "free") +
#   aes(x = n, y = cancer_abbr, fill = direction) +
#   geom_bar(data = df_stat_all, color = NA, alpha = 0.95, stat = "identity", position = position_stack()) +
#   geom_bar(data = df_stat_up, color = NA, alpha = 0.95, stat = "identity", position = position_stack()) +
#   geom_text(data = df_stat_up, aes(label = n), position = position_stack(vjust = 0.5)) +
#   geom_text(data = df_stat_psudo, aes(x = n_pos, label = n), position = position_stack(vjust = 0.5)) +
#   labs(x = "", y = "Cancer") +
#   scale_fill_manual(values = sample_color[1:2] %>% setNames(c('Up', 'Down'))) +
#   scale_x_continuous(
#     trans  = scales::pseudo_log_trans(sigma = 1),  # allows 0 baseline. sigma=5 to avoid 4+5>=10
#     breaks = c(1, 10, 30, 100, 300),
#     expand = expansion(mult = c(0, 0))
#   ) +
#   theme_classic(base_size = 8) +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(size = 10, face = "bold"))

# v4 complicated code
# Keep the palette identical to the bars
# Keep palette identical to bars
pal <- sample_color[1:2] %>% setNames(c("Up","Down"))

# One row per bar with Up/Down and total
df_lab <- df_stat %>%
  group_by(DEP.Type, cancer_abbr, direction) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(
    Up    = coalesce(Up,   0L),
    Down  = coalesce(Down, 0L),
    total = Up + Down
  )

# Build one subplot with slot-based placement (no monospace)
# Slots are in DATA UNITS and constant within each subplot:
# [start_gap] [Slot-1 for Up (width W1)] [gap G] [pipe] [gap G] [Slot-2 for Down (width W2)]
make_dep_plot_slots_du <- function(dep_type,
                                   start_gap_k = 0.00,  # 0 = flush to bar end; small (>0) = tiny gap
                                   char_k      = 0.040, # data-units per “character” of slot width
                                   min_col_du  = 0.8,   # minimum slot width (data units)
                                   gap_k       = 0.020, # gap G (data units as fraction of span) on each side of "|"
                                   pipe_size   = 2.6,
                                   label_size  = 2.6) {
  
  ds <- dplyr::filter(df_stat, DEP.Type == dep_type)
  dl <- dplyr::filter(df_lab,  DEP.Type == dep_type)
  
  # per-subplot range
  x_max <- max(dl$total, na.rm = TRUE)
  span  <- if (is.finite(x_max) && x_max > 0) x_max else 1
  
  # longest digit counts in this subplot (determine slot widths W1/W2)
  up_chars_max   <- dl %>% dplyr::filter(Up   > 0) %>% pull(Up)   %>% as.character() %>% nchar()
  down_chars_max <- dl %>% dplyr::filter(Down > 0) %>% pull(Down) %>% as.character() %>% nchar()
  up_chars_max   <- if (length(up_chars_max))   max(up_chars_max)   else 0L
  down_chars_max <- if (length(down_chars_max)) max(down_chars_max) else 0L
  
  # constants within this subplot (in DATA UNITS)
  W1 <- max(up_chars_max   * char_k * span, min_col_du)  # Slot-1 (Up) width
  W2 <- max(down_chars_max * char_k * span, min_col_du)  # Slot-2 (Down) width
  G  <- gap_k * span                                     # fixed gap on each side of "|"
  S0 <- start_gap_k * span                               # fixed gap from bar to label block
  
  # Left edge of the whole label block (must align to bar end when S0=0)
  # x0 is the bar end + optional tiny start gap
  dl_pos <- dl %>%
    dplyr::mutate(
      x0     = total + S0,        # LEFT EDGE of the label block (align this to bar end)
      x_up   = x0,                # Up starts at left edge (left-aligned)
      x_pipe = x0 + W1 + G,       # pipe center between slots
      x_dn   = x0 + W1 + 2*G      # Down starts at Slot-2 (left-aligned)
    )
  
  # right margin so labels never clip
  right_margin <- S0 + W1 + 2*G + W2 + 0.5
  
  ggplot(ds, aes(x = n, y = cancer_abbr, fill = direction)) +
    geom_bar(color = NA, alpha = 0.95, stat = "identity", position = position_stack()) +
    
    # --- LABELS (all anchored to bar end via x0) ---
    # Up only: left edge = bar end
    geom_text(
      data = dplyr::filter(dl_pos, Up > 0, Down == 0),
      aes(x = x_up, y = cancer_abbr, label = Up),
      inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = label_size, color = pal[["Up"]]
    ) +
    # Down only: left edge = bar end
    geom_text(
      data = dplyr::filter(dl_pos, Down > 0, Up == 0),
      aes(x = x0, y = cancer_abbr, label = Down),
      inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = label_size, color = pal[["Down"]]
    ) +
    # Both present: Up at x0 (left-aligned) …
    geom_text(
      data = dplyr::filter(dl_pos, Up > 0, Down > 0),
      aes(x = x_up, y = cancer_abbr, label = Up),
      inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = label_size, color = pal[["Up"]]
    ) +
    # … pipe centered at fixed offset …
    geom_text(
      data = dplyr::filter(dl_pos, Up > 0, Down > 0),
      aes(x = x_pipe, y = cancer_abbr), label = "|",
      inherit.aes = FALSE, hjust = 0.5, vjust = 0.5, size = pipe_size, color = "grey40"
    ) +
    # … Down at Slot-2 left edge
    geom_text(
      data = dplyr::filter(dl_pos, Up > 0, Down > 0),
      aes(x = x_dn, y = cancer_abbr, label = Down),
      inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = label_size, color = pal[["Down"]]
    ) +
    
    labs(x = "", y = "", title = dep_type) +
    scale_fill_manual(values = pal) +
    scale_x_continuous(limits = c(0, x_max + right_margin), expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}


# Build three subplots (NO facets), with the same rules you asked (fixed spacing per subplot)
plot_spec     <- make_dep_plot_slots_du("Tumor.specific",  start_gap_k = 0.00, char_k = 0.040, min_col_du = 0.8, gap_k = 0.020)
plot_ledep    <- make_dep_plot_slots_du("LEDEP",           start_gap_k = 0.00, char_k = 0.040, min_col_du = 0.8, gap_k = 0.020)
plot_enriched <- make_dep_plot_slots_du("Tumor.enriched",  start_gap_k = 0.00, char_k = 0.040, min_col_du = 0.8, gap_k = 0.020)

plot.bars <- ggpubr::ggarrange(plot_spec, plot_ledep, plot_enriched, nrow = 1, ncol = 3,
                               widths = c(2, 1, 1), common.legend = T, legend = 'right')
ggsave('output/DEP_subsets_bar_narrow.pdf', plot.bars, width = 4, height = 4)


df_stat %>% rename(Cancer = cancer_abbr) %>%
  select(DEP.Type, Cancer, direction, n) %>% 
  arrange(DEP.Type, Cancer) %>% 
  rio::export('output/DEP_subsets_bar.xlsx')




# stat
length(intersect(tumor_specific$protein, tumor_ledep$Protein)) # 131
tumor_ledep %>% count(cancer_abbr) %>% nrow() # 13
tumor_ledep %>% count(direction)
# direction   n
# 1      Down 177
# 2        Up  41


