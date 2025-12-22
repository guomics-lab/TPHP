pacman::p_ver("canvasXpress")
# [1] ‘1.29.6’
# the version is vital. ‘1.46.9.1’ and ‘1.57.4’ are not worked.


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../source/source_code.R")
library(dplyr)
library(magrittr)
library(writexl)
library(canvasXpress)


class2abbr <- df_abbr %>% pull(Abbr, Entire)
df1 <- rio::import('../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V5.xlsx')

# canvas -----
stat.lib <- rio::import('DDA_lib_stat.xlsx')
df3 <- df1 %>%
  filter(sample_type != "p",
         !is.na(class_abbr)) %>%
  mutate(
    sample_type = dplyr::case_when(
      sample_type == "F"  ~ -50000,
      sample_type == "N" ~ -10000,
      sample_type == "NT"  ~ 30000,
      sample_type == "T"  ~ 70000,
      TRUE                ~ NA_real_
    )
  ) %>%
  mutate(FileName = str_replace_all(FileName, '\\-', '.')) %>% 
  arrange(class_abbr, sample_type) %>%
  rename(
    Peptides = `# peptides`,
    Proteins = `# proteins`
  ) %>% 
  left_join(stat.lib)
df3 %<>% mutate(DDA_lib_type = ifelse(is.na(DDA_lib_type), class_abbr, DDA_lib_type),
                DDA_lib_type = ifelse(DDA_lib_type == 'NA', 'NAIL', DDA_lib_type),
                DDA_lib_pro = ifelse(is.na(DDA_lib_pro), 0, DDA_lib_pro))

# build Y (make sure it's numeric!)
tmpY <- df3[, c("Peptides", "Proteins", "sample_type", "DDA_lib_pro")]
tmpY[] <- lapply(tmpY, as.numeric)   # force numeric
datay <- as.data.frame(t(tmpY))
colnames(datay) <- df3$FileName

datax <- data.frame(tissue_name = df3$DDA_lib_type,
                    row.names   = df3$FileName)

dataz <- data.frame(ring = c(1, 2, 3, 4),
                    row.names = c("Peptides", "Proteins", "sample_type", "DDA_lib_pro"))

canvasXpress::canvasXpress(
  data = datay,
  smpAnnot = datax,
  varAnnot = dataz,
  circularArc = 340,
  circularType = "radar",
  graphType = "Circular",
  legendPosition = "top",
  ringGraphType = list("bar", "bar", "heatmap", "bar"),
  smpDendrogramPosition = "ins",
  arcSegmentsSeparation = 3,
  showTransition = FALSE,
  showYaxis = FALSE,
  ringGraphWeight = list(22, 22, 3, 53),
  ringsOrder = list("labels", "overlays", "dendrogram", "data"),
  segregateVariablesBy = list("ring"),
  segregateSamplesBy = list("tissue_name"),
  title = "Overview of TPHP proteome data",
  smpOverlays = list("tissue_name", "-", "sample_type"),
  colorScheme = "Tableau",
  circularAnchors2Align = "inside",
  showSampleNames = FALSE,
  
  # foreground = "rgba(0.1,0.1,0.1,0.1)",
  # showHeatmapIndicator = FALSE,
  # showLegend = FALSE,
)


## output data -----
write_xlsx(df3, "20251031TPHP_canvas_plot.xlsx")
write.csv(datay, "datay_v2_20251031.csv")
write.csv(datax, "datax_v2_20251031.csv")
write.csv(dataz, "dataz_v2_20251031.csv")
