setwd(r"(\\172.16.13.136\tphp\code\2_batch_effect_evaluation_correction)")
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(scales)

my_pvcaBatchAssess <- function (theDataMatrix, expInfo, threshold) {
  # Cite: Bushel P (2024). pvca: Principal Variance Component Analysis (PVCA). R package version 1.44.0.
  # https://doi.org/10.1002/9780470685983.ch12
  # theDataMatrix, row as probs, col as samples;
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)
  theDataMatrixCentered_transposed = apply(theDataMatrix, 1, 
                                           scale, center = TRUE, scale = FALSE)
  theDataMatrixCentered = t(theDataMatrixCentered_transposed)
  theDataCor <- cor(theDataMatrixCentered)
  eigenData <- eigen(theDataCor)
  eigenValues = eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix = eigenData$vectors
  eigenValuesSum = sum(eigenValues)
  percents_PCs = eigenValues/eigenValuesSum
  
  exp_design <- as.data.frame(expInfo)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  my_counter_2 = 0
  my_sum_2 = 1
  for (i in ev_n:1) {
    my_sum_2 = my_sum_2 - percents_PCs[i]
    if ((my_sum_2) <= threshold) {
      my_counter_2 = my_counter_2 + 1
    }
  }
  pc_n <- ifelse(my_counter_2 < 3, 3, my_counter_2)
  pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN * 
                                               pc_n), ncol = 1)
  mycounter = 0
  for (i in 1:pc_n) {
    for (j in 1:expDesignRowN) {
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter, 1] = eigenVectorsMatrix[j, 
                                                        i]
    }
  }
  AAA <- exp_design[rep(1:expDesignRowN, pc_n), ]
  Data <- cbind(AAA, pc_data_matrix)
  variables <- c(colnames(exp_design))
  for (i in 1:length(variables)) {
    Data$variables[i] <- as.factor(Data$variables[i])
  }
  op <- options(warn = (-1))
  on.exit(options(op))
  effects_n = expDesignColN + choose(expDesignColN, 2) + 1
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  model.func <- c()
  index <- 1
  for (i in 1:length(variables)) {
    mod = paste("(1|", variables[i], ")", sep = "")
    model.func[index] = mod
    index = index + 1
  }
  for (i in 1:(length(variables) - 1)) {
    for (j in (i + 1):length(variables)) {
      mod = paste("(1|", variables[i], ":", variables[j], 
                  ")", sep = "")
      model.func[index] = mod
      index = index + 1
    }
  }
  function.mods <- paste(model.func, collapse = " + ")
  for (i in 1:pc_n) {
    y = (((i - 1) * expDesignRowN) + 1)
    funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lme4::lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
                                         expDesignRowN), ], REML = TRUE, verbose = FALSE, 
                        na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(lme4::VarCorr(Rm1ML)), 
                                  resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(lme4::getME(Rm1ML, "cnms")), "resid")
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, 
                                     ncol = effects_n)
  for (i in 1:pc_n) {
    mySum = sum(randomEffectsMatrix[i, ])
    for (j in 1:effects_n) {
      randomEffectsMatrixStdze[i, j] = randomEffectsMatrix[i, j]/mySum
    }
  }
  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, 
                                      ncol = effects_n)
  for (i in 1:pc_n) {
    weight = eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n) {
      randomEffectsMatrixWtProp[i, j] = randomEffectsMatrixStdze[i, j] * weight
    }
  }
  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
  totalSum = sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, 
                                         ncol = effects_n)
  for (j in 1:effects_n) {
    randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
  }
  return(list(dat = randomEffectsMatrixWtAveProp, label = effectsNames, matrixWtProp = set_colnames(randomEffectsMatrixWtProp, effectsNames), eigenData = eigenData))
}

# input
info <- read.csv("input/batch_info_tab_batchserver.csv",stringsAsFactors = FALSE)
str(info)
metadata <- info %>%
  column_to_rownames('FileName') %>% 
  select(instrument, trans, batch, sample_type, new, tissue_type_major)

path_list <- list.files('output', pattern = '3', full.names = T)
names(path_list) <- str_extract(path_list, '(combat|limma|step2_limma)\\w+')
path <- c(imp = 'input/imputed_raw_12754proteins_2957samples.csv',
          qimp = 'input/quantiled_log2_imputated_matrix_2957_12797_imputed.csv')
path %<>% append(path_list)
pm_list <- lapply(path, function(X){
  rio::import(X) %>% column_to_rownames('V1')
})
# pm_list$combat<-read.csv("output/batch_corrected_data_combat3.csv",row.names = 1, check.names = FALSE)

sapply(pm_list, function(X) quantile(as.matrix(X)))
#            imp     qimp combat_mean_only3   combat3    limma3 step2_limma3
# 0%    2.674865  6.38579          5.921994 -1.310841  3.949413     3.161037
# 25%   2.674865  6.38579          6.407567  6.388891  6.386089     6.401242
# 50%   2.674865  6.38579          6.490974  6.971537  6.892888     7.064453
# 75%  11.688395 12.35507         12.365733 12.357216 12.365974    12.446083
# 100% 28.796569 22.61057         22.872299 33.539784 24.994231    24.906838
#         combat
# 0%   -1.310841
# 25%   6.388891
# 50%   6.971537
# 75%  12.357216
# 100% 33.539784

# # for each protein matrix in `pm_list`
# pvca_list<-function(pm_list){
#   res_pvca <- list()
#   for(i in seq_along(names(pm_list))){
#     cat(i, '...\n')
#     theDataMatrix <- pm_list[[i]] %>% column_to_rownames('V1') %>% t()
#     expInfo <- metadata
#     threshold <- 0.8
#     pvca <- my_pvcaBatchAssess(theDataMatrix, expInfo, threshold)
#     res_pvca[[i]] <- pvca
#     }
#    return(res_pvca)
# }

res_pvca <- list()
for(i in seq_along(names(pm_list))){
  message(i, '-', names(pm_list)[i], '...\n')
  theDataMatrix <- pm_list[[i]]
  expInfo <- metadata
  threshold <- 0.8
  pvca <- my_pvcaBatchAssess(theDataMatrix, expInfo, threshold)
  res_pvca[[i]] <- pvca
}
names(res_pvca) <- names(pm_list)
save(res_pvca, file = 'res_pvca_beca_V1201.RData')



# visualization
mycolors <- unique(c(ggsci::pal_npg()(10), ggsci::pal_d3()(10), ggsci::pal_cosmic()(10)))
if(!dir.exists('pvca_plot')){
  dir.create('pvca_plot')
}

for(i in seq_along(res_pvca)){
  cat(i, '...\n')
  nm <- names(res_pvca)[i]
  pvca <- res_pvca[[i]]
  pca_colors <- mycolors[1:length(pvca$label)] %>% setNames(pvca$label)
  
  # Average proportion
  df_pvca <- data.frame(t(pvca$dat)) %>% 
    cbind(pvca$label) %>% 
    setNames(c('RandomEffectWtAveProp', 'EffectName')) %>% 
    arrange(RandomEffectWtAveProp) %>% 
    mutate(EffectName = factor(EffectName, levels = EffectName)) %>% 
    arrange(desc(RandomEffectWtAveProp)) %>% 
    mutate(y_label = cumsum(RandomEffectWtAveProp) - 0.5 * RandomEffectWtAveProp)
  plot_pvca <- ggplot(df_pvca, aes(x = 2, y = RandomEffectWtAveProp, fill = EffectName)) +
    geom_bar(stat = 'identity', color = 'white') +
    coord_polar(theta = 'y', start = 0) +
    theme(legend.position = 'none') +
    geom_text(aes(y = y_label, label = round(RandomEffectWtAveProp*100, 2)), color = 'white', size = 6) +
    labs(title = str_glue('{nm}_pvca_average')) +
    scale_fill_manual(values = pca_colors) +
    theme_void() +
    xlim(0.5, 2.5)
  ggsave(str_glue('pvca_plot/V1009/{nm}_pvca_average.pdf'), plot_pvca, width = 10, height = 10)
  
  # every PC
  df_prop <- data.frame(pvca$matrixWtProp)
  df_prop_percent <- data.frame(t(apply(df_prop, 1, function(x) x / sum(x))))
  df_prop_percent$PC <- str_glue('PC{1:nrow(df_prop_percent)}: ({round(100 * pvca$eigenData$values[1:nrow(df_prop_percent)] / sum(pvca$eigenData$values[1:nrow(df_prop_percent)]), 2)} %)')
  df_prop_percent$PC %<>% factor(., levels = .)
  
  plot_pvca_individual <- plyr::dlply(df_prop_percent, 'PC', function(dfsub){
    tbl <- dfsub %>% select(-PC) %>% 
      setNames(., str_replace(names(.), '\\.', ':')) %>% 
      t() %>% data.frame() %>% 
      setNames('RandomEffectWtAveProp') %>% 
      rownames_to_column('EffectName') %>% 
      arrange(RandomEffectWtAveProp) %>% 
      mutate(EffectName = factor(EffectName, levels = EffectName)) %>% 
      arrange(desc(RandomEffectWtAveProp)) %>% 
      mutate(y_label = cumsum(RandomEffectWtAveProp) - 0.5 * RandomEffectWtAveProp)
    
    set.seed(1000)
    ggplot(tbl, aes(x = 2, y = RandomEffectWtAveProp, fill = EffectName)) +
      geom_bar(stat = 'identity', color = 'white') +
      coord_polar(theta = 'y', start = 0) +
      geom_text(aes(y = y_label, label = round(RandomEffectWtAveProp*100, 2)), color = 'white', size = 6) +
      scale_fill_manual(values = pca_colors) +
      labs(subtitle = dfsub$PC) +
      theme_void() +
      theme(legend.position = 'none') +
      xlim(0.5, 2.5)
  })
  print(length(plot_pvca_individual)) # 402 for i==5
  legend_pies <- ggpubr::get_legend(plot_pvca)
  p <- ggpubr::ggarrange(plotlist = plot_pvca_individual, nrow = 11, ncol = 40) %>%
    ggpubr::ggarrange(legend_pies, widths = c(0.95, 0.05)) %>% 
    ggpubr::annotate_figure(top = ggpubr::text_grob(str_glue('{nm}_pvca_individuals'),
                                                    color = "black",
                                                    face = "bold",
                                                    size = 14))
  # plot(cumsum(pvca$eigenData$values / sum(pvca$eigenData$values)), xlab = 'PC_n', ylab = 'Summed eigenvalues proportion')
  ggsave(str_glue('pvca_plot/V1009/{nm}_pvca_individuals.pdf'), p, width = 4*40, height = 4*11, limitsize = F)
  
  
  # top3 PCs
  pvca$eigenData$values[1:3]
  dfsub <- apply(df_prop_percent[1:3, 1:which(colnames(df_prop_percent) == 'resid')], 2, function(y) {
    y * pvca$eigenData$values[1:3]
  }) %>% as.data.frame() %>%
    mutate_all(sum) %>% slice(1)
  dfsub <- dfsub / sum(dfsub)
  PVCA_top3 <- dfsub %>%
    setNames(., str_replace(names(.), '\\.', ':')) %>% 
    t() %>% data.frame() %>% 
    setNames('RandomEffectWtAveProp') %>% 
    rownames_to_column('EffectName') %>% 
    arrange(RandomEffectWtAveProp) %>% 
    mutate(EffectName = factor(EffectName, levels = EffectName)) %>% 
    arrange(desc(RandomEffectWtAveProp)) %>% 
    mutate(y_label = cumsum(RandomEffectWtAveProp) - 0.5 * RandomEffectWtAveProp)
  
  plot_PVCA_top3 <- ggplot(PVCA_top3, aes(x = 2, y = RandomEffectWtAveProp, fill = EffectName)) +
    geom_bar(stat = 'identity', color = 'white') +
    coord_polar(theta = 'y', start = 0) +
    geom_text(aes(y = y_label, label = round(RandomEffectWtAveProp*100, 2)), color = 'white', size = 6) +
    # geom_text(aes(x = 1, y = y_label, label = round(RandomEffectWtAveProp*100, 2)), color = 'black', size = 6) +
    scale_fill_manual(values = pca_colors) +
    labs(title = str_glue('{nm}_pvca_top3PCs')) +
    theme_void() +
    xlim(0.5, 2.5)
  
  ggsave(str_glue('pvca_plot/V1009/{nm}_pvca_top3PCs.pdf'),plot_PVCA_top3, width = 10, height = 10)
  
  list(average = df_pvca,
       individuals = df_prop_percent,
       top3PCs = PVCA_top3) %>% 
    rio::export(str_glue('pvca_plot/V1009/{nm}_pvca_data.xlsx'))
}

