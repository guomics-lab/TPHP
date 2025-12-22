# set order
# pm_data<-readRDS("mat_med_na_list.rds")
setwd(r"(\\172.16.13.136\tphp\code\4_tissue_specificity_analysis)") 
library(tidyr)
library(dplyr)
library(stringr)
library(magrittr)
library(tibble)
library(stringr)
pm_data<-readRDS("//172.16.13.136/tphp/code/4_tissue_specificity_analysis/output/20251201_tissue_Specificity_naive_05minImputated_3.rds")
source('//172.16.13.136/share/members/jiangwenhao/TPHP/my_fun.R')
######
label_info<-pm_data$metadata %>% filter(sample_type == "N")%>%
  select(anatomical_classification, class_abbr)%>%
  distinct()
label_info$tissue_type<-str_match(label_info$anatomical_classification, "-(.*?)$")[,2]
label_info[is.na(label_info$tissue_type),]$tissue_type<-label_info[is.na(label_info$tissue_type),]$anatomical_classification
rownames(label_info)<-label_info$class_abbr
pm1 <- pm_data$med.matrix
pm1<-data.frame(pm1, row.names = pm1$class_abbr)[,-which(names(pm1) == "class_abbr")]
pm2 <- pm1[,colnames(pm1) %in% unique(pm_data$ts.Enrich$Protein)]
# dim(pm2)
# pm2<-pm1
rownames(pm2)<-label_info[rownames(pm1),]$tissue_type
# pm2[pm2==min]<-NA
pm2 %<>% log2()
head(pm2[,1:5])
# pm2 <- pm_data$imp
# pm2[is.na(pm2)]<-(min(pm2,na.rm=T))*0.5
set.seed(10) 
pm_heat_cor <- cor(t(pm2), method = 'pearson')
pm_heat_dist <- as.dist(1 - pm_heat_cor) # 1-Spearman's rho
pm_heat_tree <- hclust(pm_heat_dist, method="complete")
organs_ordered <- pm_heat_tree$labels[pm_heat_tree$order]
organs_ordered %<>% rev() # for mapping to barplot

### ===cluster -- tissue enriched proteins===

# set colors -->
require(RColorBrewer) 
require(dendextend)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
set.seed(10) 
# label_colors <- sample(col_vector, 50)
label_colors <- colorRampPalette(brewer.pal(12, "Set3"))(10)

pm_heat_dend <- as.dendrogram(pm_heat_tree) # create dendrogram object 
# plot(pm_heat_dend, horiz = T )#, leaflab = "none"
library(viridis)
label_colors <- viridis(10)

pdf("output/20251204_TPHP_DIA_heat_tree_tissue_enriched_only.pdf",width = 8,height = 6) 
par(mar=c(1,1,1,7))
pm_heat_dend %>%
  set("labels_cex", 0.5) %>%
  set("labels_col", value = label_colors, k=10) %>% 
  set("branches_k_color", value = label_colors, k = 10) %>% #c("skyblue", "orange", "grey")
  set("branches_lwd", 1.5) %>%  # Significantly increase line width
  plot(horiz=T, axes=FALSE)
dev.off()
# abline(v = 350, lty = 2)
# p3 <- recordPlot()
# p4<-rev(p3)  
############################################################


####### speciticity protein class protein number
organs_ordered <- pm_heat_tree$labels[pm_heat_tree$order]

# ---------  DIA identification generation ---------
### directly from classification result
X<- pm_data$ts.df
X$tissue_type<-label_info[X$tissue_type,]$tissue_type
# X <- rio::import('20231113TPHP_51tissue_ulhen_classification.xlsx')
Y <- X %>%
  pivot_longer(-tissue_type, names_to = 'prot', values_to = 'class') %>%
  plyr::ddply(c('tissue_type', 'class'), function(dfsub){
    df_ret <- dfsub %>% slice(1) %>% select(tissue_type, class)
    df_ret$n <- nrow(dfsub)
    return(df_ret)
  }) %>%
  setNames(c('Sample type', 'Protein specificity', 'protein_number'))

# only in c('tissue enriched', 'group enriched', 'tissue enhanced')
# Y %<>% filter(`Protein specificity` %in% c('tissue enriched', 'group enriched', 'tissue enhanced'))


# get abbreviation of sample types, if possible
# df_abbr <- read.delim('//172.16.13.136/share/members/jiangwenhao/TPHP/input/sample_types_abbr_20230113.txt', stringsAsFactors = F, check.names = F, na.strings = '')
# df_DIA <- get_abbr(Y, 'Sample type', df_abbr = df_abbr)
df_DIA<-Y
df_DIA$`Protein specificity` %<>% str_to_sentence()

# remove 'not detected'
df_DIA %<>% filter(`Protein specificity` != 'Not detected')
df_DIA$`Protein specificity` %<>%
  factor(levels = rev(c('Tissue enriched', 'Group enriched', 'Tissue enhanced', 'Mixed', 'Expressed in all tissues')), ordered = T)#"Not detected",

# df_DIA$`Protein specificity` %<>%
#   factor(levels = rev(c('Tissue enriched', 'Group enriched', 'Tissue enhanced')), ordered = T)

df_DIA$`Sample type` %<>% factor(levels = rev(organs_ordered), ordered = T)
df_DIA %<>% arrange(`Sample type`, `Protein specificity`)


# --------- Visualization ---------------
pacman::p_load('ggpubr', 'RColorBrewer', 'ggsci')
# function
draw_F2A_barplot <- function(df, char_title, vec_color){
  df$protein_number <- df$protein_number / 1000
  p <- ggplot(df, aes(x=`Sample type`, y=`protein_number`, fill = `Protein specificity`)) + geom_bar(stat="identity", width=0.9, position="stack")
  
  p <- p +
    scale_fill_manual(values = vec_color)+
    labs(
      y = "Protein number (×1000)",
      title = char_title
    )+
    guides(fill=guide_legend(nrow = 2))+
    theme(panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12,color = "black"),
          axis.line = element_line(colour = "black"),
          axis.line.x = element_blank(),
          plot.subtitle = element_text(size = 30, hjust = 0, color="black")
    )+
    theme(legend.text = element_text(size = 15, color = "black"),legend.position = 'bottom',
          legend.title = element_text(size = 18, color="black"),)+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, vjust = 0,color = "black", angle = 90),
          axis.ticks.x = element_blank()
    )+
    theme(axis.title.y = element_text(size = 15, hjust = 0.5, color="black", angle = 90),
          axis.text.y = element_text(size = 12,color = "black", angle = 0),
          axis.ticks.y = element_blank()
    )
  return(p)
}

draw_F2A_barplot_percen <- function(df, char_title, vec_color){
  # Calculate the percentage of each protein-specificity category within each sample type
  df_percent <- df %>%
    group_by(`Sample type`) %>%
    mutate(percentage = protein_number / sum(protein_number) * 100) %>%
    ungroup()
  
  p <- ggplot(df_percent, aes(x = `Sample type`, y = percentage, fill = `Protein specificity`)) + 
    geom_bar(stat = "identity", width = 0.9, position = "stack")
  
  p <- p +
    scale_fill_manual(values = vec_color) +
    scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) + # Set y-axis range to 0–100%
    labs(
      y = "Abundance (%)",
      title = char_title
    ) +
    guides(fill = guide_legend(nrow = 2)) +
    theme(panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12, color = "black"),
          axis.line = element_line(colour = "black"),
          axis.line.x = element_blank(),
          plot.subtitle = element_text(size = 30, hjust = 0, color = "black")
    ) +
    theme(legend.text = element_text(size = 15, color = "black"),
          legend.position = 'bottom',
          legend.title = element_text(size = 18, color = "black")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, vjust = 0, color = "black", angle = 90),
          axis.ticks.x = element_blank()
    ) +
    theme(axis.title.y = element_text(size = 15, hjust = 0.5, color = "black", angle = 90),
          axis.text.y = element_text(size = 12, color = "black", angle = 0),
          axis.ticks.y = element_blank()
    )
  return(p)
}


# parameters
ls_color <- list(c(214, 133, 122),
                 c(248, 176, 91),
                 c(240, 228, 66),
                 c(96, 217, 125),
                 c(94, 151, 181),
                 c(129, 117, 153))
vec_color <- hex_col(ls = ls_color)
my_colors <- c("#FC8D62", "#E5C494", "#984EA3")
# rev(c('Tissue enriched', 'Group enriched', 'Tissue enhanced', 'Mixed', 'Expressed in all tissues'))
# c("#A6D854", "#CCCCCC", "#FC8D62", "#E5C494", "#984EA3")
# p1 <- draw_F2A_barplot(df_DDA, 'DDA qualification', vec_color)
my_colors<-c('Mixed'="#CCCCCC","Expressed in all tissues"="#A4BED5FF","Not detected"="#A6D854",'Tissue enriched'="#984EA3",'Group enriched'="#E5C494",'Tissue enhanced'="#FC8D62")

# my_colors<-c('Tissue enriched'="#984EA3",'Group enriched'="#E5C494",'Tissue enhanced'="#FC8D62")
df_DIA1<-df_DIA%>%filter(`Protein specificity` %in% c('Tissue enriched', 'Group enriched', 'Tissue enhanced'))
p1 <- draw_F2A_barplot(df_DIA1, 'DIA qualification', my_colors)
p2 <- draw_F2A_barplot_percen(df_DIA, 'DIA Abundance', my_colors)

p <- ggarrange(p1, p2,
               ncol = 1, nrow = 2,
               labels = c("A","B"),
               font.label = list(size = 14, face = "bold"))

# output
stringr::str_c('TPHP_DIA_qualification_', as.character(Sys.Date()), '.pdf') %>%
  ggsave(p, width = 20, height = 9)


