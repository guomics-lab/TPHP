# 0.config ------------
# rm(list = ls())
# pacman::p_unload(pacman::p_loaded(), character.only = T)

wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
source('../source/source_code.R')


# 1.Read data --------
info5 <- rio::import('../0_process_DIA-NN_data/input/20251009_PUH_sample_information_3005files_V6.xlsx')
info5 %<>% filter(!Tissue_heterogeneity_abnormal, !Low_protein_IDs, !FileName %in% pancancer_low_pearson)

# Inputs
res_ts <- readRDS('../4_tissue_specificity_analysis/output/20251201_tissue_Specificity_naive_05minImputated_3.rds')

# df_enrich <- res_ts$ts.Enrich
df_tis_en <- res_ts$ts.tisEnrich
df_tis_en_ <- df_tis_en %>% inner_join(info5 %>% distinct(class_abbr, class_abbr))
ts.med <- res_ts$med.matrix

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


datall <- rio::import('../0_process_DIA-NN_data/input/mapped_pg_matrix_2856_13609.csv') %>% column_to_rownames('V1')

# 2.Intensity rank --------
# all<-read_xlsx("//172.16.13.136/share/members/yuel/2022/tables/20220719_rm_randomNA_normal.xlsx")
# info<-read_xlsx("//172.16.13.136/share/members/yuel/2022/tables/20230223_PUH_sample_information_1781files_info_edited_v7.xlsx")
# dim(info) # 1781 18
# head(colnames(all),n=35) #1:15 for sample annotation
# table(all$tissue_type)

ai <- datall %>% as.data.frame() %>% 
  rownames_to_column('FileName') %>% 
  inner_join(info5, .)
rownames(ai)<-ai$FileName
head(colnames(ai),n=40)
table(ai$class_abbr) #80
tissue<-unique(ai$class_abbr)
npm<-apply(ai[,-c(1:34)], c(1,2),as.numeric) %>% log2
imp<-min(npm,na.rm = T)+log2(0.5)
tissue_list<-list()
dotlist<-list()
boxlist<-list()
for (i in 1:length(tissue)){
  # i=1
  cat('Processing loop ', i, '...\r')
  t1<-2^npm[which(ai$class_abbr==tissue[i]),,drop=F]
  t2<-t1[,apply(t1, 2, function(v) sum(is.na(v))/nrow(t1)<0.5),drop=F]
  # dim(t2)
  t2[is.na(t2)]<-2^imp
  t3<-apply(t2,2, summary)
  tissue_list[[i]]<-as.data.frame(t2)
  boxlist[[i]]<-as.data.frame(t3)
  su<-sum(t2)
  t4<-data.frame(sum=apply(t2, 2, sum))
  t4$per<-log10(t4$sum/su)
  t4$rank<-rank(-t4$sum)
  dotlist[[i]]<-as.data.frame(t4)
}
names(tissue_list)=names(dotlist)=names(boxlist)<-tissue
# .write_excel(tissue_list,"20251208_tissue_list.xlsx")
# .write_excel(dotlist,"20251208_dot_list.xlsx")
# .write_excel(boxlist,"20251208_box_list.xlsx")

# setwd("//172.16.13.114/share/members/yuel/1project/1ref")
# ulhen<-read_xlsx("20230207TPHP_51tissue_ulhen_classification.xlsx")
ulhen <- res_ts$ts.all %>%
  pivot_wider(id_cols = class_abbr, names_from = Protein, values_from = `Tissue specificity`) %>%
  dplyr::rename(tissue_type = class_abbr)

pdf("20251208_abundance_curve.pdf", width = 18, height = 10)
par(mfrow=c(3, 4))
ra<-list()
for (i in 1:nrow(ulhen)){
  # i=1
  cat('Generating plot ', i, '...\r')
  t<-ulhen$tissue_type[i]
  dp<-dotlist[[t]]
  dp<-dp[order(dp$rank,decreasing = F),,drop=F]
  dp$pert<-cumsum(10^(dp$per))
  en<-t(ulhen[which(ulhen$tissue_type==t),
              intersect(rownames(dp),colnames(ulhen)),
              drop=F])
  plot(dp$rank, dp$pert, col = '#00000033', pch = 19,xlab = 'Abundance rank', ylab = '%. of total abundance', main = t)
  abline(h = 0.8, v = min(dp$rank[dp$pert>=0.8]), lty = 2, lwd = 1)
  ex  <- en=="tissue enriched"
  points(dp$rank[ex], dp$pert[ex], col = 1,bg = brewer.pal(9,"YlOrRd")[6], pch = 21, cex = 2)
  ea<- en=="Expressed in all tissues"
  points(dp$rank[ea], dp$pert[ea], col = 1,bg = brewer.pal(11,"RdBu")[9], pch = 21, cex = 1)
  eh<-en=="tissue enhanced"
  mx<-en=="mixed"
  ge<-en=="group enriched"
  nt<-en=="not detected"
  ra[[i]]<-c(sum(10^dp$per[ex]),sum(10^dp$per[ea]),sum(10^dp$per[eh]),sum(10^dp$per[mx]),sum(10^dp$per[ge]),sum(10^dp$per[nt]))
}
dev.off()

# per_tab<-do.call(rbind,ra)
# colnames(per_tab)<-c("tissue enriched","Expressed in all tissues","tissue enhanced","mixed","group enriched","not detected")
# rownames(per_tab)<-ulhen$tissue_type
# write.csv(per_tab,"tissue_classification_abundance_percentage.csv")

# 3.GOBP enrichment -------
##GOBP enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
bf.list <- as.list(rep(NA, length(tissue)))#<-list()
# for (i in 1:length(tissue)){
for (i in which(tissue %in% c('CE', 'LI'))){
  # i=1
  cat(i, '...\r')
  bx<-as.data.frame(t(boxlist[[i]]))
  # entrez<-bitr(rownames(bx),fromType = "UNIPROT",toType = c("ENTREZID","GENENAME"),OrgDb = "org.Hs.eg.db")
  bx$UNIPROT<-rownames(bx)
  # bf<-left_join(bx,entrez,by="UNIPROT")
  bf.list[[i]]<-enrichGO(bx$UNIPROT,OrgDb = "org.Hs.eg.db",keyType = "UNIPROT",ont = "BP",pvalueCutoff = 0.05)%>% as.data.frame()
  pro.list<-strsplit(bf.list[[i]]$geneID,split = "/")
  bf.list[[i]]$medain<-sapply(pro.list,function(v){median(bx[v,]$Median)})
  bf.list[[i]]$sum<-sapply(pro.list,function(v){sum(bx[v,]$Median)})  
}
names(bf.list)<-tissue
# .write_excel(bf.list,"tissue_GOBP.xlsx")

pdf("20230717_GOBP_abundance_curve_sum.pdf", width = 18, height = 10)
par(mfrow=c(3, 4))
ra<-list()
for (i in 1:length(tissue)){
  
  t<-bf.list[[i]]
  t$rank2<-rank(-t$sum)
  # t$rank<-rank(-t$medain)
  gene.ratio<-sapply(strsplit(t$GeneRatio,split = "/"), function(v) as.numeric(v[1])/as.numeric(v[2]))
  t1<-t[t$p.adjust<0.001 & gene.ratio> 0.02,]
  # plot(t$rank, log10(t$medain), col = '#00000033', pch = 19,xlab = 'Abundance rank', ylab = 'log10 tansformed total abundance', main = tissue[i])
  # points(t1$rank, log10(t1$medain), col = 1,bg = brewer.pal(9,"YlOrRd")[6], pch = 21, cex = 2)
  # text(label=t1$Description[which.min(t1$rank)],x=t1$rank[which.min(t1$rank)]+120, y=log10(t1$medain[which.min(t1$rank)])-0.005)
  
  
  # # gene.ratio<-sapply(strsplit(t$GeneRatio,split = "/"), function(v) as.numeric(v[1])/as.numeric(v[2]))
  # t1<-t[t$p.adjust<0.001 & gene.ratio> 0.02,]
  plot(t$rank2, log10(t$sum), col = '#00000033', pch = 19,xlab = 'Abundance rank', ylab = 'log10 tansformed total abundance', main = tissue[i])
  points(t1$rank2, log10(t1$sum), col = 1,bg = brewer.pal(9,"YlOrRd")[6], pch = 21, cex = 2)
  text(label=t1$Description[which.min(t1$rank2)],x=t1$rank[which.min(t1$rank2)]+120, y=log10(t1$sum[which.min(t1$rank2)])-0.005)
  
}
dev.off()

# pacman::p_unload(clusterProfiler, enrichplot, ReactomePA, AnnotationDbi, BiocGenerics, IRanges, S4Vectors, Biobase)

# 4.Focus on liver and brain ------------------------------------------------
# protein - gene mapping from DDA library building result
dfprot <- libprot
dfprot$`Indistinguishable Protein IDs` <- 
  sapply(dfprot$`Indistinguishable Proteins`, function(x){
    if(x != '') {
      x <- str_split(x, ', ') %>% unlist() %>% str_split('\\|') %>% sapply(function(x) x[2]) %>% str_c(collapse = ';')
    }
    return(x)
  })
dfprot$`Protein Group` <- str_c(dfprot$`Protein ID`, dfprot$`Indistinguishable Protein IDs`, sep = ';') %>% str_remove(';$')
df_pg_gn <- dfprot %>% dplyr::select(`Protein ID`, `Protein Group`, Gene)


go2gene <- toTable(org.Hs.egGO2ALLEGS)
df_gene <- go2gene %>% dplyr::filter(Ontology == 'BP')
df_idmap <- clusterProfiler::bitr(df_gene$gene_id, fromType = "ENTREZID", toType = "UNIPROT", org.Hs.eg.db)
df_gene <- df_idmap %>% full_join(df_gene, by = c(ENTREZID = 'gene_id'))
df_idmap <- clusterProfiler::bitr(df_gene$ENTREZID, fromType = "ENTREZID", toType = "SYMBOL", org.Hs.eg.db)
df_gene <- df_idmap %>% full_join(df_gene)
setdiff(df_pg_gn$`Protein ID`, df_gene$UNIPROT) %>% length() # 1533 not matched


GOBP.summary <- rio::import('N_tissue_enriched_GOBP_network_plot_v1202.xlsx', sheet = 'GOBP.summary')
enrich.network <- rio::import('N_tissue_enriched_GOBP_network_plot_v1202.xlsx', sheet = 'enrich.network')

## liver ---------
t<-'LI'
dp<-dotlist[[t]]
dp<-dp[order(dp$rank,decreasing = F),]
dp$pert<-cumsum(10^(dp$per))
en<-t(ulhen[which(ulhen$tissue_type==t),intersect(rownames(dp),colnames(ulhen))])
ex  <- en=="Tissue enriched"
ea<- en=="Expressed in all tissues"
eh<-en=="Tissue enhanced"
mx<-en=="Mixed"
ge<-en=="Group enriched"
nt<-en=="Not detected"
writeClipboard(rownames(ex)[ex]) # for Metascape v3.5.20250701
# bx<-as.data.frame(t(boxlist[[t]]))
# bx$UNIPROT<-rownames(bx)
# bxGO<-enrichGO(bx$UNIPROT,OrgDb = "org.Hs.eg.db",keyType = "UNIPROT",ont = "BP",pvalueCutoff = 0.05,
#                readable = T, minGSSize = 1, maxGSSize = 30000)
# bxGO_filter<-DOSE::gsfilter(bxGO, min=10, max=5000)

metascape_go_top_liver <- c('GO:0019752', 'GO:0044282', 'GO:0016053', 'GO:0008202', 'GO:0006753', 'GO:0072329', 'GO:0006790', 'GO:0042180', 'GO:0008652', 'GO:0009069')
df_gene_liver <- df_gene %>% dplyr::filter(go_id %in% metascape_go_top_liver)
ex_GOmap_index <- which(rownames(ex) %in% intersect(df_gene_liver$UNIPROT, rownames(ex)[ex]))

pdf('abundance_curve_liver_v1202.pdf', width = 8.1, height = 6)
plot(dp$rank, dp$pert, col = '#00000033', pch = 19,xlab = 'Abundance rank', ylab = '%. of total abundance', main = str_to_title(t))
abline(h = 0.8, v = min(dp$rank[dp$pert>=0.8]), lty = 2, lwd = 1)
points(dp$rank[ea], dp$pert[ea], col = 1,bg = brewer.pal(11,"RdBu")[9], pch = 21, cex = 1)
points(dp$rank[ex], dp$pert[ex], col = 1,bg = brewer.pal(9,"YlOrRd")[6], pch = 21, cex = 2)
text(x = nrow(dp) * 0.8, y = 0.5, str_c('# Tissue enriched proteins: ', sum(ex)), col = 'black', cex = 1)
text(x = nrow(dp) * 0.8, y = 0.6, str_c('# Quantifiable proteins: ', nrow(dp)), col = 'black', cex = 1)
legend('bottomright', inset = .02, pch = 16, title = 'Specificity class', legend = c('Tissue enriched', 'Expressed in all'), col=c("#FC4E2A", "#4393C3"), lty=1:2, cex=0.6, pt.cex = 2)

protein_labeled <- c("Q7Z4W1", "P31327", "P05062", "P42765", "Q16822")
# dp_text <- dp[ex_GOmap_index, ][c("Q7Z4W1", "P31327", "P05062", "P42765", "Q16822"), ]
# dp_text <- df_gene_liver %>% dplyr::select(SYMBOL, UNIPROT, go_id) %>% distinct() %>% 
#   dplyr::group_by(SYMBOL, UNIPROT) %>% 
#   dplyr::summarise(go_id = str_c(unique(sort(go_id)), collapse = ';'), .groups = 'drop') %>%
#   merge(dp_text, by.x = 'UNIPROT', by.y = 'row.names') %>% arrange(rank)
# dp_text$note <- c('Glycolysis', 'Glyoxylate detoxification', 'Lipid metabolism', 'Gluconeogenesis', 'Steroid metabolism')
# dp_text$label <- str_glue("{dp_text$UNIPROT} ({dp_text$SYMBOL}; {dp_text$note})")


protein_labeled <- c(
  "Q9UBR1", "P23141", "Q7Z4W1", "Q14032"
)
dp_text <- dp_text[protein_labeled, ] %>%
  mutate(note = c(
    'Amino acid biosynthesis', 'Cellular ketone metabolism',
    'Monocarboxylic acid catabolism', 'Serine metabolism'
  ),
  label = str_glue("{UniprotID} ({Gene}; {note})"))

text(x = dp_text$rank + 1800, y = dp_text$pert, dp_text$label, col = 'black', cex = 1)
segments(x0 = dp_text$rank, y0 = dp_text$pert, x1 = dp_text$rank + 800, y1 = dp_text$pert, lwd = 1)
graphics.off()
rio::export(dp_text, 'abundance_curve_liver_v1202.xlsx')

## brain ---------
t<-'CE'
dp<-dotlist[[t]]
dp<-dp[order(dp$rank,decreasing = F),]
dp$pert<-cumsum(10^(dp$per))
en<-t(ulhen[which(ulhen$tissue_type==t),intersect(rownames(dp),colnames(ulhen))])
ex  <- en=="Tissue enriched"
ea<- en=="Expressed in all tissues"
eh<-en=="Tissue enhanced"
mx<-en=="Mixed"
ge<-en=="Group enriched"
nt<-en=="Not detected"
writeClipboard(rownames(ex)[ex])
# bx<-as.data.frame(t(boxlist[[t]]))
# bx$UNIPROT<-rownames(bx)
# bxGO<-enrichGO(bx$UNIPROT,OrgDb = "org.Hs.eg.db",keyType = "UNIPROT",ont = "BP",pvalueCutoff = 0.05,
#                readable = T, minGSSize = 1, maxGSSize = 30000)
# bxGO_filter<-DOSE::gsfilter(bxGO, min=10, max=5000)

metascape_go_top_brain <- c('GO:0050804','GO:0099536','GO:0050808','GO:0042391','GO:0050803','GO:0031175','GO:0099504','GO:0031344','GO:0060078','GO:0051966','GO:0099175','GO:0007613','GO:0048168','GO:0001956','GO:0098660')
# metascape_go_top_brain <- c('GO:0050804', 'GO:0031175')

df_gene_brain <- df_gene %>% dplyr::filter(go_id %in% metascape_go_top_brain)
ex_GOmap_index <- which(rownames(ex) %in% intersect(df_gene_brain$UNIPROT, rownames(ex)[ex]))

pdf('abundance_curve_brain_v1202.pdf', width = 8.1, height = 6)
plot(dp$rank, dp$pert, col = '#00000033', pch = 19,xlab = 'Abundance rank', ylab = '%. of total abundance', main = str_to_title(t))
abline(h = 0.8, v = min(dp$rank[dp$pert>=0.8]), lty = 2, lwd = 1)
points(dp$rank[ea], dp$pert[ea], col = 1,bg = brewer.pal(11,"RdBu")[9], pch = 21, cex = 1)
points(dp$rank[ex], dp$pert[ex], col = 1,bg = brewer.pal(9,"YlOrRd")[6], pch = 21, cex = 2)
text(x = nrow(dp) * 0.8, y = 0.5, str_c('# Tissue enriched proteins: ', sum(ex)), col = 'black', cex = 1)
text(x = nrow(dp) * 0.8, y = 0.6, str_c('# Quantifiable proteins: ', nrow(dp)), col = 'black', cex = 1)
legend('bottomright', inset = .02, pch = 16, title = 'Specificity class', legend = c('Tissue enriched', 'Expressed in all'), col=c("#FC4E2A", "#4393C3"), lty=1:2, cex=0.6, pt.cex = 2)

# GOBP.summary %>% filter(ID %in% metascape_go_top_brain)
dp_text <- enrich.network %>% filter(ID %in% metascape_go_top_brain) %>%
  select(ID, Description, Protein.list) %>% 
  separate_rows(Protein.list)
dp_text %<>% count(Protein.list) %>% filter(n == 1) %>% semi_join(dp_text, .)
# dp %>% rownames_to_column('Protein.list') %>% inner_join(dp_text) %>% count(ID)
# dp %>% rownames_to_column('Protein.list') %>% inner_join(dp_text) %>% View()
# dp %>% rownames_to_column('Protein.list') %>% inner_join(dp_text) %>%
#   group_by(Description) %>%
#   arrange()
dp_text <- dp %>% rownames_to_column('Protein.list') %>%
  inner_join(dp_text) %>% 
  inner_join(dfprot %>% distinct(`Protein ID`, Gene) %>%
               rename(Protein.list = `Protein ID`)) %>%
  as.data.frame() %>% set_rownames(.$Protein.list)

# protein_labeled <- c('Q9UQM7', 'Q05586', 'Q9Y566', 'Q9UPX8', 'Q9NZ94')
# dp_text <- dp[ex_GOmap_index, ][protein_labeled, ]
# dp_text <- df_gene_brain %>% dplyr::select(SYMBOL, UNIPROT, go_id) %>% distinct() %>% 
#   dplyr::group_by(SYMBOL, UNIPROT) %>% 
#   dplyr::summarise(go_id = str_c(unique(sort(go_id)), collapse = ';'), .groups = 'drop') %>%
#   merge(dp_text, by.x = 'UNIPROT', by.y = 'row.names') %>% arrange(rank)
# dp_text$note <- c('Memory', 'Ionotropic glutamate receptor', 'Behavior', 'Behavior', 'Behavior')
# dp_text$label <- str_glue("{dp_text$UNIPROT} ({dp_text$SYMBOL}; {dp_text$note})")
protein_labeled <- c(
  "Q92752", "P17600", "P51674"#, "Q8IW52", "O15083", "Q14831"
)
dp_text <- dp_text[protein_labeled, ] %>%
  mutate(note = c(
    'Synaptic transmission', 'Synaptic vesicle', 'Synapse structure'#, 'Synapse structure', 'Synaptic vesicle', 'Synaptic transmission'
  ),
         label = str_glue("{Protein.list} ({Gene}; {note})"))

text(x = dp_text$rank + 1800, y = dp_text$pert, dp_text$label, col = 'black', cex = 1)
segments(x0 = dp_text$rank, y0 = dp_text$pert, x1 = dp_text$rank + 800, y1 = dp_text$pert, lwd = 1)
graphics.off()
rio::export(dp_text, 'abundance_curve_brain_v1202.xlsx')


# for source data -----
ulhen.class.tbl <- ulhen %>% column_to_rownames('tissue_type') %>% 
  t() %>% as.data.frame() %>% rownames_to_column('protein') %>% 
  inner_join(dfprot %>%
               distinct(`Protein ID`, Gene) %>%
               rename(protein = `Protein ID`), .)
ulhen.med.tbl <- ts.med %>% column_to_rownames('class_abbr') %>% 
  t() %>% as.data.frame() %>% rownames_to_column('protein') %>% 
  inner_join(ulhen.class.tbl %>% select(protein, Gene), .)

list(specificity.median = ulhen.med.tbl,
     specificity.class = ulhen.class.tbl) %>% 
  .write_excel('source_data_for_tableS3.xlsx')


# Radarchart obsolete ------------------------------------------------------


# 
# library(fmsb)
# pdf("TPHP_carcinoma_20organs_14functions_radarchart_nafill_20220906.pdf", width = 30, height = 10)
# # radarchart(mat,
# #            axistype = 1, pcol = curve_colors, #maxmin = F,
# #            plwd = 1, plty = 1,
# #            cglcol = "gray", cglty = 1, axislabcol = "black",
# #            caxislabels = seq(min_range, max_range, length.out = 5),
# #            cglwd = 0.8,
# #            vlcex = 1.6
# # )
# # legend(x = 2, y = 1, legend = rownames(mat)[-(1:2)], bty = "n", pch = 20, col = curve_colors, text.col = "black", cex = 2, pt.cex = 2)
# # dev.off()
# sum.list<-sapply(bf.list, function(v){
#   v$Description
# })
# all_go<-reduce(sum.list,intersect)
# # x <- tibble::rownames_to_column(sum.list,"tissue")
# bp_tab<-as.data.frame(matrix(ncol=length(sum.list),nrow=length(all_go)))
# for (i in 1:length(sum.list)){
#   # i=1
#   t<-bf.list[[i]]
#   rownames(t)<-t$Description
#   bp_tab[,i]<-t[all_go,c("sum")]
# }
# colnames(bp_tab)<-tissue
# rownames(bp_tab)<-all_go
# mat <- as.matrix(log10(bp_tab))
# 
# fivenum(mat)
# min_range <- floor(fivenum(mat)[1])
# max_range <- ceiling(fivenum(mat)[5])
# mat <- rbind(max_range = max_range, min_range = min_range, mat) 
# # mat[is.na(mat)] <- min_range
# sum_intensity<-apply(mat[-c(1,2),], 1, sum)
# rank3<-rank(-sum_intensity)
# highlight<-rownames(mat)[-c(1:2)][which(rank3<=10)]
# # set color
# qual_color_pals <- brewer.pal.info %>% filter(category == "qual")
# all_colors <- Map(brewer.pal, qual_color_pals$maxcolors, rownames(qual_color_pals)) %>%
#   unlist() %>%
#   unique()
# set.seed(1)
# curve_colors <- c(rep("grey",127),sample(all_colors,10))
# names(curve_colors)<-c(setdiff(all_go,highlight),highlight)
# mat<-mat[c("max_range","min_range",names(curve_colors)),]%>%as.data.frame()
# pdf("TPHP_carcinoma_20organs_137functions_radarchart_20220906.pdf", width = 30, height = 10)
# radarchart(mat,
#            axistype = 1, pcol = curve_colors,# maxmin = F,
#            plwd = 1, plty = 1,
#            cglcol = "gray", cglty = 1, axislabcol = "gray",
#            caxislabels = seq(min_range, max_range, length.out = 5),
#            cglwd = 0.8,
#            vlcex = 1.6
# )
# legend(x = 2, y = 1, legend = rownames(mat)[-(1:2)], bty = "n", pch = 20, col = curve_colors, text.col = "black", cex = 2, pt.cex = 2)
# dev.off()
# # 