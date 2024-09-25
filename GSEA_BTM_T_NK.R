library(tidyverse)
library(tibble)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(scales)
library(enrichplot)
library(NCmisc) #v1.1.5

# See sessionInfo() at end of script for version information.

# Modify path below. All subsequent paths are relative.

#==================================================================================================
df <- read.csv("T_NK_Gene_PsOvsHC1.csv")
df<-marker_norm
# Add z-score based on two sided null hypothesis.
df$z <- p.to.Z(df$p_val) * sign(df$avg_log2FC)
df$z.adj <- p.to.Z(df$p_val_adj) * sign(df$avg_log2FC)

# Randomize rows to reduce overplotting issues.
df <- df[sample(nrow(df)), ]

levels=c("CD4+T Naive","CD8+T Naive","CD4+TCM",
         "CD4+TCM/TEM","CD4+TEM","CD8+TCM","CD8+TEM", 
         "MAIT","Treg","gdT", "CD56+NK","CD56bright NK",
         "NKT like","Cycling_T")
df$celltype_factor <- factor(df$cellType,  levels=levels, ordered=T)
#==================================================================================================
names(mypal[1:14]) <- levels(df$celltype_factor)

# Make custom color column to facilitate grey coloring by threshold.
df$col<-NA
df$col[df$p_val_adj > 0.05] <- "#D3D3D3"
df$col[df$avg_log2FC > 0 & df$p_val_adj < 0.05] <- "#c70606"
df$col[df$avg_log2FC < 0 & df$p_val_adj < 0.05] <-"#3C5488FF"

col <- mypal[df$celltype_factor]
 # grey
df$col <- as.factor(col)

q <- ggplot(df, aes(x = celltype_factor, y = z, color = col)) +
  geom_jitter(width = 0.40, alpha = .55, size = 1) +
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 10), axis.title.x = element_blank()) +
  ggtitle("Changes in gene expression with condition") +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  labs(y = "Z-score") +
  theme(legend.position="none") +
  scale_color_manual(values = levels(df$col)) +
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1))+
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
print(q)
ggsave("TNK_celltype_DEG.pdf",height = 6, width = 8)

library(tidyverse)
library(biomaRt)
library(fgsea)
library(DT)
library(msigdbr)
library(org.Hs.eg.db)
library(msigdf)
library(NCmisc)
library(clusterProfiler)
library(GSA)

BTM <-  GSA.read.gmt("BTM_for_GSEA_20131008.gmt")
len_vec=c()           # Now create a vector for containing the length of genes at each position
len_vec[1] = 3
for(i in 1:length(BTM$genesets)){len_vec[i] <- c(length(BTM$genesets[[i]]))}
pathway_vec <- unlist(Vectorize(rep.int)(BTM$geneset.names, len_vec),use.names = FALSE) # Now create a vector for all the pathways in the data 
desired_df <- as.data.frame(cbind(pathway_vec,unlist(BTM$genesets,use.names = FALSE))) # This gives your desired dataframe
head(desired_df)



colnames(df)[7]<-"celltype"
gsea_all <- NULL
df$celltype<-gsub("/","",df$celltype)
levels<-gsub("/","",levels)
for (cell in levels[13:14]) {
  
  print(cell)
  
  # Filter for celltype, reduce to just gene and DE p value statistic.
  df2 <- df %>% dplyr::filter(celltype == cell)
  print(dim(df2))
  
  ##GSEA
  gselist <- df2$z
  names(gselist) <- df2$gene
  gselist <- sort(gselist, decreasing = T)
  
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, gene_symbol)
  #m_t2g$gs_name<- gsub("GO","",m_t2g$gs_name)
  gsea <- GSEA(gselist, TERM2GENE = desired_df)
  #write.csv(gsea@result,file=paste0("gsea_H.csv"))
  for (gs in seq(1:nrow(gsea))){
    label<- paste0(cell,gsea$ID[gs],"BTM.pdf")
    pdf(label, width=7,height=4)  
    p <- gseaplot2(gsea, geneSetID = gsea$ID[gs], title = gsea$Description[gs], ES_geom = "line",rel_heights = c(1.5, 0.5, 1), base_size = 21, subplots = 1:3, pvalue_table = T)
    print(p)
    dev.off()}

  fgseaResTidy <- gsea@result %>%
    dplyr::select(-leading_edge, -enrichmentScore) %>% 
    as_tibble() %>%
    arrange(desc(NES))
  print(head(fgseaResTidy))
  
  fgseaResTidy <- as.data.frame(fgseaResTidy)
  fgseaResTidy$`celltype` <- cell
  
  fname <- paste0(cell, "BTM_fgsea_", Sys.Date(),".txt")
  
  # Save individual tables
  write.table(fgseaResTidy, file=fname, sep = "\t")
  gsea_all <- rbind(gsea_all, fgseaResTidy)
}

# Save complete dataframe
fname <- paste0("NK_T_BTM_fgsea_", Sys.Date(),".txt")
write.table(gsea_all, file=fname, sep = "\t")

# Optional: Inspect results in browser.
gsea_all %>% 
  arrange(-NES) %>% 
  DT::datatable()

fgseaResTidy%>% 
  arrange(-NES) %>% 
  DT::datatable()
d<-gsea_all
d$celltype <- factor(d$celltype,  levels=levels, ordered=T)
d$Description <- gsub("HALLMARK_", "", d$Description)
d$Description <- gsub("_", " ", d$Description)
d$NES[d$padj > 0.05] <- NA # grey out non-significant
#d <- d %>% filter(!celltype %in% c("Neurons", "Mural_cells"))

p <- ggplot(d, aes(x=celltype, y=Description, fill=NES)) + geom_tile(colour="white",size=0.25) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), na.value="grey80") +
  theme_minimal() +
  labs(x="", y="") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
ggsave(paste0("BTM_heatmap_", Sys.Date(), ".pdf"), height=7, width=7)
sessionInfo()
# R version 3.4.3 (2017-11-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.14

# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     

# other attached packages:
#  [1] clusterProfiler_3.6.0 DOSE_3.4.0            NCmisc_1.1.5         
#  [4] msigdf_5.2            org.Mm.eg.db_3.5.0    AnnotationDbi_1.40.0 
#  [7] IRanges_2.12.0        S4Vectors_0.16.0      Biobase_2.38.0       
# [10] BiocGenerics_0.24.0   DT_0.5                fgsea_1.4.1          
# [13] Rcpp_1.0.0            biomaRt_2.34.2        bindrcpp_0.2.2       
# [16] eulerr_5.1.0          reshape2_1.4.3        forcats_0.3.0        
# [19] stringr_1.3.1         dplyr_0.7.6           purrr_0.2.5          
# [22] readr_1.1.1           tidyr_0.8.1           tibble_2.0.1         
# [25] ggplot2_3.1.0         tidyverse_1.2.1   

