library(Seurat)
library(ggplot2)
library(Hmisc)
library(clustree)
library("ggsci")
library("scales")
mypal <- c(pal_lancet("lanonc")(9)[c(2,3,4,5,7,6)], pal_futurama("planetexpress")(9)[c(1,5,7,3,8)], pal_jama()(7)[c(2,3,4,5,6,7)]);  show_col(mypal)

load("Allpbmc_k35l3.RData")
all.pbmc@meta.data$cellTypes<-Idents(all.pbmc)
T_NK <- subset(all.pbmc,  cellTypes %in% c("Memory CD4 T","Memory CD8 T","Naive CD8 T","Naive CD4 T","gdT","NK","NK & NKT","MAIT","ILC"))
dim(T_NK)
T_NK <- subset(T_NK,  celltype2 %nin% c("Platelet","Eryth","Plasmablast","Doublet","HSPC","cDC1","cDC2","pDC"))
T_NK <- subset(T_NK,  BCR %nin% c("HK","HL","HLK","L","K","H"))
# T_NK <- SCTransform(T_NK, verbose = FALSE)
T_NK <- FindVariableFeatures(T_NK,selection.method = "vst", nfeatures = 3000)
ig_genes <- grep("^IG[JHKL]", VariableFeatures(T_NK))
VariableFeatures(all.pbmc) <- VariableFeatures(T_NK)[-ig_genes]


T_NK <- RunPCA(T_NK, verbose = FALSE)
T_NK <- RunUMAP(T_NK, dims = 1:20, n.neighbors =55, verbose = FALSE) 
T_NK <- FindNeighbors(T_NK, dims = 1:20, verbose = FALSE)
T_NK <- FindClusters(T_NK, algorithm = "leiden",method = "igraph",verbose = FALSE,resolution = c(seq(0.6,1.0,0.1)))
T_NK <- FindClusters(T_NK, algorithm = "leiden",method = "igraph",verbose = FALSE,resolution = 0.7)
pdf("T_NK_clustree_RNA.pdf",height = 8, width = 13)
clustree(T_NK@meta.data, prefix = "RNA_snn_res.")
dev.off()

map(c(40,seq(from=45,by=80,length=10)) , function(x) { T_NK %>%  RunUMAP(n.neighbors = x,n.epochs=200,dims = 1:20) %>% DimPlot()}) %>% cowplot::plot_grid(plotlist = .,pt.size = 0.1)
map(c(5,seq(from=20,by=50,length=10)) , function(x) { T_NK %>%  RunUMAP(n.neighbors = 55,n.epochs=x,dims = 1:20) %>% DimPlot()}) %>% cowplot::plot_grid(plotlist = .)
map(c(2,seq(from=1,by=5,length=10)) , function(x) { T_NK %>%  RunUMAP(n.neighbors = 55,n.epochs=200,negative.sample.rate = x,dims = 1:20) %>% DimPlot()}) %>% cowplot::plot_grid(plotlist = .)

DimPlot(T_NK,group.by = c("RNA_snn_res.0.9","RNA_snn_res.1"), reduction = "umap",ncol=2,label = TRUE,pt.size = 0.1,order = T, repel = T)
DimPlot(T_NK,group.by = c("celltype2","SCT_snn_res.0.6"), reduction = "umap",ncol=2,label = TRUE,pt.size = 0.1,order = T, repel = T)
DimPlot(T_NK,group.by = c("BCR","TCR"), reduction = "umap",ncol=2,label = TRUE,pt.size = 0.5,order = T, repel = T)
pdf("T_NK_Monaco_DICE.pdf",height = 9, width = 20)
DimPlot(T_NK,group.by = c("Monaco","DICE"), reduction = "umap",ncol=2,label = TRUE,pt.size = 0.1,order = T, repel = T)
dev.off()

DimPlot(T_NK,reduction = "umap",label = TRUE,pt.size = 1,order = T, repel = T)
DimPlot(T_NK,split.by = "RNA_snn_res.1", reduction = "umap",ncol=3,label = TRUE,pt.size = 0.5,order = T, repel = T)
DimPlot(T_NK,reduction = "umap",split.by="celltype2",label = TRUE,pt.size = 0.5,ncol = 8,order = T, repel = T)
Idents(T_NK)<-T_NK[["RNA_snn_res.1"]]
#Remove Tcell-T_NK complexes 



pdf("T_NK_celltype2.pdf",height = 6, width = 7)
DimPlot(T_NK,group.by = "celltype",reduction = "umap",label = TRUE,pt.size = 0.5, repel = T,cols = mypal)
dev.off()

T_NK <- CellCycleScoring(object = T_NK, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

pdf("CellCycle_T_NK1.pdf",height = 6, width = 7)
DimPlot(T_NK ,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 0.5, cols =c("#FF8000FF","blue","#4D9221"))
dev.off()



# DE analysis
# scale data is used
merged.markers <- FindAllMarkers(T_NK, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0,slot='data')#
merged.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
Top10_markers<-merged.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(merged.markers,file = "T_NK_markerGene_cluster_noNK.csv")


T_NK_feature<-c("CD14","S100A12","S100A8","LYZ",
                "FCGR3A","TCF7L2","SIGLEC10","CDKN1C",
                "FCER1A","CD1C","CLEC10A","FLT3",
                "LILRA4","CLEC4C","IL3RA","TCF4",
                "CD34","CYTL1","IL1B","PRSS57")                
svg("heatmap_T_NK.svg",height = 3.5, width = 6)

DoHeatmap(subset(T_NK1, downsample = 100), features = T_NK_feature, size = 3,slot='scale.data', raster=F)+ scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")) ) + guides(color=FALSE)
dev.off()


T_NK$celltype<-factor(T_NK$celltype,levels=rev(levels(T_NK)))
clus_feature<-c("CD3E","CD3D", "CD4", "CCR7","TCF7", "CD8A", "CD8B", "PTPRC", 
                "CD44", "SELL", "IL7R", "CCR6", "S100A11", "GZMK", "HLA-DRA", 
                "CCL5", "GZMH", "NKG7", "TRAV1-2", "CXCR6", "FOXP3", "TIGIT",
                "TRDV2","TRGV9","TRDC","NCAM1", "GZMB", "GNLY", "PRF1",  
                "KLRB1","KLRC2","FCGR3A", "XCL1", "XCL2","STAT3","STAT4")


pdf("T_NK_clusterGene2.pdf",height = 4.5, width = 15)
DotPlot(T_NK, features = clus_feature, group.by ="celltype",scale = T,scale.by = "size") + RotatedAxis()+
  scale_colour_gradient2(low="darkblue",mid="white",high="#c70606",midpoint=0)
dev.off()
DotPlot(T_NK[,T_NK$celltype %in% c("CD4+TCM","CD4+TEM","CD4+TCM/TEM","CD8+TCM","CD8+TEM")], features = clus_feature, group.by ="sample_type",scale = T,scale.by = "size") + RotatedAxis()+
   scale_colour_gradient2(low="darkblue",mid="white",high="#c70606",midpoint=0)



#T_NK@meta.data$seurat_clusters<-Idents(T_NK)
T_NK@meta.data$seurat_clusters<-T_NK$RNA_snn_res.0.9
Idents(T_NK)<-T_NK$RNA_snn_res.0.9
# Annotate the clusters
# After remove c12 & c17, reclustered with resolution 0.9

table(Idents(T_NK),T_NK$celltypist)
T_NK$celltype<-Idents(T_NK)

T_NK$celltype<-factor(T_NK$celltype,levels=rev(c("CD4+T Naive","CD8+T Naive","CD4+TCM",
                                             "CD4+TCM/TEM","CD4+TEM","CD8+TCM","CD8+TEM", 
                                             "MAIT","Treg","gdT", "CD56+NK","CD56bright NK",
                                             "NKT like","Cycling_T")))
all.pbmc$condition<-factor(all.pbmc$condition, levels=c("PsO","HC"))
T_NK$condition<-factor(T_NK$condition, levels=c("PsO","HC"))
T_NK@meta.data %>%
  group_by(condition,orig.ident,celltype) %>%
  dplyr::summarise(n=n())%>%
  spread(celltype,n)%>%
  replace(is.na(.), 0)%>%
  mutate_if(is.numeric, funs(. + 10^-100))%>%
  gather("celltype","n", -c(1,2))%>%
  mutate(C = sum(n)) %>%
  mutate(percent = n/C*100) %>%
  ggplot(aes(x=percent, y=condition)) +
  geom_boxplot(width = 0.4,position=position_dodge(0.3))+ #绘制箱线图
  geom_jitter(aes(fill= condition),position = position_jitter(width = .01),shape = 21,size=2.5)+#aes(fill= orig.ident)
  #scale_fill_manual(values=c(mypal_HC,mypal_PsO))+
  facet_wrap(~celltype, scales="free_y", ncol=4) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("T_NK_celltype_comparison.pdf",height = 7, width = 7)  


all.pbmc@meta.data %>%
  group_by(condition,orig.ident,metaCelltype) %>%
  dplyr::summarise(n=n())%>%
  spread(metaCelltype,n)%>%
  replace(is.na(.), 0)%>%
  mutate_if(is.numeric, funs(. + 10^-100))%>%
  gather("metaCelltype","n", -c(1,2))%>%
  mutate(C = sum(n)) %>%
  mutate(percent = n/C*100) %>%
  ggplot(aes(x=percent, y=condition)) +
  geom_boxplot(width = 0.4,position=position_dodge(0.3))+ #绘制箱线图
  geom_jitter(aes(fill= condition),position = position_jitter(width = .01),shape = 21,size=2.5)+#aes(fill= orig.ident)
  #scale_fill_manual(values=c(mypal_HC,mypal_PsO))+
  facet_wrap(~metaCelltype, scales="free_y", ncol=5) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("all_celltype_comparison.pdf",height = 9, width = 9)  

#mypal <- c(pal_lancet("lanonc")(9)[c(1,2,3,4,5,6,7,9)], pal_futurama("planetexpress")(9)[c(1,3,5,7,9)], pal_jama()(7)[c(1,2,4,5,6,7)]);  show_col(mypal)

T_NK@meta.data %>% 
  group_by(orig.ident, celltype) %>% 
  dplyr::summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)%>%
  ggplot(aes(x = orig.ident, y = percent, fill = celltype))+
  scale_fill_manual(values=mypal)+
  geom_bar(stat = "identity")
ggsave("T_NK_celltype_barchart1.pdf",height = 6, width = 7)  

library(rstatix)
library(readr)
library(dplyr)
r1<-all.pbmc@meta.data %>%
  #filter(orig.ident!="M1")  %>%
  group_by(condition,orig.ident,metaCelltype)   %>%
  summarise(n=n())%>%
  spread(metaCelltype,n)%>%
  replace(is.na(.), 0)%>%
  mutate_if(is.numeric, funs(. + 10^-100))%>%
  gather("metaCelltype","n", -c(1,2))%>%
  mutate(C = sum(n)) %>%
  mutate(percent = n/C*100) %>%
  group_by(metaCelltype) %>% 
  rstatix::wilcox_test(percent ~ condition, data = .) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")



## DE analysis between condition for each cluster 
## about slot:
## slot RNA:
# counts: Stores unnormalized data such as raw counts or TPMs
# data: Normalized data matrix
# scale.data: Scaled data matrix
## slot SCT:
# counts: corrected counts
# data: log1p(counts)
# scale.data: pearson residuals


T_NK@meta.data$sample_type <- paste(T_NK@meta.data$condition, T_NK@meta.data$celltype, sep = "_")
#marker_norm<-FindMarkers(T_NK, group.by="sample_type",ident.1 = "PsO_1", ident.2 = "HC_1", logfc.threshold = 0, min.pct = 0.1,slot = "data")
#marker_norm<-FindMarkers(T_NK, group.by="sample_type",ident.1 = "PsO_CD4+T Naive", ident.2 = "HC_CD4+T Naive", logfc.threshold = 0, min.pct = 0.0,slot = "data")

#marker_norm$cellType<-"CD4+T Naive"
marker_norm<-NULL
for (i in unique(T_NK$celltype)){
  ident.1= paste0("PsO_",i)
  ident.2= paste0("HC_",i)
  markerx<-FindMarkers(T_NK, group.by="sample_type",ident.1 = ident.1, ident.2 = ident.2, logfc.threshold = 0, min.pct = 0.0,slot = "data")
  markerx$cellType<-i
  markerx$gene<-rownames(markerx)
  marker_norm<-rbind(marker_norm,markerx)
}

marker_norm<-marker_norm[!duplicated(marker_norm), ]
write.csv(marker_norm, file = "T_NK_Gene_PsOvsHC1.csv")
marker_norm1<-marker_norm[marker_norm$p_val_adj<0.05,]
marker_norm_up<-marker_norm1[marker_norm1$avg_log2FC>0,]
marker_norm_down<-marker_norm1[marker_norm1$avg_log2FC<0,]



# Functional annotation of DEGs
library(clusterProfiler)
library(org.Hs.eg.db)
for (i in unique(T_NK@meta.data$cellTypes)) {
  list<- rownames(marker_norm_up)[marker_norm_up$cellType==i]
  GO<- enrichGO(list, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
  pdf(paste0(i,"_GO_up.pdf"), width=7, height=6)
  d<-dotplot(GO,showCategory=20)
  print(d)
  dev.off()
  write.csv(GO, file = paste0(i,"_GO_up.csv"))
}
GO<- enrichGO(rownames(marker_norm)[marker_norm$avg_log2FC<0], OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
dotplot(GO,showCategory=20)
list = bitr(rownames(marker_norm)[marker_norm$avg_log2FC<0], fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Hs.eg.db, drop = TRUE)


for (i in unique(T_NK@meta.data$cellTypes)) {
  list = bitr(rownames(marker_norm_down)[marker_norm_down$cellType==i], fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Hs.eg.db, drop = TRUE)
  kegg <- enrichKEGG(list$ENTREZID, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 1,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  pdf(paste0(i,"_kegg_down.pdf"), width=7, height=6)
  d<-dotplot(kegg,title="Enrichment KEGG_dot")
  print(d)
  dev.off()
  write.csv(kegg, file = paste0(i,"_kegg_down.csv"))
}
list = bitr(rownames(marker_norm)[marker_norm$cellType=="CD14+ T_NK"], fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Hs.eg.db, drop = TRUE)
kegg <- enrichKEGG(list$ENTREZID, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 1,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
dotplot(kegg,title="Enrichment KEGG_dot")

##GSEA
library(msigdbr)
library(enrichplot)
dir.create("./gsea")
for (i in unique(T_NK@meta.data$cellTypes)) {
  gselist <- marker_norm$avg_log2FC[marker_norm$cellType==i]
  names(gselist) <- rownames(marker_norm)[marker_norm$cellType==i]
  gselist <- sort(gselist, decreasing = T)
  
  m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
    dplyr::select(gs_name, gene_symbol)
  m_t2g$gs_name<- gsub("GO","",m_t2g$gs_name)
  gsea <- GSEA(gselist, TERM2GENE = m_t2g)
  write.csv(gsea@result,file=paste0("./gsea/",i,"_gsea_c5.csv"))
  for (gs in seq(1:nrow(gsea))){
    label<- paste0("./gsea/",i,gsea$ID[gs],".pdf")
    pdf(label, width=14,height=8)  
    p <- gseaplot2(gsea, geneSetID = gsea$ID[gs], title = gsea$Description[gs], ES_geom = "line",rel_heights = c(1.5, 0.5, 1), base_size = 21, subplots = 1:3, pvalue_table = T)
    print(p)
    dev.off()
  }
}
save(T_NK, file="T_NK.RData")


