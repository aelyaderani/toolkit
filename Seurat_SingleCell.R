library(Seurat)
library(dplyr)
library(Matrix)
library(emdi)
library(scatterplot3d)
library(plotly)
library(cellrangerRkit)
#library(rgl)

Sys.setenv("plotly_username"="tgen_winnie")
Sys.setenv("plotly_api_key"="45ceI7o7qbK5kGhqaw0D")
projectName <- "5Prime_ND_NUC"

Firstpart <- "/Users/aelyaderani/Desktop/pdf_seurat/"
pathSetwd <- paste(Firstpart,projectName, sep = "")
setwd(pathSetwd)

geneData <- ".....outs/filtered_feature_bc_matrix"


pbmc.data <- Read10X(data.dir = geneData)
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200, project = "10x_human_brain")

mito.features <- grep(pattern = "^MT-", x = rownames(x = pbmc), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = pbmc, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = pbmc, slot = 'counts'))

pbmc[['percent.mito']] <- percent.mito

pdf(file=paste(projectName,"_VlnPlot.pdf",sep=""))
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

pdf(file=paste(projectName,"_MitoPlot.pdf",sep=""))
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mito")
dev.off()

pdf(file=paste(projectName,"_FeaturePlot.pdf",sep=""))
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mito < 0.05)

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(object = pbmc, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
#VlnPlot(object = pbmc, features = c("percent.mito"), ncol = 1)
#VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc), vars.to.regress = c("nCount_RNA", "percent.mito"))
#VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
pbmc <- ProjectDim(object = pbmc)

v <- pbmc@reductions$pca@stdev^2
PCA_percentage <-v/sum(v)
rm(v)
print(PCA_percentage)

pdf(file=paste(projectName,"_elbowPlot.pdf",sep=""))
ElbowPlot(object = pbmc)
dev.off()
#------------------------------------------------- JackStraw
#pbmc <- JackStraw(object = pbmc, num.replicate = 100)
#pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
#JackStrawPlot(object = pbmc, dims = 1:20)

#PCApercent <- 0
#NumPCaForWhile <- 13
#countNUM <- 1
#while (countNUM <= NumPCaForWhile) 
#  {
#  PCApercent = PCApercent + PCA_percentage[countNUM]
#  countNUM=countNUM+1
#  }
#print(PCApercent)

pbmc <- FindNeighbors(object = pbmc, dims = 1:10)
pbmc <- FindClusters(object = pbmc, resolution = 0.4)
pbmc_tsne <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
pbmc_umap <- RunUMAP(object = pbmc, reduction = "pca", assay = "RNA", dims = 1:10)

#DimPlot(pbmc, reduction.use = "umap")

pdf(file=paste(projectName,"_tsne_noName.pdf",sep=""))
DimPlot(object = pbmc_tsne, reduction = 'tsne')
dev.off()

pdf(file=paste(projectName,"_umap_noName.pdf",sep=""))
DimPlot(object = pbmc_umap, reduction = 'umap')
dev.off()

#########################-----------------------------3D TSNE
#pbmc <- RunTSNE(object = pbmc, reduction.use = "pca", dims.use = 1:6, dim.embed = 3)
pbmc <- RunUMAP(object = pbmc, reduction = "pca", assay = "RNA", dims = 1:10, n.components = 3L)

#dr <- pbmc[["tsne"]]
dr <- pbmc[["umap"]]

umap_1 <- dr@cell.embeddings[,1]
umap_2 <- dr@cell.embeddings[,2]
umap_3 <- dr@cell.embeddings[,3]

graphclust <- as.factor(pbmc@active.ident)
# get number of clusters
n_graphclust <- length(levels(graphclust))
# setup color palette 
colors_graphclust <- colorRampPalette(brewer.pal(10, "Set3"))(n_graphclust)


dat <- data.frame(umap.1 = umap_1,umap.2 = umap_2, umap.3 = umap_3)
dat$FID <- seq.int(nrow(dat))
#[,grepl("-1" ,colnames(T2_raw_data))]
p <- plot_ly(dat,x=~`umap.1`,y=~`umap.2`,z=~`umap.3`,color=graphclust, colors=colors_graphclust, marker=list(size=3.5), width = "1300", height = "1000") 
p

api_create(p, filename = projectName)


#########################----------------------------- END
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.35)

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
class_cluster_list <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


dat_final <- data.frame(gene=class_cluster_list$gene, cluster=class_cluster_list$cluster, avg_logFC=class_cluster_list$avg_logFC, p_val=class_cluster_list$p_val, p_val_adj=class_cluster_list$p_val_adj)


write.csv(dat_final,file = paste(pathSetwd,"/",projectName,".csv", sep = ""), row.names=FALSE)

#---------------------------------- PART 2!!
pdf(file=paste(projectName,"_MarkerGenes.pdf",sep=""))
FeaturePlot(object = pbmc, features = c("MT1X","MT-ATP6","GRIA2","MT-ND4","HIF3A","GAD1"))
dev.off()
#FeaturePlot(object = pbmc, features = c("THY1", "GAD1","SLC17A7", "MBP", "ALDOC", "CLDN5", "ACTA2", "RASGRF2", "RORB", "PLCXD2", "FOXP2", "NR4A2"))
#FeaturePlot(object = pbmc, features = c("CACNA1C","GRIK2","CNTNAP2","KCND2","GRIK3","ROBO2","KCNH7","GRM8","SEMA6D"))
pdf(file=paste(projectName,"_MarkerGenesRidge.pdf",sep=""))
RidgePlot(pbmc, features= c("MT1X","MT-ATP6","GRIA2","MT-ND4","HIF3A","GAD1"))
dev.off()
#---------------------------------- PART 3!!
new.cluster.ids <- c("endothelial","endothelial","Pyramidal Neuron (PFC)","endothelial (all_Dif_Mito)","astrocytes","GABAergic neurons")
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids)

pdf(file=paste(projectName,"_tsne_WithName.pdf",sep=""))
DimPlot(object = pbmc, reduction = 'tsne', label = TRUE, pt.size = 0.4, label.size = 3.8)
dev.off()
#---------------------------------- PART 4!!
markers.to.plot <- c("MT1X","MT-ATP6","GRIA2","MT-ND4","HIF3A","GAD1")
pdf(file=paste(projectName,"_MarkerGenes_Percent.pdf",sep=""))
DotPlot(pbmc, features = markers.to.plot,x.lab.rot = T, plot.legend = T, dot.scale = 20, do.return = T)
dev.off()
#---------------------------------- Testing for NUC & WC
pdf(file=paste(projectName,"_MarkersForNuc.pdf",sep=""))
FeaturePlot(object = pbmc, features = c("MALAT1","SNHG11","MEG3","SNAP25","CALM1","CALM2","RTN1"))
dev.off()

pdf(file=paste(projectName,"_MarkersForWC.pdf",sep=""))
FeaturePlot(object = pbmc, features = c("MT-ND1","MT-ND2","MT-ND4"))
dev.off()

markers.to.plot <- c("MALAT1","SNHG11","MEG3","SNAP25","CALM1","CALM2","RTN1")
pdf(file=paste(projectName,"_MarkersForNuc_Percent.pdf",sep=""))
DotPlot(pbmc, features = markers.to.plot,x.lab.rot = T, plot.legend = T, dot.scale = 20, do.return = T)
dev.off()

markers.to.plot <- c("MT-ND1","MT-ND2","MT-ND4")
pdf(file=paste(projectName,"_MarkersForWC_Percent.pdf",sep=""))
DotPlot(pbmc, features = markers.to.plot,x.lab.rot = T, plot.legend = T, dot.scale = 20, do.return = T)
dev.off()

#-------------------- Everything bellow is in testing. DO NOT USE!!!
#---------------------------------- Fix Cluster
new.cluster.ids <- c("0", "1","2", "3","0","4","5","6","7")
names(x = new.cluster.ids) <- levels(x = pbmc)
pbmc <- RenameIdents(object = pbmc, new.cluster.ids)
DimPlot(object = pbmc, reduction = 'tsne', label = TRUE, pt.size = 0.5, label.size = 4)
#---------------------------------- scater/cowplot
library(scater)
library(cowplot)
pbmc_sce <- Convert(from = pbmc, to = "sce")
p1 <- plotExpression(object = pbmc_sce, features = "MALAT1", x = "ident") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- plotPCA(object = pbmc_sce, colour_by = "ident")
plot_grid(p1, p2)
