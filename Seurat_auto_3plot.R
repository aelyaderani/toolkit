library(Seurat)
library(dplyr)
library(Matrix)
library(emdi)
library(scatterplot3d)
library(plotly)
library(cellrangerRkit)
library(profvis)

setwd("/Users/aelyaderani/Desktop/auto_seurat")
filepathlist <- read_excel("/Users/aelyaderani/Desktop/listoffiles.xlsx",col_names = FALSE)

StartPath <- "/Volumes/NOMIS/scRNA_Analysis/CellRanger_v3.0.2/"
endpath <- "/outs/filtered_feature_bc_matrix"
Rnum <- 1

while (Rnum <= nrow(filepathlist)){

filepath <- filepathlist[Rnum,1]


geneData <- paste(StartPath,filepath,endpath, sep = "")

pbmc.data <- Read10X(data.dir = geneData)
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200, project = "10x_human_brain")

mito.features <- grep(pattern = "^MT-", x = rownames(x = pbmc), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = pbmc, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = pbmc, slot = 'counts'))

pbmc[['percent.mito']] <- percent.mito

#profvis::pause(10)

pdf(file=paste(filepath,"_V.pdf",sep=""))
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

pdf(file=paste(filepath,"_M.pdf",sep=""))
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mito")
dev.off()

pdf(file=paste(filepath,"_F.pdf",sep=""))
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

Rnum = Rnum + 1
print(Rnum)
}
