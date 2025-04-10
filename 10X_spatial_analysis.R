library(gtools)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(magrittr)

Csparse_validate <- "CsparseMatrix_validate"

#Set working directory
datadir <- ('740/outs/')
setwd(datadir)
#imagedir<- ('/740/outs/spatial/')
#Image740<- Read10X_Image(image.dir = imagedir, filter.matrix = TRUE)
##Load10X visium data##
G740_1<- Load10X_Spatial(data.dir = datadir, assay = "Spatial",slice = "slice740", filter.matrix = TRUE, type = "740") (https://github.com/satijalab/seurat/issues/7530)

####Load second sample data #4 #####
datadir <- ('/4/outs/')
setwd(datadir)
G4_1<- Load10X_Spatial(data.dir = datadir, assay = "Spatial",slice = "slice4", filter.matrix = TRUE, type = "4")
###change wd to main directory
datadir <- ('/Seurat_10Xvisium_BPai/')
setwd(datadir)
##Merge #4 and #740 samples
#GBM<- merge(G4_1,G740_1)
saveRDS(G4_1, file = "4_10xVis_HuGBM_unfilt.rds")
saveRDS(G740_1, file = "740_10xVis_HuGBM_unfilt.rds") 
G4<- readRDS("4_10xVis_HuGBM_unfilt.rds")
G740<- readRDS("740_10xVis_HuGBM_unfilt.rds")

####Quality Control#####
#GBM <- PercentageFeatureSet(GBM, "^MT-", col.name = "percent_mito")
#GBM <- PercentageFeatureSet(GBM, pattern = "^RP[LS]", col.name = "percent.ribo")
rm(GBM)

G4_1 <- PercentageFeatureSet(G4_1, "^MT-", col.name = "percent_mito")
G740_1<- PercentageFeatureSet(G740_1, "^MT-", col.name = "percent_mito")

###plot QC##
pdf("4_VlnPlt_QC_10xVisium_HumanGBM.pdf")
Plot<- VlnPlot(G4_1, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"), pt.size = 0.1, ncol = 3) + NoLegend()
print(Plot)
dev.off()
rm(Plot)
pdf("740_VlnPlt_QC_10xVisium_HumanGBM.pdf")
Plot<- VlnPlot(G740_1, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"), pt.size = 0.1, ncol = 3) + NoLegend()
print(Plot)
dev.off()
rm(Plot)

pdf("740_QC_featrPlt_10xVisHuGBM.pdf", width = 25, height = 20)
plot<-SpatialFeaturePlot(object = G740_1, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"))
print(plot)
dev.off()
rm(plot)

pdf("4_QC_featrPlt_10xVisHuGBM.pdf", width = 25, height = 20)
plot<-SpatialFeaturePlot(object = G4, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"))
print(plot)
dev.off()
rm(plot)

###QC filter###

G4_filt<- subset(G4, percent_mito < 15 & 12000 > nFeature_Spatial & nFeature_Spatial> 400 & nCount_Spatial > 1000)
G740_filt<- subset(G740, percent_mito < 15 & 12000 > nFeature_Spatial & nFeature_Spatial> 400 & nCount_Spatial > 1000)

###plot QC##
pdf("4filt_VlnPlt_QC_10xVisium_HumanGBM.pdf")
Plot<- VlnPlot(G4_filt, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"), pt.size = 0.1, ncol = 3) + NoLegend()
print(Plot)
dev.off()
rm(Plot)
pdf("740filt_VlnPlt_QC_10xVisium_HumanGBM.pdf")
Plot<- VlnPlot(G740_filt, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"), pt.size = 0.1, ncol = 3) + NoLegend()
print(Plot)
dev.off()
rm(Plot)
pdf("740filt_QC_featrPlt_10xVis.pdf", width = 25, height = 20)
plot<-SpatialFeaturePlot(object = G740_filt, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"))
print(plot)
dev.off()
rm(plot)
pdf("4filt_QC_featrPlt_10xVis.pdf", width = 25, height = 20)
plot<-SpatialFeaturePlot(object = G4_filt, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"))
print(plot)
dev.off()
rm(plot)


####Plot3##
G740<- readRDS("740_10xVisium.rds")

G740[["percent.mt"]] <- PercentageFeatureSet(G740,pattern = "^mt-")
pdf("QC_VlnPlot_740_10XvisHumanGBM.pdf")
Plot<-VlnPlot(
  G740, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
print(Plot)
dev.off()
rm(Plot)

####SCtransform #740 ###
G740<-G740_filt
G740 <- SCTransform(G740, assay = "Spatial", verbose = FALSE)
G740 <- RunPCA(G740, assay = "SCT", verbose = FALSE)
G740 <- FindNeighbors(G740, reduction = "pca", dims = 1:30)
G740 <- FindClusters(G740, verbose = FALSE)
G740 <- RunUMAP(G740, reduction = "pca", dims = 1:30)
saveRDS(G740, file = "740filt_10xVis_SCT_dims30.rds")
##plot UMAP##
pdf("UMAP_740filt_SCT_dims30.pdf")
p1 <- DimPlot(G740, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(G740, label = TRUE, label.size = 3)
Plot<- p1 + p2
print(Plot)
dev.off()
rm(Plot)

####SCtransform #4 ###
G4_filt <- SCTransform(G4_filt, assay = "Spatial", verbose = FALSE)
G4_filt <- RunPCA(G4_filt, assay = "SCT", verbose = FALSE)
G4_filt <- FindNeighbors(G4_filt, reduction = "pca", dims = 1:30)
G4_filt <- FindClusters(G4_filt, verbose = FALSE)
G4_filt <- RunUMAP(G4_filt, reduction = "pca", dims = 1:30)
saveRDS(G4_filt, file = "4_filt_10xVis_SCT_dims30.rds")
##plot UMAP##
pdf("UMAP_4filt_SCT_dims30.pdf")
p1 <- DimPlot(G4_filt, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(G4_filt, label = TRUE, label.size = 3)
Plot<- p1 + p2
print(Plot)
dev.off()
rm(Plot)

##DE analysis###
de_markers <- FindMarkers(G740, ident.1 = 5, ident.2 = 6)
write.csv(de_markers, file = "740filt_deMarkers_SCT.csv")
##Plot top3 markers on spatial image##
pdf("740_Top3_DEmarkers.pdf")
plot<-SpatialFeaturePlot(object = G740, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
print(plot)
dev.off()
rm(plot)

###Find spatially variable features##
G740 <- FindSpatiallyVariableFeatures(G740, assay = "SCT", features = VariableFeatures(G740)[1:1000],
                                       selection.method = "moransi")
###Visualize##
top.features <- head(SpatiallyVariableFeatures(G740, selection.method = "moransi"), 6)
pdf("740_Top3_TopVar_SCTdims30.pdf")
plot<- SpatialFeaturePlot(G740, features = top.features, ncol = 3, alpha = c(0.1, 1))
print(plot)
dev.off()
rm(plot)

######GeneExpression Visualization#
pdf("740_FeatrPlt_Select_genes.pdf", width = 25, height = 20)
plot<- SpatialFeaturePlot(G740, features = c("MDM4", "AC007402.1", "TERT", "DENND2A", "IGF2BP3" ))
print(plot)  ########## stored in Figures folder
dev.off()
rm(plot)

pdf("G4_FeatrPlt_Select_genes.pdf", width = 25, height = 20)
plot<- SpatialFeaturePlot(G4, features = c("MDM4", "AC007402.1", "TERT", "DENND2A", "IGF2BP3" ))
print(plot)  ########## stored in Figures folder
dev.off()
rm(plot)

######GeneExpression Visualization##
pdf("740_FeatrPlt_GSCgenes.pdf", width = 25, height = 20)
plot<- SpatialFeaturePlot(G740, features = c("SOX2", "EGFR", "NES"))
print(plot)
dev.off()
rm(plot)



####Labeltransfer from scRNAseq data

HumanGBM<- readRDS('obj.rds')

anchors<- FindTransferAnchors(reference = HumanGBM, query = G740, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = HumanGBM$CellType, prediction.assay = TRUE,
                                  weight.reduction = G740[["pca"]], dims = 1:30)
G740[["predictions"]] <- predictions.assay
DefaultAssay(G740)<- predictions.assay
pdf("G740_FeaturePlot_CellType_Predictions.pdf")
SpatialFeaturePlot(G740, features= c("GSC", "MG-ML", "Oligo", "Neuron", "Astrocyte"),pt.size.factor = 1.6, ncol = 2, crop = TRUE)
dev.off()

######
#####Add Modulescore invasion signatures ANASTASSIOU_MULTICANCER_INVASIVENESS_SIGNATURE to GBM_confCI snRNAseq###
dat<- read.csv("ANASTASSIOU_MULTICANCER_INVASIVENESS_SIGNATURE.txt", header = FALSE)$V1
dat<- list(dat)
DefaultAssay(GBM_confCI)<-"SCT"
GBM_confCI<- AddModuleScore(GBM_confCI, features = dat, name = "Invasion_")
pdf("Invasion_FetrPlt_GBMconfCI_dims15.pdf", width = 15, height = 15)
FeaturePlot(GBM_confCI, features = "Invasion_1")+ggtitle("Multicancer Invasion signature score")
dev.off()
###above did not work since these are scores##
###try adding Module score to the G740 directly####
DefaultAssay(G740)<- "SCT"
G740<- AddModuleScore(G740, features = dat, name = "Invasion_")
pdf("aInvasion_FetrPlt_G740_SCT.pdf")
FeaturePlot(G740, features = "Invasion_1")+ggtitle("Multicancer Invasion signature score")
dev.off()

pdf("Invasion_spatialFetrPlt_G740_SCT.pdf")
SpatialFeaturePlot(G740, features = "Invasion_1")+ggtitle("Multi-cancer Invasion signature score")
dev.off()

####EMT signatures#####
dat<- read.csv("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt", header = FALSE)$V1
dat<- list(dat)
G740<- AddModuleScore(G740, features = dat, name = "EMT_")
pdf("EMT_FetrPlt_G740_SCT.pdf")
FeaturePlot(G740, features = "EMT_1")+ggtitle("EMT signature score")
dev.off()
pdf("aEMT_spatialFetrPlt_G740_SCT.pdf")
SpatialFeaturePlot(G740, features = "EMT_1",pt.size.factor = 3 )+ggtitle("EMT signature score")+ scale_fill_viridis_c()
dev.off()

#######Adding signature scores from MSigDB#######
library(msigdbr)
GOBP.C5<- msigdbr(species = "human", category = "C5", subcategory = "GO:BP" )
GS.migration <- GOBP.C5%>% dplyr::distinct("GOBP_CELL_MIGRATION", gene_symbol)%>% as.data.frame()
write.csv(GS.migration, file = "GOBP_Cell_Migration_Geneset.csv")
dat<- read.csv("GOBP_Cell_Migration_Geneset.csv", header = FALSE)$V1
dat<- list

DefaultAssay(G740)<- "SCT"
DefaultAssay(G4)<- "SCT"
G740<- AddModuleScore(G740, features = dat, name = "GOBPCellmigration")
G4<- AddModuleScore(G4, features = dat, name = "GOBPCellmigration")

pdf("G740_GOBPmigration_spatialFetrPlt.pdf")
SpatialFeaturePlot(G740, features = "GOBPCellmigration1",pt.size.factor = 3 )+ggtitle("GOBP_Cell_Migration")+ scale_fill_viridis_c()
dev.off()

pdf("G4_GOBPmigration_spatialFetrPlt.pdf")
SpatialFeaturePlot(G4, features = "GOBPCellmigration1",pt.size.factor = 3 )+ggtitle("GOBP_Cell_Migration")+ scale_fill_viridis_c()
dev.off()