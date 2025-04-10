library(gtools)
library(infercnv)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(magrittr)
library(DoubletFinder)


#Set working directory
datadir <- ('C://Users/')
setwd(datadir)

#######Create seurat object, specify folder with the required files
Obj<-Read10X('C://Users/')

Obj <- CreateSeuratObject(counts = Obj, project = "Core", min.cells = 3, min.features = 200)

Obj$replicate <- Obj$orig.ident

Objidents<- Obj@meta.data$orig.ident
write.csv(Objidents, file = "ObjidentsNew.csv")

######Replace Edge with SampleID number here####

Objnewidents<- read.csv(file = "ObjidentsNew.csv")
Objnewidents<- Objnewidents[,-1]

Obj@meta.data$orig.ident<- as.factor(Objnewidents)

####add feature metadata like mitochondria and ribosomal genes

Obj[["percent.mt"]] <- PercentageFeatureSet(Obj, pattern = "^MT-")


Obj[["percent.rb"]]<- PercentageFeatureSet(Obj, pattern = "^RP[SL]")

saveRDS(Obj, file = "CoreObjraw.rds")

#####Plot features to visualize

VlnPlot(Obj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))

#### QC filtering

Obj_filt<- subset(Obj, subset = percent.mt < 15 & 10000 > nFeature_RNA & nFeature_RNA> 400 & nCount_RNA > 1000)
VlnPlot(Obj_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
Obj_filt
Obj

##### SCTransform Filtered Object

Obj_filt<- SCTransform(Obj_filt, vst.flavor = "v2", method = "glmGamPoi_offset", ncells = 2000, exclude_poisson = TRUE)
Obj_filt<- RunPCA(Obj_filt, verbose = FALSE)
ElbowPlot(Obj_filt)
Obj_filt<- RunUMAP(Obj_filt, dims = 1:20, verbose = FALSE)
Obj_filt<- FindNeighbors(Obj_filt, dims = 1:20, verbose = FALSE)
Obj_filt<- FindClusters(Obj_filt, verbose = FALSE)
DimPlot(Obj_filt, reduction = "umap", repel = TRUE, label = TRUE) + ggtitle("Obj_filtered")
saveRDS(Obj_filt, file = "Obj_Filt_SCTopt3_060423.rds")

########DoubletFinder########################
Obj<- readRDS("Obj_Filt_SCTopt3_060423.rds")
sweep.res.list_Obj <- paramSweep_v3(Obj, PCs = 1:20, sct = TRUE)
sweep.stats_Obj<- summarizeSweep(sweep.res.list_Obj, GT = FALSE)
bcmvn_Obj <- find.pK(sweep.stats_Obj)

ggplot(bcmvn_Obj, aes(pK, BCmetric, group =1))+
  geom_point() +
  geom_line()

pK<- bcmvn_Obj %>% # select the pK that corresponds to the max bcmvn to optimize the doublet detection
  filter (BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))
annotation <- Obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotation)
nExp_poi <- round(0.017*nrow(Obj@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
Obj_dbltFin<- doubletFinder_v3(Obj,
                                  PCs = 1:20,
                                  pN = 0.25,
                                  pK = pK,
                                  nExp = nExp_poi.adj,
                                  reuse.pANN = FALSE, sct = TRUE)

names(Obj_dbltFin@meta.data)

DimPlot(Obj_dbltFin, reduction = "umap", group.by = "DF.classifications_0.25_0.07_33")
table(Obj_dbltFin@meta.data$DF.classifications_0.25_0.07_33)


saveRDS(Obj_dbltFin, file = "Obj_DoubletFinder.rds")

#########################Doublet Removal#############################
Obj_dbltFin<- readRDS("Obj_DoubletFinder_070423.rds")
Obj_dbRemvd<- subset(Obj_dbltFin, subset = DF.classifications_0.25_0.07_33 == "Singlet")
Obj_dbltFin
Obj_dbRemvd

Obj_dbRemvd<- SCTransform(Obj_dbRemvd, vst.flavor = "v2", method = "glmGamPoi_offset", ncells = 2000, exclude_poisson = TRUE)
Obj_dbRemvd<- RunPCA(Obj_dbRemvd, verbose = FALSE)
ElbowPlot(Obj_dbRemvd)
Obj_dbRemvd<- RunUMAP(Obj_dbRemvd, dims = 1:20, verbose = FALSE)
Obj_dbRemvd<- FindNeighbors(Obj_dbRemvd, dims = 1:20, verbose = FALSE)
Obj_dbRemvd<- FindClusters(Obj_dbRemvd, verbose = FALSE)
DimPlot(Obj_dbRemvd, reduction = "umap", repel = TRUE, label = TRUE) + ggtitle("Obj_dbRemvd")
saveRDS(Obj_dbRemvd, file = "Obj_dbRemvd_SCTopt3.rds")
