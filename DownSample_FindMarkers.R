library(gtools)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(magrittr)
Csparse_validate <- "CsparseMatrix_validate"

setwd("wd/")
#read object
GBM_Integrated_SCT<- readRDS("Obj.rds")
#get cell barcodes of each cell type

Core1<- colnames(Obj)[which(Obj$patientCellType.CI == 'GSC_Core_GBM1')]
Core2<- colnames(Obj)[which(Obj$patientCellType.CI == 'GSC_Core_GBM2')]
Core4<- colnames(Obj)[which(Obj$patientCellType.CI == 'GSC_Core_GBM4')]
Core6<- colnames(Obj)[which(Obj$patientCellType.CI == 'GSC_Core_GBM6')]

Edge1<- colnames(Obj)[which(Obj$patientCellType.CI == 'GSC_Edge_GBM1')]
Edge2<- colnames(Obj)[which(Obj$patientCellType.CI == 'GSC_Edge_GBM2')]
Edge4<- colnames(Obj)[which(Obj$patientCellType.CI == 'GSC_Edge_GBM4')]
Edge6<- colnames(Obj)[which(Obj$patientCellType.CI == 'GSC_Edge_GBM6')]


set.seed(1)

#Switch from integrated assay to RNA
DefaultAssay(Obj) <- "SCT"

#Loop "num.loop" times of differential analysis for sub-sampled obj
num.loop <- 100 # modify as required

#Create an empty object to store the FindMarkers results
pv.mtx1 <- NULL
pv.mtx2 <- NULL
pv.mtx3 <- NULL

setwd("wd/Run1/")
for (i in 1:num.loop){  
  
  #randomly sample 25 cell barcodes from GBM
  downCore<- c((sample(Core1,1576)),(sample(Core2,1888)), (sample(Core4,345)),(sample(Core6,778)))
  downEdge<- c((sample(Edge1,1576)),(sample(Edge2,1888)), (sample(Edge4,345)),(sample(Edge6,778)))
  
  #Do findMarkers with no return threshold in order to obtain the p-value from all genes
  GSC_CorevsEdge<- FindMarkers(Obj, ident.1 = downCore, ident.2 = downEdge,logfc.threshold = 0,
                               min.pct =0, latent.vars = 'patient',test.use = 'LR',verbose = FALSE)
  #mytime <- format(Sys.time(), "%b_%d_%H_%M_%S_%Y")
  myfile <- file.path(getwd(), paste0("scGSCCorevsEdge_", i, ".csv"))
  write.csv(GSC_CorevsEdge, file = paste(myfile, sep = "" ))
  
  #save the p-values
  idx.order <- match(row.names(Obj), row.names(GSC_CorevsEdge))
  pv.mtx1 <- cbind(pv.mtx1, GSC_CorevsEdge$p_val_adj[idx.order])
  pv.mtx2 <- cbind(pv.mtx2, GSC_CorevsEdge$avg_log2FC[idx.order])
  pv.mtx3 <- cbind(pv.mtx3, GSC_CorevsEdge$p_val[idx.order])

}

# 
row.names(pv.mtx1) <- row.names(Obj)
row.names(pv.mtx2) <- row.names(Obj)
row.names(pv.mtx3) <- row.names(Obj)
colnames(pv.mtx1) <- paste0("Run_",1:num.loop)
colnames(pv.mtx2) <- paste0("Run_",1:num.loop)
colnames(pv.mtx3) <- paste0("Run_",1:num.loop)
#c(1:2)) #colnames(pv.mtx[1]))

#save all the padj-values into a csv file
pv.mtx<- merge(pv.mtx1, pv.mtx2, by = 'row.names', all = TRUE)
write.csv(pv.mtx, "1_scDEGSCcorevsEdge.csv")

#calculate average padj-value of each gene and number of time pv < 0.05
avg.pvadj <- apply(pv.mtx1, 1, mean)
count <- apply(pv.mtx1, 1, function(x){length(which(x < 0.05))})

avg.fc2 <- apply(pv.mtx2, 1, mean)
avg.pv3 <- apply(pv.mtx3, 1, mean)


#create a table
summary.table <- as.data.frame(cbind(as.data.frame(cbind(as.data.frame(cbind(avg.pvadj, count)), as.data.frame(avg.fc2))), as.data.frame(avg.pv3)))

row.names(summary.table) <- names(avg.pvadj)
colnames(summary.table) <- c("avg_pvadj", "count_sig_pv", "avg_logFC", "avg_pv")

#sort the table by average padj-value
summary.table <- summary.table[order(summary.table$avg_pvadj,decreasing=FALSE),]


write.csv(summary.table, "1_scDEA_SummaryTable_GSCcorevsEdge_pADJ.csv")


