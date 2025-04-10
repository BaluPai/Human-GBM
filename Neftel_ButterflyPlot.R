library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)
library(ggpubr)

obj <- readRDS('SeuratObjectDirectory.rds')

#Read table of metamodule genes defined by Neftel et al. (2019)
gene_sets_df <- read.csv('IDHwt.GBM.MetaModules.tsv', header = T, sep = "\t")

#Combine lists of NPClike1 and NPClike2 together and MESlike1 and MESlike2 together to plot as four quadrants
gene_sets_df <- data.frame('AClike' = c(gene_sets_df$AClike), 'MESlike' = c(gene_sets_df$MESlike1, gene_sets_df$MESlike2), 'NPClike' = c(gene_sets_df$NPClike1, gene_sets_df$NPClike2), 'OPClike' = c(gene_sets_df$OPClike))

#Convert into a list (required by the AddModuleScore function)
gene_setsL <- as.list(gene_sets_df)

#Removed NAs
gene_setsL$AClike <- gene_setsL$AClike[!is.na(gene_setsL$AClike)]

#Annotate list components
names(gene_setsL) <- colnames(gene_sets_df)

#Run module score analysis
obj <- AddModuleScore(obj, features = gene_setsL, name = names(gene_setsL))

#Calculate which module each cell belongs to based on score from meta.data columns AClike_1, MESlike_2, NPClike_3, OPClike_4
colNames <- colnames(obj@meta.data[, (length((colnames(obj@meta.data))) - (length(names(gene_setsL))-1)):length((colnames(obj@meta.data)))])
obj$Neftel.ident <- colNames[max.col(obj@meta.data[, (length((colnames(obj@meta.data))) - (length(names(gene_setsL))-1)):length((colnames(obj@meta.data)))], ties.method = "first")]

#Change OPClike_4, NPClike_3, AClike_1, and MESlike_2 to whatever the meta.data column name is (in my case, AddModuleScore adds "_" and a number corresponding to the list components, e.g. AClike is the first on the list so the meta.data column is calld "AClike_1" rather than "AClike")
obj$y_axis <- sapply(1:nrow(obj@meta.data), function(x) max(c(obj$OPClike_4[x], obj$NPClike_3[x])) - max(c(obj$AClike_1[x], obj$MESlike_2[x])))
obj$x_axis <- 0
obj$x_axis[obj$y_axis > 0] <- (obj$OPClike_4 - obj$NPClike_3)[obj$y_axis > 0] 
obj$x_axis[obj$y_axis < 0] <- (obj$AClike_1 - obj$MESlike_2)[obj$y_axis < 0]

#Generate butterfly plot
p1 <- ggplot(obj@meta.data, aes(x=x_axis, y=y_axis, color=Neftel.ident)) + 
  geom_point(alpha=0.5, shape=19, size=.8) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(xintercept = 0,linetype = 'dashed') +
  geom_hline(yintercept = 0,linetype = 'dashed') +
  annotate("text", x = -1, y  = 1.5, label = "NPC-like") + 
  annotate("text", x = 1, y = 1.5, label = "OPC-like") +
  annotate("text", x = -1, y = -1.5, label = "MES-like") + 
  annotate("text", x = 1, y = -1.5, label = "AC-like") +
  scale_color_manual(breaks = c("AClike_1", "NPClike_3", "MESlike_2", "OPClike_4"), values = c("#BD9E60", "#6B9FC7", "#9F4F5B", "#5AA272")) 

pdf("NeftelPlot.pdf")
p1
dev.off()
