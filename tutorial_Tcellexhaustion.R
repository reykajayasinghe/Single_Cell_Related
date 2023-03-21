##scRNA-seq tutorial for evaluating Exhausted T cell Markers in scRNA-seq Data 

#Load in required packages for running script - these may need to be installed
library(Seurat)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dittoSeq)

#Load in rds object (name might need to be changed)
seurat_obj<-readRDS("myeloma_object.rds")

#Gene List of Exhausted T Cell Markers
exhausted_T_cell <- c("VSIR", "TIGIT", "ICOS", "EOMES", "HAVCR2", "PDCD1", "BTLA", "CD244", "LAG3", "CD160", "CTLA4", "CD96")

#Visualize Cell Type Distribution on UMAP
#celltype label in "group.by" make need to be changed based on metadata for your object
DimPlot(seurat_obj,group.by=c("celltype"),pt.size=1,label=TRUE,raster.dpi = c(1000, 1000))&coord_equal()&theme(legend.position = 'right')
#Visualize Marker Expression on UMAP
FeaturePlot(seurat_obj,features=exhausted_T_cell,group.by=c("celltype"),cols=c("paleturquoise1", "magenta4"),order=TRUE)&coord_equal()
#Visualize Marker Expression using DotPlot
DotPlot(seurat_obj,features=exhausted_T_cell,group.by=c("celltype"),cols=c("paleturquoise1", "magenta4"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5),legend.position = "right")

#Subset only the T cell population from your seurat object
Idents(seurat_obj)<-seurat_obj@meta.data$celltype
T_cell_only_object<-subset(seurat_obj,idents=c("T_Cells"))

#Calculate module score for exhausted T cell features
T_cell_only_object<-AddModuleScore(T_cell_only_object,features=list(c("VSIR", "TIGIT", "ICOS", "EOMES", "HAVCR2", "PDCD1", "BTLA", "CD244", "LAG3", "CD160", "CTLA4", "CD96")),name="Exhausted_Score")

#Visualize exhausted score on UMAP and using dotplot by progression status
#group.by may need to be changed based on your metadata
FeaturePlot(T_cell_only_object,features="Exhausted_Score1",group.by=c("progression_status"),cols=c("paleturquoise1", "magenta4"),order=TRUE)&coord_equal()
DotPlot(T_cell_only_object,features="Exhausted_Score1",group.by=c("progression_status"),cols=c("paleturquoise1", "magenta4"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5),legend.position = "right")
#visualize exhausted score using violin plot
multi_dittoPlot(T_cell_only_object, vars = c("Exhausted_Score1"), 
                group.by = "progression_status", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Expression", 
                theme = theme_classic())