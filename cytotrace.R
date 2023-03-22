###Cytotrace Usage
#Last updated 03/22/2023
#reyka@wustl.edu
##

##Install
#install.packages("devtools")
library(devtools)
#download cytotrace from here: https://cytotrace.stanford.edu/
devtools::install_local("/path/to/directory/CytoTRACE_0.3.3.tar.gz")

#Load libraries of interest
library(Seurat)
library(CytoTRACE)
library(tidyr)

#Load in seurat object
obj<-readRDS("myobject.rds")

#Extract matrix for running in cytotrace
exp.rawdata <- as.matrix(obj$RNA@counts)
write.table(exp.rawdata, file=paste0("rawdata.txt"), sep="\t", quote = FALSE, row.names = TRUE) #Save for posterity

#Run Cytotrace
raw_counts<-read.table(file=paste0("rawdata.txt"),sep="\t")
results <- CytoTRACE(raw_counts, ncores = 40)

#Extract barcode and cytotrace score to add to seurat object
cytoout<-as.data.frame(results$CytoTRACE)
colnames(cytoout)<-c("cytotrace")
cytoout$rownames<-rownames(cytoout)
#cytotrace replaces barcode names -1 with .1 - need to convert back to add to seurat object
cytoout$rownames<-gsub("\\.","-",as.character(cytoout$rownames))

rownames(cytoout)<-cytoout$rownames
cytoout$rownames<-NULL
write.table(cytoout, file=paste0("cytotraceoutput.txt"), sep="\t", quote = FALSE, row.names = TRUE)

##Add to metadata of seurat object
cytotrace<-read.table(file=paste0("cytotraceoutput.txt"), sep="\t", header = TRUE)
obj <- AddMetaData(
  object = obj,
  metadata = cytotrace)