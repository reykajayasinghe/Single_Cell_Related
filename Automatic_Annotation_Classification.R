###Single-cell type annotation 
###reyka@wustl.edu
###Last modified: March 11, 2023
#This script helps with annotation of single-cell data using three tools:
#[1] SingleR: https://github.com/dviraran/SingleR
#[2] sc-type: https://github.com/IanevskiAleksandr/sc-type/
#[3] scCATCH: https://cran.r-project.org/web/packages/scCATCH/vignettes/tutorial.html

####Install before starting####
#prior to installling these github packages you may need to run the following to 
#set up github path
#create_github_token() #Take link from this command and paste in browser and follow details to generate login
#Sys.setenv(GITHUB_PAT = 'TOKEN HERE') #copy token thats created to this path

devtools::install_github('dviraran/SingleR')
install.packages(pkgs = 'devtools')
devtools::install_github('ZJUFanLab/scCATCH')
remotes::install_github("LTLA/celldex")

#load packages
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)

myobject=readRDS("object.rds") #Load in your seurat_object of interest

##########################################################
###Add Single R Annotation for Cell Type Classification###
##########################################################

library(SingleR)

#Convert seurat object to single-cell format for singleR 
#https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html
Seurat_Object_Diet <- DietSeurat(myobject, graphs = "umap") #https://github.com/satijalab/seurat/issues/4633
obj.sce <- as.SingleCellExperiment(Seurat_Object_Diet)

#Reference SingleR
#https://bioconductor.org/packages/3.13/data/experiment/vignettes/celldex/inst/doc/userguide.html

ref <- MonacoImmuneData()
#singleR usage
#http://bioconductor.org/books/devel/SingleRBook/using-the-classic-mode.html#annotating-the-test-dataset
pred <- SingleR(method = "single",sc_data =  obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
ref_pred<-as.data.frame(pred)
#Extract columns of interest to add to seurat object
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "MonacoImmuneData"
myobject <- AddMetaData(
  object = myobject,
  metadata = ref_predictions)

ref <- NovershternHematopoieticData()
pred <- SingleR(method = "single",sc_data = obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
ref_pred<-as.data.frame(pred)
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "NovershternHematopoieticData"
myobject <- AddMetaData(
  object = myobject,
  metadata = ref_predictions)

ref <- DatabaseImmuneCellExpressionData()
pred <- SingleR(method = "single",sc_data = obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
ref_pred<-as.data.frame(pred)
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "DatabaseImmuneCellExpressionData"
myobject <- AddMetaData(
  object = myobject,
  metadata = ref_predictions)

ref <- HumanPrimaryCellAtlasData()
pred <- SingleR(method = "single",sc_data = obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
ref_pred<-as.data.frame(pred)
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "HumanPrimaryCellAtlasData"
myobject <- AddMetaData(
  object = myobject,
  metadata = ref_predictions)

ref <- BlueprintEncodeData()
pred <- SingleR(method = "single",sc_data = obj.sce, ref = ref, 
    labels = ref$label.fine, assay.type.test=1)
ref_pred<-as.data.frame(pred)
ref_predictions<- ref_pred %>%
  select(pruned.labels)
colnames(ref_predictions) <- "BlueprintEncodeData"
myobject <- AddMetaData(
  object = myobject,
  metadata = ref_predictions)

saveRDS(myobject,"object_celltypeannotation_singleR.rds")

##########################################################
###Add sc-type Annotation for Cell Type Classification###
##########################################################
#https://github.com/IanevskiAleksandr/sc-type/
# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)

# get cell-type-specific gene sets from our in-built database (DB)
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = myobject[["SCT"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(myobject@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(myobject@meta.data[myobject@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(myobject@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

myobject@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  myobject@meta.data$customclassif[myobject@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

saveRDS(myobject,"object_celltypeannotation_singleR_sctype.rds")


##########################################################
###Add scCATCH Annotation for Cell Type Classification###
##########################################################

##scCATCH cell type annotation
#https://cran.r-project.org/web/packages/scCATCH/vignettes/tutorial.html
####### THIS ONE TAKES THE LONGEST TIME - if you filter the gene list by tissue type like below it will run a lot faster
library(scCATCH)

#create scratch object
obj <- createscCATCH(data = myobject[['RNA']]@data, cluster = as.character(Idents(myobject)))
#demo_marker()
#Filter marker list of interest
# filter cellmatch
###Modify this tissue list based on your tissue type of interest
#I dont reccomend doing the cancer type filter - gets rid of a lot of celltypes.
#https://github.com/ZJUFanLab/scCATCH/wiki/human_tissues
cellmatch_new <- cellmatch[cellmatch$species == "Human" & cellmatch$tissue %in% c("Blood", "Peripheral blood", "Plasma", "Serum", "Umbilical cord blood","Bone marrow"), ]
# you can change the content of tissue as "interested"
cellmatch_new$tissue <- "interested"

# use custom cellmatch by setting tissue as "interested"
obj <- findmarkergene(object = obj,
                      species = "Human",
                      marker = cellmatch_new,
                      tissue = "interested")

# Evidence-based score and annotation for each cluster with findcelltype()
obj <- findcelltype(object = obj)
#Grab cell type annotation and add to Seurat Object
hea#Grab cell type annotation and add to Seurat Object
sccatch_anno<-obj@celltype %>% as.data.frame()

write.table(sccatch_anno,"sccatch_celltype_annotation.tsv",row.names=FALSE,quote=FALSE)

###Extract columns of interest for annotation
clusteridents<-sccatch_anno %>% select(cluster,cell_type)
#   cluster                       cell_type
#1        7                   Memory B Cell
#2        4             Natural Killer Cell
#3       12 Monocyte Derived Dendritic Cell
#4        5     Natural Killer T (NKT) Cell

myobject@meta.data$sccatch_anno <- clusteridents$cell_type[match(myobject@meta.data$seurat_clusters, clusteridents$cluster)]

saveRDS(myobject,"object_celltypeannotation_all.rds")

