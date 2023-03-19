#########################
#CASE CONTROL COMPARISON#
#########################
#reyka@wustl.edu#########
#Last updated March 10th, 2023
###############################

##This script is doing pairwise DEG analysis between a set of controls and patients using a Seurat object

##packages of intereset
library(Seurat)
library(pheatmap)
library(tidyverse)

 controlsd28<-c("None_Patient_17_Day_28-10","None_Patient_27_Day_28-16","None_Patient_33_Day_28-24","None_Patient_36_Day_28-27","None_Patient-009-Day-28-1","None_Patient-026-Day-28-10","None_Patient_11_Day_28-7")
 patientsd28<-c("GvHD_Patient_24_Day_28-13","GvHD_Patient_29_Day_28-19","GvHD_Patient_31_Day_28-22","GvHD_Patient_38_Day_28-30","GvHD_Patient-016-Day-28-4")

  #controlsd28<-c("Control-1-Day-28-13","Control-2-Day-28-15","Control_5_Day_28-1","Control_6_Day_28-3")
  #patientsd28<-c("Patient_17_Day_28-10","Patient_11_Day_28-7","Patient_24_Day_28-13","Patient_27_Day_28-16","Patient_29_Day_28-19","Patient_31_Day_28-22","Patient_33_Day_28-24","Patient_36_Day_28-27","Patient_38_Day_28-30","Patient-009-Day-28-1","Patient-016-Day-28-4","Patient-026-Day-28-10")

###To name output file
file="D28"

seurat_obj=patient_monomacro_d28
metadataforgroup="acute_GVHD_patient"

  datalist = list()
  Idents(seurat_obj)<-seurat_obj@meta.data[[metadataforgroup]]
  for (control in controlsd28){
    for (patient in patientsd28){
      compgroup="D28" ###Might want to change this but doesnt matter. 
      group1=control
      group2=patient
      i=paste0(group1,"_",group2)
      print(group1)
      print(group2)
      #Calculate differentially expressed genes between two groups of interest
      df <- FindMarkers(object = seurat_obj,only.pos=FALSE,min.pct=0.3,ident.1=c(paste0(group1)),ident.2=paste0(group2))
      dfup = df[df$avg_log2FC  >= 0.5 & df$p_val_adj <= 0.05,] ###Change this if you want to make less stringent with respect to avg_log2FC
      dfdown = df[df$avg_log2FC  <= -0.5 & df$p_val_adj <= 0.05,] ###Change this if you want to make less stringent with respect to avg_log2FC
      upreggenes<-nrow(dfup)
      downreggenes<-nrow(dfdown)
      print(paste0(i,"Upregulated genes:",upreggenes," Downregulated genes:",downreggenes))
      if ( (upreggenes >= 1) && (downreggenes >= 1) ){
        dfup$ident1<-paste0(group1)
        dfup$ident2<-paste0(group2)
        dfup$gene<-rownames(dfup)
        dfdown$ident1<-paste0(group1)
        dfdown$ident2<-paste0(group2)
        dfdown$gene<-rownames(dfdown)
        ###If this pairing has dfup and dfdown as non empty
        dat<-rbind(dfup,dfdown)
        rownames(dat)<-NULL
        dat$i <- i 
        dat$group<-compgroup
        datalist[[i]] <- dat
      } else if ( (upreggenes == 0) && (downreggenes >= 1) ) {
        dfdown$ident1<-paste0(group1)
        dfdown$ident2<-paste0(group2)
        dfdown$gene<-rownames(dfdown)
        dat<-dfdown
        rownames(dat)<-NULL
        dat$i <- i 
        dat$group<-compgroup
        datalist[[i]] <- dat
      } else if ( (upreggenes >= 1) && (downreggenes == 0) ){
        dfup$ident1<-paste0(group1)
        dfup$ident2<-paste0(group2)
        dfup$gene<-rownames(dfup)
        dat<-dfup
        rownames(dat)<-NULL
        dat$i <- i 
        dat$group<-compgroup
        datalist[[i]] <- dat
      } else{
        print("SHIIIIT")
      }
    }
  }

  #Combine all data into one big thing
  big_data_d28 = do.call(rbind, datalist)
  big_data_d28$comparison<-paste0(big_data_d28$ident1,"_",big_data_d28$ident2)
  #        p_val avg_log2FC pct.1 pct.2 p_val_adj        ident1   ident2 i
  #H1FX        0   1.541079 0.990 0.963         0 Plerixafor(H) G-CSF(H) 1
  metadata <- unique(data.frame(big_data_d28$comparison,big_data_d28$group,big_data_d28$ident1,big_data_d28$ident2))
  colnames(metadata) <- c("comparison","group","ident1","ident2") 
  row.names(metadata)<-metadata$comparison

  big_datamatrix_d28=big_data_d28 %>% select(comparison,gene,avg_log2FC) %>%
    pivot_wider(names_from = comparison, values_from = avg_log2FC, values_fill = 0) %>% 
    as.data.frame() %>%
    column_to_rownames("gene")
  #set color breaks
  paletteLength=50
  myColor <- colorRampPalette(c("blue","white","darkred"))(paletteLength)

  myBreaksd28 <- c(seq(min(big_data_d28$avg_log2FC), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(big_data_d28$avg_log2FC)/paletteLength, max(big_data_d28$avg_log2FC), length.out=floor(paletteLength/2)))

  #Change annotation colors based on your data - or no need to annotate if you dont care
  #annotation_colors = list(
  #  group = c("Day28"="#e76f51","Day60"="#2a9d8f"),
  #  ident2 = c("Control-1-Day-28-13"="#f47068","Control-1-Day-60-14"="#ffb3ae","Control-2-Day-28-15"="#0e606b","Control-2-Day-60-16"="#1697a6"),
  #  ident1=c("Patient-009-Day-28-1"="#f4f1de","Patient-009-Day-60-2"="#eab69f","Patient-016-Day-28-4"="#e07a5f","Patient-016-Day-60-5"="#3d405b","Patient-022-Day-60-8"="#81b29a","Patient-026-Day-28-10"="#f2cc8f")
  #)

  #subsetted heatmaps
  pdf(paste("heatmap_DEG_",file,".pdf", sep=""), width=15, height=35)
      pheatmap(big_datamatrix_d28, display_numbers = F,color = myColor,breaks= myBreaksd28,cluster_rows = F, cluster_cols = F, annotation_colors = annotation_colors)
  dev.off()