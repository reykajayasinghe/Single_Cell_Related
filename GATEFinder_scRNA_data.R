###Run GATEfinder on sc-RNAseq Data
##Modified on April 2, 2023
##reyka@wustl.edu


#Need to determine target population
#Collect cell.ids for cell group of interest
Idents(sobj)<-sobj@meta.data$seurat_clusters
cellsofinterest<-WhichCells(sobj, idents = c("0","2","3"))
#Extract data for specific genes of interest
PM_prediction_DEG<-c("A2M","ABCB1","ACTB","ACTN1","ADAM10","ADAM19","ADAM23","ADD3","ADRB2","AHNAK","AKAP11","ALCAM","ANK3","ANXA1","ANXA5","APBB1IP","APMAP","AQP3","AREG","ARHGAP19","ARHGAP33","ARHGAP9","ARHGEF1","ARHGEF6","ARL4C","ARPC5","ASPM","ATP5PB","B2M","B4GALT1","BAIAP2L1","BCL2","BIN2","BIRC6","BORCS5","BTN3A2","C1orf43","C1QBP","C4orf3","CALR","CANX","CAT","CCR4","CCR6","CCR7","CCR9","CCT4","CD2","CD247","CD27","CD28","CD300A","CD38","CD3D","CD3E","CD3G","CD4","CD40LG","CD44","CD48","CD5","CD52","CD55","CD59","CD63","CD69","CD7","CD74","CD79A","CD81","CD83","CD84","CD8A","CD8B","CD8B2","CD96","CD99","CDC14A","CDC20","CDC25A","CDC42EP3","CDC42SE1","CDC42SE2","CDC45","CDCA2","CDCA3","CDCA4","CDCA5","CDCA7","CDCA7L","CDCA8","CDK1","CDK17","CDK2","CDK6","CDKN1A","CDKN1B","CDKN2C","CDKN2D","CDKN3","CDR2","CDV3","CEP55","CFL1","CFLAR","CKAP5","CLDN1","CLDND1","CLEC2B","CLEC2D","CLIC1","COLQ","COX4I1","COX6C","COX7C","COX8A","CR1","CRLF3","CRTAM","CSGALNACT1","CTLA4","CXCR4","CYB5B","CYB5R3","CYTH1","DDX58","DIAPH3","DMTF1","DNAJA4","DNAJC1","DOCK8","ECE1","EEF1A1","EEF1D","EEF2","EMB","EMP3","ERBIN","ERN1","ERP29","EVA1C","EVI2B","EZR","FAAH2","FCGRT","FCMR","FGR","FHIT","FLNA","FURIN","FXYD5","FYB1","FYN","GAPDH","GAS7","GLG1","GNAI2","GNG2","GPR155","GPR183","GYPA","GZMH","HCST","HERPUD1","HNRNPK","HSP90AA1","HSP90AB1","HSP90B1","HSPA5","HSPA8","HSPD1","ICOS","IDH2","IFFO2","IFNGR1","IL10RA","IL17RA","IL18R1","IL18RAP","IL2RA","IL2RB","IL2RG","IL4R","IL6ST","IL7R","IQGAP1","IQGAP2","ITGAL","ITGB1","ITGB2","ITK","ITM2A","ITM2B","ITM2C","ITPRIP","KLRB1","KLRD1","KLRF1","KLRG1","LAG3","LAT","LCP1","LDLR","LDLRAD4","LGALS1","LIPA","LMBR1","LMNA","LNPEP","LRIG1","LRRC59","LRRN3","LTB","LY9","LYPD3","LZTFL1","MAN1A2","MAN1C1","MBP","MCL1","MCUR1","ME1","MGAT4A","MPHOSPH9","MPZL3","MSN","MYADM","MYH9","NCR3","NDUFA4","NDUFB3","NEU1","NHSL2","NKTR","NT5E","ODC1","OSBPL8","OST4","PAM","PCBP2","PDE4B","PDE4D","PDIA4","PDZD8","PEBP1","PECAM1","PFN1","PHB2","PHF20L1","PIGX","PIK3IP1","PKM","PLCB1","PLP2","PLXDC1","PLXNA4","PPIA","PPP1R16B","PPP1R18","PRF1","PRKCA","PRNP","PTGER2","PTGER4","PTK2B","PTPRC","PTPRE","PTPRJ","RAC2","RCCD1","RELT","RHOA","ROBO1","RPL18","RPL19","RPL21","RPL27A","RPL28","RPL3","RPL32","RPL35","RPL35A","RPL37A","RPL9","RPLP2","RPN2","RPS15A","RPS23","RPS3","RPSA","RSL1D1","S100A4","S100A6","SARAF","SELENOK","SELL","SELPLG","SERINC5","SERP1","SIRPG","SKAP1","SLAMF1","SLAMF6","SLAMF7","SLC1A4","SLC25A3","SLC29A1","SLC2A3","SLC37A4","SLC38A2","SLC4A10","SLCO3A1","SLFN12L","SLMAP","SMDT1","SMIM24","SNRPD2","SOD1","SORL1","SPN","SSR2","ST6GAL1","STK17B","STOM","STX16","SYNE1","SYNE2","TAGLN2","TGFB1","TGFBR2","TGFBR3","TGOLN2","TMEM106C","TMEM123","TMPO","TMX4","TNF","TNFRSF1B","TOMM20","TOMM7","TOR1AIP1","TPI1","TRABD2A","TRAC","TRAF3IP3","TRAT1","TRBC1","TRBC2","TRGC1","TRGC2","TRIM59","TSPAN14","TUBB","TUBB4B","TXNIP","TYROBP","UBC","UBE2C","UBE2D3","UQCRFS1","UTRN","VAMP2","VIM","VSIR","XBP1","YBX1","YWHAQ","YWHAZ","ZC3HAV1")
SCT_assay<-FetchData(object = sobj, vars = PM_prediction_DEG)
#True false list of target population of interest to feed into gatefinder tool
targetpop<-rownames(SCT_assay) %in% cellsofinterest
results_targeted_1=GateFinder(SCT_assay, targetpop)
#Need to determine target population
#Collect cell.ids for cell group of interest
cellsofinterest_2<-WhichCells(sobj, idents = c("1","4","6","5","7"))
#True false list of target population of interest to feed into gatefinder tool
targetpop_2<-rownames(SCT_assay) %in% cellsofinterest_2
results_targeted_2=GateFinder(SCT_assay, targetpop_2)

#Plot results using below command
pdf(paste("Gatefinder_testrun_threegroups.pdf", sep=""))
  plot(results)
  #plot(SCT_assay, results, c(2,3), targetpop)&ggtitle("Pops:2,6")
  DotPlot(tissue_harmony_version1,features=c("CD3D","KLRB1","CDC42SE1","CXCR4"),cols=c("paleturquoise1","magenta4"))&ggtitle("Pops:2,6")
  FeaturePlot(tissue_harmony_version1,features=c("CD3D","KLRB1","CDC42SE1","CXCR4"),order=TRUE,cols=c("paleturquoise1","magenta4"))&coord_equal()&ggtitle("Pops:2,6")
  FeatureScatter(object = tissue_harmony_version1, feature1 = 'CD3D', feature2 = 'KLRB1',slot="scale.data")&ggtitle("Pops:2,6")
  FeatureScatter(object = tissue_harmony_version1, feature1 = 'CDC42SE1', feature2 = 'CXCR4',slot="scale.data")&ggtitle("Pops:2,6")

  plot(results_2)
  DotPlot(tissue_harmony_version1,features=c("CD2","FYB1","CXCR4","EZR"),cols=c("paleturquoise1","magenta4"))&ggtitle("Pops:0,3,5,8")
  FeaturePlot(tissue_harmony_version1,features=c("CD2","FYB1","CXCR4","EZR"),order=TRUE,cols=c("paleturquoise1","magenta4"))&coord_equal()&ggtitle("Pops:0,3,5,8")
  FeatureScatter(object = tissue_harmony_version1, feature1 = 'CD2', feature2 = 'FYB1',slot="scale.data")&ggtitle("Pops:0,3,5,8")
  FeatureScatter(object = tissue_harmony_version1, feature1 = 'CXCR4', feature2 = 'EZR',slot="scale.data")&ggtitle("Pops:0,3,5,8")

  plot(results_3)
  DotPlot(tissue_harmony_version1,features=c("EZR","NHSL2","CXCR4","KLRB1"),cols=c("paleturquoise1","magenta4"))&ggtitle("Pops:1,4,7")
  FeaturePlot(tissue_harmony_version1,features=c("EZR","NHSL2","CXCR4","KLRB1"),order=TRUE,cols=c("paleturquoise1","magenta4"))&coord_equal()&ggtitle("Pops:1,4,7")
    FeatureScatter(object = tissue_harmony_version1, feature1 = 'EZR', feature2 = 'NHSL2',slot="scale.data")&ggtitle("Pops:1,4,7")
  FeatureScatter(object = tissue_harmony_version1, feature1 = 'CXCR4', feature2 = 'KLRB1',slot="scale.data")&ggtitle("Pops:1,4,7")

dev.off()