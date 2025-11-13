score_subtypes <- function(expression_mat, norm = TRUE, method = "vst", plot = TRUE){
  
  require('GSVA')
  require('edgeR')
  require('DESeq2')
  
  ## Collison geneset (https://pubmed.ncbi.nlm.nih.gov/21460848/)
  collison_gs <- list(
    Collison_Exocrine_Like = c('REG1B','REG3A','REG1A','PNLIPRP2','CEL','PNLIP','PLA2G1B','CELA3A','CPB1','CELA3B','CTRB2','CLPS','CELA2B','PRSS2','PRSS1','GP2','SLC3A1','CFTR','SLC4A4','SPINK1'),
    Collison_Basal = c('AIM2','FAM26F','GPM6B','S100A2','KRT14','CAV1','LOX','SLC2A3','TWIST1','PAPPA','NT5E','CKS2','HMMR','SLC5A3','PMAIP1','PHLDA1','SLC16A1','FERMT1','HK2','AHNAK2'),
    Collison_Classical = c('TMEM45B','SDR16C5','GPRC5A','AGR2','S100P','FXYD3','ST6GALNAC1','CEACAM5','CEACAM6','TFF1','TFF3','CAPN8','FOXQ1','ELF3','ERBB3','TSPAN8','TOX3','LGALS4','PLS1','GPX2','ATP10B','MUC13'))
  
  ## Moffitt geneset (https://pubmed.ncbi.nlm.nih.gov/26343385/)
  moffit_gs <- list(
    Moffit_Basal = c("VGLL1","UCA1","S100A2","LY6D","SPRR3","SPRR1B","LEMD1","KRT15","CTSL2","DHRS9","AREG","CST6","SERPINB3",
                     "SERPINB4","KRT6C","KRT6A","FAM83A","SCEL","FGFBP1","KRT7","KRT17","GPR87","TNS4","SLC2A1","ANXA8L2"),
    Moffit_Classical = c("BTNL8","FAM3D","ATAD4","AGR3","CTSE","LOC400573","LYZ","TFF2","TFF1","ANXA10","LGALS4","PLA2G10",
                         "CEACAM6","VSIG2","TSPAN8","ST6GALNAC1","AGR2","TFF3","CYP3A7","MYO1A","CLRN3","KRT20","CDH17","SPINK4","REG4")
  )
  
  ## Baily geneset (https://www.nature.com/articles/nature16965)
  baily_gs <- list()
  baily_gs[['Baily_Basal']] <- c('PGAM1','ANGPTL4','ARPC2','ITGA3','PTPN1','GSDMC','DRAP1','ANXA8','LDHA','ARHGAP23','HCAR2','ULBP2','KIAA1609','PTHLH','CEBPB','SLC16A3','EGFR','PLXNA1','ZBED2','IL2A','TPD5212','DUSP14','FSCN1','RND3','GPR87','PANX1','S100A2','ANXA1','KRT6A','ADAM17','PTRF','EMLIN1','FBXL7','CERCAM','DDR2','GXYLT2','CHSY3','VSTM4','COL8A1','PRRX1','COL5A2','SPARC','FBN1','FAP','COL6A3','COL1A2','BNC2','EFEMP2','FSTL1','COL3A1','GSG2','KIF18A','CDCA5','CCNA2','MCM10','DEPDC1','UBE2C','KIF4A','CCNB1','KIF2C','SPAG5','TPX2','HIST1H2B0','BUB1','CKAP2L','PTTG1','CDC20','ERCC6L','NCAPG','NCAPH','PSMA7','CCT6A','CCT7','DKC1','HSPE1','EIF2S2','MRPL11','MRPL17','LYAR','PSMA1','EIF4A3','PPAT','HAT1','PAICS','BRIX1','RUVBL1','CKS1B','TOMM40','CCT5','ALYREF')
  baily_gs[['Baily_Classical']] <- c('LRRC66','GJB1','LGALS4','PLA2G10','PLS1','BCAS1','CAPN5','HNF4G','ARHGEF38','CAPN8','AP001187.9','ABP1','ERBB3','C9orf152','SLC44A4','IHH','HNF4A','FAM83E','PRR15L','EPS8L3')
  baily_gs[['Baily_ADEX']] <- c('REG3G','SYCN','REG1P','SERPINI2','CPB1','RP11-331F4.4','AQP8','RBPJL','CUZD1','PLA2G1B','CLPS','CEL','CELA3B','CELA3A','CTRB1','CTRC','CTRB2','PNLIPRP1','CPA1','CPA2','UNC79','GCK','GJD2','SYT4','KCNK16','NOL4','SCGN','INS','CABP7','CHGB','BEX1','SVOP','MR7-3HG','ABCC8','HMGCLL1','SLC30A8','SST','CELF3','PCSK2','SCG3')
  # baily_gs[['Baily_Immunogenic']] <-  c('IGHV3-15','IGHV3-7','IGHV1-46','IGKV2-28','IGHV3-23','IGHV3-53','IGHV5-51','IGHA1','IGKV4-1','IGLV2-23','IGLV1-40','IGKV3-11','IGLL5','IGJ','IGLC3','IGKVLD-39','IGLV2-14','AC096579.7','IGKV3-20','IGKC','CSF1R','TYROBP','FCER1G','TLR8','C3AR1','CD86','WAS','SPI1','HAVCR2','SASH3','CYBB','MS4A4A','BTK','LAPTM5','PTPRC','MARCH1','CD53','AIF1','DOCK2','NCKAP1L','SIT1','GIMAP4','GPR174','SEPT1','CXCR6','ITK','TBC1D10C','GIMAP5','CD96','ZNF831','GZMK','PTPRCAP','SLAMF1','P2RY10','CD48','THEMIS','CD3E','CD3D','SH2D1A','CD2')
  
  subtypes_geneSets <- list("Collison_Exocrine_Like" = collison_gs$Collison_Exocrine_Like, "Baily_ADEX" = baily_gs$Baily_ADEX,
                            "Collison_Classical" = collison_gs$Collison_Classical, "Baily_Classical" = baily_gs$Baily_Classical,"Moffit_Classical" = moffit_gs$Moffit_Classical,
                            "Collison_Basal" = collison_gs$Collison_Basal, "Baily_Basal" = baily_gs$Baily_Basal, "Moffit_Basal" = moffit_gs$Moffit_Basal)
 
  if(norm){
    
    if(method == "vst"){
    
    message("normalizing using VST")
    expression_mat <- expression_mat + 1
    expression_mat <- vst(as.matrix(expression_mat))
    
    } else if (method == "logcpm") {
    
    message("normalizing using logCPM")
    expression_mat <- cpm(as.matrix(expression_mat), log = T)
    
    } else if (method == 'rlog') {
      
    message("normalizing using rlog")
    expression_mat <- expression_mat + 1
    expression_mat <- rlog(as.matrix(expression_mat))
    
    }
    
    kcdf = "Gaussian"
    
  } else {
    
    message("Skipping normalization")
    kcdf = "Poisson"
  }
  
  pdac_subtype_gvsa <- gsva(expr = expression_mat,
                            gset.idx.list = subtypes_geneSets,
                            kcdf = kcdf)
  
  if(plot){
    
  sign_type <- data.frame(c("NA","NA","Classical","Classical","Classical","Basal","Basal","Basal"), row.names = names(subtypes_geneSets)) %>%
    `colnames<-`("Class")
  ann_colors <- list(Class = c("NA" = 'black', "Classical" = "blue", "Basal" = "red"),
                     cell_type = c("Acinar" = 'yellow', "ADM" = "orange", "PanIN" = "darkred", "Duct" = "darkgreen", "Glandular_Tumor" = "darkblue", "PoorlyDiffTumor" = "purple"))
  
  pheatmap(t(pdac_subtype_gvsa),
           annotation_col = sign_type, annotation_colors = ann_colors,
           gaps_col = c(2,5), scale = "none", fontsize_row = 6,
           cluster_cols = F, cluster_rows = T)
  }
  
  return(pdac_subtype_gvsa)
  

}

