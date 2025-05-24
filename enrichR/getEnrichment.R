getEnrichment <- function(df,tag){
  ###########
  dirEnrichDirection <- paste0(dirEnrich,tag,"/")
  
  if (!dir.exists(dirEnrichDirection)){
    dir.create(dirEnrichDirection)
  }else{
    print(paste("The directory",dirEnrichDirection,"already exists"))
  }
  ##############
  
  DisGeNET <- df$DisGeNET
  BP <- df$GO_Biological_Process_2021
  MF <- df$GO_Molecular_Function_2021
  
  KEGG <- df$KEGG_2021_Human
  
  TF <- df$TRANSFAC_and_JASPAR_PWMs
  
  getEnrichmentPlot(DisGeNET,"DisGeNET",top_term,thr_pval,dirEnrichDirection)
  getEnrichmentPlot(BP,"GO_BP",top_term,thr_pval,dirEnrichDirection)
  getEnrichmentPlot(MF,"GO_MF",top_term,thr_pval,dirEnrichDirection)
  getEnrichmentPlot(KEGG,"KEGG",top_term,thr_pval,dirEnrichDirection)
 # getEnrichmentPlot(TF,"TF",top_term,thr_pval,dirEnrichDirection)
  
}

