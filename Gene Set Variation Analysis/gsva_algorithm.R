####### Dated Feb 24th 2023 ############
####### Author: Prathyusha Bachali, Email - prathyusha.bachali@ampelbiosolutions.com###
######### Function to run GSVA on the log2 expression values or the Expression Set Objects (Eset)##############
##### The input object for the IQR function is Expression Set Object (Eset) or the matrix of the log2 expression values ####
#### If the matrix of log2 expression values are provided as input, load the gene annotations and the patient data or metadata separately ####
#### The IQR function is applied to filter low expressed genes prior to run GSVA. The results will be different if one opt out the IQR filteration technqiue ###



library(limma)
library(DESeq2)
library(GSEABase)
library(GSVA)
library(dplyr)





do_iqr <- function(eset){
  
  IQR_df <- data.frame()
  index <- numeric()
  exprs_df <- data.frame()
  annot_df <- data.frame()
  exprs_df_filt <- data.frame()
  exprs_df_genes <- data.frame()
  exprs_df_entrez <- data.frame()
  
  for(i in 1:nrow(exprs(eset))) {
    row <- as.vector(exprs(eset)[i,])
    #print(row)
    IQR_df[i,1] <- IQR(row)
    colnames(IQR_df) <- "IQR"
  }
  annot_df <- fData(eset)
  exprs_df <- as.data.frame(exprs(eset))
  if("probe"%in%names(annot_df)){
    IQR_df$probe <- rownames(exprs(eset))
    exprs_df$probe <- rownames(exprs_df)
    exprs_df <- merge(exprs_df,annot_df,by="probe")
    IQR_df <- IQR_df[(IQR_df$IQR>0.00),] #20,251(Could change the 0.05 to 0)
    exprs_df_filt <- merge(exprs_df,IQR_df,by="probe")
  }
  else{
    IQR_df$ensembl <- rownames(exprs(eset))
    exprs_df$ensembl <- rownames(exprs_df)
    exprs_df <- merge(exprs_df,annot_df,by="ensembl")
    IQR_df <- IQR_df[(IQR_df$IQR>0.00),] #20,251(Could change the 0.05 to 0)
    exprs_df_filt <- merge(exprs_df,IQR_df,by="ensembl")
  }
  exprs_df_filt = exprs_df_filt[ order(-exprs_df_filt$IQR), ]
  exprs_df_genes <- exprs_df_filt[ !duplicated(exprs_df_filt$geneSymbol), ]
  exprs_df_genes <- exprs_df_genes[!is.na(exprs_df_genes$geneSymbol), ]
  rownames(exprs_df_genes) <- exprs_df_genes$geneSymbol
  exprs_df_entrez <- exprs_df_filt[ !duplicated(exprs_df_filt$geneEntrezID), ]
  exprs_df_entrez <- exprs_df_entrez[!is.na(exprs_df_entrez$geneEntrezID), ]
  rownames(exprs_df_entrez) <- exprs_df_entrez$geneEntrezID
  return(list("IQR_BA"=IQR_df,"exprs_df_filtered"=exprs_df_filt,
              "exprs_filtered_genes"=exprs_df_genes,
              "exprs_filtered_entrez"=exprs_df_entrez))
  
}


####### Input parameters - Expression Set Object (Eset), Predefined geneset modules - list object, var1 - geneSymbol, var2 - Category or column stored the module ####### 

run_gsva <- function(eset,geneSig,var1,var2){
  results <- do_iqr(eset)
  eset_genes <- results$exprs_filtered_genes
  gsva_input <- eset_genes[,-which(names(eset_genes) %in% c("ensembl","geneSymbol","geneEntrezID","IQR","probe","geneName","geneEnsembl"))]
  cols_of_interest <- c(var1,var2)
  geneSig <- geneSig[,which(colnames(geneSig)%in%cols_of_interest)]
  testList <- tapply(geneSig[[var1]],geneSig[[var2]],c)
  
  n <- names(testList)
  uniqueList <- lapply(testList, unique)
  
  makeSet <- function(geneIds, n) {
    GeneSet(geneIds, geneIdType=SymbolIdentifier(), setName=n)
  }
  
  gsList <- gsc <- mapply(makeSet, uniqueList[], n)
  newlist = gsList %>% lapply(function(x) trimws(x@geneIds))
  # 
  gsva_output <- gsva(as.matrix(gsva_input), newlist, method = c("gsva"),
                      min.sz=1, max.sz=Inf, verbose=TRUE)
  
  overlap_genes <- merge(eset_genes,geneSig,by=var1)
  overlap_genes <- overlap_genes[,-which(names(overlap_genes) %in% c("ensembl","geneEntrezID","IQR","probe","geneName","geneEnsembl"))]
  
  #
  scores_df <- as.data.frame(gsva_output)
  scores_df$category <- rownames(scores_df)
  mods <- rownames(scores_df)
  scores_t <- as.data.frame(t(as.data.frame(gsva_output)))
  #scores_t$friendlyName <- rownames(scores_t)
  #return(eset_genes)
  
  return(list("scores" = scores_df,
              "overlap_genes" = overlap_genes,
              "scores_t" = scores_t))
  
}

