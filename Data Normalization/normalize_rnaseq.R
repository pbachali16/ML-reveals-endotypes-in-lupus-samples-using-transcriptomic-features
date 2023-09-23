##### The script to normalize the RNAseq counts using DeSeq method #######


library(dplyr)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(genefilter)
library(vsn)
library(biomaRt)
library(Biobase)
library(limma)
library(HTSFilter)



normalize_rnaseq <- function(filepath,metadata){
  setwd(filepath)
  txt_files_ls = list.files(pattern="*.txt")
  txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x, header = T, sep ="\t")})
  # Combine them
  counts_df <- do.call("cbind", lapply(txt_files_df, as.data.frame))
  rownames(counts_df) <- counts_df$Geneid
  c <- 1:ncol(counts_df)
  c%%7 #Give columns numeric designations by groups of three
  counts_df <- counts_df[, (c%%7==0)] #Drop all of the columns with probe names
  colnames(counts_df) <- gsub("X.mnt.disks.ssd1.working.RD_","",colnames(counts_df))
  colnames(counts_df) <- gsub("_R1_001.fastq.gz_STARAligned.out.sorted.bam","",colnames(counts_df))
  metadata <- metadata
  e <- identical(rownames(metadata),colnames(counts_df))
  dds <- DESeqDataSetFromMatrix(countData = counts_df,
                                colData = metadata,
                                design = ~cohort) ### The design matrix is modified based on the experimental design of the dataset ####
  ####perform a minimal pre-filtering to keep only rows that have at least 10 reads total####
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  gene.list <- as.data.frame(rownames(counts_df))
  colnames(gene.list) <- "ensembl"
  gene.list <- substr(gene.list$ensembl,1,15)
  gene.list <- as.data.frame(gene.list)
  colnames(gene.list) <- c("ensembl")
  ensembl <- useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl",host = "http://www.ensembl.org")
  #listMarts(host="www.ensembl.org")
  genemap <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id","hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values = gene.list$ensembl,
                   mart = ensembl,useCache = FALSE)
  idx <- match(gene.list$ensembl, genemap$ensembl_gene_id)
  gene.list$geneEntrezID <- genemap$entrezgene[idx]
  gene.list$geneSymbol <- genemap$hgnc_symbol[idx]
  rownames(dds) <- substr(rownames(dds),1,15)
  annotation <- gene.list[match(rownames(dds), gene.list$ensembl),]
  all(rownames(dds) == annotation$ensembl)
  mcols(dds) <- cbind(mcols(dds),annotation)
  
  dds$cohort <- as.factor(dds$cohort)
  idx <- as.character(dds$cohort) #Make a vector of WHO variables for HTS filter
  pdf("HTS-Filter_Plot.pdf", width = 16, height = 9)
  filtered <- HTSFilter(dds,
                        conds = idx,
                        pAdjustMethod = "BH",
                        s.len = 100,
                        s.min = 1,
                        s.max = 200,
                        normalization = "DESeq",
                        plot = T)
  dev.off()
  filtered <- filtered$filteredData #Get the data frame of the filtered transcript reads
  dim(filtered) #13049 remaining features
  filtered@rowRanges@elementMetadata@listData$ensembl <- droplevels(filtered@rowRanges@elementMetadata@listData$ensembl)
  
  pdf(file="ALL_HTS_Dispersion_Estimates.pdf",height = 9,width = 16)
  plotDispEsts( dds, ylim = c(1e-8, 1e1) ) #Check dispersion estimates 
  dev.off()
  
  dds_norm <- DESeq(filtered)
  
  ##### Log2 transformation ######
  log2 <- normTransform(dds_norm)
  log2.counts <- assay(log2)
  log2.counts <- as.data.frame(log2.counts)
  gene.ids <- annotation[match(rownames(log2.counts), annotation$ensembl),]
  log2.counts$geneSymbol <- gene.ids$geneSymbol
  log2.counts$geneEntrezID <- gene.ids$geneEntrezID
  
  ### VST normalized data #####
  vsd <- varianceStabilizingTransformation(dds_norm, blind=FALSE)
  vsd.counts <- assay(vsd)
  vsd.counts <- as.data.frame(vsd.counts)
  gene.ids <- annotation[match(rownames(vsd.counts), annotation$ensembl),]
  vsd.counts$geneSymbol <- gene.ids$geneSymbol
  vsd.counts$geneEntrezID <- gene.ids$geneEntrezID
  
  pData <- metadata
  cols_to_remove <- c("ensembl","geneSymbol","geneEntrezID")
  exprs <- log.counts[,-which(colnames(log2.counts)%in%cols_to_remove)]
  fData <- as.data.frame(log2.counts[,c("ensembl","geneSymbol","geneEntrezID")])
  fData$ensembl <- rownames(fData)
  exprs <- exprs[order(rownames(exprs)), ]
  fData <- fData[(order(fData$ensembl)),]
  e2 <- identical(rownames(fData),rownames(exprs))
  if (e2==TRUE){
    exprs.eset <- as.matrix(exprs)
    rownames(exprs.eset) <- rownames(exprs)
    
    pData.meta <- data.frame(labelDescription = colnames(pData), row.names = colnames(pData) )
    pData.eset <- new("AnnotatedDataFrame", data = data.frame(pData), varMetadata = pData.meta)
    
    fData.meta <- data.frame(labelDescription = colnames(fData), row.names = colnames(fData) )
    fData.eset <- new("AnnotatedDataFrame", data = data.frame(fData), varMetadata = fData.meta)
    
    expData <- new("MIAME", name = "PB",
                   lab = "ResearchDx TO4", 
                   contact = "PB",title = "m-rna",
                   abstract = "blah,blah,blah",
                   url = "ampel.com",
                   other = list(notes = "stuff and things"));
    
    eset <- ExpressionSet(assayData=exprs.eset,
                          phenoData = pData.eset,
                          featureData = fData.eset,
                          experimentData = expData,
                          annotation = "Lupus_RNA-seq")
    
    
    save(eset, file="ESET_DESeq-Log2-counts_WithHTS-Filter.RData")
  }
  print ("The order doesn't match")
}