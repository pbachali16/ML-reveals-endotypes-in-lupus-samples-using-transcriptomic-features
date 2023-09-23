####### Dated Feb 24th 2023 ############
####### Author: Prathyusha Bachali, Email - prathyusha.bachali@ampelbiosolutions.com###
#######Function to create Expression Set Objects(Eset) for Illumina platform#######

##########Illumina Eset############
create_illumina_eset <- function(path,num,gpl,metadata){
  file_list <- list.files(path)
  nrow(read.delim(file_list[1],header = FALSE,stringsAsFactors = FALSE))
  dataset <- data.frame(matrix("", ncol = 1, nrow = 48803))
  for (file in file_list){
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
      dataset <- read.table(file, header = FALSE, sep="\t",stringsAsFactors = FALSE)
    }
    # if the merged dataset does exist, append to it
    if (exists("dataset")){
      temp_dataset <-read.table(file, header = FALSE, sep="\t",stringsAsFactors = FALSE)
      dataset<-cbind(dataset, temp_dataset)
      rm(temp_dataset)
    }
  }
  dataset <- dataset[,-1] #Get rid of the first empty column
  rownames(dataset) <- dataset$V1 #Make the first probe column the rownames
  c <- 1:ncol(dataset) #Make an object counting the number of columns
  c%%3 #Give columns numeric designations by groups of three
  dataset <- dataset[, !(c%%3==1)] #Drop all of the columns with probe names
  file_list <- substr(file_list,1,nchar(file_list)-10) #Remove file appendix from file list names
  file_list <- paste(file_list,"Sample",sep=".") #Append ".Sample" to the file list names
  index <- (2*(1:num))-1 #Make an index of numbers representing every other column beginning with column 1
  samples <- dataset[,index] #Collect just the non-normalized sample probe intensity values
  colnames(samples) <- file_list #Make the column names of the non-normalized probe values that of file list sample names
  file_list <- gsub("Sample","Detection",file_list) #Append "Detection" to the file list names
  index <- 2*(1:num) #Make an index of numbers representing every other column beginning with column 2
  detection <- dataset[,index] #Collect just the non-normalized sample detection p-values
  colnames(detection) <- file_list #Make the column names of the detection p-values that of file list sample names
  temp <- cbind(samples,detection) #Join the sample probe values and their detection p-values
  n <- ncol(samples) #Get the number of columns
  dataset <- temp[,kronecker(1:n, c(0, n), "+") ] #Insert the detection p-value columns behind each of their corresponding probe intensity value columns
  dataset <- cbind(ID_Ref = rownames(dataset),dataset) #Add a separate column of the probe IDs onto the dataset
  exprs.raw <- dataset #Rename dataset as exprs.raw
  rm(c,n,file,file_list,temp,samples,detection,index,dataset)
  
  fData <- gpl
  fData <- fData[1:48834,]
  index <- which(fData$ID%in%rownames(exprs.raw))
  fData <- fData[index,]
  rownames(fData) <- fData$ID
  fData <- fData[order(rownames(fData)), ]
  exprs.raw <- exprs.raw[ order(rownames(exprs.raw)), ]
  x <- identical(rownames(fData),rownames(exprs.raw))
  # 
  # # Append annotation columns to the exprs.raw frame and write out a new file for read.illmn
  exprs.annot <- cbind(exprs.raw, fData)
  write.table(exprs.annot, file="annotated.txt", sep="\t", quote=FALSE, row.names=FALSE)
  # 
  annotation.cols <- colnames(fData)
  elist <- read.ilmn(files="annotated.txt",expr="Sample",other.columns="Detection",
                               probeid="ID_Ref", annotation=annotation.cols, verbose=TRUE)
  eset.neqc <- neqc(elist)
  
  eset.neqc.exprs <- eset.neqc$E
  rownames(eset.neqc.exprs) <- rownames(exprs.raw)
  
  fData <- fData[ order(fData$ID), ]
  eset.neqc.exprs <- eset.neqc.exprs[ order(rownames(eset.neqc.exprs)), ]
  y <- identical(rownames(fData),rownames(eset.neqc.exprs))
  
  pData <- metadata
  rownames(pData) <- pData$geo_accession
  #index <- 2*(1:6)
  #colnames(eset.neqc.exprs) <- colnames(exprs.raw[,index]) #Skip if using the family based data file
  eset.neqc.exprs <- eset.neqc.exprs[ , order(colnames(eset.neqc.exprs))]
  rownames(pData) <- pData$geo_accession
  colnames(eset.neqc.exprs) <- rownames(pData)
  x <- identical( rownames(pData),colnames(eset.neqc.exprs) )
  if(x & y ==TRUE){
    exprs.eset <- as.matrix(eset.neqc.exprs)
    rownames(exprs.eset) <- rownames(eset.neqc.exprs)
    fData.meta <- data.frame(labelDescription = colnames(fData), row.names = colnames(fData) )
    fData.eset <- new("AnnotatedDataFrame", data = data.frame(fData), varMetadata = fData.meta)
    
    pData.meta <- data.frame(labelDescription = colnames(pData), row.names = colnames(pData) )
    pData.eset <- new("AnnotatedDataFrame", data = data.frame(pData), varMetadata = pData.meta)
    
    expData <- new("MIAME", name = "Prat",
                   lab = "AMPEL",
                   contact = "AMPEL BX Lead",title = "Illumina",
                   abstract = "blah,blah,blah",
                   url = "ampel.com",
                   other = list(notes = "stuff and things"));
    
    eset <- ExpressionSet(assayData=exprs.eset,
                          phenoData = pData.eset,
                          featureData = fData.eset,
                          experimentData = expData,
                          annotation = "Illumina Eset compiled by AMPEL Biosolutions")
    
  }
  
  else{
    print("The order of rows/columns in Eset and featuredata/metadata do not match")
  }
  
  return(eset)
  
}





