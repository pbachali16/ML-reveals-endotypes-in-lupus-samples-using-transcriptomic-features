####### Dated Feb 24th 2023 ############
####### Author: Prathyusha Bachali, Email - prathyusha.bachali@ampelbiosolutions.com###
#######Function to create Expr3ssion object to pre-process  and  normalize the raw microarray data#######

library(limma)
library(GEOquery)
library(affy)
library(simpleaffy) # see http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/
library(affycoretools) # contains plotPCA


##### Provide the GSE number in order to down the raw data, pheno data from NCBI GEO ########
gse <- getGEO("", GSEMatrix = TRUE)
show(gse)
###### Download the pheno or patient data ########
metadata <- gse$_series_matrix.txt.gz
metadata <- pData(metadata)
######Access the rawdata from GEO####
filePaths = getGEOSuppFiles("")
filePaths
####Decompress the CEL files#####
untar("GSE13917_RAW.tar")
celfiles <- list.celfiles()
for(i in 1:length(celfiles))
{
  gunzip(celfiles[i])
}
celfilenames <- as.data.frame(list.celfiles())
colnames(celfilenames)[1] <- "CELFILES"
order(celfilenames$CELFILES)
rownames(metadata) <- celfilenames$CELFILES

pheno <- metadata
#### Downlaod the right annotation packages ######
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezg.db)

cdf="HGU133Plus2_Hs_ENTREZG" # this from the BrainArray website where the CDFs were obtained. Look at the "Custom CDF Name"
# column of http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/19.0.0/entrezg.asp

abatch_brainarray<- ReadAffy(verbose=TRUE, cdfname=cdf, phenoData=pheno)

# BRAINARRAY ENTREZ GCRMA
eset_brainarray_gcrma_ampel <- gcrma(abatch_brainarray, cdfname=cdf)
eset_brainarray_rma_ampel <- rma(abatch_brainarray, cdfname=cdf)

###################################################################
# Annotate BRAINARRAY esets using BrainArray probe annotations DB
columns(hgu133plus2hsentrezg.db) # list available gene annotations
# prepare new fData frame for BrainArray sets
fData = data.frame(probe=rownames(eset_brainarray_gcrma_ampel))
rownames(fData) = fData$probe
keys=rownames(fData(eset_brainarray_gcrma_ampel))
geneSymbol <- data.frame(geneSymbol=convertIDs( keys, "PROBEID", "SYMBOL", hgu133plus2hsentrezg.db ))
geneName <- data.frame(geneName=convertIDs( keys, "PROBEID", "GENENAME", hgu133plus2hsentrezg.db ))
geneEntrezID <- data.frame(geneEntrezID=convertIDs( keys, "PROBEID", "ENTREZID", hgu133plus2hsentrezg.db ))
geneEnsembl <- data.frame(geneEnsembl=convertIDs( keys, "PROBEID", "ENSEMBL", hgu133plus2hsentrezg.db ))
fDataComplete = cbind( fData, geneSymbol, geneName, geneEntrezID, geneEnsembl  )

fData(eset_brainarray_gcrma_ampel) = fDataComplete
pData(eset_brainarray_gcrma_ampel) = pheno


save(eset_brainarray_gcrma_ampel,file="eset_entrez_gcrma.RData")



