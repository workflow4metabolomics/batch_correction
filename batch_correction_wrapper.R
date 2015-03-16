#!/usr/bin/Rscript --vanilla --slave --no-site-file

################################################################################################
# batch_correction_wrapper                                                                     #
#                                                                                              #
# Author: Marion LANDI / Jean-Francois MARTIN / Melanie Petera                                 #
# User: Galaxy                                                                                 #
# Original data: --                                                                            #
# Starting date: 22-07-2014                                                                    #
# Version 1: 22-07-2014                                                                        #
# Version 2: 08-12-2014                                                                        #
# Version 2.1: 09-01-2015 modification in Error message of sample matching                     #
#                                                                                              #
#                                                                                              #
# Input files: dataMatrix.txt ; sampleMetadata.txt ; variableMetadata.txt (for DBC)            #
# Output files: graph_output.pdf ; corrected table ; diagnostic table                          #
#                                                                                              #
################################################################################################


library(batch) #necessary for parseCommandArgs function
args = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects

source_local <- function(fname){
	argv <- commandArgs(trailingOnly = FALSE)
	base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
	source(paste(base_dir, fname, sep="/"))
}

#Import the different functions
source_local("Normalisation_QCpool.r")
source_local("RcheckLibrary.R")


## Reading of input files
idsample=read.table(args$sampleMetadata,header=T,sep='\t')
iddata=read.table(args$dataMatrix,header=T,sep='\t')

### Table match check 
#if(length(which(colnames(iddata)[-1]%in%idsample[,1]))!=(dim(iddata)[2]-1) ||
#     length(which(idsample[,1]%in%colnames(iddata)[-1]))!=dim(idsample)[1]){
#  stop("\nData matrix and sample metadata do not match regarding sample identifiers.\n",
#       "Please check your data (including check of special characters in identifiers).\n",
#	   "Note: data must not contain duplication in samples' identifiers.")
#}

### Table match check 
table.check <- match2(iddata,idsample,"sample")



### Formating
idsample[[1]]=make.names(idsample[[1]])
dimnames(iddata)[[1]]=iddata[[1]]

### Transposition of ions data
idTdata=t(iddata[,2:dim(iddata)[2]])
idTdata=data.frame(dimnames(idTdata)[[1]],idTdata)
	
### Merge of 2 files (ok even if the two dataframe are not sorted on the same key)
id=merge(idsample, idTdata, by.x=1, by.y=1)

id$batch=as.factor(id$batch)
ids=id[id$sampleType == 'pool' | id$sampleType == 'sample',]
nbid=dim(idsample)[2]
	
### Checking the number of sample and pool
	
# least 2 samples
if(length(which(ids$sampleType == "sample"))<2){
	table.check <- c(table.check,"\nError: less than 2 samples specified in Sample meta-data.",
	       "\nMake sure this is not due to errors in sampleType coding.\n")
}
	
# least 2 pools per batch for all batchs
B <- rep(0,length(levels(ids$batch)))
for(nbB in length(levels(ids$batch))){
	B[nbB]<-length(which(ids[which(ids$batch==(levels(ids$batch)[nbB])),]$sampleType == "pool"))
}
if(length(which(B>1))==0){
	table.check <- c(table.check,"\nError: less than 2 pools specified in each batch in Sample meta-data.",
	       "\nMake sure this is not due to errors in sampleType coding.\n")
}
	
### Factor of interest 
factbio=args$ref_factor


if(args$analyse == "batch_correction") {
	## Reading of Metadata Ions file
	metaion=read.table(args$variableMetadata,header=T,sep='\t')
	## Table match check 
	table.check <- c(table.check,match2(iddata,metaion,"variable"))
	check.err(table.check)
	
	## variables
	detail=args$detail
	method=args$method
	
	## outputs
	outfic=args$variable_for_simca
	outlog=args$graph_output
	
	## Launch
	res = norm_QCpool(ids,nbid,outfic,outlog,factbio,metaion,detail,F,F,method,args$span)
	save(res, file=args$rdata_output)
	write.table(res[[1]], file=args$dataMatrix_out, sep = '\t', row.names=F, quote=F)
	write.table(res[[2]], file=args$variableMetadata_out, sep = '\t', row.names=F, quote=F)
}else{
	## error check
	check.err(table.check)
	
	## outputs
	out_graph_pdf=args$out_graph_pdf
	out_preNormSummary=args$out_preNormSummary
	
	## Launch
	plotsituation(ids,nbid,out_graph_pdf,out_preNormSummary,factbio,args$span)
}

rm(args)