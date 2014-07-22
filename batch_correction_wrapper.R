#!/usr/bin/Rscript --vanilla --slave --no-site-file

################################################################################################
# batch_correction_wrapper                                                                     #
#                                                                                              #
# Author : Marion LANDI                                                                        #
# User : Galaxy                                                                                #
# Original data : --                                                                           #
# Starting date : 22-07-2014                                                                   #
# Version 1 : 22-07-2014                                                                       #
# Version 2 :                                                                        #
#                                                                                              #
#                                                                                              #
# Input files : dataMatrix.txt ; sampleMetadata.txt                                            #
# Output files : graph_output.pdf                                                              #
#                                                                                              #
################################################################################################


library(batch) #necessary for parseCommandArgs function
args = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects

source_local <- function(fname){
	argv <- commandArgs(trailingOnly = FALSE)
	base_dir <- dirname(substring(argv[grep(--file=, argv)], 8))
	source(paste(base_dir, fname, sep=/))
}

#Import the different functions
source_local(Normalisation_QCpool.r)


## Reading of Metadata Samples file
idsample=read.table(args$sampleMetadata,header=T,sep='\t')
idsample[[1]]=make.names(idsample[[1]])
	
### Reading the data file XCMS size (Ions in lines)
iddata=read.table(args$dataMatrix,header=T,sep='\t')
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
#No perfect match
wrng <- ""
if((nrow(id)!=nrow(idsample)) || (nrow(id)!=nrow(idTdata))) {
	wrng <- paste("Warning: Sample meta-data table and Data matrix are not a perfect match.",
	              "\nMake sure this is not due to errors in sample identifiers.\n\n\n")
}
	
# least 2 samples
if(length(which(ids$sampleType == "sample"))<2){
	stop(c("\n\nError: less than 2 samples specified in Sample meta-data.",
	       "\nMake sure this is not due to errors in sampleType coding.\n\n",wrng))
}
	
# least 2 pools per batch for all batchs
B <- rep(0,length(levels(ids$batch)))
for(nbB in length(levels(ids$batch))){
	B[nbB]<-length(which(ids[which(ids$batch==(levels(ids$batch)[nbB])),]$sampleType == "pool"))
}
if(length(which(B>1))==0){
	stop(c("\n\nError: less than 2 pools specified in each batch in Sample meta-data.",
	       "\nMake sure this is not due to errors in sampleType coding.\n\n",wrng))
}
	
### Factor of interest 
factbio=args$ref_factor

if(args$analyse=='batch_correction') {
	## Reading of Metadata Ions file
	metaion=read.table(args$variableMetadata,header=T,sep='\t')
	
	## variables
	detail=args$detail
	method=args$method
	
	## outputs
	outfic=args$variable_for_simca
	outlog=args$graph_output
	
	## Launch
	res = norm_QCpool(ids,nbid,outfic,outlog,factbio,metaion,detail,F,F,method)`);
	save(res, file=args$rdata_output)
	write.table(res[[1]], file=args$dataMatrix_out, sep = '\t', row.names=F, quote=F)
	write.table(res[[2]], file=args$variableMetadata_out, sep = '\t', row.names=F, quote=F)
}else{
	## outputs
	out_graph_pdf=args$out_graph_pdf
	out_preNormSummary=args$out_preNormSummary
	
	## Launch
	plotsituation(ids,nbid,out_graph_pdf,out_preNormSummary,factbio)
}

rm(args)