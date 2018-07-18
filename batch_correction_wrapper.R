#!/usr/bin/env Rscript

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
# Version 2.2: 16-03-2015 inclusion of miniTools' functions for special characters             #
# Version 2.90: 18-08-2015 new parameter valnull                                               #
# Version 2.91: 25-08-2016 error message improvment                                            #
#                                                                                              #
#                                                                                              #
# Input files: dataMatrix.txt ; sampleMetadata.txt ; variableMetadata.txt (for DBC)            #
# Output files: graph_output.pdf ; corrected table ; diagnostic table                          #
#                                                                                              #
################################################################################################


library(batch) #necessary for parseCommandArgs function

##------------------------------
## test help option
##------------------------------

# Prog. constants
argv.help <- commandArgs(trailingOnly = FALSE)
script.path <- sub("--file=", "", argv.help[grep("--file=", argv.help)])
prog.name <- basename(script.path)

# Test Help
if (length(grep('-h', argv.help)) > 0) {
  cat("Usage: Rscript ", 
    prog.name,
    "{args} \n",
    "parameters: \n",
    "\tanalyse {val}: must be set to \"batch_correction\"",
    "\tdataMatrix {file}: set the input data matrix file (mandatory) \n",
    "\tsampleMetadata {file}: set the input sample metadata file (mandatory) \n",
    "\tvariableMetadata {file}: set the input variable metadata file (mandatory) \n",
    "\tmethod {opt}: set the method; can set to \"linear\", \"lowess\" or \"loess\" (mandatory) \n",
    "\tspan {condition}: set the span condition; set to \"none\" if method is set to \"linear\" (mandatory) \n", 
    "\tref_factor {value}: set the ref_factor value; (if span value is set to NULL, optional) \n",
    "\tdetail {value}: set the detail value; (if span value is set to NULL, optional) \n",
    "\tdataMatrix_out {file}: set the output data matrix file (mandatory) \n",
    "\tvariableMetadata_out {file}: set the output variable metadata file (mandatory) \n",
    "\tgraph_output {file}: set the output graph file (mandatory) \n",
    "\trdata_output {file}: set the output Rdata file (mandatory) \n",
    "\tbatch_col_name {val}: the column name for batch. Default value is \"batch\".\n",
    "\tinjection_order_col_name {val}: the column name for the injection order. Default value is \"injectionOrder\".\n",
    "\tsample_type_col_name {val}: the column name for the sample types. Default value is \"sampleType\".\n",
    "\tsample_type_tags {val}: the tags used inside the sample type column, defined as key/value pairs separated by commas (example: blank=blank,pool=pool,sample=sample).\n",
    "\n")
  quit(status = 0)
}

##------------------------------
## init. params
##------------------------------

args = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects

# Set default col names
if ( ! 'batch_col_name' %in% names(args))
	args[['batch_col_name']] <- 'batch'
if ( ! 'injection_order_col_name' %in% names(args))
	args[['injection_order_col_name']] <- 'injectionOrder'
if ( ! 'sample_type_col_name' %in% names(args))
	args[['sample_type_col_name']] <- 'sampleType'
if ( ! 'sample_type_tags' %in% names(args))
	args[['sample_type_tags']] <- 'blank=blank,pool=pool,sample=sample'

# Parse sample type tags
sample.type.tags <- list()
for (kv in strsplit(strsplit(args$sample_type_tags, ',')[[1]], '='))
	sample.type.tags[[kv[[1]]]] <- kv[[2]]
if ( ! all(c('pool', 'blank', 'sample') %in% names(sample.type.tags)))
	stop("All tags pool, blank and sample must be defined in option sampleTypeTags.")
args$sample_type_tags <- sample.type.tags

##------------------------------
## init. functions
##------------------------------

source_local <- function(...){
	argv <- commandArgs(trailingOnly = FALSE)
	base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
	for(i in 1:length(list(...))){source(paste(base_dir, list(...)[[i]], sep="/"))}
}
#Import the different functions
source_local("Normalisation_QCpool.r","easyrlibrary-lib/RcheckLibrary.R","easyrlibrary-lib/miniTools.R")


## Reading of input files
idsample=read.table(args$sampleMetadata,header=T,sep='\t',check.names=FALSE,comment.char = '')
iddata=read.table(args$dataMatrix,header=T,sep='\t',check.names=FALSE,comment.char = '')

### Table match check 
table.check <- match2(iddata,idsample,"sample")
if(length(table.check)>1){check.err(table.check)}

### StockID
samp.id <- stockID(iddata,idsample,"sample")
iddata<-samp.id$dataMatrix ; idsample<-samp.id$Metadata ; samp.id<-samp.id$id.match

### Checking mandatory variables
mand.check <- ""
for(mandcol in c(args$sample_type_col_name, args$injection_order_col_name, args$batch_col_name)){
  if(!(mandcol%in%colnames(idsample))){
    mand.check <- c(mand.check,"\nError: no '",mandcol,"' column in sample metadata.\n",
                    "Note: table must include this exact column name (it is case-sensitive).\n")
  }
}
if(length(mand.check)>1){
  mand.check <- c(mand.check,"\nFor more information, see the help section or:",
                  "\n http://workflow4metabolomics.org/sites/",
                  "workflow4metabolomics.org/files/files/w4e-2016-data_processing.pdf\n")
  check.err(mand.check)
}

### Formating
idsample[[1]]=make.names(idsample[[1]])
dimnames(iddata)[[1]]=iddata[[1]]

### Transposition of ions data
idTdata=t(iddata[,2:dim(iddata)[2]])
idTdata=data.frame(dimnames(idTdata)[[1]],idTdata)
	
### Merge of 2 files (ok even if the two dataframe are not sorted on the same key)
id=merge(idsample, idTdata, by.x=1, by.y=1)

id[[args$batch_col_name]]=as.factor(id[[args$batch_col_name]])
ids=id[id[[args$sample_type_col_name]] == args$sample_type_tags$pool | id[[args$sample_type_col_name]] == args$sample_type_tags$sample,]
nbid=dim(idsample)[2]
	
### Checking the number of sample and pool
	
# least 2 samples
if(length(which(ids[[args$sample_type_col_name]] == args$sample_type_tags$sample))<2){
	table.check <- c(table.check,"\nError: less than 2 samples specified in sample metadata.",
	       "\nMake sure this is not due to errors in sampleType coding.\n")
}
	
# least 2 pools per batch for all batchs
B <- rep(0,length(levels(ids[[args$batch_col_name]])))
for(nbB in length(levels(ids[[args$batch_col_name]]))){
	B[nbB]<-length(which(ids[which(ids[[args$batch_col_name]]==(levels(ids[[args$batch_col_name]])[nbB])),][[args$sample_type_col_name]] == args$sample_type_tags$pool))
}
if(length(which(B>1))==0){
	table.check <- c(table.check,"\nError: less than 2 pools specified in each batch in sample metadata.",
	       "\nMake sure this is not due to errors in sampleType coding.\n")
}
	
### Factor of interest 
factbio=args$ref_factor


if(args$analyse == "batch_correction") {
	## Reading of Metadata Ions file
	metaion=read.table(args$variableMetadata,header=T,sep='\t',check.names=FALSE,comment.char = '')
	## Table match check 
	table.check <- c(table.check,match2(iddata,metaion,"variable"))
	check.err(table.check)
	
	## variables
	detail=args$detail
	method=args$method
	
	## outputs
	outlog=args$graph_output
	
	## Launch
	res = norm_QCpool(ids,nbid,outlog,factbio,metaion,detail,F,F,method,args$span,args$valnull)
	save(res, file=args$rdata_output)
	write.table(reproduceID(res[[1]],res[[3]],"sample",samp.id)$dataMatrix, file=args$dataMatrix_out, sep = '\t', row.names=F, quote=F)
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
