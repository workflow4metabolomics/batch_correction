###############################################################################################################
# batch_correction_3Llauncher                                                                                 #
#                                                                                                             #
# Authors: Jean-Francois MARTIN / Melanie Petera                                                              #
# Starting date: 04-08-2020                                                                                   #
# Based on batch_correction_wrapper.R version 2.91                                                            #
# Version 1: xx-xx-2020                                                                                       #
#            - split of tool-linked code and script-linked one                                                #
#                                                                                                             #
# Input files: dataMatrix.txt, sampleMetadata.txt, variableMetadata.txt (BC only)                             #
# Output files: graph.pdf, corrected table (BC only), diagnostic table (DBC only), variableMetadata (BC only) #
#                                                                                                             #
###############################################################################################################

meth3L <- function(idsample,iddata,sample_type_col_name,injection_order_col_name,batch_col_name,sample_type_tags,
                   factbio,analyse,metaion,detail,method,outlog,span,valnull,
                   rdata_output,dataMatrix_out,variableMetadata_out,out_graph_pdf,out_preNormSummary){

## Reading of input files
idsample=read.table(idsample,header=TRUE,sep='\t',check.names=FALSE,comment.char = '')
iddata=read.table(iddata,header=TRUE,sep='\t',check.names=FALSE,comment.char = '')

### Table match check 
table.check <- match2(iddata,idsample,"sample")
if(length(table.check)>1){check.err(table.check)}

### StockID
samp.id <- stockID(iddata,idsample,"sample")
iddata<-samp.id$dataMatrix ; idsample<-samp.id$Metadata ; samp.id<-samp.id$id.match

### Checking mandatory variables
mand.check <- ""
for(mandcol in c(sample_type_col_name, injection_order_col_name, batch_col_name)){
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

### Formating ########## to update: check if make.names needed, change colnames for mandatory var
idsample[[1]]=make.names(idsample[[1]])
dimnames(iddata)[[1]]=iddata[[1]]

### Transposition of ions data
idTdata=t(iddata[,2:dim(iddata)[2]])
idTdata=data.frame(dimnames(idTdata)[[1]],idTdata)
	
### Merge of 2 files (ok even if the two dataframe are not sorted on the same key)
id=merge(idsample, idTdata, by.x=1, by.y=1)

id[[batch_col_name]]=as.factor(id[[batch_col_name]])
ids=id[id[[sample_type_col_name]] == sample_type_tags$pool | id[[sample_type_col_name]] == sample_type_tags$sample,]
nbid=dim(idsample)[2]
	
### Checking the number of sample and pool
	
# least 2 samples
if(length(which(ids[[sample_type_col_name]] == sample_type_tags$sample))<2){
	table.check <- c(table.check,"\nError: less than 2 samples specified in sample metadata.",
	       "\nMake sure this is not due to errors in sampleType coding.\n")
}
	
# least 2 pools per batch for all batchs
B <- rep(0,length(levels(ids[[batch_col_name]])))
for(nbB in length(levels(ids[[batch_col_name]]))){
	B[nbB]<-length(which(ids[which(ids[[batch_col_name]]==(levels(ids[[batch_col_name]])[nbB])),][[sample_type_col_name]] == sample_type_tags$pool))
}
if(length(which(B>1))==0){
	table.check <- c(table.check,"\nError: less than 2 pools specified in each batch in sample metadata.",
	       "\nMake sure this is not due to errors in sampleType coding.\n")
}

### BC/DBC-specific processing

if(analyse == "batch_correction") {
	## Reading of Metadata Ions file
	metaion=read.table(metaion,header=T,sep='\t',check.names=FALSE,comment.char = '')
	## Table match check 
	table.check <- c(table.check,match2(iddata,metaion,"variable"))
	check.err(table.check)
	
	## Launch
	res = norm_QCpool(ids,nbid,outlog,factbio,metaion,detail,F,F,method,span,valnull)
	save(res, file=rdata_output)
	write.table(reproduceID(res[[1]],res[[3]],"sample",samp.id)$dataMatrix, file=dataMatrix_out, sep = '\t', row.names=FALSE, quote=FALSE)
	write.table(res[[2]], file=variableMetadata_out, sep = '\t', row.names=FALSE, quote=FALSE)
}else{
	## error check
	check.err(table.check)
	
	## Launch
	plotsituation(ids,nbid,out_graph_pdf,out_preNormSummary,factbio,span)
}

}#end of meth3L