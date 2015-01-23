###############################################
# R check library
# Coded by : M.Petera, 
# - -
# R functions to use in R scripts  
# (management of various generic subroutines)
# - -
# V0 : script structure + first functions
###############################################


# Generic function to return an error if problems have been encountered - - - -

check.err <- function(err.stock){
	
	# err.stock = vector of results returned by check functions
	
	if(length(err.stock)!=0){ stop("\n- - - - - - - - -\n",err.stock,"\n- - - - - - - - -\n") }
	
}




# Table match check functions - - - - - - - - - - - - - - - - - - - - - - - - -

# To check if dataMatrix and (variable or sample)Metadata match regarding identifiers
match2 <- function(dataMatrix, Metadata, Mtype){

	# dataMatrix = data.frame containing dataMatrix
	# Metadata = data.frame containing sampleMetadata or variableMetadata
	# Mtype = "sample" or "variable" depending on Metadata content
	
	err.stock <- NULL # error vector

	id2 <- Metadata[,1]
	if(Mtype=="sample"){ id1 <- colnames(dataMatrix)[-1] }
	if(Mtype=="variable"){ id1 <- dataMatrix[,1] }
	
	if( length(which(id1%in%id2))!=length(id1) || length(which(id2%in%id1))!=length(id2) ){
		err.stock <- c("\nData matrix and ",Mtype," metadata do not match regarding ",Mtype," identifiers.\n",
					"Please check your data.\n")
	}
	
	return(err.stock)
	
}

# To check if the 3 standard tables match regarding identifiers
match3 <- function(dataMatrix, sampleMetadata, variableMetadata){
	
	# dataMatrix = data.frame containing dataMatrix
	# sampleMetadata = data.frame containing sampleMetadata
	# variableMetadata = data.frame containing variableMetadata
	
	err.stock <- NULL # error vector
	
	id1 <- colnames(dataMatrix)[-1]
	id2 <- sampleMetadata[,1]
	id3 <- dataMatrix[,1]
	id4 <- variableMetadata[,1]
	
	if( length(which(id1%in%id2))!=length(id1) || length(which(id2%in%id1))!=length(id2) ){
		err.stock <- c(err.stock,"\nData matrix and sample metadata do not match regarding sample identifiers.")
	}
	
	if( length(which(id3%in%id4))!=length(id3) || length(which(id4%in%id3))!=length(id4) ){
		err.stock <- c(err.stock,"\nData matrix and variable metadata do not match regarding variable identifiers.")
	}
	
	if(length(err.stock)!=0){ err.stock <- c(err.stock,"\nPlease check your data.\n") }
	
	return(err.stock)
	
}