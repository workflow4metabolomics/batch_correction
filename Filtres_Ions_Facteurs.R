################################################################################################
# IONS AND FACTORS FILTERING                                                                   #
#                                                                                              #
# User : Galaxy                                                                                #
# Original data : --                                                                           #
# Starting date : 17-10-2013                                                                   #
# Version 1 : 22-10-2013                                                                       #
# Version 1.1 : 22-10-2013                                                                     #
# Version 2 : 28-11-2013 (Filter according to factors' levels added)                           #
# Version 3 : 13-03-2014 (CV ratio added)                                                      #
# Version 4 : 14-04-2014 (Voc adjustment, extension of KW to all pvalues, new CV filter added) #
# Version 4.1 : 24-04-2014 (handling of double pool deletion)                                  #
#                                                                                              #
#                                                                                              #
# Input files : file_IONS.txt ; file_PROTOC.txt ; file_MetaDataIon.txt                         #
# Output files : file_IONS_fl.txt ; file_PROTOC_fl.txt ; file_MetaDataIon_fl.txt               #
#                                                                                              #
################################################################################################

# Parameters (for dev)
if(FALSE){
  
  ion.file.in <- "test/ressources/inputs/ex_data_IONS.txt"  #tab file
  meta.samp.file.in <- "test/ressources/inputs/ex_data_PROTOCOLE1.txt"  #tab file
  meta.ion.file.in <- "test/ressources/inputs/ex_data_METAION.txt"  #tab file
  
  sep1 <- "\t" # ion.file.in separator
  sep2 <- "\t" # meta.samp.file.in separator
  sep3 <- "\t" # meta.ion.file.in separator
  
  ion.file.out <- "test/ressources/outputs/ex_data_IONS_fl.txt"  #tab file
  meta.samp.file.out <- "test/ressources/outputs/ex_data_PROTOCOLE1_fl.txt"  #tab file
  meta.ion.file.out <- "test/ressources/outputs/ex_data_METAION_fl.txt"  #tab file
  
  CV <- TRUE ; if(CV){Compa<-TRUE;seuil<-1.25}else{Compa<-NULL;seuil<-NULL}
  testp <- FALSE ; if(testp){pval<-0.05;pvar<-"5"}else{pval<-NULL;pvar<-NULL}
  POOL <- FALSE
  FACT <- FALSE ; if(FACT){ls.fact <- list(c("3","A"),c("6","s1"))}else{ls.fact <- NULL}
  
}

filter <- function(ion.file.in, meta.samp.file.in, meta.ion.file.in, sep1, sep2, sep3,
                   CV, Compa, seuil, testp, pval, pvar, POOL, FACT, ls.fact,
                   ion.file.out, meta.samp.file.out, meta.ion.file.out){
  # This function allows to filter ions, pools and samples according to differents parameters. 
  # It needs 3 datasets : the ions' data, the ions' metadata, the samples' metadata. 
  # It generates 3 new datasets corresponding to the 3 inputs filtered (if required). 
  #
  # Parameters :
  # - xxx.in : input files' names
  # - sepX : file separator
  # - xxx.out : output files' names
  # - CV : CV filter yes/no
  # | > Compa : filter comparing pool and sample CVs (TRUE) or according to pool CV (FALSE)
  # | > seuil : maximum ratio tolerated between pools' and samples' CV
  # - testp : filtering according to p-values yes/no
  # | > pval : p-value limit for testp filter
  # | > pvar : p-values column number
  # - POOL : pool deletion yes/no
  # - FACT : filter according to factors yes/no
  # | > ls.fact : factors' list for filter
  
  
# Input --------------------------------------------

ion.data <- read.table(ion.file.in,sep=sep1,header=TRUE)
meta.samp.data <- read.table(meta.samp.file.in,sep=sep2,header=TRUE)
meta.ion.data <- read.table(meta.ion.file.in,sep=sep3,header=TRUE)

# Error vector
err.stock <- "\n"


# Function 1 : CV filter ---------------------------
# Allows to filter ions according to the Coefficient of Variation (CV) :
# Compa=TRUE :
# 	CV of pools and CV of samples are compared ; if the ration between pools' one
# 	and samples' one is higher than a given ration, corresponding ion is deleted. 
# Compa=FALSE :
# 	only CV of pools are considered ; when the CV is higher than a given threshold,
# 	corresponding ion is deleted. 
if(CV){
  
  # Checking the sampleType variable
  if(is.null(meta.samp.data$sampleType)){
    err.stock <- c(err.stock,"\n-------",
                   "\nWarning : no 'sampleType' variable detected in sample meta-data !",
                   "\nCV will not be tested.\n-------\n")
  }else{
    if(!("pool"%in%levels(factor(meta.samp.data$sampleType)))){
      err.stock <- c(err.stock,"\n-------",
                     "\nWarning : no 'pool' detected in 'sampleType' variable (sample meta-data) !",
                     "\nCV will not be tested.\n-------\n")
    }else{
      if((!("sample"%in%levels(factor(meta.samp.data$sampleType))))&(Compa)){
        err.stock <- c(err.stock,"\n-------",
                       "\nWarning : no 'sample' detected in 'sampleType' variable (sample meta-data) !",
                       "\nCV will not be tested.\n-------\n")
      }else{
  
  # Statement
  tmp.ion <- data.frame(CV.ind=rep(NA,nrow(ion.data)),CV.samp=rep(NA,nrow(ion.data)),
                        CV.pool=rep(NA,nrow(ion.data)),ion.data,stringsAsFactors=FALSE)
  # CV samples
  tmp.samp <- which(colnames(tmp.ion)%in%meta.samp.data[which(meta.samp.data$sampleType=="sample"),1])
  tmp.ion$CV.samp <- apply(tmp.ion[,tmp.samp],1,sd) / rowMeans(tmp.ion[,tmp.samp])
  # CV pools
  tmp.samp <- which(colnames(tmp.ion)%in%meta.samp.data[which(meta.samp.data$sampleType=="pool"),1])
  tmp.ion$CV.pool <- apply(tmp.ion[,tmp.samp],1,sd) / rowMeans(tmp.ion[,tmp.samp])
  # CV indicator
  if(Compa){tmp.ion$CV.ind <- ifelse((tmp.ion$CV.pool)/(tmp.ion$CV.samp)>seuil,0,1)
  }else{tmp.ion$CV.ind <- ifelse((tmp.ion$CV.pool)>seuil,0,1)}
  # filter and storage ion.data
  tmp.ion <- tmp.ion[which(tmp.ion$CV.ind==1),]
  ion.data <- tmp.ion[,-c(1:3)]
  rownames(ion.data) <- NULL
  # filter and storage meta.ion.data
  meta.ion.data <- meta.ion.data[which(meta.ion.data[,1]%in%ion.data[,1]),]
  rownames(meta.ion.data) <- NULL
  
  rm(tmp.ion,tmp.samp)
  
      }}}
  
} # end if(CV)



# Function 2 : filtering according to p-values --------------
# Allows to exclude ions non-significant according to prior test 
# (according to a p-value column given in variable meta-data).
if(testp){
  
  # Checking the pvar variable
  if(!(as.numeric(pvar)%in%(1:ncol(meta.ion.data)))) {
    err.stock <- c(err.stock,"\n-------",
                   "\nWarning : no column ",pvar," detected in variable meta-data!",
                   "\nP-value filter will not be executed.\n-------\n")
  }else{
  if(!(is.numeric(meta.ion.data[,as.numeric(pvar)]))) {
    err.stock <- c(err.stock,"\n-------",
                   "\nWarning : column ",pvar," in variable meta-data is not numeric!",
                   "\nP-value filter will not be executed.\n-------\n")
  }else{
  
  # Statement
  tmp.ion <- merge(x=meta.ion.data[,c(1,as.numeric(pvar))], 
                   y=ion.data,by.x=1,by.y=1)
  # Filter
  tmp.ion <- tmp.ion[which(tmp.ion[,2]<pval),]
  # Storage
  ion.data <- tmp.ion[,-2]
  rownames(ion.data) <- NULL
  meta.ion.data <- meta.ion.data[which(meta.ion.data[,1]%in%ion.data[,1]),]
  rownames(meta.ion.data) <- NULL
  
  rm(tmp.ion)
  
  }}
  
} # end if(testp)



# Function 3 : Pools deletion ----------------------
# All lines containing 'pool' (for pools) in column "sampleType" (column defining the type of sample)
# of sample meta-data will be deleted.
if(POOL){
  
  # Checking the sampleType variable
  if(is.null(meta.samp.data$sampleType)){
    err.stock <- c(err.stock,"\n-------",
                   "\nWarning : no 'sampleType' variable detected in sample meta-data !",
                   "\nPools will not be removed.\n-------\n")
  }else{
    if(!("pool"%in%levels(factor(meta.samp.data$sampleType)))){
      err.stock <- c(err.stock,"\n-------",
                     "\nWarning : no 'pool' detected in 'sampleType' variable (sample meta-data) !",
                     "\nPools will not be removed.\n-------\n")
    }else{
  
  # Pool deletion in ion.data and meta.samp.data
  tmp.samp <- which(colnames(ion.data)%in%meta.samp.data[which(meta.samp.data$sampleType=="pool"),1])
  ion.data <- ion.data[,-c(tmp.samp)]
  meta.samp.data <- meta.samp.data[which(meta.samp.data$sampleType!="pool"),]
  
  rm(tmp.samp)
  
  }}
  
} # end if(POOL)



# Function 4 : Filter according to factors ---------
# Allows to delete all samples corresponding to selected value of designated factor.
if(FACT){

  # For each factor to filter
  for(i in 1:length(ls.fact)){
    
    # Checking the columns and factors variables
    numcol <- as.numeric(ls.fact[[i]][1])
    if(!(numcol%in%(1:ncol(meta.samp.data)))) {
    err.stock <- c(err.stock,"\n-------",
                   "\nWarning : no column ",ls.fact[[i]][1]," detected in Sample meta-data !",
                   "\nFiltering impossible for this factor.\n-------\n") 
    }else{
    if(!(ls.fact[[i]][2]%in%levels(as.factor(meta.samp.data[,numcol])))){
      err.stock <- c(err.stock,"\n-------",
                     "\nWarning : no ",ls.fact[[i]][2]," level detected in column ",numcol,
                     " (sample meta-data) !\nFiltering impossible for this factor.\n-------\n")
    }else{
      
    # Filtering
    if(length(which(meta.samp.data[,numcol]==ls.fact[[i]][2]))!=0){
      meta.samp.data <- meta.samp.data[-c(which(meta.samp.data[,numcol]==ls.fact[[i]][2])),]
      ion.data <- ion.data[,c(1,which(colnames(ion.data)%in%meta.samp.data[,1]))]
    }
  }}}

} # end if(FACT)



# Output -------------------------------------------

# Error checking
if(length(err.stock)>1){
  stop(err.stock)
}else{

write.table(ion.data, ion.file.out, sep=sep1, row.names=FALSE, quote=FALSE)
write.table(meta.samp.data, meta.samp.file.out, sep=sep2, row.names=FALSE, quote=FALSE)
write.table(meta.ion.data, meta.ion.file.out, sep=sep3, row.names=FALSE, quote=FALSE)

}


} # end of filter function


# Typical function call
#filter(ion.file.in, meta.samp.file.in, meta.ion.file.in, sep1, sep2, sep3,
#       CV, Compa, seuil, testp, pval, pvar, POOL, FACT, ls.fact,
#       ion.file.out, meta.samp.file.out, meta.ion.file.out)

