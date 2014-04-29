#################################################################################################
# NORMALISATION ACCORDING TO A GIVEN FACTOR                                                     #
#                                                                                               #
# User : Galaxy                                                                                 #
# Original data : --                                                                            #
# Starting date : 18-10-2013                                                                    #
# Version 1 : 23-10-2013                                                                        #
# Version 2 : 06-11-2013 (plots before/after added)                                             #
#                                                                                               #
#                                                                                               #
# Input files : file_IONS.txt ; file_PROTOC.txt                                                 #
# Output files : file_IONS_norm.txt                                                             #
#                                                                                               #
# Called libraries : ade4                                                                       #
#                                                                                               #
#################################################################################################

if(FALSE){
	# For dev
	rm(list=ls())
	setwd("J:/WorkSpace/ScriptR_divers/Galaxy")
	ion.file.in <- "DataExemple/NormVDK_mat_201_mz_IONS.txt"  #tab file
	meta.samp.file.in <- "DataExemple/NormVDK_mat_201_mz_PROTOCOLE.txt"  #tab file
	ion.file.out <- "DataExemple/NormVDK_mat_201_mz_IONS_norm.txt"  #tab file
	norm.type <- "weighted.means"
	norm.factor <- "centre"
	weight.factor <- "repondeur"
	graph.pdf.out <- paste("Norm_by_",norm.factor,"_plots.pdf",sep="")
	acp.pdf.out <- paste("Norm_by_",norm.factor,"_acp.pdf",sep="")
}

# Parameters
if(FALSE){
  
  ion.file.in <- "DataExemple/ex_data_IONS.txt"  #tab file
  meta.samp.file.in <- "DataExemple/ex_data_PROTOCOLE.txt"  #tab file
  ion.file.out <- "DataExemple/ex_data_IONS_norm.txt"  #tab file
  
  norm.type <- "weighted.means"
  
  norm.factor <- "centre"
  weight.factor <- "repondeur" 
  
  graph.pdf.out <- paste("Norm_by_",norm.factor,"_plots.pdf",sep="") #pdf file name
  acp.pdf.out <- paste("Norm_by_",norm.factor,"_acp.pdf",sep="") #pdf2 file name
  
  
}


norm.by.fact <- function(
  ion.file.in, meta.samp.file.in, norm.type, norm.factor, weigth.factor,ion.file.out){
  # This function allows to normalise ions' intensity according to a given factor, taking 
  # into account a weight factor corresponding to a variable potentially important that may 
  # not be well distributed relatively to the factor of normalisation. 
  # 
  # Parameters :
  # - ion.file.in : path to the ions' intensity file
  # - meta.samp.file.in : path to the samples' metadata file
  # - ion.file.out : path for the normalised ions' intensity file
  # - norm.factor : label of the factor of interest (FI) used for normalisation
  # - weight.factor : label of the weight factor (WF) 
  #
  # Calculation : 
  # For each FI's modality, means of corresponding samples are calculated from the means of 
  # each WF's modality, weighted by these modalities' global frequencies.
  # From these FI's modality means, a normalisation is done by dividing each intensity by 
  # its corresponding FI's modality mean, and then multiplying by the global mean of original
  # intensities. 
  #
  # Warning :
  # This function is meant to be used for ions' intensity data without pools. 
  # Please make sure that the ions' intensity file and the samples' metadata file do NOT
  # contain pools. 


# Input -----------------------------------------------------------------------------------------

ion.data <- read.table(ion.file.in,sep="\t",header=TRUE,stringsAsFactors=FALSE)
meta.samp.data <- read.table(meta.samp.file.in,sep="\t",header=TRUE,stringsAsFactors=FALSE)


# Preliminary processing ------------------------------------------------------------------------

# Merging the tables
tmp.data <- t(ion.data) 
colnames(tmp.data) <- tmp.data[1,]
tmp.data <- tmp.data[-1,]
all.data <- merge(x=meta.samp.data,y=tmp.data,by.x=1,by.y=0)
all.data[,-c(1:ncol(meta.samp.data))] <- 
  apply(apply(all.data[,-c(1:ncol(meta.samp.data))],2,as.character),2,as.numeric)

# Ion intensity means
ion.mean <-rowMeans(ion.data[,-1])

# Columns numbers
numcol.norm <- which(colnames(all.data)==norm.factor)
numcol.weight <- which(colnames(all.data)==weight.factor)

# Modality counts
norm.count <- as.data.frame(table(all.data[,numcol.norm]))
weight.count <- as.data.frame(table(all.data[,numcol.weight]))


#################################################################################################
###### Normalisation of ions' intensity according to weighted means of modality #################

if(norm.type=="weighted.means"){

# Weighted means of each modality for the factor of interest ------------------------------------

# Statement
adj.means <- matrix(0,nrow=nrow(norm.count),ncol=nrow(ion.data))
rownames(adj.means) <- norm.count[,1] ; colnames(adj.means) <- ion.data[,1]

# Matrix filling
for(modality in as.character(norm.count[,1])){
  # Means by modalities for the weight factor 
  tmp.vect <-  matrix(NA,nrow=nrow(weight.count),ncol=nrow(ion.data))
  rownames(tmp.vect) <- weight.count[,1] ; colnames(tmp.vect) <- ion.data[,1]
  for(wgt in as.character(weight.count[,1])){
    tmp.data <- all.data[which(all.data[,numcol.norm]==modality),]
    tmp.data <- tmp.data[which(tmp.data[,numcol.weight]==wgt),]
    if(nrow(tmp.data)!=0){#if wgt exists for the modality
    tmp.vect[which(rownames(tmp.vect)==wgt),] <- colMeans(tmp.data[,-c(1:ncol(meta.samp.data))])
    }else{tmp.vect[which(rownames(tmp.vect)==wgt),] <- NA}
  }
  # Mean for a modality of the factor of interest
  # (weighted by global frequencies of weight factor)
  glob.count <- 0
  for(wgt in as.character(weight.count[,1])){
    if(length(which(is.na(tmp.vect[which(rownames(tmp.vect)==wgt),])))==0){#if wgt exists
    adj.means[which(rownames(adj.means)==modality),] <-
      adj.means[which(rownames(adj.means)==modality),] + 
      tmp.vect[which(rownames(tmp.vect)==wgt),] * weight.count[which(weight.count[,1]==wgt),2]
    glob.count <- glob.count + weight.count[which(weight.count[,1]==wgt),2]
    }
  }
  adj.means[which(rownames(adj.means)==modality),] <-
    adj.means[which(rownames(adj.means)==modality),] / glob.count
}
adj.means[adj.means==0] <- 1


# Normalisation according to weighted means of modality -----------------------------------------

for (samp in 1:nrow(all.data)) {
  all.data[samp,-c(1:ncol(meta.samp.data))] <- 
    all.data[samp,-c(1:ncol(meta.samp.data))] * ion.mean / 
    adj.means[which(rownames(adj.means)==all.data[samp,numcol.norm]),]
  if(length(which(is.na(all.data)))!=0){print(samp)} # dev : print le num d'echant si pb de NA
}



} # End if(norm.type)


#################################################################################################
###### Plots "Before/After" #####################################################################

# "Before" table
col.plot <- as.factor(meta.samp.data[,numcol.norm])
levels(col.plot) <- rep(c("purple","green","orange","blue","red"),5)[1:length(levels(col.plot))]
col.plot <- as.character(col.plot)
tempo1 <- data.frame(col.plot,t(ion.data[,-1]),stringsAsFactors=FALSE)
tempo1 <- tempo1[order(tempo1[,1],rownames(tempo1)),]

# "After" table
tempo2 <- data.frame(all.data[,numcol.norm],all.data[,-c(1:ncol(meta.samp.data))])
rownames(tempo2) <- all.data[,1]
levels(tempo2[,1]) <- rep(c("purple","green","orange","blue","red"),5)[1:length(levels(tempo2[,1]))]
tempo2[,1] <- as.character(tempo2[,1])
tempo2 <- tempo2[order(tempo2[,1],rownames(tempo2)),]



pdf(graph.pdf.out,height=1.5,width=12)

for(i in 1:nrow(ion.data)) {
  par(mfrow=c(1,4),mai=c(0.5,0.3,0.2,0.2))
  
# Plots before/after colored by norm.factor
plot(tempo1[,i+1],col=tempo1[,1],ylab="",xlab="",cex.main=0.7,cex.axis=0.7,cex=0.5,lwd=0.5,
     main=paste(ion.data[i,1],"Before normalisation",sep="  -  "))
plot(tempo2[,i+1],col=tempo2[,1],ylab="",xlab="",cex.main=0.7,cex.axis=0.7,cex=0.5,lwd=0.5,
     main=paste(ion.data[i,1],"After normalisation",sep="  -  "))

# plot.design before/after
plot.design(t(ion.data[i,-1])~as.factor(meta.samp.data[,numcol.norm])+
              as.factor(meta.samp.data[,numcol.weight]),xlab="",ylab="",xaxt="n",
            main="Before",cex.main=0.7,cex.axis=0.7,cex=0.5)
plot.design(all.data[,i+ncol(meta.samp.data)]~as.factor(all.data[,numcol.norm])+
              as.factor(all.data[,numcol.weight]),xlab="",ylab="",xaxt="n",
            main="After",cex.main=0.7,cex.axis=0.7,cex=0.5)

}

dev.off()



# ACP plots : factorial axes 1-2 & 3-4 before/after
pdf(acp.pdf.out)
par(mfrow=c(2,2))

require(ade4)

# "Before" plots
res.pca <- dudi.pca(t(ion.data[,-1]),scale=TRUE,center=TRUE,scannf=FALSE,nf=4)
s.class(dfxy=res.pca$li,col=levels(as.factor(col.plot)),xax=1,yax=2,possub="topright",cpoint=0.8,
    fac=as.factor(meta.samp.data[,numcol.norm]),sub="Before  \nFactorial axes 1&2  ",grid=FALSE)
s.class(dfxy=res.pca$li,col=levels(as.factor(col.plot)),xax=3,yax=4,possub="topright",cpoint=0.8,
    fac=as.factor(meta.samp.data[,numcol.norm]),sub="Before  \nFactorial axes 3&4  ",grid=FALSE)

# "After" plots
res.pca <- dudi.pca(all.data[,-c(1:ncol(meta.samp.data))],scale=TRUE,center=TRUE,scannf=FALSE,nf=4)
s.class(dfxy=res.pca$li,col=levels(as.factor(col.plot)),xax=1,yax=2,possub="topright",cpoint=0.8,
    fac=as.factor(meta.samp.data[,numcol.norm]),sub="After  \nFactorial axes 1&2  ",grid=FALSE)
s.class(dfxy=res.pca$li,col=levels(as.factor(col.plot)),xax=3,yax=4,possub="topright",cpoint=0.8,
    fac=as.factor(all.data[,numcol.norm]),sub="After  \nFactorial axes 3&4  ",grid=FALSE)

dev.off()


#################################################################################################
#################################################################################################

# Final processing and output -------------------------------------------------------------------

# Formating all.data into ion.data format
rownames(all.data) <- all.data[,1]
all.data <- all.data[,-c(1:ncol(meta.samp.data))]
all.data <- t(all.data)
all.data <- data.frame(ions=rownames(all.data),all.data)
colnames(all.data)[1] <- colnames(ion.data)[1]

# Export
write.table(all.data, ion.file.out, sep="\t", row.names=FALSE, quote=FALSE)



} # end of norm.factor function



# Typical function call
norm.by.fact(ion.file.in, meta.samp.file.in, norm.type, norm.factor, weigth.factor,ion.file.out)

