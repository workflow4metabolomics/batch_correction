# Author: jfmartin
# Modified by : mpetera
###############################################################################
# Correction of analytical effects inter and intra batch on intensities using quality control pooled samples (QC-pools)
# according to the algorithm mentioned by Van der Kloet (J Prot Res 2009).
# Parameters : a dataframe of Ions intensities and an other of samples? metadata which must contains at least the three following columns :
#   "batch" to identify the batches of analyses ; need at least 3 QC-pools for linear adjustment and 8 for lo(w)ess adjustment
#   "injectionOrder" integer defining the injection order of all samples : QC-pools and analysed samples
#   "sampleType" indicates if defining a sample with "sample" or a QC-pool with "pool"
# NO MISSING DATA are allowed
# Version 0.91 insertion of ok_norm function to assess correction feasibility 
# Version 0.92 insertion of slope test in ok_norm
# Version 0.93 name of log file define as a parameter of the correction function
# Version 0.94 Within a batch, test if all QCpools or samples values = 0. Definition of an error code in ok_norm function (see function for details)
# Version 0.99 include non linear lowess correction. 
# Version 1.00 the corrected result matrix is return transposed in Galaxy
# Version 1.01 standard deviation=0 instead of sum of value=0 is used to assess constant data in ok_norm function. Negative values in corrected matrix are converted to 0. 
# Version 1.02 plotsituation create a result file with the error code of non execution of correction set by function ok_norm
# Version 1.03 fix bug in plot with "reg" option. suppression of ok_norm=4 condition if ok_norm function
# Version 2.00 Addition of loess function, correction indicator, plots ; modification of returned objects' format, some plots' displays and ok_norm ifelse format
# Version 2.01 Correction for pools negative values earlier in norm_QCpool
# Version 2.10 Script refreshing ; vocabulary adjustment ; span in parameters for lo(w)ess regression ; conditionning for third line ACP display ; order in loess display
# Version 2.11 ok1 and ok2 permutation (ok_norm) ; conditional display of regression (plotsituation) ; grouping of linked lignes + conditioning (normX) ; conditioning for CVplot
# Version 2.20 acplight function added from previous toolBox.R [# Version 1.01 "NA"-coding possibility added in acplight function]
# Version 2.30 addition of suppressWarnings() for known and controlled warnings ; suppression of one useless "cat" message ; change in Rdata names ; 'batch(es)' in cat
# Version 2.90 change in handling of generated negative and Inf values
# Version 2.91 Plot improvement

ok_norm=function(qcp,qci,spl,spi,method) {
  # Function used for one ion within one batch to determine whether or not batch correction is possible
  # ok_norm values :
  #   0 : no preliminary-condition problem
  #   1 : standard deviation of QC-pools or samples = 0
  #   2 : insufficient number of QC-pools within a batch (n=3 for linear, n=8 for lowess or loess)
  #   3 : significant difference between QC-pools' and samples' means
  #   4 : denominator =0 when on 1 pool per batch <> 0
  #   5 : (linear regression only) the slopes ratio ?QC-pools/samples? is lower than -0.2
  
  ok=0
  if (method=="linear") {minQC=3} else {minQC=8}
  if (length(qcp)<minQC) { ok=2 
  } else {
    if (sd(qcp)==0 | sd(spl)==0) { ok=1 
    } else {
    cvp= sd(qcp)/mean(qcp); cvs=sd(spl)/mean(spl)
    rttest=t.test(qcp,y=spl)
    reslsfit=lsfit(qci, qcp)
    reslsfitSample=lsfit(spl, spi)
    ordori=reslsfit$coefficients[1]
    penteB=reslsfit$coefficients[2]
    penteS=reslsfitSample$coefficients[2]
    # Significant difference between samples and pools
    if (rttest$p.value < 0.01) { ok=3 
    } else {
      # to avoid denominator =0 when on 1 pool per batch <> 0
      if (method=="linear" & length(which(((penteB*qci)+ordori)==0))>0 ){ ok=6 
      } else {
      # different sloop between samples and pools
      if (method=="linear" & penteB/penteS < -0.20) { ok=5 }
  }}}}
  ok_norm=ok
}

plotsituation <- function (x, nbid,outfic="plot_regression.pdf", outres="PreNormSummary.txt",fact=args$batch_col_name,span="none") {
    #	Check for all ions in every batch if linear or lo(w)ess correction is possible.
    #	Use ok_norm function and create a file (PreNormSummary.txt) with the error code.
    #	Also create a pdf file with plots of linear and lo(w)ess regression lines.
    #	x: dataframe with ions in columns and samples in rows ; x is the result of concatenation of sample metadata file and ions file
    #	nbid: number of samples description columns (id and factors) with at least : "batch","injectionOrder","sampleType"
    #	outfic: name of regression plots pdf file
    #	fact: factor to be used as categorical variable for plots and PCA.  
    indfact	=which(dimnames(x)[[2]]==fact)
    indtypsamp	=which(dimnames(x)[[2]]==args$sample_type_col_name) 
    indbatch	=which(dimnames(x)[[2]]==args$batch_col_name)
    indinject	=which(dimnames(x)[[2]]==args$injection_order_col_name)
    lastIon=dim(x)[2]
    nbi=lastIon-nbid # Number of ions = total number of columns - number of identifying columns
    nbb=length(levels(x[[args$batch_col_name]])) # Number of batch = number of levels of "batch" comlumn (factor)
    nbs=length(x[[args$sample_type_col_name]][x[[args$sample_type_col_name]]=="sample"])# Number of samples = number of rows with "sample" value in sampleType
    pdf(outfic,width=27,height=7*ceiling((nbb+2)/3))
    cat(nbi," ions ",nbb," batch(es) \n")
    cv=data.frame(matrix(0,nrow=nbi,ncol=2))# initialisation de la dataset qui contiendra les CV
    pre_bilan=matrix(0,nrow=nbi,ncol=3*nbb) # dataset of ok_norm function results	    
    for (p in 1:nbi) {# for each ion
        par (mfrow=c(ceiling((nbb+2)/3),3),ask=F,cex=1.2)
        labion=dimnames(x)[[2]][p+nbid]
        indpool=which(x[[args$sample_type_col_name]]=="pool") # QCpools subscripts in x 
        pools1=x[indpool,p+nbid]; cv[p,1]=sd(pools1)/mean(pools1)# CV before correction
        for (b in 1:nbb) {# for each batch...
            xb=data.frame(x[(x[[args$batch_col_name]]==levels(x[[args$batch_col_name]])[b]),c(indtypsamp,indinject,p+nbid)])
            indpb = which(xb[[args$sample_type_col_name]]=="pool")# QCpools subscripts in the current batch
            indsp = which(xb[[args$sample_type_col_name]]=="sample")# samples subscripts in the current batch
            indbt = which(xb[[args$sample_type_col_name]]=="sample" | xb[[args$sample_type_col_name]]=="pool")# indices de tous les samples d'un batch pools+samples  
            normLinearTest=ok_norm(xb[indpb,3],xb[indpb,2], xb[indsp,3],xb[indsp,2],method="linear")
            normLoessTest=ok_norm(xb[indpb,3],xb[indpb,2], xb[indsp,3],xb[indsp,2],method="loess")
            normLowessTest=ok_norm(xb[indpb,3],xb[indpb,2], xb[indsp,3],xb[indsp,2],method="lowess")
            #cat(dimnames(x)[[2]][p+nbid]," batch ",b," loess ",normLoessTest," linear ",normLinearTest,"\n")
            pre_bilan[ p,3*b-2]=normLinearTest
            pre_bilan[ p,3*b-1]=normLoessTest
            pre_bilan[ p,3*b]=normLowessTest
          if(length(indpb)>1){
            if(span=="none"){span1<-1 ; span2<-2*length(indpool)/nbs}else{span1<-span ; span2<-span}
            resloess=loess(xb[indpb,3]~xb[indpb,2],span=span1,degree=2,family="gaussian",iterations=4,surface="direct") 
            resloessSample=loess(xb[indsp,3]~xb[indsp,2],span=2*length(indpool)/nbs,degree=2,family="gaussian",iterations=4,surface="direct") 
            reslowess=lowess(xb[indpb,2],xb[indpb,3],f=span2)
            reslowessSample=lowess(xb[indsp,2],xb[indsp,3])
            liminf=min(xb[indbt,3]);limsup=max(xb[indbt,3])
            plot(xb[indsp,2],xb[indsp,3],pch=16, main=paste(labion,"batch ",b),ylab="intensity",xlab="injection order",ylim=c(liminf,limsup))
            points(xb[indpb,2], xb[indpb,3],pch=5)
            points(cbind(resloess$x,resloess$fitted)[order(resloess$x),],type="l",col="green3")
            points(cbind(resloessSample$x,resloessSample$fitted)[order(resloessSample$x),],type="l",col="green3",lty=2)
            points(reslowess,type="l",col="red"); points(reslowessSample,type="l",col="red",lty=2)
            abline(lsfit(xb[indpb,2],xb[indpb,3]),col="blue")
            abline(lsfit(xb[indsp,2],xb[indsp,3]),lty=2,col="blue")
            legend("topleft",c("pools","samples"),lty=c(1,2),bty="n")
            legend("topright",c("linear","lowess","loess"),lty=1,col=c("blue","red","green3"),bty="n")
          }
        }
# series de plot avant et apres correction
minval=min(x[p+nbid]);maxval=max(x[p+nbid])
plot( x[[args$injection_order_col_name]], x[,p+nbid],col=x[[args$batch_col_name]],ylim=c(minval,maxval),ylab=labion,
      main=paste0("before correction (CV for pools = ",round(cv[p,1],2),")"))
suppressWarnings(plot.design( x[c(indtypsamp,indbatch,indfact,p+nbid)],main="factors effect before correction"))
    }
dev.off()
pre_bilan=data.frame(pre_bilan)
labion=dimnames(x)[[2]][nbid+1:nbi]
for (i in 1:nbb) {
  dimnames(pre_bilan)[[2]][3*i-2]=paste("batch",i,"linear")
  dimnames(pre_bilan)[[2]][3*i-1]=paste("batch",i,"loess")
  dimnames(pre_bilan)[[2]][3*i]=paste("batch",i,"lowess")
}
bilan=data.frame(labion,pre_bilan)
write.table(bilan,file=outres,sep="\t",row.names=F,quote=F)
}


normlowess=function (xb,detail="no",vref=1,b,span=NULL) {
  #	Correction function applied to 1 ion in 1 batch. Use a lowess regression computed on QC-pools in order to correct samples intensity values
  #	xb : dataframe for 1 ion in columns and samples in rows.
  # vref : reference value (average of ion)
  # b : batch subscript
  # nbid: number of samples description columns (id and factors) with at least : "batch","injectionOrder","sampleType"
  indpb = which(xb[[args$sample_type_col_name]]=="pool") # pools subscripts of current batch
  indsp = which(xb[[args$sample_type_col_name]]=="sample") # samples of current batch subscripts
  indbt = which(xb[[args$sample_type_col_name]]=="sample" | xb[[args$sample_type_col_name]]=="pool");# batch subscripts of all samples and QC-pools
  labion=dimnames(xb)[[2]][3]
  newval=xb[[3]] # initialisation of corrected values = intial values
  ind <- 0 # initialisation of correction indicator
  normTodo=ok_norm(xb[indpb,3],xb[indpb,2], xb[indsp,3],xb[indsp,2],method="lowess")
  #cat("batch:",b," dim xb=",dim(xb)," ok=",normTodo,"\n")
  if (normTodo==0) {
    if(length(span)==0){span2<-2*length(indpb)/length(indsp)}else{span2<-span}
    reslowess=lowess(xb[indpb,2],xb[indpb,3],f=span2) # lowess regression with QC-pools 
    px=xb[indsp,2];  # vector of injectionOrder values only for samples  
    for(j in 1:length(indbt)) {	 
      if (xb[[args$sample_type_col_name]][j]=="pool") {
        if (reslowess$y[which(indpb==j)]==0) reslowess$y[which(indpb==j)] <- 1
        newval[j]=(vref*xb[j,3]) / (reslowess$y[which(indpb==j)])} 
      else { # for samples, the correction value cor correspond to the nearest QCpools 
        cor= reslowess$y[which(abs(reslowess$x-px[which(indsp==j)])==min(abs(reslowess$x - px[which(indsp==j)])))]
        if (length(cor)>1) {cor=cor[1]}
        if (cor <= 0) {cor=vref} # no modification of initial value
        newval[j]=(vref*xb[j,3]) / cor
      }
    }
    if (detail=="reg") {
      liminf=min(xb[indbt,3]);limsup=max(xb[indbt,3])
      plot(xb[indsp,2],xb[indsp,3],pch=16,main=paste(labion,"batch ",b),ylab="intensity",xlab="injection order",ylim=c(liminf,limsup))
      points(xb[indpb,2], xb[indpb,3],pch=5)
      points(reslowess,type="l",col="red")
    }
    ind <- 1
  } else {# if ok_norm <> 0 , we perform a correction based on batch samples average
    moySample=mean(xb[indsp,3]);if (moySample==0) moySample=1
    newval[indsp] = (vref*xb[indsp,3])/moySample
    if(length(indpb)>0){
      moypool=mean(xb[indpb,3]) ; if (moypool==0) moypool=1
      newval[indpb] = (vref*xb[indpb,3])/moypool
    }
  }
  newval <- list(norm.ion=newval,norm.ind=ind)
  return(newval)
}

normlinear <- function (xb,detail="no",vref=1,b,valneg=0) {
  # Correction function applied to 1 ion in 1 batch.
  # Use a linear regression computed on QC-pools in order to correct samples intensity values
  # xb: dataframe with ions in columns and samples in rows; x is a result of concatenation of sample metadata file and ion file 
  # nbid: number of sample description columns (id and factors) with at least "batch", "injectionOrder" and "sampleType"
  # b: which batch it is
  # valneg: to determine what to do with generated negative and Inf values
  indpb = which(xb[[args$sample_type_col_name]]=="pool")# pools subscripts of current batch
  indsp = which(xb[[args$sample_type_col_name]]=="sample")# samples of current batch subscripts
  indbt = which(xb[[args$sample_type_col_name]]=="sample" | xb[[args$sample_type_col_name]]=="pool") # QCpools and samples of current batch subscripts
  labion=dimnames(xb)[[2]][3]
  newval=xb[[3]] # initialisation of corrected values = intial values	
  ind <- 0 # initialisation of correction indicator
  normTodo=ok_norm(xb[indpb,3],xb[indpb,2], xb[indsp,3],xb[indsp,2],method="linear")
  if (normTodo==0) {
    ind <- 1
    reslsfit=lsfit(xb[indpb,2],xb[indpb,3])       # linear regression for QCpools 
    reslsfitSample=lsfit(xb[indsp,2],xb[indsp,3]) # linear regression for samples
    ordori=reslsfit$coefficients[1]
    pente=reslsfit$coefficients[2]
    if (detail=="reg") {
      liminf=min(xb[indbt,3]);limsup=max(xb[indbt,3])
      plot(xb[indsp,2],xb[indsp,3],pch=16,
           main=paste(labion,"batch ",b),ylab="intensity",xlab="injection order",ylim=c(liminf,limsup))
      points(xb[indpb,2], xb[indpb,3],pch=5)
      abline(reslsfit)
      abline(reslsfitSample,lty=2)
    }
    # correction with rescaling of ion global intensity (vref)
    newval = (vref*xb[indbt,3]) / (pente * (xb[indbt,2]) + ordori)
	newval[which((pente * (xb[indbt,2]) + ordori)<1)] <- -1 # to handle cases where 0<denominator<1
	# handling if any negative values (or null denominators)
	if(length(which((newval==Inf)|(newval<0)))!=0){
	  toajust <- which((newval==Inf)|(newval<0))
	  if(valneg=="NA"){
	    newval[toajust] <- NA
	  } else {
	    newval[toajust] <- vref * (xb[indbt,3][toajust]) / mean(xb[indbt,3])
    ### Other possibility
	##  if(pente>0){ # slope orientation
	##    newval[toajust]<-(vref*(xb[indbt,3][toajust]))/(pente*ceiling(-ordori/pente+1.00001)+ordori)
	##  }else{
	##    newval[toajust]<-(vref*(xb[indbt,3][toajust]))/(pente*floor(-ordori/pente-1.00001)+ordori)
	##  }
	  }
	}
  } else {# if ok_norm!=0 , we perform a correction based on batch samples average.
    moySample=mean(xb[indsp,3]); if (moySample==0) moySample=1
    newval[indsp] = (vref*xb[indsp,3])/moySample
    if(length(indpb)>0){
      moypool=mean(xb[indpb,3]) ; if (moypool==0) moypool=1
      newval[indpb] = (vref*xb[indpb,3])/moypool
    }
  }
  newval <- list(norm.ion=newval,norm.ind=ind)
  return(newval)
}


normloess <- function (xb,detail="no",vref=1,b,span=NULL) {
    #	Correction function applied to 1 ion in 1 batch. 
    #   Use a loess regression computed on QC-pools in order to correct samples intensity values.
    #	xb : dataframe for 1 ion in columns and samples in rows.
    #   detail : level of detail in the outlog file.
    #   vref : reference value (average of ion)
    #   b : batch subscript
    indpb = which(xb[[args$sample_type_col_name]]=="pool") # pools subscripts of current batch
    indsp = which(xb[[args$sample_type_col_name]]=="sample") # samples of current batch subscripts
    indbt = which(xb[[args$sample_type_col_name]]=="sample" | xb[[args$sample_type_col_name]]=="pool");# batch subscripts of all samples and QCpools
    labion=dimnames(xb)[[2]][3]
    newval=xb[[3]] # initialisation of corrected values = intial values
    ind <- 0 # initialisation of correction indicator
    normTodo=ok_norm(xb[indpb,3],xb[indpb,2], xb[indsp,3],xb[indsp,2],method="loess")
    #cat("batch:",b," dim xb=",dim(xb)," ok=",normTodo,"\n")
    if (normTodo==0) { 
        if(length(span)==0){span1<-1}else{span1<-span}
        resloess=loess(xb[indpb,3]~xb[indpb,2],span=span1,degree=2,family="gaussian",iterations=4,surface="direct") # loess regression with QCpools 
        cor=predict(resloess,newdata=xb[,2])
		    cor[cor<=1] <- 1 
        newval=(vref*xb[,3]) / cor 
        if(length(which(newval>3*(quantile(newval)[4])))>0){ # in this case no modification of initial value
			newval <- xb[,3]} else {ind <- 1} # confirmation of correction
        if ((detail=="reg")&(ind==1)) { # plot
            liminf=min(xb[indbt,3]);limsup=max(xb[indbt,3])
            plot(xb[indsp,2],xb[indsp,3],pch=16,main=paste(labion,"batch ",b),ylab="intensity",xlab="injection order",ylim=c(liminf,limsup))
            points(xb[indpb,2], xb[indpb,3],pch=5)
            points(cbind(resloess$x,resloess$fitted)[order(resloess$x),],type="l",col="red")
        }
    } 
    if (ind==0) {# if ok_norm != 0 or if correction creates outliers, we perform a correction based on batch samples average
      moySample=mean(xb[indsp,3]);if (moySample==0) moySample=1
      newval[indsp] = (vref*xb[indsp,3])/moySample
      if(length(indpb)>0){
        moypool=mean(xb[indpb,3]) ; if (moypool==0) moypool=1
        newval[indpb] = (vref*xb[indpb,3])/moypool
      }
    }
    newval <- list(norm.ion=newval,norm.ind=ind)
    return(newval)
}



norm_QCpool <- function (x, nbid, outlog, fact, metaion, detail="no", NormMoyPool=F, NormInt=F, method="linear",span="none",valNull="0")
{
  ### Correction applying linear or lo(w)ess correction function on all ions for every batch of a dataframe.
  # x: dataframe with ions in column and samples' metadata
  # nbid: number of sample description columns (id and factors) with at least "batch", "injectionOrder", "sampleType"
  # outlog: name of regression plots and PCA pdf file
  # fact: factor to be used as categorical variable for plots  
  # metaion: dataframe of ions' metadata
  # detail: level of detail in the outlog file. detail="no" ACP + boxplot of CV before and after correction.
  #          detail="plot" with plot for all batch before and after correction.
  #          detail="reg" with added plots with regression lines for all batches.
  # NormMoyPool: not used
  # NormInt: not used
  # method: regression method to be used to correct : "linear" or "lowess" or "loess"
  # valNull: to determine what to do with negatively estimated intensities
	indfact		=which(dimnames(x)[[2]]==fact)
	indtypsamp=which(dimnames(x)[[2]]==args$sample_type_col_name)
	indbatch	=which(dimnames(x)[[2]]==args$batch_col_name)
	indinject	=which(dimnames(x)[[2]]==args$injection_order_col_name)
	lastIon=dim(x)[2]
	valref=apply(as.matrix(x[,(nbid+1):(lastIon)]),2,mean) # reference value for each ion used to still have the same rought size of values
	nbi=lastIon-nbid # number of ions
	nbb=length(levels(x[[args$batch_col_name]])) # Number of batch(es) = number of levels of factor "batch" (can be =1)
	nbs=length(x[[args$sample_type_col_name]][x[[args$sample_type_col_name]]=="sample"])# Number of samples 
	nbp=length(x[[args$sample_type_col_name]][x[[args$sample_type_col_name]]=="pool"])# Number of QCpools
	Xn=data.frame(x[,c(1:nbid)],matrix(0,nrow=nbp+nbs,ncol=nbi))# initialisation of the corrected dataframe (=initial dataframe)
  dimnames(Xn)=dimnames(x)
	cv=data.frame(matrix(0,nrow=nbi,ncol=2))# initialisation of dataframe containing CV before and after correction
	dimnames(cv)[[2]]=c("avant","apres")
	if (detail!="reg" && detail!="plot" && detail!="no") {detail="no"}
	pdf(outlog,width=27,height=20)
	cat(nbi," ions ",nbb," batch(es) \n")
	if (detail=="plot") {if(nbb<6){par(mfrow=c(3,3),ask=F,cex=1.5)}else{par(mfrow=c(4,4),ask=F,cex=1.5)}}
  res.ind <- matrix(NA,ncol=nbb,nrow=nbi,dimnames=list(dimnames(x)[[2]][-c(1:nbid)],paste("norm.b",1:nbb,sep="")))
	for (p in 1:nbi) {# for each ion
	  labion=dimnames(x)[[2]][p+nbid]
        if (detail == "reg") {if(nbb<6){par(mfrow=c(3,3),ask=F,cex=1.5)}else{par(mfrow=c(4,4),ask=F,cex=1.5)}}
		indpool=which(x[[args$sample_type_col_name]]=="pool")# QCpools subscripts in all batches
		pools1=x[indpool,p+nbid]; cv[p,1]=sd(pools1)/mean(pools1)# CV before correction
		for (b in 1:nbb) {# for every batch
			indpb = which(x[[args$batch_col_name]]==levels(x[[args$batch_col_name]])[b] & x[[args$sample_type_col_name]]=="pool")# QCpools subscripts of the current batch	
			indsp = which(x[[args$batch_col_name]]==levels(x[[args$batch_col_name]])[b] & x[[args$sample_type_col_name]]=="sample")# samples subscripts of the current batch
      indbt = which(x[[args$batch_col_name]]==levels(x[[args$batch_col_name]])[b] & (x[[args$sample_type_col_name]]=="pool" | x[[args$sample_type_col_name]]=="sample")) # subscripts of all samples
      # cat(dimnames(x)[[2]][p+nbid]," indsp:",length(indsp)," indpb=",length(indpb)," indbt=",length(indbt)," ")
			sub=data.frame(x[(x[[args$batch_col_name]]==levels(x[[args$batch_col_name]])[b]),c(indtypsamp,indinject,p+nbid)])
			if (method=="linear") { res.norm = normlinear(sub,detail,valref[p],b,valNull)
			} else { if (method=="loess"){ res.norm <- normloess(sub,detail,valref[p],b,span)
			  } else { if (method=="lowess"){ res.norm <- normlowess(sub,detail,valref[p],b,span)
          } else {stop("\n--\nNo valid 'method' argument supplied.\nMust be 'linear','loess' or 'lowess'.\n--\n")}
			}}
			Xn[indbt,p+nbid] = res.norm[[1]]
			res.ind[p,b] <- res.norm[[2]]
			# CV batch test : if after normaliszation, CV before < CV after initial values are kept
#  			moypoolRaw=mean(x[indpb,p+nbid])  ; if (moypoolRaw==0) moypoolRaw=1
#			moySampleRaw=mean(x[indsp,p+nbid]); if (moySampleRaw==0) moySampleRaw=1
#			moypool=mean(Xn[indpb,p+nbid]) ; if (moypool==0) moypool=1
#			#moySample=mean(Xn[indsp,p+nbid]); if (moySample==0) moySample=1
#			if (sd( Xn[indpb,p+nbid])/moypool>sd(x[indpb,p+nbid])/moypoolRaw) { 
#					Xn[indpb,p+nbid] = (valref[p]*x[indpb,p+nbid])/moypoolRaw
#					Xn[indsp,p+nbid] = (valref[p]*x[indsp,p+nbid])/moySampleRaw
#			}
		}	
		Xn[indpool,p+nbid][Xn[indpool,p+nbid]<0] <- 0
		pools2=Xn[indpool,p+nbid]; cv[p,2]=sd(pools2,na.rm=TRUE)/mean(pools2,na.rm=TRUE)# CV apres correction
		if (detail=="reg" || detail=="plot" ) {
		  	# plot before and after correction
		  	minval=min(cbind(x[p+nbid],Xn[p+nbid]),na.rm=TRUE);maxval=max(cbind(x[p+nbid],Xn[p+nbid]),na.rm=TRUE)
		  	plot( x[[args$injection_order_col_name]], x[,p+nbid],col=x[[args$batch_col_name]],ylab=labion,ylim=c(minval,maxval),
              main=paste0("before correction (CV for pools = ",round(cv[p,1],2),")"))
              points(x[[args$injection_order_col_name]][indpool],x[indpool,p+nbid],col="maroon",pch=".",cex=2)
		  	plot(Xn[[args$injection_order_col_name]],Xn[,p+nbid],col=x[[args$batch_col_name]],ylab="",ylim=c(minval,maxval),
             main=paste0("after correction (CV for pools = ",round(cv[p,2],2),")"))
              points(Xn[[args$injection_order_col_name]][indpool],Xn[indpool,p+nbid],col="maroon",pch=".",cex=2)
		  	suppressWarnings(plot.design( x[c(indtypsamp,indbatch,indfact,p+nbid)],main="factors effect before correction"))
		  	suppressWarnings(plot.design(Xn[c(indtypsamp,indbatch,indfact,p+nbid)],main="factors effect after correction"))
		}
	}
  ### Replacement of post correction negative values by chosen value
	Xnn=Xn
	for (i in c((nbid+1):dim(Xn)[2])) {
	  cneg=which(Xn[[i]]<0)
	  Xnn[[i]]=replace(Xn[[i]],cneg,as.numeric(valNull))
	}
  Xn=Xnn

	if (detail=="reg" || detail=="plot" || detail=="no") {
		if (nbi > 3) {
			par(mfrow=c(3,4),ask=F,cex=1.2) # PCA Plot before/after, normed only and ions plot
			acplight(x[,c(indtypsamp,indbatch,indtypsamp,indfact,(nbid+1):lastIon)],"uv",TRUE)
			norm.ion <- which(colnames(Xn)%in%(rownames(res.ind)[which(rowSums(res.ind)>=1)]))
			acplight(Xn[,c(indtypsamp,indbatch,indtypsamp,indfact,(nbid+1):lastIon)],"uv",TRUE,norm.ion)
			if(length(norm.ion)>0){acplight(Xn[,c(indtypsamp,indbatch,indtypsamp,indfact,norm.ion)],"uv",TRUE)}
			par(mfrow=c(1,2),ask=F,cex=1.2) # Before/after boxplot
			cvplot=cv[!is.na(cv[[1]])&!is.na(cv[[2]]),]
      if(nrow(cvplot)>0){
			  boxplot(cvplot[[1]],ylim=c(min(cvplot),max(cvplot)),main="CV before correction")
        boxplot(cvplot[[2]],ylim=c(min(cvplot),max(cvplot)),main="CV after correction")
      }
      dev.off()
		}
	}
  if (nbi<=3) {dev.off()}
  # transposed matrix is return  (format of the initial matrix with ions in rows)
  Xr=Xn[,-c(1:nbid)]; dimnames(Xr)[[1]]=Xn[[1]]
  Xr=t(Xr) ; Xr <- data.frame(ions=rownames(Xr),Xr)
  
  res.norm[[1]] <- Xr ; res.norm[[2]] <- data.frame(metaion,res.ind) ; res.norm[[3]] <- x[,c(1:nbid)]
  names(res.norm) <- c("dataMatrix","variableMetadata","sampleMetadata")
  return(res.norm)
}





acplight <- function(ids, scaling="uv", indiv=FALSE,indcol=NULL) {
	suppressPackageStartupMessages(library(ade4))
	suppressPackageStartupMessages(library(pcaMethods))
  # Make a PCA and plot scores and loadings.
  # First column must contain samples' identifiers.
  # Columns 2 to 4 contain factors to colour the plots. 
  for (i in 1:3) {
    idss=ids[which(ids[,i+1]!="NA"),]
    idss=data.frame(idss[idss[,i+1]!="",])
    classe=as.factor(idss[[i+1]])
    idsample=as.character(idss[[1]])
    colour=1:length(levels(classe))
    ions=as.matrix(idss[,5:dim(idss)[2]])
	# Removing ions containing NA (not compatible with standard PCA)
	ions=t(na.omit(t(ions)))
	if(i==1){if(ncol(ions)!=(ncol(idss)-4)){cat("Note:",(ncol(idss)-4)-ncol(ions),"ions were ignored for PCA display due to NA in intensities.\n")}}
    # Scaling choice: "uv","none","pareto"
    object=suppressWarnings(prep(ions, scale=scaling, center=TRUE))
	if(i==1){if(length(which(apply(ions,2,var)==0))>0){cat("Warning: there are",length(which(apply(ions,2,var)==0)),"constant ions.\n")}}
    # ALGO: nipals,svdImpute, Bayesian, svd, probalistic=F
    result <- pca(object, center=F, method="svd", nPcs=2)
    # ADE4 : to plot samples' ellipsoid for each class
    s.class(result@scores, classe, cpoint = 1,xax=1,yax=2,col=colour,sub=sprintf("Scores - PCs %sx%s",1,2), possub="bottomright")
    #s.label(result@loadings,label = ions, cpoint = 0, clabel=0.4, xax=1,yax=2,sub="Loadings",possub="bottomright")
    if(i==1){resulti <- result}
  }
  if(indiv) {
    colour <- rep("darkblue",length(resulti@loadings)) ; if(!is.null(indcol)) {colour[-c(indcol)] <- "red"}
    plot(resulti@loadings,col=colour,main="Loadings",xaxt="n",yaxt="n",pch=20,
         xlab=bquote(PC1-R^2==.(resulti@R2[1])),ylab=bquote(PC2 - R^2 == .(resulti@R2[2])))
    abline(h=0,v=0)}
}


