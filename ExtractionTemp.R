####### Module de remplacement temporaire de diffreport sous Galaxy



if(FALSE){
	nbSamp <- 36
	#rdatafile <- "J:/WorkSpace/Galaxy12-[xset.group.retcor.group.RData].rdata"
	rdatafile <- "J:/WorkSpace/Galaxy17-[xset.group.retcor.group.fillPeaks.RData].rdata"
	varmd.out <- "J:/WorkSpace/VariableMetadata.txt"
	dm.out <- "J:/WorkSpace/DataMatrix.txt"
	Perfwhm <- 0.6
	Polarity <- "positive"
	
}

extractionTemp <- function(rdatafile,nbSamp,Perfwhm,Polarity,varmd.out,dm.out){
# rdatafile = chemin et nom du Rdata
# nbSamp = nb de samples (contient pool & sample & blank & tout autre type d'echantillon)
# varmd.out = chemin et nom du fichier de sortie des variable metadata
# dm.out = chemin et nom du fichier de sortie de la data matrix

load(rdatafile)

library(xcms)
library(CAMERA)

xsetPnofill <- xset

##...xcmsSet, group, retcorr,...puis la derniere étape fillPeaks
xsetP <- fillPeaks(xsetPnofill)
## puis CAMERA
an <-xsAnnotate(xsetP)
an <-groupFWHM(an, perfwhm = Perfwhm)
an <-findIsotopes(an, mzabs = 0.01)
an <-groupCorr(an, cor_eic_th = 0.75)
anP <-findAdducts(an, polarity=Polarity)

thelist <- getPeaklist(anP)

fin.var <- ncol(thelist) - 3 - nbSamp

varmd.tb <- data.frame(thelist[,1:fin.var],thelist[,(ncol(thelist)-2):(ncol(thelist))])
write.table(varmd.tb, file=varmd.out,sep="\t", row.names=F)

dm.tb <- data.frame(thelist[,1,drop=FALSE],thelist[,(fin.var+1):(ncol(thelist)-3)])
write.table(dm.tb, file=dm.out,sep="\t", row.names=F)


}

# Typical function call
#extractionTemp(rdatafile,nbSamp,Perfwhm,Polarity,varmd.out,dm.out)
