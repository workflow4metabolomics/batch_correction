####### Module de remplacement temporaire de diffreport sous Galaxy



if(FALSE){
  nbSamp <- 36
  rdatafile <- "test/ressources/inputs/Galaxy15-[xset.group.retcor.group.fillpeaks.RData].rdata"
  dm.out <- "test/ressources/outputs/DataMatrix.txt"
  varmd.out <- "test/ressources/outputs/VariableMetadata.txt"
  Perfwhm <- 0.6
  #Cor_eic_th <- 0.75
  Polarity <- "positive"
  intval <- "maxo"
}

extractionTemp <- function(rdatafile,nbSamp,Perfwhm,Polarity,varmd.out,dm.out,intval="maxo"){
# rdatafile = chemin et nom du Rdata
# nbSamp = nb de samples (contient pool & sample & blank & tout autre type d'echantillon)
# varmd.out = chemin et nom du fichier de sortie des variable metadata
# dm.out = chemin et nom du fichier de sortie de la data matrix

library(xcms)
library(CAMERA)

load(rdatafile)

xsetP <- xset

## puis CAMERA
an <-xsAnnotate(xsetP)
an <-groupFWHM(an, perfwhm = Perfwhm)
an <-findIsotopes(an, mzabs = 0.01)
#an <-groupCorr(an, cor_eic_th = Cor_eic_th)
anP <-findAdducts(an, polarity=Polarity)

thelist <- getPeaklist(anP)

# Addition of RT in minutes and Variable ID 
thelist <- thelist[,c(4,1,4,(2:ncol(thelist)))]
colnames(thelist)[c(1,3,6)] <- c("ions","rt","rtsec")
thelist[,3] <- round(thelist[,6]/60,digits=1)
if(Polarity=="positive"){imode <- "p"}else{imode <- "n"} 
thelist[,1] <- paste(imode,round(thelist[,2],digits=2),"T",thelist[,3],sep="")

fin.var <- ncol(thelist) - 3 - nbSamp

# Etape de groupval
gint=groupval(xsetP,method="medret",value=intval,intensity=intval)
thelist[,(fin.var+1):(ncol(thelist)-3)] <- gint

varmd.tb <- data.frame(thelist[,1:fin.var],thelist[,(ncol(thelist)-2):(ncol(thelist))])
write.table(varmd.tb, file=varmd.out,sep="\t", row.names=F, quote=FALSE)

dm.tb <- data.frame(thelist[,1,drop=FALSE],thelist[,(fin.var+1):(ncol(thelist)-3)])
write.table(dm.tb, file=dm.out,sep="\t", row.names=F, quote=FALSE)


}

# Typical function call
#extractionTemp(rdatafile,nbSamp,Perfwhm,Polarity,varmd.out,dm.out)
