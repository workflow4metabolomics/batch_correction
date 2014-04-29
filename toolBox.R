# Author: jfmartin
# Modified by : mpetera
###############################################################################
# Version 1.01 "NA"-coding possibility added in acplight function
###############################################################################



###### Fonction qui plotte les meilleures correlations entre les var de 2 matrices A et B
# Elimine les data manquantes par couples
plotrsup = function(A, B, rmin)
{
        for(i in 1:dim(A)[2]) {
                if(is.numeric(A[,i]))
                        for(j in 1:dim(B)[2]) {
                                if(is.numeric(B[,j])) 
                                {
                                  AnNa=A[!is.na(A[,i]) & !is.na(B[,j]),i]
                                  BnNa=B[!is.na(A[,i]) & !is.na(B[,j]),j]
                                  r = cor(AnNa,BnNa)
                                  if(abs(r) >= rmin && r<0.99999)
                                  {
                                    plot(AnNa, BnNa, xlab = dimnames(A)[[2]][i], ylab = dimnames(B)[[2]][j])
                                    abline(lsfit(AnNa, BnNa))
                                    title(paste("r=", as.character(round(r, 2))))
                                  }
                                }
                        }
        }
}

##### Plot A,B avec symboles differents en fonction du facteur s
plotby = function(s, A, B, lx ="", ly ="", reg=T ,regby=F, leg=T ,rmin)
{
if(!is.factor(s))  stop("Le premier argument doit etre un FACTEUR")
if(!is.numeric(A)) stop("Le deuxieme argument doit etre NUMERIQUE")
if(!is.numeric(B)) stop("Le troisieme argument doit etre NUMERIQUE")

   AnNa=A[!is.na(A)& !is.na(B)]
   BnNa=B[!is.na(A)& !is.na(B)]
   snNa=s[!is.na(A)& !is.na(B)]
   D <- data.frame(f=snNa, x=AnNa, y=BnNa)
   r <- cor (D$x, y=D$y)
   if(abs(r) >= rmin && r<0.99999) {

      plot (D$x , D$y , type="n", xlab=lx , ylab=ly )
      if ( reg==T ) abline(lsfit(D$x , D$y), lty=1)
      title (paste("r=", as.character(round(r, 2))))
      nl <- length(levels(D$f))
      cl=c(1,17,16,2,3,4); sy=c(1,16,17,3,5,12)
      for (i in 1:nl) {
        sD <- D[ D$f==levels(D$f)[i],]
        points (sD$x , sD$y, pch=sy[i], col=cl[i])
        if ( regby==T ) abline(lsfit( sD$x , sD$y), lty=i+1)
      }
      if (leg==T) {
         py <- max(D$y)
         if (r >= 0) px <- min(D$x)
         if (r <  0) px <- min(D$x)+(max(D$x)-min(D$x))/2
         legend(px, py , levels(D$f),pch=sy[1:nl])
      }
   }
}


##### Plot les boxplot des var d un dataframe qui comporte des obs a +/- nsd ecartypes
plot_outlier = function (ids, by, nsd)
{
for (i in 1:dim(ids)[2])
  if (is.numeric(ids[,i])) 
  {
    mini=min(ids[!is.na(ids[,i]),i])
    maxi=max(ids[!is.na(ids[,i]),i])
    moy=mean(ids[!is.na(ids[,i]),i])
    ect=sd(ids[!is.na(ids[,i]),i])  
    if (mini < moy- nsd*ect | maxi > moy+nsd*ect) boxplot (ids[,i]~ by, ylab=dimnames(ids)[[2]][i]) 
  }
}

library(ade4); library(pcaMethods)

acpEllipse=function(ids,fact,firstvar,lastvar=dim(ids)[2],scaling="uv",meth="svd",plotloading=T,tc=0.66) {
	# fait une ACP sur ids sachant que la colonne 1 contient l'identificateur d'individu
	# tc  est la taille du caractere utilise pour les plots
	# la colonne 2 contient le facteur d?finissant la couleur des individus
	# Les colonnes suivantes 3:dim(ids)[2] sont les donnees quantitatives utilisees pour l'ACP
	# L'ACP est realise avec le scaling "uv","none","pareto"
	# L'Algo d'ACP est : nipals (si donnees manquantes) ou svdImpute ou Bayesian ou svd ou probalistic
	# Le loading de l'ACP est plotte si plotloading==T (valeur par defaut)
	classe=factor(ids[[fact]])
	idions=dimnames(ids[,firstvar:lastvar])[2]
	colour=1:length(levels(classe))
	ions=as.matrix(ids[,c(firstvar:lastvar)])
	# choix du scaling : "uv","none","pareto"
	object=prep(ions, scale=scaling, center=TRUE)
	# ALGO: nipals,svdImpute, Bayesian, svd, probalistic=F
	result <- pca(object, center=F, method=meth, nPcs=2)
	# ADE4 : representation des ellipsoides des individus de chaque classe
	s.class(result@scores, classe, cpoint = 1,xax=1,yax=2,col=colour,sub=sprintf("Scores - PCs %sx%s",1,2), possub="bottomright")
	if (plotloading==T) 
		{s.label(result@loadings,cpoint = 0,boxes=F,clabel=tc, xax=1,yax=2,sub="Loadings",possub="bottomright")}
}

acplight <- function(ids, scaling="uv", indiv=FALSE,indcol=NULL) {
  # fait une ACP sur ids sachant que la colonne 1 contient l'identificateur d'individu
  # la colonne 2:nf contient les facteurs definissant la couleur des individus
  for (i in 1:3) {
    idss=ids[which(ids[,i+1]!="NA"),]
    idss=data.frame(idss[idss[,i+1]!="",])
    classe=as.factor(idss[[i+1]])
    idsample=as.character(idss[[1]])
    colour=1:length(levels(classe))
    ions=as.matrix(idss[,5:dim(idss)[2]])
    # choix du scaling : "uv","none","pareto"
    object=prep(ions, scale=scaling, center=TRUE)
    # ALGO: nipals,svdImpute, Bayesian, svd, probalistic=F
    result <- pca(object, center=F, method="svd", nPcs=2)
    # ADE4 : representation des ellipsoides des individus de chaque classe
    s.class(result@scores, classe, cpoint = 1,xax=1,yax=2,col=colour,sub=sprintf("Scores - PCs %sx%s",1,2), possub="bottomright")
    #s.label(result@loadings,label = ions, cpoint = 0, clabel=0.4, xax=1,yax=2,sub="Loadings",possub="bottomright")
    if(i==1){resulti <- result}
  }
  if(indiv) {
    colour <- rep("darkblue",length(resulti@loadings)) ; if(!is.null(indcol)) {colour[-c(indcol)] <- "red"}
    plot(resulti@loadings,col=colour,main="Loadings",xaxt="n",yaxt="n",pch=20)
    abline(h=0,v=0)}
}

# Programme de plot des ions avec mis en evidence des facteurs

library(gplots)

plotInterCI = function (ids,indf1,indf2,firstvar) {
# Function qui plotte une variable (ions) avec 2 facteurs crois?s. Si cinetique l'utiliser en facteur 2 
# En entree, ids input dataframe, indf1 et indf2 indice dans la dataframe des 2 facteurs
# ids=x;i=25;indf1=6;indf2=10;indSubject=4
lastvar=dim(ids)[2]
for (i in firstvar:lastvar) { 

	par(las=2)
	xp=ids[!is.na(ids[[i]]),c(indf1,indf2,i)];
	xp=data.frame(xp,interaction(xp[[1]],xp[[2]],sep="."));
	dimnames(xp)[[2]][4]="Inter";
	levInter=c(1:length(levels(xp[[4]])));nl1=length(levels(xp[[1]]));nl2=length(levels(xp[[2]]));
	levInter=array(data=levInter,dim=c(nl1,nl2))
	connection=vector("list",nl1)
	for (i in 1:nl1) {connection[[i]]=levInter[i,]}
	plotmeans(xp[[3]]~ xp[[4]],
			connect=connection,
			ccol=c(1:nl1), barwidth=1.5, pch=16,
			main=dimnames(xp)[[2]][3],xlab=" ",	
			ylab="Intensity"
			);
#	plot.design(xp[[3]]~xp[[1]]+xp[[2]]+xp[[4]],fun=mean,main=dimnames(xp)[[2]][3],xlab=" ",ylab=" ");
# 	boxplot(xp[[3]]~xp[[4]],ylab="Intensity");
}
}


BoxHistDesign = function (ids, fact, firstvar,lastvar=dim(ids)[2]) {
# ids input dataset; fact type vecteur contenant les indices des facteurs etudies; firstvar 1ere var quantitative
# Function qui pour une dataset ids plotte un histogramme de chaque var quantitative AVEC 2 FACTEURS Minimum
# plotte un boxplot by facteur fact[i]
# plotte le plo.design des facteurs
# S'il y a une cinetique, mettre l'indice en premier dans le vecteur d'indice "fact"
#	lastvar=dim(ids)[2]
	nbfact=length(fact); nbinter=0.5*(nbfact^2)-0.5*nbfact; gco=c(1,2,3,3,4);gli=c(3,3,4,6,6)
	for (i in firstvar:lastvar) {
		if (nbfact < 6) { par(mfrow=c(gli[nbfact],gco[nbfact]))} else { par(mfrow=c(nbfact-1,nbinter-1))}
		# cat("var ",i," -> ",names(ids)[i],"\n")
		# eliminantion des donn?es manquantes 
		ods=ids[!is.na(ids[[i]]),]; lab=names(ods)[i]
		xpf=ods[,fact]
		xpv=ods[,i]
#		xpi=array(data="inter",nrow=dim(ods)[1],ncol=nbinter)	
		xb=data.frame(xpf,xpv); indvar=dim(xb)[2];lab=names(ods)[i];names(xb)[indvar]=lab
		hist(xb[[indvar]],main=lab,cex=2.5,xlab=" ",ylab=" ")
		plot.design(xb,fun=mean,main=" ",xlab=" ",ylab=" ")
		for (b in 1:nbfact) {
			boxplot(xb[[indvar]]~xb[[b]],ylab="Intensity")
		}
		## Plot des interactions en fonction du nombre de facteurs
		if (nbfact>1) {
			y=0
			for (j in 1:(nbfact-1)){ for (u in (j+1):nbfact){
				y=y+1
				xpi=interaction(ods[[fact[j]]],ods[[fact[u]]],sep="*")	
				xp=data.frame(xpf,xpi,xpv);
				indinter=length(fact)+1; names(xp)[indinter]="interaction";		
				indvar=dim(xp)[2]; names(xp)[indvar]=names(ods)[i];
				indinter=length(fact)+1; names(xp)[indinter]="interaction";
				boxplot(xp[[indvar]]~xp[[indinter]],ylab="Intensity")
				interaction.plot(xp[[j]],xp[[u]],xp[[indvar]])	
			}}
		}
	}

}

##### Fonction best_Q2 (mixOmics) qui recherche pour la meilleure combinaison de nombre de composante (<=maxcomp)
##### et de loading (vecteur spk) pour avoir le meilleur Q2 pour une spls (sparse PLS mixOmics) avec CVM cross validation
best_Q2 = function (labvar,X,Y,maxcomp,spk,CVM) {

	i=0 ; nr=maxcomp*length(spk)
	resu=matrix(NA,nr,maxcomp+3)
	for (nbc in 1:maxcomp) {
     	for (k in 1:length(spk)) {

	i=i+1
        resu[i,1]=labvar
	resu[i,2]=nbc
	resu[i,3]=spk[k]
cat("i=",i,"Nb comp = ",nbc," Nb loadings kept= ",spk[k]," -> ")
	mtprf.valid=valid(X,Y,
                    ncomp=nbc,
                    keepX=rep(spk[k],nbc),
            #        scaleY=T,
                    mode="regression",
            #        criterion="q2",
                    method="spls",max.iter=500,
                    validation="Mfold",
                    M=CVM)
	for (q in 1:maxcomp) {resu[i,3+q]=mtprf.valid$Q2$q2[q]}
cat(" Q2=", mtprf.valid$Q2$q2,"\n")
      }
      }
resu
}

### exemples de source R divers et varies utilisation de prcomp et pca
bidon = function () {

BoxHistDesign(Q,1,2)

par(mfrow=c(1,2),ask=T)
pdf("ACPconso.pdf",width=27,height=20)
acpEllipse(Q,1,2,tc=1.5)

setwd("C:/Users/jfmartin/Documents/PROJETS_METABO/ANR/PR0092_Phenomenep/ASM0136_NM_jul11_UrinePos/Bioinfo/xcms2/Norm_VDK/ResGLM3wNormConsoSexBatchBH/ANOVAbyAliment")
Qi=read.table(file="UrinePos_932mz_filtrePresenceAvecConso.txt", header=T, sep="\t")
sQi=Qi[,c(4,993:1014)]
acpEllipse(sQi,1,2,tc=1.5)
dev.off()

result <- pca(Qi[,c(993:1014)], center=T, scale="uv", method="svd", nPcs=5)
plot(result)
# definition d'un vecteur "color" issus d'un facteur de dataframe pour definir la couleur des labels en fonction des niveaux du facteur
qcol=as.character(Qi$conso)
qcol[which(qcol=="Q1")]="red"
qcol[which(qcol=="Q4")]="blue"
slplot(result,pcs=c(1,2), sl=Qi$conso,lcex=0.66,scex=0.66,scol=qcol)

rpca <- prcomp(Qi[,c(993:1014)],scale=T,nPcs=5)
plot(rpca)

pdf("bilanConsoCorrel.pdf",width=27,height=20)
corrgram(x=Qcons,lower.panel=panel.pie, upper.panel=panel.pts,order=T,
		text.panel=panel.txt, 
		diag.panel=panel.minmax  )
dev.off()



}

