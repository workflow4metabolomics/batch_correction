#! perl

use strict ;
use Carp ;
use warnings ;

use Getopt::Long ;
use Statistics::R ;

use FindBin ; ## Permet de localisez le repertoire du script perl d'origine

###########################################################################################################
### norma_VDK_or_Lowess.pl
### Author : Franck GIACOMONI, Jean-Francois MARTIN and Marion LANDI
### Email : fgiacomoni\@clermont.inra.fr (for perl) or jean-francois.martin\@clermont.inra.fr (for R)
### Version : 1.2
### Created : 27/07/2012
###########################################################################################################

my ( $help, $samplefile, $ionfile, $ionmeta) = ( undef, undef, undef, undef, undef );
my ( $ref_factor, $detail, $method ) = ( undef, undef, undef, undef );
my ( $graph_output, $csv_output_for_simca, $csv_output, $csv_output_metaion, $rdata_output ) = ( undef, undef, undef, undef, undef, undef );
my ( $module_Norma , $module_pfemR, $R_bin, $Confs) = ( undef, undef, undef, undef ) ;

my ( $sep1, $sep2, $sep3) = ( "tab", "tab", "tab" );

my $Path = $FindBin::Bin ;

## commande d'execution :
# $ perl norma_VDK_or_Lowess.pl -samplefile MetadataSample_phenomenep_xcmsMF.txt -sep1 tab -ionfile phenomenep_xcmsMF.txt -sep2 tab -ref_factor conso -detail reg -method linear -graph_output logVDK.pdf -csv_output NormVDK_phenomenep_xcmsMF.txt -rdata_output NormVDK_phenomenep_xcmsMF.Rdata

#commande des options gere en ligne de commande
GetOptions(
        "help|h:s"		=>\$help,  			## only help !
        "samplefile:s"	=>\$samplefile,		## input data matrix with samples in rows and ions in columns
        "ionfile:s"		=>\$ionfile,	 	## input with ions data (extract from xcms)
        "ionmeta:s"		=>\$ionmeta,		## liste and information on IONS (metadata ions)
        "sep1:s"		=>\$sep1,			## caractere sep in input samplefile
        "sep2:s"		=>\$sep2,			## caractere sep in input ionfile
        "sep3:s"		=>\$sep3,			## caractere sep in input ionfile
        "ref_factor:s"	=>\$ref_factor,		## Ref factor for plotting
        "detail:s"		=>\$detail,			## used algo for norma [no, plot,reg]
        "method:s"		=>\$method,			## model used for normalization : lowess (none linear), loess (none linear) or vdk (linear) [lowess, loess, linear]
        "graph_output:s"=>\$graph_output,	## plot output
        "csv_output_for_simca:s"	=>\$csv_output_for_simca,		## csv output
        "csv_output:s"	=>\$csv_output,		## csv output
        "csv_output_metaion:s"	=>\$csv_output_metaion,		## csv output with metadata ions
        "rdata_output:s"=>\$rdata_output		## rdata output
);

## si on met l'option -help ou -h on lance la fonction d'aide
if ( defined($help) ){ help(); }

## get all params for diff envt using:
foreach my $conf ( <$Path/*.cfg> ) {
	$Confs = as_conf($conf) ;
}

if ( ( defined $Confs->{'R_NORMA_TB'} ) and ( defined $Confs->{'R_PFEM_TB'} ) and ( defined $Confs->{'R_BIN'} ) ) {
	$module_Norma = $Confs->{'R_NORMA_TB'} ;
	$module_pfemR = $Confs->{'R_PFEM_TB'} ;
	$R_bin = $Confs->{'R_BIN'} ;
}
else { 	 die "Problem with R envt : $!\n" ; }

if ( ( defined $module_Norma ) and ( -e $module_Norma ) and ( defined $module_pfemR ) and ( -e $module_pfemR ) ) {
	## Envt checking 
	if 		($sep1 eq "tab"){  	$sep1 = '\t'; }
	elsif 	($sep1 eq "tabular"){	$sep1 = '\t';	 }
	elsif 	($sep1 eq "space"){	$sep1 = ' ';	 }
	if 		($sep2 eq "tab"){  	$sep2 = '\t'; }
	elsif 	($sep2 eq "tabular"){	$sep2 = '\t';	 }
	elsif 	($sep2 eq "space"){	$sep2 = ' ';	 }
	if 		($sep3 eq "tab"){  	$sep3 = '\t'; }
	elsif 	($sep3 eq "tabular"){	$sep3 = '\t';	 }
	elsif 	($sep3 eq "space"){	$sep3 = ' ';	 }
	
	# Declaration de l'objet
	my $R = Statistics::R->new(
	    "r_bin"   => $R_bin, 		# path exe R
	    "log_dir" => $FindBin::Bin, # espace de travail, repertoire du pont entre R et Perl
	) or die "Problem with R : $!\n";
	
	# Ouverture du pont
	$R->startR;
	$R->send(qq`source("$module_pfemR")`) ;
	$R->send(qq`source("$module_Norma")`) ;

	## Lecture du fichier metadata sample 
	$R->send(qq`samplefile="$samplefile"`) ;
	$R->send(qq`idsample=read.table(samplefile,header=T,sep="$sep1")`) ;
	$R->send(qq`idsample[[1]]=make.names(idsample[[1]])`) ;
	
	### Lecture du fichier de donn�es format xcms (ions en lignes)
	$R->send(qq`infic="$ionfile"`) ;
	$R->send(qq`iddata=read.table(infic,header=T,sep="$sep2")`) ;
	$R->send(qq`dimnames(iddata)[[1]]=iddata[[1]]`) ;
	
	### Lecture du fichier metadata ions
	$R->send(qq`ionmeta="$ionmeta"`) ;
	$R->send(qq`metaion=read.table(ionmeta,header=T,sep="$sep3")`) ;
	
	### Transposition des donnees ions
	$R->send(qq`idTdata=t(iddata[,2:dim(iddata)[2]])`) ;
	$R->send(qq`idTdata=data.frame(dimnames(idTdata)[[1]],idTdata)`) ;
	
	### Merge des 2 fichiers (ok m�me si les 2 dataframe ne sont pas tri�s sur la m�me clef
	$R->send(qq`id=merge(idsample, idTdata, by.x=1, by.y=1)`) ;
	$R->send(qq`id\$batch=as.factor(id\$batch)`) ;
	
	$R->send(qq`ids=id[id\$typsample == "p" | id\$typsample == "s",]`) ;
	$R->send(qq`nbid=dim(idsample)[2]`) ;
	
	### Facteur biologique 
	$R->send(qq`factbio="$ref_factor"`) ;
	
	## variables
	$R->send(qq`detail="$detail"`) ;
	$R->send(qq`method="$method"`) ;
	
	## outputs
	$R->send(qq`outfic="$csv_output_for_simca"`) ;
	$R->send(qq`outlog="$graph_output"`) ;
	
	## Launch
	$R->send(qq`res = norm_QCpool(ids,nbid,outfic,outlog,factbio,metaion,detail,F,F,method)`);
	$R->send(qq`save(res, file="$rdata_output")`) ;
	$R->send(qq`write.table(res[[1]], file="$csv_output", sep = "\t", row.names=F, quote=F)`) ;
	$R->send(qq`write.table(res[[2]], file="$csv_output_metaion", sep = "\t", row.names=F, quote=F)`) ;
		
	# Fermeture du pont
	$R->stopR();
	
}
else {
	die "Absence of needed R modules ($module_Norma, $module_pfemR)\n";
}



=head2 METHOD as_conf

	## Description : permet de cr�er l'object conf � partir d'un fichier de conf de type KEY=VALUE
	## Input : $file
	## Ouput : $oConf (a hash)
	## Usage : my ( $oConf ) = as_conf( $file ) ;
	
=cut
## START of SUB
sub as_conf {
	## Retrieve Values
    my ( $file, $separator ) = @_ ;
    
    if (!defined $separator) { $separator = ';' } ## set separator to ;
    
    if ( !defined $file )  {  croak "Can't create object with an none defined file\n" ; }
    
    my %Conf = () ; ## Hash devant contenir l'ensemble des parametres locaux
	
	if (-e $file) {
		open (CFG, "<$file") or die "Can't open $file\n" ;
		while (<CFG>) {
			chomp $_ ;
			if ( $_ =~ /^#(.*)/)  {	next ; }
			elsif ($_ =~/^(.*)=(.*)/) {
				
				my ($key, $value) = ($1, $2) ;
				
				if ( $value=~/$separator/ ) { ## is a list to split
					my @tmp = split(/$separator/ , $value) ;
					$Conf{$key} = \@tmp ;
				}
				else {
					$Conf{$key} = $value ;
				}
			}
		}
		close(CFG) ;
	}
	else { 
		croak "Can't create object with an none existing file\n" ;
	}
	
    return ( \%Conf ) ;
}
## END of SUB


#====================================================================================
# Help subroutine called with -h option
# number of arguments : 0
# Argument(s)        :
# Return           : 1
#====================================================================================
sub help {
    print STDERR "
	$0
	
	# norma_VDK_or_Lowess.pl
	# Input : 
	# Author : Franck GIACOMONI, Jean-Francois MARTIN and Marion LANDI
	# Email : fgiacomoni\@clermont.inra.fr or marion.landi\@clermont.inra.fr (for perl) or jean-francois.martin\@clermont.inra.fr (for R)
	# Version : 1.2
	# Created : 27/07/2012 -- release 05/06/2013
	USAGE :
	        norma_VDK_or_Lowess.pl -help
	        example : norma_VDK_or_Lowess.pl -samplefile MetadataSample_phenomenep_xcmsMF.txt -sep1 tab -ionfile phenomenep_xcmsMF.txt -sep2 tab -ref_factor conso -detail reg -method linear -graph_output logVDK.pdf -csv_output NormVDK__phenomenep_xcmsMF.txt -rdata_output NormVDK_phenomenep_xcmsMF.Rdata
	!!! ATTENTION !!! ce logiciel necessite :
		# fonction de normalisation necessite en entree les 2 fichiers : ion frame avec les valeurs des ions et aussi samples frame avec
		# 3 colonnes : 
		# 	'batch' pour identifier les series d'analyses (entre 2 calibrations ou lavage de source) Il  faut at least 2 pools par batch
		#	'injection' type integer contenant les ordres d'injection de tous les echantillons 
		# 	'typsample' contenant le type d'echantillon: 'p' pour pool ou 's' pour sample
		# les signaux (ions) doivent etre dans le fichier ionfile
		# detail : definis le niveau de detail dans les resultats : detail='no' pour juste ACP+histogramme des CV avant et apres
		# normalisation. detail='plot' pour les plots par batch des points toujours avant apres. detail='reg' pour faire apparaitre
		# ATENTION : Pas de donnees manquantes admises pour les intensites des ions !!!!!!!!!!!!!!!!
		# detail='no' : ACP+histogramme des CV avant et apres normalisation
		# detail='plot' : idem 'no' + plots par batch des points avant apres
		# detail='reg' : idem 'plot' + plots des regressions par batch des QC pools et des samples
		\n";
		
    exit(1);
}
