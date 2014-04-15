#! perl

use strict ;
use Carp ;
use warnings ;

use Getopt::Long ;
use Statistics::R ;

use FindBin ; ## Permet de localisez le repertoire du script perl d'origine

###########################################################################################################
### Determine_VdK.pl
### Author : Franck GIACOMONI, Jean-Francois MARTIN and Marion LANDI
### Email : fgiacomoni\@clermont.inra.fr (for perl) or jean-francois.martin\@clermont.inra.fr (for R)
### Version : 1.0
### Created : 26/09/2013
###########################################################################################################

my ( $help, $samplefile, $ionfile, $ref_factor, $out_graph_pdf, $out_preNormSummary ) = ( undef, undef, undef, undef, undef, undef, undef, undef );
my ( $module_Norma , $module_pfemR, $R_bin, $Confs) = ( undef, undef, undef, undef ) ;

my ( $sep1, $sep2 ) = ( "tab", "tab" ) ;

my $Path = $FindBin::Bin ;

#commande des options gere en ligne de commande
GetOptions(
        "help|h:s"		=>\$help,  			## only help !
        "samplefile:s"	=>\$samplefile,		## input data matrix with samples in rows and ions in columns
        "ionfile:s"		=>\$ionfile,	 	## input with ions data (extract from xcms)
        "sep1:s"		=>\$sep1,			## caractere sep in input samplefile
        "sep2:s"		=>\$sep2,			## caractere sep in input ionfile
        "ref_factor:s"	=>\$ref_factor,		## Ref factor for plotting
        "out_graph_pdf:s"=>\$out_graph_pdf,	## plot output
		"out_preNormSummary:s"=>\$out_preNormSummary	## tabular bilan
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
	if 		($sep1 eq "tabulation"){  	$sep1 = '\t'; }
	elsif 	($sep1 eq "tabular"){	$sep1 = '\t';	 }
	elsif 	($sep1 eq "space"){	$sep1 = ' ';	 }
	if 		($sep2 eq "tabulation"){  	$sep2 = '\t'; }
	elsif 	($sep2 eq "tabular"){	$sep2 = '\t';	 }
	elsif 	($sep2 eq "space"){	$sep2 = ' ';	 }
	
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
	
	### Lecture du fichier de données format xcms (ions en lignes)
	$R->send(qq`infic="$ionfile"`) ;
	$R->send(qq`iddata=read.table(infic,header=T,sep="$sep2")`) ;
	$R->send(qq`dimnames(iddata)[[1]]=iddata[[1]]`) ;
	
	### Transposition des donnees ions
	$R->send(qq`idTdata=t(iddata[,2:dim(iddata)[2]])`) ;
	$R->send(qq`idTdata=data.frame(dimnames(idTdata)[[1]],idTdata)`) ;
	
	### Merge des 2 fichiers (ok même si les 2 dataframe ne sont pas triés sur la même clef
	$R->send(qq`id=merge(idsample, idTdata, by.x=1, by.y=1)`) ;
	$R->send(qq`id\$batch=as.factor(id\$batch)`) ;
	
	$R->send(qq`ids=id[id\$typsample == "p" | id\$typsample == "s",]`) ;
	$R->send(qq`nbid=dim(idsample)[2]`) ;
	
	### Facteur biologique 
	$R->send(qq`factbio="$ref_factor"`) ;
	
	## outputs
	$R->send(qq`out_graph_pdf="$out_graph_pdf"`) ;
	$R->send(qq`out_preNormSummary="$out_preNormSummary"`) ;
	
	## Launch
	$R->send(qq`plotsituation(ids,nbid,out_graph_pdf,out_preNormSummary,factbio)`) ;
	### TODO
	
	# Fermeture du pont
	$R->stopR();
	
}
else {
	die "Absence of needed R modules ($module_Norma, $module_pfemR)\n";
}



=head2 METHOD as_conf

	## Description : permet de créer l'object conf à partir d'un fichier de conf de type KEY=VALUE
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
	
	# Determine_VdK.pl
	# Input : 
	# Author : Franck GIACOMONI, Jean-Francois MARTIN and Marion LANDI
	# Email : fgiacomoni\@clermont.inra.fr or marion.landi\@clermont.inra.fr (for perl) or jean-francois.martin\@clermont.inra.fr (for R)
	# Version : 1.0
	# Created : 26/09/2013 -- release XX/XX/20XX
	USAGE :
	        Determine_VdK.pl -help
	        example : Determine_VdK.pl -samplefile MetadataSample_phenomenep_xcmsMF.txt -sep1 tabULATION -ionfile phenomenep_xcmsMF.txt -sep2 tabulation -ref_factor batch -out_graph_pdf out_graph_pdf.pdf -out_preNormSummary out_preNormSummary.txt 
	!!! ATTENTION !!! ce logiciel necessite :
		# fonction de normalisation necessite en entree les 2 fichiers : ion frame avec les valeurs des ions et aussi samples frame avec
		# 3 colonnes : 
		# 	'batch' pour identifier les series d'analyses (entre 2 calibrations ou lavage de source) Il  faut at least 2 pools par batch
		#	'injection' type integer contenant les ordres d'injection de tous les echantillons 
		# 	'typsample' contenant le type d'echantillon: 'p' pour pool ou 's' pour sample
		# les signaux (ions) doivent etre dans le fichier ionfile
		\n";
		
    exit(1);
}

