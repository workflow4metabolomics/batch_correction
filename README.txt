## ****** Norma_VdK_Lowess + Determine_VdK environnemnt : ****** ##
# version jul 2014 M Landi / F Giacomoni

## --- PERL compilator / libraries : --- ##
$ perl -v
This is perl, v5.10.1 (*) built for x86_64-linux-thread-multi

# libs CORE PERL : 
use strict ;
use warnings ;
use Carp qw (cluck croak carp) ;
use Data::Dumper ;
use Getopt::Long ;
use FindBin ;

# libs CPAN PERL : 
$ perl -e 'use Statistics::R'
$ sudo perl -MCPAN -e shell
cpan> install Statistics::R

# libs pfem PERL : 
NA
--

## --- R bin and Packages : --- ##
$ R --version
R version 3.0.1 (2013-05-16) -- "Good Sport"
Platform: x86_64-redhat-linux-gnu (64-bit)

The dependent libs are :
>install.packages("batch", dep=TRUE)
-- 

## --- Binary dependencies --- ##
NA
--

## --- Config : --- ##
NA
--

## --- XML HELP PART --- ##
two tools - two images :
Normalization_Van-der-Kloet.png
Determine_Van-der-Kloet.png
--

## --- DATASETS --- ##
No data set ! waiting for galaxy pages
--

## --- ??? COMMENTS ??? --- ##
NA
--
