## ****** Norma_VdK_Lowess environnemnt : ****** ##
# version avr 2014 M Landi / F Giacomoni

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
No dependency with pfem lib
--

## --- R bin and Packages : --- ##
$ R --version
R version 3.0.1 (2013-05-16) -- "Good Sport"
Platform: x86_64-redhat-linux-gnu (64-bit)

This script has two R dependencies availables in the "Tool Dependency Packages" category.
The PFEM_R package is available on the ABIMS toolshed : /Tool Dependency Packages/pfem_r
This package can be deployed in /usr/local/share/R

The Norma_VdK_Lowess tool need in the PFEM_R "tool" only the following files :
normalization/Normalisation_QCpool.r
toolbox/toolBox.R
-- 

## --- Binary dependencies --- ##
No external binary.
--

## --- Config : --- ##
Edit the config file : ~/Norma_VdK_Lowess/normalization.cfg
R_BIN=/your/R/bin
R_NORMA_TB=/your/PFEM_R/package/path/normalization/Normalisation_QCpool.r
R_PFEM_TB=/your/PFEM_R/package/path/toolbox/toolBox.R
--

## --- XML HELP PART --- ##
Copy the following images in ~/static/images/metabolomics/Workflow_position :
Normalization_Van-der-Kloet.png
--

## --- DATASETS --- ##
No data set ! waiting for galaxy pages
--

## --- ??? COMMENTS ??? --- ##
--