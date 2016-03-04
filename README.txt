## ****** Determine_bc + batch_correction environnemnt : ****** ##
# version March 2016 M Landi / F Giacomoni / M Petera / E Thevenot

## --- PERL compilator / libraries : --- ##
NA
--

## --- R bin and Packages : --- ##
$ R --version
R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
Platform: x86_64-redhat-linux-gnu (64-bit)

The dependent libs are :
>install.packages("batch", dep=TRUE)
>install.packages("ade4", dep=TRUE)

>source("http://www.bioconductor.org/biocLite.R")
>biocLite("pcaMethods")
>biocLite("ropls")
-- 

## --- Binary dependencies --- ##
NA
--

## --- Config : --- ##
NA
--

## --- XML HELP PART --- ##
two tools - two images :
batch_correction.png
determine_batch_correction.png
--

## --- DATASETS --- ##
No data set ! waiting for galaxy pages
--

## --- ??? COMMENTS ??? --- ##
!WARNING! : Two tools in the same directory !
--

## --- Changelog/News --- ##
Version 2.0.3:
- Update: plot layout improvement
- Update: suppression of combined table for simca
- Improvement: addition of a new parameter "Null values"
--