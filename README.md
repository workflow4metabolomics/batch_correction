Signal drift and batch-effect correction  
========================================

A Galaxy module from the [Workflow4metabolomics](http://workflow4metabolomics.org) infrastructure  

Status: [![Build Status](https://travis-ci.org/workflow4metabolomics/batch_correction.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/batch_correction).

### Description

**Version:** 3.0.0  
**Date:** 2020-10-02    
**Author:** Jean-Francois Martin (INRAE, AXIOM), Melanie Petera (INRAE, PFEM), Marion Landi (PFEM), Franck Giacomoni (INRAE, PFEM), and Etienne A. Thevenot (CEA, LIST)  
**Email:** [jean-francois.martin(at)inrae.fr](mailto:jean-francois.martin@inrae.fr), [melanie.petera(at)inrae.fr](mailto:melanie.petera@inrae.fr), [etienne.thevenot(at)cea.fr](mailto:etienne.thevenot@cea.fr)  
**Citation:**  
**Licence:** CeCILL  
**Reference history:** [W4M00001_Sacurine-statistics; DOI:10.15454/1.4811121736910142E12](http://dx.doi.org/10.15454/1.4811121736910142E12) and [W4M00002_Sacurine-comprehensive; DOI:10.15454/1.481114233733302E12](http://dx.doi.org/10.15454/1.481114233733302E12)  
**Funding:** Agence Nationale de la Recherche ([MetaboHUB](http://www.metabohub.fr/index.php?lang=en&Itemid=473) national infrastructure for metabolomics and fluxomics, ANR-11-INBS-0010 grant)  

### Installation

* Configuration files:
    + `batch_correction.xml` (for the "Batch Correction" module)
    + `determine_bc.xml` (for the "Determine Batch Correction" module)  
* Image files: 
    + `static/images/batch_correction.png`    
    + `static/images/determine_batch_correction.png`          
* Wrapper files:
    + `batch_correction_3Lwrapper.R`  
    + `batch_correction_all_loess_wrapper.R`  
* Script files:
    + `batch_correction_3Lfct.R`  
    + `batch_correction_3Llauncher.R`
    + `batch_correction_all_loess_script.R`  
* R packages
  + **batch** from CRAN  
  
    ```r
    install.packages("batch", dep=TRUE)  
    ```

  + **ade4** from CRAN  
  
    ```r
    install.packages("ade4", dep=TRUE)  
    ```

 + **pcaMethods** from Bioconductor  
  
    ```r
    source("http://www.bioconductor.org/biocLite.R")  
    biocLite("pcaMethods")      
    ```

 + **ropls** from Bioconductor  
  
    ```r
    source("http://www.bioconductor.org/biocLite.R")  
    biocLite("ropls")      
    ```

### Tests

Some functional tests are provided in Galaxy format (see XML wrapper and test-data folder).

WIP: The code in the batch_correction_wrapper can be tested by running the `test/batchcorrection_runtests.R` R file  
You will need to install **RUnit** package in order to make it run:
```r
install.packages('RUnit', dependencies = TRUE)
```

### Working example  

See the reference histories [W4M00001_Sacurine-statistics; DOI:10.15454/1.4811121736910142E12](http://dx.doi.org/10.15454/1.4811121736910142E12) and [W4M00002_Sacurine-comprehensive; DOI:10.15454/1.481114233733302E12](http://dx.doi.org/10.15454/1.481114233733302E12)  
 

### News

###### CHANGES IN VERSION 3.0.0  

NEW FEATURES 

 * Specific names for the 'sampleType', 'injectionOrder', and 'batch' from sampleMetadata are now available in a dedicated parameter section
 * Addition of a sum of ions before/after plot for linear/lowess/loess methods
 * Addition of a third option in "Null values" parameter (renamed "unconsistant values") in linear/lowess/loess methods
 * linear/lowess/loess methods now handle NA in intensities and allow "blank" samples in the dataset

INTERNAL MODIFICATIONS

 * XML optimisation using macros
 * Output name changes
 * linear/lowess/loess methods: disabling of RData output
 * linear/lowess/loess methods: split of tool-linked code and script-linked one
 * linear/lowess/loess methods: adjustments in the normalisation process to match matters linked to NA acceptance
 * linear/lowess/loess methods: better handling of special characters in IDs and column names

###### CHANGES IN VERSION 2.2.4  

INTERNAL MODIFICATIONS  

Fixed bug for pool selection ("all_loess" methods)

###### CHANGES IN VERSION 2.2.2  

INTERNAL MODIFICATIONS  

Fixed bug for color plot ("all_loess" methods)  

###### CHANGES IN VERSION 2.2.0  

NEW FEATURE  

Specific names for the 'sampleType', 'injectionOrder', and 'batch' from sampleMetadata can be selected by the user (for compatibility with the MTBLS downloader)  

##### CHANGES IN VERSION 2.1.2  

INTERNAL MODIFICATIONS  

 * Minor modifications in config file  

##### CHANGES IN VERSION 2.1.0  

INTERNAL MODIFICATIONS  

 * For PCA figure display only ("all_loess_" options): missing values are set to the minimum value before PCA computation is performed (with svd)
 
 * Additional running and installation tests added with planemo, conda, and travis  

 * Modification of the 'all_loess_wrapper.R' file to handle 'ropls' package versions of 1.3.15 and above (i.e. after switching to S4 classes) 
