Signal drift and batch-effect correction  
========================================

A Galaxy module from the [Workflow4metabolomics](http://workflow4metabolomics.org) infrastructure  

Status: [![Build Status](https://travis-ci.org/workflow4metabolomics/batchcorrection.svg?branch=master)](https://travis-ci.org/workflow4metabolomics/batchcorrection).

### Description

**Version:** 2.1.0  
**Date:** 2017-03-20  
**Author:** Jean-François Martin (INRA, AXIOM), Mélanie Pétéra (INRA, PFEM), Marion Landi (INRA, PFEM), Franck Giacomoni (INRA, PFEM), and Etienne A. Thévenot (CEA, LIST)  
**Email:** [jean-francois.martin(at)toulouse.inra.fr](mailto:jean-francois.martin@toulouse.inra.frr), [melanie.petera(at)clermont.inra.fr](mailto:melanie.petera@clermont.inra.fr), [etienne.thevenot(at)cea.fr](mailto:etienne.thevenot@cea.fr)  
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
    + `static/images/pdf_plotsituation.png`    
    + `static/images/Vdk_pdf1.png`    
    + `static/images/Vdk_pdf2.png`        
* Wrapper files:
    + `batch_correction_wrapper.R`  
    + `batch_correction_all_loess_wrapper.R`  
* Script files:
    + `Normalisation_QCpool.r`  
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

The code in the batch_correction_wrapper can be tested by running the `runit/batchcorrection_runtests.R` R file  

You will need to install **RUnit** package in order to make it run:
```r
install.packages('RUnit', dependencies = TRUE)
```

### Working example  

See the reference histories [W4M00001_Sacurine-statistics; DOI:10.15454/1.4811121736910142E12](http://dx.doi.org/10.15454/1.4811121736910142E12) and [W4M00002_Sacurine-comprehensive; DOI:10.15454/1.481114233733302E12](http://dx.doi.org/10.15454/1.481114233733302E12)  
 

### News

##### CHANGES IN VERSION 2.1.0  

INTERNAL MODIFICATIONS  

 * For PCA figure display only ("all_loess_" options): missing values are set to the minimum value before PCA computation is performed (with svd)
 
 * Additional running and installation tests added with planemo, conda, and travis  

 * Modification of the 'all_loess_wrapper.R' file to handle 'ropls' package versions of 1.3.15 and above (i.e. after switching to S4 classes) 