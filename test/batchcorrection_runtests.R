#!/usr/bin/env Rscript

## Package
##--------

library(RUnit)

## Constants
##----------

testOutDirC <- "output"
argVc <- commandArgs(trailingOnly = FALSE)
scriptPathC <- sub("--file=", "", argVc[grep("--file=", argVc)])


## Functions
##-----------

## Reading tables (matrix or data frame)
readTableF <- function(fileC, typeC = c("matrix", "dataframe")[1]) {

    	file.exists(fileC) || stop(paste0("No output file \"", fileC ,"\"."))

        switch(typeC,
               matrix = return(t(as.matrix(read.table(file = fileC,
                   header = TRUE,
                   row.names = 1,
                   sep = "\t",
                   stringsAsFactors = FALSE)))),
               dataframe = return(read.table(file = fileC,
                   header = TRUE,
                   row.names = 1,
                   sep = "\t",
                   stringsAsFactors = FALSE)))

}

## Call wrapper
wrapperCallF <- function(paramLs, allLoessL) {

    ## Set program path
    wrapperPathC <- file.path(dirname(scriptPathC), "..",
                              ifelse(allLoessL,
                                     "batch_correction_all_loess_wrapper.R",
                                     "batch_correction_wrapper.R"))

    ## Set arguments
    argLs <- NULL
    for (parC in names(paramLs))
        argLs <- c(argLs, parC, paramLs[[parC]])

    ## Call
    wrapperCallC <- paste(c(wrapperPathC, argLs), collapse = " ")

    if(.Platform$OS.type == "windows")
        wrapperCallC <- paste("Rscript", wrapperCallC)


    print(wrapperCallC)
    

    wrapperCodeN <- system(wrapperCallC)

    if (wrapperCodeN != 0)
        stop(paste0("Error when running 'batch_correction_",
                    ifelse(allLoessL, "all_loess_", ""),
                    "wrapper.R'"))

    ## Get output
    outLs <- list()
    if ("dataMatrix_out" %in% names(paramLs))
        outLs[["datMN"]] <- readTableF(paramLs[["dataMatrix_out"]], "matrix")
    if ("sampleMetadata_out" %in% names(paramLs))
        outLs[["samDF"]] <- readTableF(paramLs[["sampleMetadata_out"]], "dataframe")
    if ("variableMetadata_out" %in% names(paramLs))
        outLs[["varDF"]] <- readTableF(paramLs[["variableMetadata_out"]], "dataframe")
    if("information" %in% names(paramLs))
        outLs[["infVc"]] <- readLines(paramLs[["information"]])

    if("out_preNormSummary" %in% names(paramLs))
        outLs[["sumDF"]] <- readTableF(paramLs[["out_preNormSummary"]], "dataframe")

    return(outLs)
    
}

## Setting default parameters
defaultArgF <- function(testInDirC, determineL) {

    defaultArgLs <- list()

    if(file.exists(file.path(dirname(scriptPathC), testInDirC, "dataMatrix.tsv")))
        defaultArgLs[["dataMatrix"]] <- file.path(dirname(scriptPathC), testInDirC, "dataMatrix.tsv")
    if(file.exists(file.path(dirname(scriptPathC), testInDirC, "sampleMetadata.tsv")))
        defaultArgLs[["sampleMetadata"]] <- file.path(dirname(scriptPathC), testInDirC, "sampleMetadata.tsv")

    if(!determineL)
        if(file.exists(file.path(dirname(scriptPathC), testInDirC, "variableMetadata.tsv")))
            defaultArgLs[["variableMetadata"]] <- file.path(dirname(scriptPathC), testInDirC, "variableMetadata.tsv")

    if(determineL) { ## determinebc

        defaultArgLs[["out_graph_pdf"]] <- file.path(dirname(scriptPathC), testOutDirC, "out_graph.pdf")
        defaultArgLs[["out_preNormSummary"]] <- file.path(dirname(scriptPathC), testOutDirC, "preNormSummary.txt")

    } else { ## batchcorrection
        
        defaultArgLs[["dataMatrix_out"]] <- file.path(dirname(scriptPathC), testOutDirC, "dataMatrix.tsv")
        defaultArgLs[["variableMetadata_out"]] <- file.path(dirname(scriptPathC), testOutDirC, "variableMetadata.tsv")
        defaultArgLs[["graph_output"]] <- file.path(dirname(scriptPathC), testOutDirC, "graph_output.pdf")
        defaultArgLs[["rdata_output"]] <- file.path(dirname(scriptPathC), testOutDirC, "rdata_output.rdata")
        
       }

    defaultArgLs

}

## Main
##-----

## Create output folder
file.exists(testOutDirC) || dir.create(testOutDirC)

## Run tests
test.suite <- defineTestSuite('tests', dirname(scriptPathC), testFileRegexp = paste0('^.*_tests\\.R$'), testFuncRegexp = '^.*$')
isValidTestSuite(test.suite)
test.results <- runTestSuite(test.suite)
print(test.results)

