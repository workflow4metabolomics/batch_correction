## Etienne Thevenot
## CEA, MetaboHUB Paris
## etienne.thevenot@cea.fr
## April 8th, 2015

library(batch) ## necessary for parseCommandArgs function
args = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects

source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}

## Import the different functions
source_local("batch_correction_all_loess_script.R")



argLs <- args

##  argLs[["method"]] is either 'all_loess_pool' or 'all_loess_sample'
##  alternative version developped by CEA
##  all variables are treated with loess
##  the reference observations for loess are either 'pool'
## ('all_loess_pool') or 'sample' ('all_loess_sample')

rssVsGalL <- FALSE

if(rssVsGalL) { ## for running with R outside the Galaxy environment during development of the script

    ## 'example' input dir
    exaDirInpC <- "example/input"

    argLs <- list(dataMatrix = file.path(exaDirInpC, "dataMatrix.tsv"), #tab file
                  sampleMetadata = file.path(exaDirInpC, "sampleMetadata.tsv"), #tab file
                  variableMetadata = file.path(exaDirInpC, "variableMetadata.tsv"), # tab file
                  ref_factor = "batch",
                  method = c("all_loess_pool", "all_loess_sample")[1],
                  span = 1)

    ## 'example' output dir
    exaDirOutC <- gsub("input", "output", exaDirInpC)

    argLs <- c(argLs,
               list(dataMatrix_out = file.path(exaDirOutC, paste(argLs[["method"]], "dataMatrix.tsv", sep="_")),

                    variableMetadata_out = file.path(exaDirOutC, paste(argLs[["method"]], "variableMetadata.tsv", sep="_")),
                    variable_for_simca = file.path(exaDirOutC, paste(argLs[["method"]], "combined_results.tsv", sep="_")),
                    graph_output = file.path(exaDirOutC, paste(argLs[["method"]],"graph.pdf", sep="_")),
                    rdata_output = file.path(exaDirOutC, paste(argLs[["method"]], "rdata.rdata", sep="_"))))

    stopifnot(file.exists(exaDirOutC))

}


##------------------------------
## Initializing
##------------------------------

## options
##--------

strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)

## libraries
##----------

library(ropls)

## constants
##----------

modNamC <- "Batch correction" ## module name

## functions
##----------

testF <- function() { ## unit tests of the functions

    options(stringsAsFactors = FALSE)
    epsN <- .Machine[["double.eps"]]

    filDatC.in <- "example/input/dataMatrix.tsv"
    filSamC.in <- "example/input/sampleMetadata.tsv"

    refC <- "pool"

    rawMN <- t(as.matrix(read.table(filDatC.in, header = TRUE, row.names = 1, sep = "\t")))
    samDF <- read.table(filSamC.in, header = TRUE, row.names = 1, sep = "\t")
    samDF[, "ordIniVi"] <- 1:nrow(rawMN)

    ordBatInjVi <- order(samDF[, "batch"], samDF[, "injectionOrder"])
    rawMN <- rawMN[ordBatInjVi, ]
    samDF <- samDF[ordBatInjVi, ]

    ## loessF

    cat("\n'loessF'...")

    sumVn <- rowSums(rawMN, na.rm = TRUE)

    pooVi <- grep("pool", samDF[, "sampleType"])
    samVi <- grep("sample", samDF[, "sampleType"])

    loeVi <- loessF(sumVn, pooVi, samVi, spnN=spnN)

    stopifnot(abs(loeVi[1] - 824401480) < 1e-1)

    cat("OK")

    ## plotBatchF

    cat("\n'plotBatchF'...")

    resLs <- plotBatchF(rawMN, samDF)

    stopifnot(abs(resLs[["sumVn"]][1] - 804524011) < epsN)
    stopifnot(abs(resLs[["tcsMN"]][1, 1] - 0.1835649) < 4.8e-08)

    cat("OK")

    ## shiftBatchCorrectF

    cat("\n'shiftBatchCorrectF'...")

    nrmMN <- shiftBatchCorrectF(rawMN, samDF, refC)

    stopifnot(abs(nrmMN[1, 1] - 17737142) < 2.7e-1)

    cat("OK")

    ## end

    cat("\nEnd\n")

} ## testF


## log file
##---------

## sink(argLs[["information"]]) ## not implemented

cat("\nStart of the '", modNamC, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")

## loading
##--------

rawMN <- t(as.matrix(read.table(argLs[["dataMatrix"]],
                                header = TRUE,
                                row.names = 1,
                                sep = "\t")))

samDF <- read.table(argLs[["sampleMetadata"]],
                    header = TRUE,
                    row.names = 1,
                    sep = "\t")

varDF <- read.table(argLs[["variableMetadata"]],
                    check.names = FALSE,
                    header = TRUE,
                    row.names = 1,
                    sep = "\t") ## not used; for compatibility only

refC <- tolower(gsub("all_loess_", "", argLs[["method"]]))

spnN <- as.numeric(argLs[["span"]])

## checking
##---------

if(packageVersion("ropls") < "0.10.16")
    cat("\nWarning: new version of the 'ropls' package is available\n", sep="")

stopifnot(refC %in% c("pool", "sample"))

if(refC == "pool" &&
   !any("pool" %in% samDF[, "sampleType"]))
    stop("No 'pool' found in the 'sampleType' column; use the samples as normalization reference instead")


##------------------------------
## Computation
##------------------------------


## ordering (batch and injection order)
##-------------------------------------

samDF[, "ordIniVi"] <- 1:nrow(rawMN)
ordBatInjVi <- order(samDF[, "batch"], samDF[, "injectionOrder"])
rawMN <- rawMN[ordBatInjVi, ]
samDF <- samDF[ordBatInjVi, ]

## signal drift and batch-effect correction
##-----------------------------------------

nrmMN <- shiftBatchCorrectF(rawMN,
                            samDF,
                            refC)

## figure
##-------

cat("\nPlotting\n")

pdf(argLs[["graph_output"]], onefile = TRUE, width = 11, height = 7)
plotBatchF(rawMN, samDF)
plotBatchF(nrmMN, samDF)
dev.off()

## returning to initial order
##---------------------------

ordIniVi <- order(samDF[, "ordIniVi"])
nrmMN <- nrmMN[ordIniVi, ]
samDF <- samDF[ordIniVi, ]
samDF <- samDF[, colnames(samDF) != "ordIniVi", drop=FALSE]


##------------------------------
## Ending
##------------------------------


## saving
##-------

datMN <- nrmMN

datDF <- cbind.data.frame(dataMatrix = colnames(datMN),
                          as.data.frame(t(datMN)))
write.table(datDF,
            file = argLs[["dataMatrix_out"]],
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

varDF <- cbind.data.frame(variableMetadata = rownames(varDF),
                          varDF) ## not modified; for compatibility only
write.table(varDF,
            file = argLs[["variableMetadata_out"]],
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

simDF <- cbind.data.frame(samDF, as.data.frame(datMN)) ## for compatibility
simDF <- cbind.data.frame(names = rownames(simDF),
                          simDF)
write.table(simDF,
            file = argLs[["variable_for_simca"]],
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

res <- list(dataMatrix_raw = rawMN,
            dataMatrix_normalized = nrmMN,
            sampleMetadata = samDF)
save(res,
     file = argLs[["rdata_output"]]) ## for compatibility

## closing
##--------

cat("\nEnd of '", modNamC, "' Galaxy module call: ",
    as.character(Sys.time()), "\n", sep = "")

## sink()

options(stringsAsFactors = strAsFacL)

rm(argLs)

