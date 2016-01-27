library(RUnit)

wrapperF <- function(argVc) {

    source("../batch_correction_all_loess_script.R")


#### Start_of_testing_code <- function() {}


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

    ## log file
    ##---------

    ## sink(argVc["information"]) ## not implemented

    cat("\nStart of the '", modNamC, "' Galaxy module call: ",
        format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")

    ## loading
    ##--------

    rawMN <- t(as.matrix(read.table(argVc["dataMatrix"],
                                    header = TRUE,
                                    row.names = 1,
                                    sep = "\t")))

    samDF <- read.table(argVc["sampleMetadata"],
                        header = TRUE,
                        row.names = 1,
                        sep = "\t")

    varDF <- read.table(argVc["variableMetadata"],
                        check.names = FALSE,
                        header = TRUE,
                        row.names = 1,
                        sep = "\t") ## not used; for compatibility only

    refC <- tolower(gsub("all_loess_", "", argVc["method"]))

    spnN <- as.numeric(argVc["span"])

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
                                refC,
                                spnN)

    ## figure
    ##-------

    cat("\nPlotting\n")

    pdf(argVc["graph_output"], onefile = TRUE, width = 11, height = 7)
    plotBatchF(rawMN, samDF, spnN)
    plotBatchF(nrmMN, samDF, spnN)
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
                file = argVc["dataMatrix_out"],
                quote = FALSE,
                row.names = FALSE,
                sep = "\t")

    varDF <- cbind.data.frame(variableMetadata = rownames(varDF),
                              varDF) ## not modified; for compatibility only
    write.table(varDF,
                file = argVc["variableMetadata_out"],
                quote = FALSE,
                row.names = FALSE,
                sep = "\t")

    simDF <- cbind.data.frame(samDF, as.data.frame(datMN)) ## for compatibility
    simDF <- cbind.data.frame(names = rownames(simDF),
                              simDF)
    write.table(simDF,
                file = argVc["variable_for_simca"],
                quote = FALSE,
                row.names = FALSE,
                sep = "\t")

    res <- list(dataMatrix_raw = rawMN,
                dataMatrix_normalized = nrmMN,
                sampleMetadata = samDF)
    save(res,
         file = argVc["rdata_output"]) ## for compatibility

    ## closing
    ##--------

    cat("\nEnd of '", modNamC, "' Galaxy module call: ",
        as.character(Sys.time()), "\n", sep = "")

    ## sink()

    options(stringsAsFactors = strAsFacL)


#### End_of_testing_code <- function() {}


    return(list(datMN = datMN))

    rm(list = ls())

}

exaDirOutC <- "output"
stopifnot(file.exists(exaDirOutC))

tesArgLs <- list(input_allLoessPool = c(method = "all_loess_pool",
                     span = "1",
                     .chkC = "checkEqualsNumeric(outLs[['datMN']][1, 1], 25803076, tolerance = 1e-3)"),
                 input_allLoessSample = c(method = "all_loess_sample",
                     span = "1",
                     .chkC = "checkEqualsNumeric(outLs[['datMN']][1, 1], 23402048, tolerance = 1e-3)"),
                 mrey_allLoessSample = c(method = "all_loess_sample",
                     span = "1",
                     .chkC = "checkEqualsNumeric(outLs[['datMN']][1, 1], 21732604, tolerance = 1e-3)"),
                 mrey_allLoessSampleSpan06 = c(method = "all_loess_sample",
                     span = "0.6",
                     .chkC = "checkEqualsNumeric(outLs[['datMN']][1, 1], 134619170, tolerance = 1e-3)"))

for(tesC in names(tesArgLs))
    tesArgLs[[tesC]] <- c(tesArgLs[[tesC]],
                          dataMatrix = file.path(unlist(strsplit(tesC, "_"))[1], "dataMatrix.tsv"),
                          sampleMetadata = file.path(unlist(strsplit(tesC, "_"))[1], "sampleMetadata.tsv"),
                          variableMetadata = file.path(unlist(strsplit(tesC, "_"))[1], "variableMetadata.tsv"),
                          dataMatrix_out = file.path(exaDirOutC, "dataMatrix.tsv"),
                          variableMetadata_out = file.path(exaDirOutC, "variableMetadata.tsv"),
                          variable_for_simca = file.path(exaDirOutC, "variable_for_simca.tsv"),
                          graph_output = file.path(exaDirOutC, "graph_output.pdf"),
                          rdata_output = file.path(exaDirOutC, "rdata_output.rdata"))

for(tesC in names(tesArgLs)) {
    print(tesC)
    outLs <- wrapperF(tesArgLs[[tesC]])
    if(".chkC" %in% names(tesArgLs[[tesC]]))
        stopifnot(eval(parse(text = tesArgLs[[tesC]][[".chkC"]])))
}

message("Checks successfully completed")
