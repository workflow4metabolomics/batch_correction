test_input_allLoessPool <- function() {

    testDirC <- "input"
    determineL <- FALSE
    allLoessL <- TRUE
    argLs <- list(method = "all_loess_pool",
                  span = "1")

    if(!allLoessL)
        argLs[["analyse"]] <- ifelse(determineL, "determine_bc", "batch_correction")
    
    argLs <- c(defaultArgF(testDirC, determineL = determineL), argLs)
    outLs <- wrapperCallF(argLs, allLoessL = allLoessL)  

    checkEqualsNumeric(outLs[['datMN']][1, 1], 25803076, tolerance = 1e-3)

}

test_input_allLoessSample <- function() {

    testDirC <- "input"
    determineL <- FALSE
    allLoessL <- TRUE
    argLs <- list(method = "all_loess_sample",
                  span = "1")

    if(!allLoessL)
        argLs[["analyse"]] <- ifelse(determineL, "determine_bc", "batch_correction")
    
    argLs <- c(defaultArgF(testDirC, determineL = determineL), argLs)
    outLs <- wrapperCallF(argLs, allLoessL = allLoessL)
    
    checkEqualsNumeric(outLs[['datMN']][1, 1], 23402048, tolerance = 1e-3)

}

test_example1_allLoessSample <- function() {

    testDirC <- "example1"
    determineL <- FALSE
    allLoessL <- TRUE
    argLs <- list(method = "all_loess_sample",
                  span = "1")

    if(!allLoessL)
        argLs[["analyse"]] <- ifelse(determineL, "determine_bc", "batch_correction")
    
    argLs <- c(defaultArgF(testDirC, determineL = determineL), argLs)
    outLs <- wrapperCallF(argLs, allLoessL = allLoessL)    

    checkEqualsNumeric(outLs[['datMN']][1, 1], 21732604, tolerance = 1e-3)

}

test_example1_allLoessSampleSpan06 <- function() {

    testDirC <- "example1"
    determineL <- FALSE
    allLoessL <- TRUE
    argLs <- list(method = "all_loess_sample",
                  span = "0.6")

    if(!allLoessL)
        argLs[["analyse"]] <- ifelse(determineL, "determine_bc", "batch_correction")
    
    argLs <- c(defaultArgF(testDirC, determineL = determineL), argLs)
    outLs <- wrapperCallF(argLs, allLoessL = allLoessL)

    checkEqualsNumeric(outLs[['datMN']][1, 1], 134619170, tolerance = 1e-3)

}

test_sacurine_allLoessPool <- function() {

    testDirC <- "sacurine"
    determineL <- FALSE
    allLoessL <- TRUE
    argLs <- list(method = "all_loess_pool",
                  span = "1")

    if(!allLoessL)
        argLs[["analyse"]] <- ifelse(determineL, "determine_bc", "batch_correction")
    
    argLs <- c(defaultArgF(testDirC, determineL = determineL), argLs)
    outLs <- wrapperCallF(argLs, allLoessL = allLoessL)
 
    checkEqualsNumeric(outLs[['datMN']]["HU_neg_017", "M53T345"], 7902.366, tolerance = 1e-3)


}

test_sacurine_determinebc <- function() {

    testDirC <- "sacurine"
    determineL <- TRUE
    allLoessL <- FALSE
    argLs <- list(ref_factor = "batch",
                  span = "none")

    if(!allLoessL)
        argLs[["analyse"]] <- ifelse(determineL, "determine_bc", "batch_correction")

    argLs <- c(defaultArgF(testDirC, determineL = determineL), argLs)
    outLs <- wrapperCallF(argLs, allLoessL = allLoessL)

    checkEqualsNumeric(outLs[['sumDF']]["M59T62", "batch.2.linear"], 3, tolerance = 1e-3)


}
