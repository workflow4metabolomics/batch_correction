#!/usr/bin/Rscript --vanilla --slave --no-site-file

################################################################################################
# batch_correction_main_wrapper                                                                #
#                                                                                              #
# Author: Nils Paulhe                                                                          #
# User: Galaxy                                                                                 #
# Original data: --                                                                            #
# Starting date: 2017-12-11                                                                    #
# Version 1: 2017-12-11                                                                        #
#                                                                                              #
#                                                                                              #
#                                                                                              #
################################################################################################

library(batch) #necessary for parseCommandArgs function

##------------------------------
## init. prog. constants
##------------------------------

argv.wrapper <- commandArgs(trailingOnly = FALSE)
script.path <- sub("--file=", "", argv.wrapper[grep("--file=", argv.wrapper)])
prog.name <- basename(script.path)

##------------------------------
## init. functions
##------------------------------

script_bypass <- function(other.script.name) {
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- dirname(script.name)
  other.script.fullpath <- paste(sep="/", script.basename, other.script.name)
  # print(paste("calling", other.script.fullpath, "from", script.name))
  other.script.cmd <- paste(sep=" ", "Rscript", other.script.fullpath, "-h")
  system(other.script.cmd, wait=TRUE)
}

source_wrapper <- function(other.script.name){
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- dirname(script.name)
  other.script.fullpath <- paste(sep="/", script.basename, other.script.name)
  #print(paste("calling", other.script.fullpath, "from", script.name))
  source(other.script.fullpath)
}

##------------------------------
## Test Help
##------------------------------

if (length(grep('-h', argv.wrapper)) > 0) {
  cat("Usage: Rscript ", 
    prog.name,
    "{args} \n",
    "parameters: \n",
    "\t-h: display this help message, call all scripts with the same option and exit (optional) \n",
    "\t--loess \"TRUE\": call the script as \"batch_correction_all_loess_wrapper.R\"; otherwise call it as \"batch_correction_wrapper.R\" one (optional) \n",
    "for other parameters, please refer to each script specific options and parameters. \n",
    "\n")
  script_bypass("batch_correction_all_loess_wrapper.R")
  script_bypass("batch_correction_wrapper.R")
  quit(status = 0)
}

##------------------------------
## check if loess or normal
##------------------------------

if (length(grep('--loess', argv.wrapper)) > 0) {
  source_wrapper("batch_correction_all_loess_wrapper.R")
} else {
  source_wrapper("batch_correction_wrapper.R")
}

