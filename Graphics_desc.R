###########################################################################
#            Quality controls and descriptive analysis plots              #
###########################################################################
# Authors : Melanie Petera                                                #
#                                                                         #
# Description : This script allows various displays of data for quality   #
#      control and descriptive analysis. The input data is a matrix of    #
#      quantitative variables, and it returns chosen plots in png format  #
#      and a table with chosen statistics.                                #
#                                                                         #
# Version 1 (XX-05-2014) : display scatterplot, boxplot, histogram,       #
#      MA plot, density plot, pairs plot, and return a table of chosen    #
#      statistics (quantiles, mean, variance, standard error of the mean) #
###########################################################################

desc_fct <- function(file.in, nacode, file.out, graph.file, stat, ploting){
  # Parameters :
  # - file.in : count matrix input (tab-separated) [file name]
  # - nacode : missing value coding character
  # - file.out : results table of chosen statistics [file name]
  # - graph.file : graphical outputs [file name]
  # - stat : should statistics be calculated ? (TRUE/FALSE)
  # - ploting : should graphics be displayed ? (TRUE/FALSE)

# Data import - - - - - - - - - - - - - - - - - 

Dataset <- read.table(file.in,header=TRUE,na.strings=nacode)


# Statistics table computation - - - - - - - - -

if(stat){
  
  
  
  
  
  
  
  
} # end if(stat)



# Graphics generation - - - - - - - - - - - - - 

if(ploting){
  
  
  
  
  
  
  
  
  
  
  
} # end if(ploting)



} # end of function

