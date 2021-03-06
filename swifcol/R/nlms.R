##----------------------------------------------------------------------------------------------------------##
##  Cedric Scherer (cedricphilippscherer@gmail.com)                                                         ##
##  Creating NLMs as .txt files                                                                             ##
##  2017-02-28                                                                                              ##
##----------------------------------------------------------------------------------------------------------b##

library(NLMR)
library(landscapetools)
source("./R/parameters.R")

## three landscape setups: 
##   * homogeneous
##   * mosaictess with low fragmentation
##   * mosaictess landscapes with high fragmentation
##
## mean breeding capacity = 4.5 females
##   -> 10 classes between 0 and 10 for heterogeneous setups

## paramters
n_col <- 25
n_row <- 50

## homogeneous setup (only one needed)
hom <- rep(4.5, n_col * n_row)
write(hom, file = paste0("./model/nlms/hom.txt"), sep = " ")

## set seed for reproducible neutral landscape model generation
set.seed(2018)

for (i in 1:n) {  ## n = number of replications -> assigned in ".R/parameters.R"
  ## low number of patches: germs = 15
  low <- nlm_mosaictess(ncol = n_col, nrow = n_row, resolution = 1, germs = 15)
  low <- util_classify(low, weighting = rep(0.1, 10))  ## 10 classes
  write(low@data@values, file = paste0("./model/nlms/frag_low_", i, ".txt"), sep = " ")
  
  ## high number of patches: germs = 150 
  high <- nlm_mosaictess(ncol = n_col, nrow = n_row, resolution = 1, germs = 150)
  high <- util_classify(high, weighting = rep(0.1, 10))  ## 10 classes
  write(high@data@values, file = paste0("./model/nlms/frag_high_", i, ".txt"), sep = " ")
}
