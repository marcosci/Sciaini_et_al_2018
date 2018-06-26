##----------------------------------------------------------------------------------------------------------##
##  Cedric Scherer (cedricphilippscherer@gmail.com)                                                         ##
##  Input and output parameters for simulation-fun.R, nlms.R, and 1_Simulation.R                            ##
##  2018-02-28                                                                                              ##
##----------------------------------------------------------------------------------------------------------##


#### INPUT ---------------------------------------------------------------------------------------------------

### PreSim()
nl_path <- "C:/Program Files/NetLogo 5.3/app"  ## alternative path: NetLogo5.3.1/app
model_path <- paste0(getwd(), "/model/SwiFCoIBM_basic.nlogo")  ## NLStart() and NLLoadModel() need absolute paths
gui <- FALSE

### Sim()

## number of replications
n <- 100 

## landscapes
nlms <- c(rep('\"hom\"', n), 
          paste0('\"frag_low_', seq(n), '\"'), 
          paste0('\"frag_high_', seq(n), '\"'))
           