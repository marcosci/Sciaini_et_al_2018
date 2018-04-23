##----------------------------------------------------------------------------------------------------------##
##  Cedric Scherer (cedricphilippscherer@gmail.com)                                                         ##
##  RNetLogo Script for setting up a parallel session                                                       ##
##  2018-02-28                                                                                              ##
##----------------------------------------------------------------------------------------------------------##

#### PREPARATION ---------------------------------------------------------------------------------------------

library(parallel)

source("./R/parameters.R")
source("./R/nlms.R")
source("./R/simulation-fun.R")


#### CLUSTER -------------------------------------------------------------------------------------------------

(proc <- detectCores() - 6)  ## 18
cl <- makeCluster(proc, type = "PSOCK", methods = FALSE, outfile = "")
showConnections()


#### PARAMETERS ----------------------------------------------------------------------------------------------

### Sim()
runs <- 1:length(nlms)  ## number of landscapes

### Seed()
s <- 1:proc  ## different seeds for each node

### Output
results <- as.list(runs)


############################################################################################################-
#### CLUSTER SIMULATIONS -------------------------------------------------------------------------------------

## export variables to nodes
clusterExport(cl, c("nlms", "runs"))

## start NetLogo in each node
invisible(parLapply(cl, 1:proc, PreSim, gui = gui, nl_path = nl_path, model_path = model_path))

## set seed for each node
parSapply(cl, 1:proc, Seed)

## start parallel simulations by calling parSapply
(start <- Sys.time()); results <- parSapply(cl, runs, SimWeekly); (time <- Sys.time() - start)

## save results
saveRDS(results, "./output/results_100x52.Rds")

## stop NetLogo + cluster
invisible(parLapply(cl, 1:proc, PostSim))
closeAllConnections()


############################################################################################################-
#### SINGLE SIMULATION --------------------------------------------------------------------------------------

library(RNetLogo)
NLStart(nl_path, gui = T)
NLLoadModel(model_path)
results <- sapply(1, SimWeekly)
PostSim()
