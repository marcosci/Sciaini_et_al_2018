##----------------------------------------------------------------------------------------------------------##
##  Cedric Scherer (cedricphilippscherer@gmail.com)                                                         ##
##  Functions for running RNetLogo on a cluster                                                             ##
##  2018-02-28                                                                                              ##
##----------------------------------------------------------------------------------------------------------##

PreSim <- function(dummy, gui, nl_path, model_path) {
  ## Load the RNetLogo package, afterwards start NetLogo and the model on each node before starting a simulation
  ## Args:
  ##   gui:        if TRUE, NetLogo will be started with GUI, FALSE will start NetLogo in headless mode
  ##   nl_path:    absolute path to NetLogo.jar
  ##   model_path: absolute path to model.nlogo 
  library(RNetLogo)
  NLStart(nl_path, gui = gui)
  NLLoadModel(model_path)
}

SimWeekly <- function(row) {
  ## Sets variables and runs the model on each node: sets fixed parameters and input variables, calls setup and runs the model
  ## Args:
  ##   row: vector from 1 to number of rows in 'input' (= number of parameter combinations)
  ## Returns:
  ##   Output variables
  library(tidyverse)
  ret <- tibble(Run = numeric(0), Week = numeric(0), IndSucs = numeric(0), IndTran = numeric(0), IndLeth = numeric(0), 
                IndImmu = numeric(0), NewTran = numeric(0), NewLeth = numeric(0), P_infected = numeric(0), P_infecitous = numeric(0), 
                F_infected  = numeric(0), F_infectious = numeric(0), DistInf = numeric(0), WeekMaxDistInf = numeric(0), WeekRelease = numeric(0))
  NLCommand('set Filename', nlms[row])
  NLCommand('setup')
  while (as.numeric(NLReport('DONE')) == 0) {
      NLCommand('go')
      ret_help <- NLReport(c(row,
                             'ticks',
                             'count turtles with [EpiStat = "esSusc"]',
                             'count turtles with [EpiStat = "esTrans"]',
                             'count turtles with [EpiStat = "esLeth"]',
                             'count turtles with [EpiStat = "esImmu"]',
                             'NewTrans',
                             'NewLeth',
                             'count patches with [infected? = 1]',
                             'count patches with [infectious? = 1]',
                             'F_infected',
                             'F_infectious',
                             'InfDist',
                             'MaxDist?',
                             'WeekRelease'))
      ret <- rbind(ret, ret_help)
    }
    names(ret) <- c("Run",
                    "Weeks",
                    "IndSusc",
                    "IndTran",
                    "IndLeth",
                    "IndImmu",
                    "NewTran",
                    "NewLeth",
                    "P_infected",
                    "P_infectious",
                    "F_infected",
                    "F_infectious",
                    "DistInf",
                    "WeekMaxDistInf",
                    "WeekRelease")

    return(ret)
}

PostSim <- function(x){
  ## Stops NetLogo application on each node
  NLQuit()
}

Seed <- function(s) {
  ## Sets a unique seed on each node
  ## Args:
  ##   s: vector from 1 to number of created nodes
  NLCommand('random-seed', s)
  NLReport('random 10')
}
