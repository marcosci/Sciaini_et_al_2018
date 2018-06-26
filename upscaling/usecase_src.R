# required packages

library(NLMR)
library(landscapetools)
library(raster)

library(tidyverse)
library(viridis)

#devtools::install_github("hadley/multidplyr")
library(multidplyr)

sim_landscape <- function(pe, ha, sz, ag, keepsize = TRUE, dropMaps = TRUE) {
  
  # make a list with the initial and 'scaled' maps
  # generate base map
  mapraw <- util_classify(nlm_mpd(sz, sz, 25, roughness = ha, verbose = FALSE), c(1 - pe, pe))
  
  # plain aggregation with mean
  # fix size differencies from midpointdisplacement
  maplist_mean <- list(crop(mapraw, extent(mapraw, 1, sz, 1, sz)))
  # scaling (aggregation) loop over given factors
  for (i in seq_along(ag)) {
    sclmp <- aggregate(maplist_mean[[1]], ag[i], fun = mean)
    # cut back into same amount of tiles if true
    if (keepsize) {
      sclmp <- disaggregate(sclmp, ag[i])
    }
    maplist_mean <- append(maplist_mean, sclmp)
  }
  maplist_mean <- set_names(maplist_mean, c("1", paste0(ag)))
  
  # plain aggregation with modal
  # fix size differencies from midpointdisplacement
  maplist_moda <- list(crop(mapraw, extent(mapraw, 1, sz, 1, sz)))
  # scaling (aggregation) loop over given factors
  for (i in seq_along(ag)) {
    sclmp <- aggregate(maplist_moda[[1]], ag[i], fun = modal)
    # cut back into same amount of tiles if true
    if (keepsize) {
      sclmp <- disaggregate(sclmp, ag[i])
    }
    maplist_moda <- append(maplist_moda, sclmp)
  }
  maplist_moda <- set_names(maplist_moda, c("1", paste0(ag)))
  
  # return statistics if maps are dropped
  r <- tibble(agf = as.integer(c(1, ag))) %>%
    # calculate desired metrics
    mutate(
      # habitat percentage (calculated)
      pc_mean = map_dbl(maplist_mean, function(x) {sum(values(x) != 0) / length(values(x))}),
      # habitat percentage increase
      pi_mean = pc_mean - pc_mean[1],
      pc_moda = map_dbl(maplist_moda, function(x) {sum(values(x) != 0) / length(values(x))}),
      pi_moda = pc_moda - pc_moda[1]
    )
  
  if (dropMaps) {
    return(list(r))
  }else{
    return(list(append(maplist_mean, maplist_moda)))
  }
}
