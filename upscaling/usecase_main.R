################################################################################
# Bocedi 2012 Paper
# 2049 x 1025 (res 25m x 25m), midpoint displacement (Saupe 1988, With 1997)
# p = .1 .3 .5 .7, H = .1 .5 .9
# 50 maps each
# upscaling 100m, 250m, 500m, 1000m (sum up then prop for coarser cell)


# Simulation Parameter
fragGra <- c(0.1, 0.5, 0.9)
habPerc <- c(0.1, 0.3, 0.5, 0.7)
size <- 1024
agfac <- as.integer(c(2, 4, 8, 16, 32))
reps <- 70

#clu <- get_default_cluster()
clu <- create_cluster(12)
clu %>%
  cluster_library(c("NLMR", "landscapetools", "raster", "tidyverse")) %>%
  cluster_assign_value("sim_landscape", sim_landscape) %>%
  cluster_assign_value("size", size) %>%
  cluster_assign_value("agfac", agfac)

system.time({
  lanSca <-
    expand.grid(H = fragGra, p = habPerc, repl = 1:reps) %>%
    as_tibble() %>%
    rowid_to_column("id") %>%
    partition(id, cluster = clu) %>%
    mutate(dta = sim_landscape(p, H, size, agfac, dropMaps = TRUE)) %>%
    collect() %>% arrange(id) %>% ungroup()
})

# output statistics
stats <- lanSca %>% 
  unnest() %>% group_by(id) %>%
  select(-starts_with("pc")) %>% gather('scm', 'pinc', starts_with("pi"))
saveRDS(stats, "stats.RDS")

# for map example visualisation rerun with dropMaps = FALSE
# reduce reps to one to reduce memory requirements
map_example <- lanSca %>% filter(p == .7, H == .9) %>% select(dta) %>% unlist() %>% .[1:12] %>% 
  set_names(
    c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10",  "11", "12")
  )
saveRDS(map_example, "map_example.RDS")
