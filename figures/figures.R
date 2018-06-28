# load packages ----
library(NLMR)
library(landscapetools)
library(tidyverse)
library(magrittr)
library(patchwork)       
library(scales)       
library(pals)       
library(readr)          

### FIG 1 ----
#### create landscapes ----
set.seed(5)
perco_lation <- nlm_percolation(128, 128, 0.4)
ran_dom <- nlm_random(128, 128)
dista_grad <-
  nlm_distancegradient(128,
    128,
    origin = c(10, 15, 10, 15)
  )
ed_ge <- nlm_edgegradient(128, 128, 90)
plan_ary <- nlm_planargradient(128, 128)
random_reccluster <- nlm_randomrectangularcluster(128, 128,
                                                   minl = 8,
                                                   maxl = 16,
                                                   rescale = TRUE
                                                  )

random_clustery <- nlm_randomcluster(128, 128,
                                     p = 0.5,
                                     ai = c(0.3, 0.6, 0.1),
                                     rescale = TRUE
                                    )

random_mpd <- nlm_mpd(100, 100, roughness = 0.7, resolution = 7.751940)
gau_ssian <- nlm_gaussianfield(128, 128, autocorr_range = 35, resolution = 10)
f_bm <- nlm_fbm(128, 128, fract_dim = 1)
mos_aic <-nlm_mosaicfield(128, 128, 20)
neigh_boorcluster <- nlm_neigh(128,
                               128,
                               p_neigh = 0.75,
                               categories = 5,
                               p_empty = 0.01,
                               neighbourhood = 4,
                               proportions = c(0.25, 0.05, 0.1, 0.3, 0.3))
cur_ds <- util_rescale(nlm_curds(curds = c(0.9, 0.63, 0.33),
                    recursion_steps = c(2, 4, 16)))
whey_s <- util_rescale(nlm_curds(curds = c(0.33, 0.33, 0.33),
                    recursion_steps = c(16, 4, 2),
                    wheyes = c(0.22, 0.22, 0.22)))
poly_scapes_1 <- nlm_mosaictess(128, 128, germs = 20)
poly_scapes_2 <- nlm_mosaicgibbs(128,128,
                                 germs = 20,
                                 R = 0.02,
                                 patch_classes = 12,
                                 rescale = TRUE)


#### collect landscapes in list ----
beastiary <- list(
"a) Random Curdling" = cur_ds,
"b) Distance Gradient" = dista_grad,
"c) Edge Gradient" = ed_ge,
"d) Fractional Brownian Motion" = f_bm,
"e) Gaussian Random Field" = gau_ssian,
"f) Mosaic Random Field" = mos_aic,
"g) Random Neighborhood" = neigh_boorcluster,
"h) Percolation" = perco_lation,
"i) Planar Gradient" = plan_ary,
"j) Polygonal Landscapes (Tesselation)" = poly_scapes_1,
"k) Polygonal Landscapes (Gibbs)" = poly_scapes_2,
"l) Random" = ran_dom,
"m) Random Cluster" = random_clustery,
"n) Midpoint Displacement" = random_mpd,
"o) Random Rectangular Cluster" = random_reccluster,
"p) Wheyed Random Curdling" = whey_s)

#### rescale and plot ----
beastiary <- purrr::map(beastiary, util_rescale)
util_facetplot(beastiary) + scale_fill_gradientn(colours=pals::parula(100)) + 
    ggplot2::theme(strip.text = ggplot2::element_text(size = 18))
ggsave("bestiary.eps", height = 20, width = 16, device=cairo_ps)

### FIG 2 ----
#### binarize, classify and merge fbm ----
binarized_raster <- util_binarize(f_bm, breaks = 0.31415)
classified_raster <- util_classify(f_bm,
                                   c(0.33, 0.33, 0.33))

merge_vis <- list(
    "a) Binarized Landscape" = binarized_raster,
    "b) Classified Landscape" = util_rescale(classified_raster),
    "c) Merged Landscape" =  util_merge(classified_raster, neigh_boorcluster)
)

#### plot landscapes ----
util_facetplot(merge_vis) +
    ggplot2::scale_fill_gradientn(colours=pals::parula(100), guide = "colourbar")+ 
    ggplot2::theme(strip.text = ggplot2::element_text(size = 18))
ggsave("landscapetools.eps", height = 4, width = 12, device=cairo_ps)

### Fig 3 ----
#### load data ----
results <- readRDS("../swiftcol/output/results_100x50.Rds")
source("../swiftcol/R/parameters.R")

#### prepare for plotting ----
scen <- nlms %>% 
    as_tibble() %>% 
    transmute(Run = 1:n(),
              Landscape = gsub("[^[a-z]", "", value)) %>%
    mutate(Landscape = factor(Landscape, levels = c("hom", "fraglow", "fraghigh")))

results_proc <- results %>% 
    t() %>% 
    tbl_df() %>% 
    unnest() %>% 
    full_join(scen, ., by = "Run") %>% 
    mutate(Run = factor(Run)) %>% 
    group_by(Run) %>% 
    mutate(WeekOutbreak = Weeks - WeekRelease,
           WeekLast = max(Weeks),
           IndAll = IndSusc + IndTran + IndLeth + IndImmu,
           IndInf = IndTran + IndLeth)

inf <- results_proc %>% 
    group_by(Landscape) %>% 
    summarize(sum_inf = sum(IndInf),
              sum_pop = sum(IndAll)) %>% 
    mutate(prop_inf = sum_inf / sum_pop)

p_inf <- results_proc %>% 
    group_by(Run, Landscape) %>%
    summarize(sum_inf = sum(IndInf)) %>% 
    #summarize(sum_inf = sum(IndInf) / sum(IndAll)) %>%  ## optional proportion
    ungroup() %>% 
    ggplot(aes(Landscape, sum_inf, col = Landscape)) +
    geom_boxplot(width = 0.75, notch = TRUE, outlier.shape = NA) + 
    geom_jitter(width = 0.25, alpha = 0.25) +
    labs(x = "Landscape type", y = "Number of infected individuals", tag = "A") + 
    theme_nlm_discrete(legend.position = "none", axis_title_size = 13, axis_title_just = 0.5, 
                       axis_text_size = 11) +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
    scale_x_discrete(limit = c("hom", "fraglow", "fraghigh"),
                     labels = c("homogeneous","low\nfragmentation","high\nfragmentation")) +
    scale_y_continuous(labels = scales::comma, breaks = seq(0, 1250000, by = 250000))

pers_year <- map_df((1:50) * 52, ~results_proc %>% 
                        filter(WeekOutbreak == .x) %>% 
                        group_by(Landscape) %>% 
                        summarise(prob = n() / n,
                                  year = .x / 52))

## add zero probability to final year (if belowmaximum simulation time of 50 years)
last <- pers_year %>% 
    group_by(Landscape) %>% 
    summarize(last = max(year)) %>% 
    filter(last < 50)

pers_year <- add_row(pers_year, Landscape = pull(last[1]), prob = 0, year = pull(last[2])) %>% 
    arrange(Landscape, year)

p_pers <- ggplot(pers_year, aes(year, prob, col = Landscape)) + 
    geom_line(size = 1) + 
    ylim(0, 1) +
    scale_x_continuous(breaks = seq(0, 50, 5)) + 
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), name = "Landscape type: ", 
                     labels = c("homogeneous", "low fragmentation", "high fragmentation")) +
    labs(x = "Years since outbreak", y = "Probability of persistence", tag = "B") +
    theme_nlm(legend.position = "right", legend_text_size = 11, legend_title_size = 11,
              axis_title_size = 13, axis_title_just = 0.5, axis_text_size = 11, panel.grid.minor = element_blank())


#### final plot ----
p <- p_inf + p_pers 
ggsave(plot = p, "swiftcol.eps", height = 4, width = 12, device=cairo_ps)

### FIG 4 ----

# land example landscapes
mapex <- readRDS('figures/map_example.RDS')
saveRDS(mapex, "figures/map_example.RDS")
names(mapex) <- c("1, unscaled",
                  "2",
                  "4",
                  "8",
                  "16",
                  "32",
                  "1, unscaled",
                  "2",
                  "4",
                  "8",
                  "16",
                  "32")

maptibb <- tibble::enframe(mapex, "id", "maps") %>% dplyr::mutate(maps = purrr::map(.$maps, 
                                                                                    util_raster2tibble)) %>% tidyr::unnest()


maptibb$id <- factor(maptibb$id, levels = c("1, unscaled",
                                            "2",
                                            "4",
                                            "8",
                                            "16",
                                            "32"))


maptibb$id2 <- NA
maptibb$id2[1:(12582902/2)] <- "Average Rule"
maptibb$id2[(12582902/2 + 1):12582902] <- "Majority Rule"
maptibb$id2[which(is.na(maptibb$id2))] <-  "Majority Rule"


plt <- ggplot2::ggplot(maptibb, ggplot2::aes_string("x", "y")) + 
    ggplot2::coord_fixed() + ggplot2::geom_raster(ggplot2::aes_string(fill = "z")) + 
    ggplot2::facet_grid(id~id2) + 
    ggplot2::scale_x_continuous(expand = c(0, 0),
                                limits = c(0, 
                                           max(maptibb$x))) + ggplot2::scale_y_continuous(expand = c(0, 
                                                                                                     0), limits = c(0, max(maptibb$y))) + ggplot2::guides(fill = FALSE) + 
    ggplot2::labs(titel = NULL, x = NULL, y = NULL) + 
    theme_facetplot() +  
    scale_fill_gradientn(colours=pals::parula(100))+ 
    ggplot2::theme(strip.text = ggplot2::element_text(size = 18))

ggsave(plot = plt, "figures/scaling.eps", height = 12, width = 4, device=cairo_ps)

### FIG 5 ----
stats <- readRDS('figures/stats.RDS')
stats$scm[stats$scm == "average_rule"] <- "Average Rule"
stats$scm[stats$scm == "majority_rule"] <- "Majority Rule"

plt2 <- ggplot(stats, aes(agf, pinc, colour = as.character(H))) +
    geom_jitter(alpha = 1/3, shape = 20) +
    stat_summary(fun.y = mean, geom = 'line', aes(group=factor(H))) +
    facet_grid(scm~p, scales = "free_y", labeller = label_value) +
    labs(y = 'Habitat cells increment [%]', colour = 'Roughness') + 
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
    theme_grey(base_size = 18) +
    theme(strip.background = ggplot2::element_rect(fill = "grey80"),
          strip.text = ggplot2::element_text(hjust  = 0,
                                             size   = 13.5,
                                             family = "Roboto Condensed"),
          panel.background = element_rect(fill = NA),
          panel.grid = element_line(colour = "#e4e4e4"),
          axis.text.x = element_text(angle=90, 
                                     size = 10,
                                     vjust= 0.5)) + 
    scale_x_continuous('Landscape resolution [map units]',
                       breaks = c(1, 2, 4, 8, 16, 32),
                       labels = c(1, 2, 4, 8, 16, 32),
                       limits = c(1, 32))


# ggsave("suitable_cells.eps", height = 4, width = 12, device=cairo_ps, fallback_resolution = 600)
ggsave(plot = plt2, "figures/suitable_cells.eps", height = 4, width = 12, device=cairo_ps)


## FIG 6 ----
mee_softwareusage <- read_csv("mee_softwareusage.csv")

mee_softwareusage %<>%
    mutate(Software=replace(Software, Software == "-", NA)) 

mee_softwareusage$Software <- factor(mee_softwareusage$Software, levels = c("R", "MATLAB", "NA", "C++", "SPSS",
                                                                            "MAXENT", "Python", "C",
                                                                            "Colony", "MARK", "Mathematica",
                                                                            "SPOT", "pinpoint", "Genstat",
                                                                            "JAVA", "STSim"
))


ggplot(mee_softwareusage, aes(Software)) +
    geom_bar() + 
    theme_nlm_grey_discrete(axis_title_size = 14) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylab("Number of publications")

ggsave("mee_softwareusage.eps")


