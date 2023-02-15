##############
# R script to fit the species richness models
# and run model checks
#######################


# set wd
setwd("~/Thuenen/ASE_fallow/ASE_fallow/")
# source helper functions
source("script/rscript/XX_helper_function.R")

# load libraries
library(tidyverse)
library(DHARMa)
library(brms)

# load data
bird <- read.csv("data/preprocessed/bird_fallow_v14.csv")

# some transformation / scaling of covariates
bird$fallow_sqrt <- sqrt(bird$fallow_per)
bird$fallow_std <- scale(bird$fallow_sqrt)
bird$edge_std <- scale(bird$edge_mha)
bird$agri_std <- scale(bird$agriculture)


# fit the models and save them in the subfolder model_output/stanmodelrich
for(yr in c(2007, 2010, 2016)){
  subdat <- subset(bird, year == yr)

  # rich model
  m_rich <- fit_stan_rich(subdat, "rich", "model_output/stanmodelrich/", paste0("rich_tn", yr))

  print(yr)
}

# grab the list of saved model object names
ff <- list.files("model_output/stanmodelrich//")
# read the model files
m_list <- lapply(ff, function(x) readRDS(paste0("model_output/stanmodelrich/", x)))

# run convergence, divergence, ESS, distribution, dispersion and spatial checks
for(i in 1:12){
  yrr <- substr(ff[i], nchar(ff[i])-8+1, nchar(ff[i])-4)
  coln <- strsplit(ff[i], split = "_")[[1]][1]
  id <- paste0(coln, "_", yrr)
  dat <- subset(bird, year == yrr)
  oo <- stan_check(m_list[[i]], dat, colname = coln, write_file = TRUE, file_id = id)
  print(oo)
}

# run posterior predictive checks
# plot density of simulated data vs original data
# first for species richness
ppd_rich <- lapply(m_list[10:12], function(x) brms::pp_check(x, type = "dens_overlay", ndraws = 100))
gg_rich <- gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs=ppd_rich, nrow=1, ncol=3))
                                                 
# save the plot
ggsave("figures/XX_ppd_rich_tn.png", gg_rich)



# make plot of spatial residuals
ss <- NULL
for(i in 1:12){
  yrr <- substr(ff[i], nchar(ff[i])-8+1, nchar(ff[i])-4)
  coln <- strsplit(ff[i], split = "_")[[1]][1]
  id <- paste0(coln, "_", yrr)
  dat <- subset(bird, year == yrr)
  ss <- rbind(ss, stan_spatres(m_list[[i]], dat, coln, id))
}

# separate id cols
ss <- separate(ss, id, c("var", "year"))
ss$var <- factor(ss$var, levels = c("rich", "edge", "field", "low"),
                 labels = c("Species richness", "Edge-breeders", "Field-breeders",
                            "Foraging visitors"))

# get the Moran I value
mm <- read.csv("model_output/model_checks.csv")
mm <- separate(mm, id, c("var", "year"))
mm$var <- factor(mm$var, levels = c("rich", "edge", "field", "low"),
                 labels = c("Species richness", "Edge-breeders", "Field-breeders",
                            "Foraging visitors"))
mm$label <- paste0("Moran's I: ", round(mm$spat, 2))
# x/y coords to write the Moran's I index
mm$x <- 3550000
mm$y <- 6087002

# a DE shape
de <- sf::st_read("data/geodata/DE_mainland.shp")
de <- sf::st_transform(de, 31467)

gg_spat <- ggplot() +
  geom_point(data = ss, aes(x=x, y=y, color=res), size =0.5) +
  geom_sf(data = de, fill = NA) +
  geom_label(data = mm, aes(x=x, y=y, label = label), size = 2.1) +
  facet_grid(rows = vars(var),
             cols = vars(year)) +
  scale_color_gradient2(midpoint = 0.5, 
                        name = "Scaled residuals") +
  theme_classic() +
  labs(x="", y="") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave("figures/XX_spatial_residuals.png", gg_spat)

# make correlogram, code can take some time to run
cc <- list()
for(i in 1:3){
  yrr <- substr(ff[i], nchar(ff[i])-8+1, nchar(ff[i])-4)
  coln <- strsplit(ff[i], split = "_")[[1]][1]
  id <- paste0(coln, "_", yrr)
  dat <- subset(bird, year == yrr)
  cc[[i]] <- spline_correlog(m_list[[i]], dat, coln)
}

# load the correlogram
# cc <- readRDS("model_output/spline_correlo.rds")

# put in df format
cc_d <- plyr::ldply(cc, function(x) data.frame(x = x$real$predicted$x[1,],
                                               corr = x$real$predicted$y[1,],
                                               lci = x$boot$boot.summary$predicted$y[2,],
                                               hci = x$boot$boot.summary$predicted$y[10,]))

cc_d$year <- rep(c(2007, 2010, 2016), each = 300)

cc_d %>%
  filter(x < 1e5) %>%
  mutate(x = x / 1e3) -> cc_f


gg_spat <- ggplot(cc_f, aes(x=x, y=corr, ymin=lci, ymax=hci)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(alpha = 0.1) +
  geom_path() +
  facet_grid(cols = vars(year)) +
  labs(x = "Distance between sampling plots (km)",
       y = "Correlation in model residuals (with 95% confidence bands)")

ggsave("figures/XX_correlogram.png", gg_spat, height = 5)
