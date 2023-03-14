############
# R script to fit species-level
# abundance models
# and meta-analysis models on
# fitted coefficients
#############################################

# set wd
setwd("~/Thuenen/ASE_fallow/gitrepo/ase_fallow/")
# source helper functions
source("script/XX_helper_function.R")

# load libraries
library(brms)
options(mc.cores = parallel::detectCores())
library(tidyverse)

# load data
bird <- read.csv("data/bird_abundance.csv")
group <- read.csv("data/bird_info.csv")

# a species x year data.frame
species_year <- expand.grid(species = unique(bird$species),
                            year = c(2007, 2010, 2016))

# loop through this
for(i in 1:nrow(species_year)){
  id <- paste(species_year[i, "species"],
              species_year[i, "year"], sep = "_")
  # subset the data
  data <- subset(bird, species == species_year[i, "species"] & 
                   year == species_year[i, "year"])
  
  # fit the model
  model <- fit_stan_abund(data, "Abundance", "model_output/standmodelabund/", file_id = id)

  # extract beta parameters of fallow effect and interaction
  betaa <- extract_beta(model, file_name = id, file_id = id)

  # compute expected bird abundance along the fallow land gradient for
  # different values of edge density
  
  # create new data frame ranging from 0% to 6% (90% quantile) fallow 
  # for 10%, mean and 90% quantile of edge density
  newdat <- expand.grid(fallow_std = seq(-1.26, 1.26, length.out = 50),
                        edge_std = c(-1.06, 0, 1.28),
                        agri_std = 0)
  # model predictions
  pp <- posterior_epred(model, newdata = newdat, re.form = NA)
  ppred <- rbind(posterior_summary(apply(pp[,1:50], 2, function(x) x / pp[,25])),
       posterior_summary(apply(pp[,51:100], 2, function(x) x / pp[,75])),
       posterior_summary(apply(pp[,101:150], 2, function(x) x / pp[,125])))
  # the output dataframe
  mm <-cbind(newdat[,1:2], ppred)
  mm$species <- species_year[i, "species"]
  mm$year <- species_year[i, "year"]
  
  write.csv(mm,
            paste0("model_output/cond_species/", id, ".csv"))

  gc()
  rm(model)
  print(paste0(id, " is done!\n"))
  
}

# meta-analysis on species-level effects
ll <- list.files("model_output/coef_species/", pattern = "csv")
spp <- plyr::ldply(ll, function(x) read.csv(paste0("model_output/coef_species/", x)))
spp$Species.name <- substr(spp$id, 1, nchar(spp$id) - 5)
spp <- inner_join(spp, group, by = "Species.name")

spp$year <- factor(substr(spp$id, nchar(spp$id) - 3, nchar(spp$id)))
spp$group <- fct_relevel(spp$group, c("low", "field", "edge")) # put low as ref level

mcoef <- brm(beta | se(beta_se, sigma = TRUE) ~ year + group,
             data = subset(spp, param == "fallow"))

saveRDS(mcoef, "model_output/stanmodelmeta/meta_species.rds")

mcoef2 <- brm(beta | se(beta_se, sigma = TRUE) ~ year + group,
              data = subset(spp, param == "fallow:edge"))

saveRDS(mcoef2, "model_output/stanmodelmeta/meta_species_int.rds")
