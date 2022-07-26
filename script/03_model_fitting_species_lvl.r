############
# R script to fit species-level
# abundance models, run model checks
# and extract slope coefficients
#############################################

# set wd
setwd("~/Dokumente/Bird-stuff/monvia_voegel/ro2_ase_fallow/")
source("script/rscript/XX_helper_function.R")

# load libraries
library(brms)
options(mc.cores = parallel::detectCores())
library(tidyverse)

# load data
bird <- read.csv("data/preprocessed/bird_fallow_species_2.csv")

# some transformation
bird$fallow_sqrt <- sqrt(bird$fallow_per)
bird$fallow_std <- scale(bird$fallow_sqrt)
bird$edge_std <- scale(bird$edge_mha)
bird$agri_std <- scale(bird$agriculture)

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
  model <- fit_stan_abund(data, "Abundance", "model_output/stanmodel/", file_id = id)
  
  # perform model checks
  checks <- stan_check(model, data, "Abundance", write_file = TRUE, file_id = id)
  
  # extract beta parameters of fallow effect and interaction
  betaa <- extract_beta(model, file_id = id)

  print(paste0(id, " is done!\n"))
  
}