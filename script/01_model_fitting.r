##############
# R script to fit the species richness
# and group-level abundance models
#######################


# set wd
setwd("~/PostDoc_Thunen/ASE-stuff/")
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


# fit the models and save them in the subfolder model_output/stanmodel
for(yr in c(2007, 2010, 2016)){
  subdat <- subset(bird, year == yr)

  # rich model
  m_rich <- fit_stan_rich(subdat, "rich", "model_output/stanmodel/", paste0("rich_tn", yr))
  # edge model
  m_edge <- fit_stan_abund(subdat, "edge", "model_output/stanmodel/", paste0("edge_zinb", yr))
  extract_beta(m_edge, "model_beta_grp", paste0("edge_", yr))
  # field model
  m_field <- fit_stan_abund(subdat, "field", "model_output/stanmodel/", paste0("field_zinb", yr))
  extract_beta(m_field, "model_beta_grp", paste0("field_", yr))
  # low model
  m_low <- fit_stan_abund(subdat, "low", "model_output/stanmodel/", paste0("low_zinb", yr))
  extract_beta(m_low, "model_beta_grp", paste0("low_", yr))
  
  print(yr)
}