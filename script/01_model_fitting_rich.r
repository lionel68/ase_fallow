##############
# R script to fit the species richness models
#######################


# set wd
setwd("~/Thuenen/ASE_fallow/gitrepo/ase_fallow/")
# source helper functions
source("script/XX_helper_function.R")

# load libraries
library(tidyverse)
library(DHARMa)
library(brms)

# load data
bird <- read.csv("data/bird_rich.csv")


# fit the models and save them in the subfolder model_output/stanmodelrich
for(yr in c(2007, 2010, 2016)){
  subdat <- subset(bird, year == yr)

  # rich model
  m_rich <- fit_stan_rich(subdat, "rich", "model_output/stanmodelrich/", paste0("rich_tn", yr))

  print(yr)
}

