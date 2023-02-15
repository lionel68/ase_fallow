############
# R script to fit species-level
# abundance models, run model checks
# and extract slope coefficients
#############################################

# set wd
setwd("~/Thuenen/ASE_fallow/ASE_fallow/")
# source helper functions
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
  model <- fit_stan_abund(data, "Abundance", "model_output/stanmodelabund/", file_id = id)
  
  # perform model checks
  checks <- stan_check(model, data, "Abundance", write_file = TRUE, file_id = id)
  
  # extract beta parameters of fallow effect and interaction
  betaa <- extract_beta(model, file_name = id)

  # compute expected bird abundance along the fallow land gradient for
  # different values of edge density
  
  # create new data frame ranging from 0% to 6% (90% quantile) fallow 
  # for 10%, mean and 90% quantile of edge density
  newdat <- expand.grid(fallow_std = seq(-1.26, 1.26, length.out = 50),
                        edge_std = (c(17, 47, 83) - mean(bird$edge_mha)) / sd(bird$edge_mha),
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

  print(paste0(id, " is done!\n"))
  
}

# run posterior predictive checks
rr <- list.files("model_output/stanmodelabund/", pattern = "20")
field_sp <- group$species[group$group == "field"]
field_sp <- gsub(" ", "_", field_sp)
edge_sp <- group$species[group$group == "edge"]
edge_sp <- gsub(" ", "_", edge_sp)
low_sp <- group$species[group$group == "low"]
low_sp <- gsub(" ", "_", low_sp)

foo <- function(x){
  m <- readRDS(paste0("model_output/stanmodelabund/", x))
  sp_lat <- paste(strsplit(x, split = "_")[[1]][1:2], collapse = " ")
  yr <- substr(x, nchar(x)-7, nchar(x)-4)
  spp <- sp_inf$English.name[sp_inf$Species.name == sp_lat]
  pp <- brms::pp_check(m, ndraws = 100) + scale_x_log10() + labs(title = paste(spp, yr, sep = " "))
  return(pp)
}

field <- plyr::llply(rr[grep(paste(field_sp, collapse = "|"), rr)], function(x) foo(x))
ff <- patchwork::wrap_plots(field, ncol = 3)
ggsave("figures/pp_check_fields_new.png", ff)

edge <- plyr::llply(rr[grep(paste(edge_sp, collapse = "|"), rr)], function(x) foo(x))
ee <- patchwork::wrap_plots(edge, ncol = 3)
ggsave("figures/pp_check_edge_new.png", ee, height = 12)

low <- plyr::llply(rr[grep(paste(low_sp, collapse = "|"), rr)], function(x) foo(x))
ll <- patchwork::wrap_plots(low[1:18], ncol = 3)
ggsave("figures/pp_check_low1_new.png", ll, height = 12)
ll2 <- patchwork::wrap_plots(low[19:36], ncol = 3)
ggsave("figures/pp_check_low2_new.png", ll2, height = 12)


# grab R2
foo_r2 <- function(x){
  m <- readRDS(paste0("model_output/stanmodel/", x))
  r2_m <- bayes_R2(m, re.form = NA)[1,1]
  r2_c <- bayes_R2(m)[1,1]
  out <- data.frame(species = substr(x, 1, nchar(x) - 9),
                    year = substr(x, nchar(x)-7, nchar(x)-4),
                    r2m = r2_m,
                    r2c = r2_c)
  return(out)
}

r2_sp <- plyr::ldply(rr, function(x) foo_r2(x))


# meta-analysis on species-level effects
ll <- list.files("model_output/coef_species/", pattern = "csv")
spp <- plyr::ldply(ll, function(x) read.csv(paste0("model_output/coef_species/", x)))

spp$group <- rep(group$group, each = 6)
spp$year <- factor(substr(spp$id, nchar(spp$id) - 3, nchar(spp$id)))
spp$group <- fct_relevel(spp$group, c("low", "field", "edge")) # put low as ref level

mcoef <- brm(beta | se(beta_se, sigma = TRUE) ~ year + group,
             data = subset(spp, param == "fallow"))

saveRDS(mcoef, "model_output/stanmodel/meta_species.rds")

mcoef2 <- brm(beta | se(beta_se, sigma = TRUE) ~ year + group,
              data = subset(spp, param == "fallow:edge"))

saveRDS(mcoef2, "model_output/stanmodel/meta_species_int.rds")

# proba of differences between groups
## p(edge > field) = 0.99
sum(posterior_samples(mcoef, pars = "groupfield") < 0) / 4000
## p(edge > low) = 0.97
sum(posterior_samples(mcoef, pars = "grouplow") < 0) / 4000
## p(field > low) = 0.15
sum(posterior_samples(mcoef, pars = "groupfield") > 
      posterior_samples(mcoef, pars = "grouplow")) / 4000

# proba of differences between years
## p(2007 > 2010) = 0.75
sum(posterior_samples(mcoef, pars = "year2010") < 0) / 4000
## p(2007 > 2016) = 0.62
sum(posterior_samples(mcoef, pars = "year2016") < 0) / 4000
## p(2010 > 2016) = 0.36
sum(posterior_samples(mcoef, pars = "year2010") > 
      posterior_samples(mcoef, pars = "year2016")) / 4000
