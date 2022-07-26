# sensitivity analysis script dropping species one at a time
# and checking the impact on the slopes and interaction

setwd("/home/hertzog/Dokumente/Bird-stuff/monvia_voegel/ro2_ase_fallow/")
#setwd("PostDoc_Thunen/ASE-stuff/")

library(brms)
options(mc.cores = parallel::detectCores())
library(dplyr)

# helper function to compare posterior of fallow parameters
comp_post <- function(x, y){
  
  eff1 <- sum(as_draws_matrix(x, variable = "b_fallow_std") <
        as_draws_matrix(y, variable = "b_fallow_std")) / 4000
  
  eff2 <- sum(as_draws_matrix(x, variable = "b_edge_std:fallow_std") 
              < as_draws_matrix(y, variable = "b_edge_std:fallow_std")) / 4000
  
  return(c(eff1, eff2))
  
}

# load bird data
bird <- read.csv("data/preprocessed/bird_fallow_species_2.csv",
                 stringsAsFactors = FALSE)
# list of original models
ff <- list.files("model_output/stanmodel/", pattern = "zinb")
m_list <- lapply(ff, function(x) readRDS(paste0("model_output/stanmodel/", x)))

# a data frame to keep everything
bird_group <- read.csv("data/birddata/bird_fallow_group.csv")
df_out <- bird_group[rep(1:nrow(bird_group), each = 3),c(3,1)]
df_out$species <- gsub(" ", "_", df_out$species)
df_out$year <- rep(c(2007, 2010, 2016), times = nrow(bird_group))
df_out$var_fallow <- NA
df_out$var_inter <- NA


# loop through the species and remove them from the computation
for(sp in unique(bird$species)){
  
  # find which group to re-compute
  gr <- bird$group[bird$species == sp][1]
  
  # load original models
  m_orig <- m_list[grep(gr, ff)]
  
  # re-compute that group
  bird %>%
    filter(species != sp) %>% 
    filter(group == gr) %>%
    group_by(species) %>%
    mutate(abundance_scale = scales::rescale(Abundance, c(0, 100))) %>%
    ungroup() %>%
    group_by(routcode, year, edge_mha, agriculture,
             NUTS_code, NUTS_year, fallow_per) %>%
    summarise(Abundance = sum(abundance_scale)) -> bird_dd2
  
  # add missing edge infos
  bird_dd2$edge_mha[is.na(bird_dd2$edge_mha)] <- 0
  
  # some transformation
  bird_dd2$fallow_sqrt <- sqrt(bird_dd2$fallow_per)
  bird_dd2$fallow_std <- scale(bird_dd2$fallow_sqrt)
  bird_dd2$edge_std <- scale(bird_dd2$edge_mha)
  bird_dd2$agri_std <- scale(bird_dd2$agriculture)
  bird_dd2$Abundance <- round(bird_dd2$Abundance, 0)
  
  # the formula
  form <- "Abundance ~ agri_std + edge_std * fallow_std +
                         (1 | NUTS_code) + (0 + edge_std | NUTS_code)"
  
  # the prior
  bprior <- c(prior(student_t(3, 0, 2.5), class = Intercept),
              prior(normal(0, 2.5), class = b),
              prior("exponential(1)", class = sd),
              prior(gamma(0.01, 0.01), class = shape))

  
  # fit the yearly models
  m_sens <- list()
  for(yr in c(2007, 2010, 2016)){
    dat <- subset(bird_dd2, year == yr)
    m_sens[[length(m_sens) + 1]] <- brm(bf(form,
                                         zi ~ agri_std + edge_std),
                                      data = dat,
                                      family = "zero_inflated_negbinomial",
                                      prior = bprior,
                                      control = list(adapt_delta = 0.9, max_treedepth = 25))
  }

  
  # compare posterior samples to original
  df_out[df_out$species == sp, 4:5] <- t(mapply(comp_post, m_orig, m_sens))
  
  print(sp)

}

write.csv(df_out, "model_output/sensitivity_onespeciesout_new.csv", row.names = FALSE)