# sensitivity analysis script dropping species one at a time
# and checking the impact on the slopes and interaction

setwd("")

library(brms)
options(mc.cores = parallel::detectCores())
library(tidyverse)

# load bird data
bird <- read.csv("data/preprocessed/bird_fallow_species.csv",
                 stringsAsFactors = FALSE)
# load covariates
# fallow data
fallow <- read.csv("data/landusedata/mhb_1kmbuffer_fallow.csv")

# bkr data
bkr <- read.csv("data/landusedata/mhb_bkr.csv")
# edge length data combining ATKIS and SWF
edge <- read.csv("data/landusedata/mhb_1kmbuffer_edge_atkisswf.csv")
edge$edge_m[is.na(edge$edge_m)] <- 0 # replace NAs by 0

# atkis data
atkis <- read.csv("data/landusedata/mhb_1kmbuffer_atkis.csv")
pro_atki <- function(data, col){
  
  #col2 <- enquo(col)
  data$type <- ""
  data$type[data$obb %in% c("Siedlung", "Verkehr")] <- "urban"
  data$type[data$obb == "Gewässer"] <- "water"
  data$type[data$obabez %in% c("Wald, Forst ", "Gehölz ")] <- "forest"
  data$type[data$obabez %in% c("Sumpf, Ried ", "Moor, Moos ")] <- "water"
  data$type[data$obabez %in% c("Gartenland ", "Sonderkultur ", "Ackerland ")] <- "agriculture"
  data$type[data$obabez %in% c("Grünland ", "Heide ")] <- "grassland"
  
  data %>%
    group_by(routcode, type) %>%
    summarise(area_ha = sum({{ col }} * 1e-4)) -> out
  
  return(out)
}
atkis %>%
  pro_atki(atkis_sqm) -> atkis_dd

# make this in wide format for merging
atkis_dd %>%
  filter(type != "") %>%
  pivot_wider(names_from = type, 
              values_from = area_ha, 
              values_fill = 0) -> atkis_dd2


# load original models
m_list <- readRDS("model_output/fitted_model.rds")

# a data frame to keep everything
bird_group <- read.csv("data/birddata/bird_fallow_group.csv")
df_out <- bird_group[,c(3,1)]
df_out$species <- gsub(" ", "_", df_out$species)
df_out$var_fallow <- NA
df_out$var_inter <- NA


# loop through the species and remove them from the computation
for(sp in unique(bird$species)){
  # find which group to re-compute
  gr <- bird$group[bird$species == sp][1]
  
  # re-compute that group
  bird %>%
    filter(species != sp) %>% 
    filter(group == gr) %>%
    group_by(species) %>%
    mutate(abundance_scale = round(scales::rescale(Abundance, c(0, 100)), 0)) %>%
    ungroup() %>%
    group_by(routcode, year) %>%
    summarise(Abundance = sum(abundance_scale)) -> bird_abu
  
  # adding the covariates
  # put together
  bird_abu %>%
    filter(year %in% c(2007, 2010, 2016)) %>%
    left_join(edge[,-1], by = "routcode") %>%
    left_join(atkis_dd2, by = "routcode") %>%
    left_join(fallow[,-1], by = c("routcode", "year")) %>%
    left_join(bkr[,-1]) %>%
    filter(!is.na(urban)) %>%
    filter(!is.na(fallow_ha)) -> bird_dd2
  
  # add missing bkr infos
  bird_dd2$bkr[bird_dd2$routcode == "ni130"] <- unique(bird_dd2$bkr[bird_dd2$routcode == "ni129"])
  bird_dd2$bkr[bird_dd2$routcode == "sh42"] <- unique(bird_dd2$bkr[bird_dd2$routcode == "sh43"])
  # add missing edge infos
  bird_dd2$edge_m[is.na(bird_dd2$edge_m)] <- 0
  
  # some transformation
  bird_dd2$year_cat <- factor(bird_dd2$year)
  bird_dd2$fallow_sqrt <- sqrt(bird_dd2$fallow_ha)
  bird_dd2$fallow_std <- scale(bird_dd2$fallow_sqrt)
  bird_dd2$edge_std <- scale(bird_dd2$edge_m)
  bird_dd2$urban_std <- scale(bird_dd2$urban)
  bird_dd2$forest_std <- scale(bird_dd2$forest)
  bird_dd2$bkr <- factor(bird_dd2$bkr)
  
  # fit the model
  m <- brm(bf(Abundance ~ year_cat + edge_std * fallow_std + (1 | bkr),
              zi ~ urban_std + forest_std),
           data = bird_dd2, 
           family = "zero_inflated_negbinomial")
  
  # compare posterior samples to original
  df_out$var_fallow[df_out$species == sp] <-
  sum(posterior_samples(m_list[[gr]], pars = "b_fallow_std", fixed = TRUE) <
        posterior_samples(m, pars = "b_fallow_std", fixed = TRUE)) / 4000
  
  df_out$var_inter[df_out$species == sp] <-
    sum(posterior_samples(m_list[[gr]], pars = "edge_std:fallow_std") 
        < posterior_samples(m, pars = "edge_std:fallow_std")) / 4000
  
  print(sp)

}

write.csv(df_out, "model_output/sensitivity_onespeciesout.csv", row.names = FALSE)