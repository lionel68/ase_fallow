############
# R script to make figures and table
# of the manuscript
#############################################

# set wd
setwd("~/Thuenen/ASE_fallow/ASE_fallow/")

# load library
library(tidyverse)
library(sf)

# load data
bird <- read.csv("data/preprocessed/bird_fallow_v14.csv")
# some transformation / scaling of covariates
bird$fallow_sqrt <- sqrt(bird$fallow_per)
bird$fallow_std <- scale(bird$fallow_sqrt)
bird$edge_std <- scale(bird$edge_mha)
bird$agri_std <- scale(bird$agriculture)

# load fitted models
ff_rich <- list.files("model_output/stanmodelrich/")
m_rich <- lapply(ff_rich, function(x) readRDS(paste0("model_output/stanmodelrich/", x)))

ff_abund <- list.files("model_output/stanmodelabund/")
m_abund <- lapply(ff_abund, function(x) readRDS(paste0("model_output/stanmodelabund/", x)))

## figure 1 fallow land over district 
# load data
district <- read.csv("data/preprocessed/district_fallow_v2.csv")

# some descriptive stats
district %>%
  mutate(fallow_pera = (SEFA / ARAB) * 100) %>%
  select(year, fallow_per, fallow_pera) %>%
  pivot_longer(fallow_per:fallow_pera) %>%
  group_by(year, name) %>%
  summarise(Mi = min(value, na.rm=TRUE),
            Ma = max(value, na.rm=TRUE),
            Med = median(value, na.rm=TRUE),
            Nb4 = sum(value > 4, na.rm=TRUE),
            Nb10 = sum(value > 10, na.rm=TRUE))

dshape <- st_read("data/geodata/250_NUTS3.shp")
# plot the fallow land values and changes between the time points

# put the district with multiple lines in one
dshape %>%
  rename(NUTS_code = NUTS_CODE) %>%
  group_by(NUTS_code) %>%
  summarise() -> dshape_m

district %>%
  arrange(NUTS_code, year) %>%
  group_by(NUTS_code) %>%
  mutate(fallow_change = fallow_per - lag(fallow_per),
         year_cat = factor(year)) %>%
  full_join(dshape_m, by = "NUTS_code") %>%
  ungroup() %>%
  mutate(fallow_cut = cut(fallow_per, c(-1, 0, 2, 5, 10, 25, 75), 
                          include.lowest = TRUE,
                          labels = c("0%", ">0 - 2%", ">2 - 5%",
                                     ">5 - 10%", ">10 - 25%", ">25%")),
         fallow_d_cut = cut(fallow_change, c(-20, -10, -2, 0, 2, 10, 45),
                            include.lowest = TRUE,
                            labels = c("-20; -10%", "-10; -2%",
                                       "-2; 0%", "0; 2%", "2; 10%",
                                       "10; 45%")))-> dd

dd <- st_as_sf(dd)
dd$is_bird <- factor(dd$is_bird, labels = c("Absent", "Present"))
dd$year_cat <- factor(dd$year, labels = c("Census year:\n2007", "Census year:\n2010", "Census year:\n2016"))
dd$year_cat2 <- factor(dd$year, labels = c("Census year:\n2007", "Abolishment\nof mandatory\nset-aside:\n2007 - 2010",
                                           "Establishment\nof EFAs:\n2010 - 2016"))


# first plot the conditions in 2007, 2010 and 2016
gg_1 <- ggplot(subset(dd, year %in% c(2007, 2010, 2016))) +
  geom_sf(aes(fill = fallow_cut), size = 0.05) +
  facet_wrap(~year_cat, ncol = 3) +
  scale_fill_brewer(palette = "Greens", 
                    name = "Proportion of fallow land\non agricultural land (%)") +
  theme(legend.position = "right",
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# now plot the changes in fallow land percentage
gg_2 <- ggplot(subset(dd, year %in% c(2010, 2016))) +
  geom_sf(aes(fill = fallow_d_cut), size = 0.05) +
  facet_wrap(~year_cat2, ncol = 2) +
  scale_fill_brewer(palette = "PRGn", 
                    name = "Changes in proportion\n of fallow land on\nagricultural land (p. p.)") +
  theme(legend.position = "right",
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# put the two together
gg_a <- gridExtra::grid.arrange(gg_1, gg_2, nrow = 2)
ggsave("figures/01_fallow_figure.png", gg_a)

## fig 2 fallow land effect on richness
# load the richness models
names(m_rich) <-c(2007, 2010, 2016)
# a data frame with the values of edge density to use
edge_df <- expand.grid(year = c(2007, 2010, 2016),
                       edge_std = seq(quantile(bird$edge_std, probs = 0.025), 
                                      quantile(bird$edge_std, probs = 0.975), 
                                      length.out = 25),
                       fallow_eff_med = 0,
                       fallow_eff_q5 = 0,
                       fallow_eff_q95 = 0)
# run the computation for each years
for(yr in c(2007, 2010, 2016)){
  edge_val <- edge_df[edge_df$year == yr, "edge_std"]
  post <- posterior::as_draws_matrix(m_rich[[as.character(yr)]], variable = "fallow", regex = TRUE)
  fallow_eff <- posterior::summarise_draws(sapply(edge_val, 
                                                  function(x) post[,1] + 
                                                    x * post[,2] + 
                                                    (x ** 2) * post[,3]))
  edge_df[edge_df$year == yr, "fallow_eff_med"] <- fallow_eff$median
  edge_df[edge_df$year == yr, "fallow_eff_q5"] <- fallow_eff$q5
  edge_df[edge_df$year == yr, "fallow_eff_q95"] <- fallow_eff$q95
}


edge_df$year_cat <- factor(edge_df$year)
# the plot
gg_fig2 <- ggplot(edge_df, aes(x=(edge_std * sd(bird$edge_mha)) + mean(bird$edge_mha), y=fallow_eff_med, 
                               ymin=fallow_eff_q5, ymax=fallow_eff_q95)) +
  geom_ribbon(alpha=0.15, aes(fill=year_cat)) +
  geom_path(aes(color=year_cat)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             size = 1.25) +
  scale_fill_brewer(name = "Survey year:", palette = "Accent") +
  scale_color_brewer(name = "Survey year:", palette = "Accent") +
  labs(x = "Edge density (meter per hectares)",
       y = "Effect of fallow land") +
  theme_classic()

ggsave("figures/02_fallow_rich_year.png", gg_fig2)

## fig 3 fallow land effect on pop size of species and species group
# two figures: for main effect of fallow and for interaction with edge density
# load coef from meta-analysis model
mcoef <- readRDS("model_output/stanmodelmeta/meta_species.rds")
group_dat <- expand.grid(year = c("2007", "2010", "2016"),
                         group = c("edge", "field", "low"),
                         beta_se = 0.14) # median of se

group_dat <- cbind(group_dat,
                   brms::posterior_summary(brms::posterior_epred(mcoef, newdata = group_dat)))
group_dat$group <- factor(group_dat$group,
                          labels = c("Edge-breeding", "Field-breeding","Foraging visitors"))
group_dat$group <- fct_relevel(group_dat$group, c("Field-breeding", "Edge-breeding",
                                                  "Foraging visitors"))
group_dat$faa <- "Group-level"

# now for interaction term
mcoef_int <- readRDS("model_output/stanmodelmeta/meta_species_int.rds")
group_dat <- expand.grid(year = c("2007", "2010", "2016"),
                         group = c("edge", "field", "low"),
                         beta_se = 0.10) # median of se

group_int <- cbind(group_dat,
                   brms::posterior_summary(brms::posterior_epred(mcoef_int, newdata = group_dat)))
group_int$group <- factor(group_int$group,
                          labels = c("Edge-breeding", "Field-breeding","Foraging visitors"))
group_int$group <- fct_relevel(group_int$group, c("Field-breeding", "Edge-breeding",
                                                  "Foraging visitors"))
group_int$faa <- "Group-level"

# now get species-level slopes
ll2 <- list.files("model_output/coef_species/", pattern = "csv")

# load sp info
sp_inf <- read.csv2("data/birddata/bird_info.csv")

sp_beta <- plyr::ldply(ll2, function(x) read.csv(paste0("model_output/coef_species/", x)))
sp_beta %>%
  separate(id, into = c("species", "year"),
           sep = -5) %>%
  mutate(year = gsub("_", "", year),
         species = gsub("_", " ", species)) %>%
  rename(Species.name = species) %>%
  arrange(param, beta) %>%
  inner_join(sp_inf, by = "Species.name") %>%
  mutate(English.name = factor(English.name, levels = unique(English.name)),
         Category = factor(Category, levels = c("Field-breeders", "Edge-breeders",
                                                "Foraging visitors"),
                           labels = c("Field-breeding", "Edge-breeding", "Forgaging visitors")),
         year = factor(year, levels = c("2007", "2010", "2016"))) -> sp_betad


gg_fallow <- ggplot(subset(sp_betad, param == "fallow"), aes(color = year)) +
  geom_hline(yintercept = 0) +
  #geom_hline(data = ssd, aes(yintercept = fallow_m, color = year), linetype = "dashed", size = 0.75) +
  geom_linerange(aes(x = English.name,
                     ymin = beta - 2 * beta_se,
                     ymax = beta + 2 * beta_se),
                 position = position_dodge2(width = 0.5)) +
  geom_point(aes(x = English.name, y = beta),
             position = position_dodge2(width = 0.5), size = 0.95) +
  facet_wrap(~Category, scales = "free") +
  coord_flip() +
  scale_y_continuous(minor_breaks = NULL, breaks = c(-1, 0, 1, 2)) +
  scale_color_brewer(name = "Survey year:", palette = "Accent", breaks = c("2016", "2010", "2007")) +
  labs(x = "",
       y = "Effect of fallow land")

gg_fallowgrp <- ggplot(group_dat, aes(x=group,y=Estimate,ymin=Q2.5, ymax=Q97.5 , color = year)) +
  geom_hline(yintercept = 0) +
  geom_linerange(position = position_dodge2(width = 0.5)) +
  geom_point(position = position_dodge2(width = 0.5), size = 0.95) +
  facet_wrap(~faa) +
  coord_flip() +
  scale_color_brewer(name = "Survey year:", palette = "Accent", breaks = c("2016", "2010", "2007")) +
  labs(x = "",
       y = "")

gg_f <- ggpubr::ggarrange(gg_fallowgrp, gg_fallow, nrow = 2, heights = c(1, 3),
                  common.legend = TRUE, legend = "right")

ggsave("figures/fallow_effect_species_new.png", gg_f)

gg_int <- ggplot(subset(sp_betad, param == "fallow:edge"), aes(color = year)) +
  geom_hline(yintercept = 0) +
  #geom_hline(data = ssd, aes(yintercept = int_m, color = year), linetype = "dashed", size = 0.75) +
  geom_linerange(aes(x = English.name,
                     ymin = beta - 2 * beta_se,
                     ymax = beta + 2 * beta_se),
                 position = position_dodge2(width = 0.5)) +
  geom_point(aes(x = English.name, y = beta),
             position = position_dodge2(width = 0.5), size = 0.95) +
  facet_wrap(~Category, scales = "free") +
  coord_flip() +
  scale_color_brewer(name = "Survey year:", palette = "Accent", breaks = c("2016", "2010", "2007")) +
  labs(x = "",
       y = "Interaction effect of fallow land: edge density")

gg_intgrp <- ggplot(group_int, aes(x=group,y=Estimate,ymin=Q2.5, ymax=Q97.5 , color = year)) +
  geom_hline(yintercept = 0) +
  geom_linerange(position = position_dodge2(width = 0.5)) +
  geom_point(position = position_dodge2(width = 0.5), size = 0.95) +
  facet_wrap(~faa) +
  coord_flip() +
  scale_color_brewer(name = "Survey year:", palette = "Accent", breaks = c("2016", "2010", "2007")) +
  labs(x = "",
       y = "")

gg_i <- ggpubr::ggarrange(gg_intgrp, gg_int, nrow = 2, heights = c(1, 3),
                          common.legend = TRUE, legend = "right")

ggsave("figures/interaction_effect_species_new.png", gg_i)

## dev: test hypothesis of fallow effect per species
tt <- plyr::ldply(m_list, function(x) hypothesis(x, "fallow_std > 0")$hypothesis)
tt$id <- ff

tt %>%
  separate(id, c("genus", "species", "year", "foo")) %>%
  unite(genus_sp, genus:species) %>%
  group_by(genus_sp) %>%
  summarise(n95 = sum(Post.Prob > 0.95),
            n90 = sum(Post.Prob > 0.90)) -> dd


# looking at the hump-shaped function for the richness models

## load models
ff <- list.files("model_output/stanmodelrich/")
m_rich <- lapply(ff, function(x) readRDS(paste0("model_output/stanmodelrich/", x)))

## get the edge value at the peak of the hump
foo <- function(x){
  psamp <- as_draws_df(x, variable = "fallow_std", regex = TRUE)
  ## identify peak effect, -b1 / 2*b2
  peak_samp <- - (psamp[,2] / (2 * psamp[,3]))
  # back-transform
  out <- posterior_summary(peak_samp * sd(bird$edge_mha) + mean(bird$edge_mha))
  return(out)
}

plyr::ldply(m_rich, function(x) foo(x))

## get the posterior proba that the funciton is hump-shaped
foo <- function(x){
  out <- hypothesis(x, "fallow_std:Iedge_stdE2 < 0")
  return(out)
}

plyr::llply(m_rich, function(x) foo(x))

## fig 4 compute geometric mean for changes in fallows at the group-level

# first get group composition
group <- read.csv("data/bird_fallow_group.csv")

# load files
ll <- list.files("model_output/cond_species/")
sp_all <- plyr::ldply(ll, function(x) read.csv(paste0("model_output/cond_species/", x)))

# add back id infos
sp_all$group <- rep(group$group, each = 50 * 9)

# compute relative changes from no fallows per species, year and landscape context
# then compute geometric mean per group, year and landscape context
sp_all %>%
  mutate(edge = factor(edge_std, labels = c("low SWF", "average SWF", "high SWF")),
         fallow_per = (fallow_std * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)) ** 2) %>%
  group_by(group, year, edge, fallow_per) %>%
  summarise(r_mean = exp(sum(log(Estimate)) / length(unique(species))),
            r_low = exp(sum(log(Q2.5)) / length(unique(species))),
            r_high = exp(sum(log(Q97.5)) / length(unique(species)))) -> dd

gg_e <- ggplot(subset(dd, group == "edge"), aes(x = fallow_per, color = factor(year))) +
  geom_ribbon(aes(ymin = r_low, ymax = r_high), alpha = 0.05) +
  geom_path(aes(y = r_low), size = 0.5) +
  geom_path(aes(y = r_high), size = 0.5) +
  geom_path(aes(y = r_mean), size = 1.5) +
  facet_grid(cols = vars(edge)) +
  scale_color_brewer(name = "Survey year:", palette = "Accent") +
  scale_fill_brewer(name = "Survey year:", palette = "Accent") +
  scale_y_log10() +
  labs(x = "", y = "Group-level bird abundance", title = "Edge-breeding")


gg_f <- ggplot(subset(dd, group == "field"), aes(x = fallow_per, color = factor(year))) +
  geom_ribbon(aes(ymin = r_low, ymax = r_high), alpha = 0.05) +
  geom_path(aes(y = r_low), size = 0.5) +
  geom_path(aes(y = r_high), size = 0.5) +
  geom_path(aes(y = r_mean), size = 1.5) +
  facet_grid(cols = vars(edge)) +
  scale_color_brewer(name = "Survey year:", palette = "Accent") +
  scale_fill_brewer(name = "Survey year:", palette = "Accent") +
  scale_y_log10() +
  labs(x = "", y = "", title = "Field-breeding")

gg_l <- ggplot(subset(dd, group == "low"), aes(x = fallow_per, color = factor(year))) +
  geom_ribbon(aes(ymin = r_low, ymax = r_high), alpha = 0.05) +
  geom_path(aes(y = r_low), size = 0.5) +
  geom_path(aes(y = r_high), size = 0.5) +
  geom_path(aes(y = r_mean), size = 1.5) +
  facet_grid(cols = vars(edge)) +
  scale_color_brewer(name = "Survey year:", palette = "Accent") +
  scale_fill_brewer(name = "Survey year:", palette = "Accent") +
  scale_y_log10() +
  labs(x = "Percent fallow land", y = "", title = "Foraging visitors")


gga <- gg_f / gg_e / gg_l + plot_layout(guides = "collect")

ggsave("figures/geometric_mean_fallow.png", gga, width = 8, height = 8)

# now the same for average SWF
dd$group <- factor(dd$group, labels = c("Edge-breeding", "Field-breeding", "Foraging visitor"))
dd$group <- fct_relevel(dd$group, c("Field-breeding", "Edge-breeding", "Foraging visitor"))

gg_e <- ggplot(subset(dd, edge == "average SWF"), aes(x = fallow_per, color = factor(year),
                                                      fill = factor(year))) +
  geom_ribbon(aes(ymin = r_low, ymax = r_high), alpha = 0.1) +
  geom_path(aes(y = r_low)) +
  geom_path(aes(y = r_high)) +
  geom_path(aes(y = r_mean), size = 1.5) +
  facet_grid(cols = vars(group)) +
  scale_color_brewer(name = "Survey year:", palette = "Accent") +
  scale_fill_brewer(name = "Survey year:", palette = "Accent") +
  scale_y_log10() +
  labs(x = "Percent fallow land", y = "Group-level bird abundance")


ggsave("figures/geometric_mean_fallow_meanswf.png", gg_e)
