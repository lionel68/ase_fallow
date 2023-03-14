############
# R script to make figures 
# of the main part of the manuscript
#############################################

# set wd
setwd("~/Thuenen/ASE_fallow/gitrepo/ase_fallow/")

# load library
library(tidyverse)
library(sf)

## figure 1 fallow land over district 
# load data
dfallow <- st_read("data/district_fallow.gpkg")
dfallow$fallow_d_cut <- fct_relevel(dfallow$fallow_d_cut, c("-20; -10%","-10; -2%", "-2; 0%", "0; 2%", "2; 10%", "10; 45%"))
dfallow$fallow_cut <- fct_relevel(dfallow$fallow_cut, c("0%", ">0 - 2%", ">2 - 5%", ">5 - 10%", ">10 - 25%", ">25%"))

# first plot the conditions in 2007, 2010 and 2016
gg_1 <- ggplot(subset(dfallow, year %in% c(2007, 2010, 2016))) +
  geom_sf(aes(fill = fallow_cut), size = 0.05) +
  facet_wrap(~year_cat, ncol = 3) +
  scale_fill_brewer(palette = "Greens", 
                    name = "Proportion of fallow land\non agricultural land (%)") +
  theme(legend.position = "right",
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# now plot the changes in fallow land percentage
gg_2 <- ggplot(subset(dfallow, year %in% c(2010, 2016))) +
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
ll_rich <- list.files("model_output/stanmodelrich/")
m_rich <- lapply(ll_rich, function(x) readRDS(paste0("model_output/stanmodelrich/", x)))
names(m_rich) <-c(2007, 2010, 2016)
# a data frame with the values of edge density to use
edge_df <- expand.grid(year = c(2007, 2010, 2016),
                       edge_std = seq(-1.46, 2.39, length.out = 25),
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
# put back the edge values in their original scale
edge_df$edge_mha <- edge_df$edge_std * 28 + 46.7
# the plot
gg_fig2 <- ggplot(edge_df, aes(x= edge_mha, y=fallow_eff_med, 
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

# now get species-level slopes
ll2 <- list.files("model_output/coef_species/", pattern = "csv")

# load sp info
sp_inf <- read.csv("data/bird_info.csv")

sp_beta <- plyr::ldply(ll2, function(x) read.csv(paste0("model_output/coef_species/", x)))
sp_beta %>%
  separate(id, into = c("species", "year"),
           sep = -5) %>%
  mutate(year = gsub("_", "", year)) %>%
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

ggsave("figures/03_fallow_effect_species.png", gg_f)


## fig 4 compute geometric mean for changes in fallows at the group-level

# load files
ll <- list.files("model_output/cond_species/")
sp_all <- plyr::ldply(ll, function(x) read.csv(paste0("model_output/cond_species/", x)))

# add back id infos
sp_all <- inner_join(sp_all, sp_inf, by = c("species" = "Species.name"))

# compute relative changes from no fallows per species, year and landscape context
# then compute geometric mean per group, year and landscape context
sp_all %>%
  mutate(edge = factor(edge_std, labels = c("low SWF", "average SWF", "high SWF")),
         fallow_per = (fallow_std * 0.61 + 1.59) ** 2) %>%
  group_by(group, year, edge, fallow_per) %>%
  summarise(r_mean = exp(sum(log(Estimate)) / length(unique(species))),
            r_low = exp(sum(log(Q2.5)) / length(unique(species))),
            r_high = exp(sum(log(Q97.5)) / length(unique(species)))) -> dd

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


ggsave("figures/04_geometric_mean_fallow_meanswf.png", gg_e)
