############
# R script to make figures and table
# of the manuscript
#############################################

# set wd
setwd("~/PostDoc_Thunen/ASE-stuff/")

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
ff <- list.files("model_output/stanmodel/", pattern = "zinb|tn")
m_list <- lapply(ff, function(x) readRDS(paste0("model_output/stanmodel/", x)))

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
ff_rich <- list.files("model_output/stanmodel/", pattern = "rich")
m_rich <- lapply(ff_rich, function(x) readRDS(paste0("model_output/stanmodel/", x)))
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
  post <- as_draws_matrix(m_rich[[as.character(yr)]], variable = "fallow", regex = TRUE)
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
# load species betas
sp_beta <- read.csv("model_output/stanbeta/model_beta_zinb.csv")

sp_beta %>%
  separate(id, into = c("species", "year"),
           sep = -5) %>%
  mutate(year = gsub("_", "", year)) %>%
  arrange(param, beta) -> sp_beta

sp_beta$species <- factor(sp_beta$species, levels = unique(sp_beta$species))


# load group betas
group_beta <- read.csv("model_output/stanbeta/model_beta_grp.csv")

group_beta %>%
  separate(id, into = c("species", "year"),
           sep = -5) %>%
  mutate(year = gsub("_", "", year)) %>%
  arrange(param, beta) -> group_beta

group_beta$species <- factor(group_beta$species, labels = c("Edge-breeding",
                                                            "Field-breeding",
                                                            "Foraging visitors"))

# the plots
gg_figsp <- ggplot(sp_beta, aes(x=species, y=beta, ymin=beta-2*beta_se,
                    ymax=beta+2*beta_se, color = factor(year))) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.01) +
  geom_linerange(position = position_dodge2(width = 0.5)) +
  geom_point(position = position_dodge2(width = 0.5), size = 0.95) +
  facet_grid(~param) +
  coord_flip() +
  scale_color_brewer(name = "Survey year:", palette = "Accent") +
  labs(x = "",
       y = "Effect of fallow land")

gg_figgrp <- ggplot(group_beta, aes(x=species, y=beta, ymin=beta-2*beta_se,
                    ymax=beta+2*beta_se, color = factor(year))) +
  geom_linerange(position = position_dodge2(width = 0.5)) +
  geom_point(position = position_dodge2(width = 0.5), size = 0.95) +
  facet_grid(~param) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.05) +
  coord_flip() +
  scale_color_brewer(name = "Survey year:", palette = "Accent") +
  labs(x = "",
       y = "")

gg_fig3 <- ggpubr::ggarrange(gg_figgrp, gg_figsp, nrow = 2, 
                             common.legend = TRUE, heights = c(1, 3))
ggsave("figures/03_fallow_abund.png", gg_fig3, width = 7.45,
       height = 6, units = "in")

# SI figure effect size of fallow on species group for given values
# of edge densities

edge_val <- c(quantile(bird$edge_std, probs = 0.1),
              mean(bird$edge_std),
              quantile(bird$edge_std, probs = 0.9))
ii <- NULL

for(i in 1:9){
  yrr <- substr(ff[i], nchar(ff[i])-8+1, nchar(ff[i])-4)
  coln <- strsplit(ff[i], split = "_")[[1]][1]
  id <- paste0(coln, "_", yrr)
  
  ii <- rbind(ii, extract_interaction(m_list[[i]], edge_val, id))
  
}

# turn into nicer format
ii %>%
  separate(is, c("species", "year")) %>%
  mutate(edge_cat = factor(rep(c("low", "med", "high"), 9),
                           levels = c("low", "med", "high"),
                           labels = c("Low edge\ndensity",
                                      "Medium edge\ndensity",
                                      "High edge\ndensity")),
         species = factor(species, labels = c("Edge-\nbreeding",
                                              "Field-\nbreeding",
                                              "Foraging\nvisitor"))) %>%
  arrange(species, year) -> ii_dd

# a plot
gg_grp <- ggplot(ii_dd, aes(x=species, y=med, ymin=lci, ymax=hci, color = year)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.01) +
  geom_point(position = position_dodge2(width = 0.4)) +
  geom_linerange(position = position_dodge2(width = 0.4)) +
  facet_grid(cols = vars(edge_cat)) +
  labs(x = "",
       y = "Effect of fallow land") +
  scale_color_brewer(name = "Survey year:", palette = "Accent")

ggsave("figures/XX_fallow_eff_grp.png", gg_grp, width = 9, height = 5, units = "in")


## SI figure, effect of edge density on response
# define fallow land values, 0, 4 and 10%
df <- data.frame(fallow_std = (sqrt(c(0, 4, 10)) - mean(bird$fallow_sqrt)) / sd(bird$fallow_sqrt))
cc <- make_conditions(df, vars = "fallow_std")
m_rich <- m_list[10:12]
eff <- NULL
for(i in 1:3){
  yrr <- c(2007, 2010, 2016)[i]
  tmp <- conditional_effects(m_rich[[i]], effects = "edge_std",
                             conditions = cc)
  tmp$edge_std$year <- yrr
  eff <- rbind(eff, tmp$edge_std[,c("year", "fallow_std","edge_std", "estimate__",
                                    "lower__", "upper__")])
}

eff$fallow_std <- factor(eff$fallow_std, labels = c("0%", "4%", "10%"))
eff <- filter(eff, edge_std < 2.4 & edge_std > -1.46)
eff$edge <- (eff$edge_std * sd(bird$edge_mha)) + mean(bird$edge_mha)
eff$year <- factor(eff$year)

gg_edge_rich <- ggplot(eff, aes(x=edge, y=estimate__, ymin=lower__, ymax=upper__)) +
  geom_ribbon(aes(fill=year), alpha = 0.1) +
  geom_path(aes(color=year)) +
  facet_wrap(~fallow_std, labeller = labeller(fallow_std = c("0%" = "Fallow proportion: 0%",
                                                             "4%" = "Fallow proportion: 4%",
                                                             "10%" = "Fallow proportion: 10%"))) +
  scale_fill_brewer(name = "Survey year:", palette = "Accent") +
  scale_color_brewer(name = "Survey year:", palette = "Accent") +
  labs(x = "Edge density (meter per hectares)",
       y = "Predicted bird richness (with 95% credible interval)")

ggsave("figures/XX_landscape_richness.png", gg_edge_rich)

# now for abundance groups
ff_abund <- list.files("model_output/stanmodel/", pattern = "zinb")
m_abund <- lapply(ff_abund, function(x) readRDS(paste0("model_output/stanmodel/", x)))
effa <- NULL
for(i in 1:length(m_abund)){
  yrr <- substr(ff_abund[i], nchar(ff_abund[i])-8+1, nchar(ff_abund[i])-4)
  coln <- strsplit(ff_abund[i], split = "_")[[1]][1]
  tmp <- conditional_effects(m_abund[[i]], effects = "edge_std",
                             conditions = cc)
  tmp$edge_std$year <- yrr
  tmp$edge_std$var <- coln
  effa <- rbind(effa, tmp$edge_std[,c("var", "year", "fallow_std","edge_std", "estimate__",
                                    "lower__", "upper__")])
}

effa$fallow_std <- factor(effa$fallow_std, labels = c("0%", "4%", "10%"))
effa <- filter(effa, edge_std < 2.4 & edge_std > -1.46)
effa$edge <- (effa$edge_std * sd(bird$edge_mha)) + mean(bird$edge_mha)
effa$year <- factor(effa$year)

gg_edge_abunde <- ggplot(subset(effa, var=="edge"), aes(x=edge, y=estimate__, ymin=lower__, ymax=upper__)) +
  geom_ribbon(aes(fill=year), alpha = 0.1) +
  geom_path(aes(color=year)) +
  facet_wrap(~fallow_std, labeller = labeller(fallow_std = c("0%" = "Fallow proportion: 0%",
                                                             "4%" = "Fallow proportion: 4%",
                                                             "10%" = "Fallow proportion: 10%"))) +
  scale_fill_brewer(name = "Survey year:", palette = "Accent") +
  scale_color_brewer(name = "Survey year:", palette = "Accent") +
  labs(x = "",
       y = "",
       title = "Edge-breeders")

gg_edge_abundl <- ggplot(subset(effa, var=="low"), aes(x=edge, y=estimate__, ymin=lower__, ymax=upper__)) +
  geom_ribbon(aes(fill=year), alpha = 0.1) +
  geom_path(aes(color=year)) +
  facet_wrap(~fallow_std, labeller = labeller(fallow_std = c("0%" = "Fallow proportion: 0%",
                                                             "4%" = "Fallow proportion: 4%",
                                                             "10%" = "Fallow proportion: 10%"))) +
  scale_fill_brewer(name = "Survey year:", palette = "Accent") +
  scale_color_brewer(name = "Survey year:", palette = "Accent") +
  labs(x = "",
       y = "",
       title = "Foraging visitors")

gg_edge_abundf <- ggplot(subset(effa, var=="field"), aes(x=edge, y=estimate__, ymin=lower__, ymax=upper__)) +
  geom_ribbon(aes(fill=year), alpha = 0.1) +
  geom_path(aes(color=year)) +
  facet_wrap(~fallow_std, labeller = labeller(fallow_std = c("0%" = "Fallow proportion: 0%",
                                                             "4%" = "Fallow proportion: 4%",
                                                             "10%" = "Fallow proportion: 10%"))) +
  scale_fill_brewer(name = "Survey year:", palette = "Accent") +
  scale_color_brewer(name = "Survey year:", palette = "Accent") +
  labs(x = "",
       y = "",
       title = "Field-breeders")

gg_edge_abund <- ggpubr::ggarrange(gg_edge_abunde, gg_edge_abundf, gg_edge_abundl,
                  nrow=3, ncol=1, common.legend = TRUE)

gg2 <- ggpubr::annotate_figure(gg_edge_abund, bottom = "Edge density (meter per hectares)",
                        left = "Predicted bird abundance (with 95% credible intervals)")

ggsave("figures/XX_landscape_abund.png", gg2, width = 7.45,
       height = 4 * 3, units = "in")

