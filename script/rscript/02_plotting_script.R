###################
# R script to derive the main plots
###############################

setwd("~/PostDoc_Thunen/ASE-stuff/")

# load libraries
library(tidyverse)
library(sf)
library(brms)

# load data
bird <- read.csv("data/preprocessed/bird_fallow_v10.csv")
district <- read.csv("data/preprocessed/district_fallow.csv")
dshape <- st_read("data/geodata/250_NUTS3.shp")

# some transformation
bird$year_cat <- factor(bird$year)
bird$fallow_sqrt <- sqrt(bird$fallow_per)
bird$fallow_std <- scale(bird$fallow_sqrt)
bird$agri_std <- scale(bird$agriculture)
bird$edge_std <- scale(bird$edge_mha)

# load models
m_list <- readRDS("model_output/fitted_model.rds")

######### Figure 1 #########
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
                    name = "Share of fallow land\non agricultural land (%)") +
  theme(legend.position = "right",
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# now plot the changes in fallow land percentage
gg_2 <- ggplot(subset(dd, year %in% c(2010, 2016))) +
  geom_sf(aes(fill = fallow_d_cut), size = 0.05) +
  facet_wrap(~year_cat2, ncol = 2) +
  scale_fill_brewer(palette = "PRGn", 
                    name = "Changes in share\n of fallow land on\nagricultural land (p. p.)") +
  theme(legend.position = "right",
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# put the two together
gg_a <- gridExtra::grid.arrange(gg_1, gg_2, nrow = 2)
ggsave("figures/01_fallow_figure.png", gg_a)

# NOT SHOWN: a table of fallow land percent across the years
bird %>%
  group_by(NUTS_code, NUTS_year) %>%
  summarise(fallow_per = mean(fallow_per)) %>%
  group_by(NUTS_code) %>%
  mutate(fallow_d = fallow_per - lag(fallow_per)) %>%
  mutate(year = substr(NUTS_year, 7, 10)) %>%
  pivot_longer(fallow_per:fallow_d) %>%
  group_by(year, name) %>%
  summarise(Mean = round(mean(value, na.rm = TRUE), 2),
            Median = round(median(value, na.rm = TRUE), 2),
            Q05 = quantile(value, probs = 0.05, na.rm = TRUE),
            Q95 = quantile(value, probs = 0.95, na.rm = TRUE),
            Min = min(value, na.rm = TRUE),
            Max = max(value, na.rm=TRUE)) -> dd_tab

dd_tab <- dd_tab[-1,]
dd_tab$name <- c("2007", "2010 - 2007", "2010", "2016 - 2010", "2016")
dd_tab$Q05_Q95 <- paste0(round(dd_tab$Q05, 2), " - ", round(dd_tab$Q95, 2))
dd_tab$min_max <- paste0(round(dd_tab$Min, 2), " - ", round(dd_tab$Max, 2))


write.csv(dd_tab[,c(2:4, 9, 10)], "fallow_table.csv", row.names = FALSE)

########## Figure 2 #######
# plot the fallow land coefficients values
# across different landscape complexity and
# for the different responses

# first make an helper function
plot_param <- function(model, edge_val, type = "field"){
  tmp <- posterior_samples(model, pars = "fallow")
  # take into account quadratic effect for richness
  if(type == "Species richness"){
    tmp2 <- sapply(edge_val, function(x) tmp[,1] + x * tmp[,2] + (x ** 2) * tmp[,3])
  } else {
    tmp2 <- sapply(edge_val, function(x) tmp[,1] + x * tmp[,2])
  }
  df_out <- as.data.frame(posterior_summary(tmp2, probs = c(0.025, 0.25, 0.75, 0.975)))
  
  df_out$type <- type
  df_out$edge <- edge_val
  df_out$param <- "avg"
 
  return(df_out)
}

# the gradient of edge values, 95% of the data
edge_val <- seq(quantile(bird$edge_std, probs = 0.025), 
                quantile(bird$edge_std, probs = 0.975), 
                length.out = 25)

# compute the values across the variables
df_param <- rbind(plot_param(m_list$rich, edge_val, "Species richness"),
                  plot_param(m_list$field, edge_val, "Field-breeders"),
                  plot_param(m_list$edge, edge_val, "Edge-breeders"),
                  plot_param(m_list$low, edge_val, "Foraging visitors"))

df_param$type <- factor(df_param$type, levels = c("Species richness", "Field-breeders",
                                                  "Edge-breeders", "Foraging visitors"))

# a second helper function to make the plots
plot_1 <- function(btype){
  dat <- subset(df_param, type == btype)
  
  if(btype == "Species richness"){
    ggplot(dat, aes(x=edge, y=Estimate, ymin=Q2.5, ymax=Q97.5)) +
      geom_hline(yintercept = 0, size = 0.5) +
      geom_ribbon(color="grey50", fill = NA, linetype = "dashed") +
      geom_line(size = 1.5) +
      #coord_flip() +
      scale_x_continuous(breaks = (c(25, 50, 75, 100)  - mean(bird$edge_mha)) / sd(bird$edge_mha),
                         labels = c(25, 50, 75, 100))+
      labs(x = "",
           y = "",
           title = btype)
  }
  else{
    ggplot(dat, aes(x=edge, y=Estimate, ymin=Q2.5, ymax=Q97.5)) +
      geom_hline(yintercept = 0, size = 0.5) +
      geom_ribbon(color="grey50", fill = NA, linetype = "dashed") +
      geom_line(size = 1.5) +
      ylim(c(-0.05, 0.45)) +
      #coord_flip() +
      scale_x_continuous(breaks = (c(25, 50, 75, 100)  - mean(bird$edge_mha)) / sd(bird$edge_mha),
                         labels = c(25, 50, 75, 100)) +
      labs(x = "",
           y = "",
           title = paste0(btype, " abundance"))
  }
}

# make the plots
dd_1 <- gridExtra::grid.arrange(plot_1("Species richness"),
                        plot_1("Field-breeders"),
                        plot_1("Edge-breeders"),
                        plot_1("Foraging visitors"),
                        bottom = "Landscape structural complexity (m edges / ha)",
                        left = "Effect of fallow land")


ggsave("figures/02_slope_fallow_land.png", dd_1, width = 8, height = 8)


########### Figure 3 ########
# predicted absolute and percent changes in birds
# given landscape complexity and some changes in fallows

# the fallow values considered:
# a. 4% corresponding to average situation in 2007
# b. 2% average value in 2016
# c. 4% target value in the post-2020 CAP
# d. 10% target value in the EU biodiversity strategy

# vector of fallow land values
fval <- (sqrt(c(4, 2, 4, 10)) - mean(bird$fallow_sqrt)) / sd(bird$fallow_sqrt)
# expand along the landscape complexity gradient
newdf <- expand.grid(fallow_std = fval,
                     edge_std = seq(quantile(bird$edge_std, probs = 0.025),
                                    quantile(bird$edge_std, probs = 0.975),
                                    length.out = 25),
                     agri_std = 0)

# an helper function computing predicted percent changes
plot_pch <- function(model, newdf, type = ""){
  pp <- posterior_epred(model, newdf, re.form = NA)
  tmp <- posterior_summary(((pp[,-c(seq(1, 100, by = 4))] - 
                              pp[,-c(seq(4, 100, by = 4))]) /
                              pp[,-c(seq(4, 100, by = 4))]) * 100)
  tmp2 <- cbind(newdf[-c(seq(1, 100, by = 4)),], tmp)
  tmp2$fallow_cat <- factor(tmp2$fallow_std, labels = c("Past changes", "Future changes:\nminimum requirements",
                                                        "Future changes:\nvoluntary additions"))
  tmp2$type <- type
  
  return(tmp2)
}
# an helper function computing predicted absolute changes
plot_ch <- function(model, newdf, type = ""){
  pp <- posterior_epred(model, newdf, re.form = NA)
  tmp <- posterior_summary(pp[,-c(seq(1, 100, by = 4))] - 
                             pp[,-c(seq(4, 100, by = 4))])
  tmp2 <- cbind(newdf[-c(seq(1, 100, by = 4)),], tmp)
  tmp2$fallow_cat <- factor(tmp2$fallow_std, 
                            labels = c("Past changes", 
                                       "Future changes:\nminimum requirements",
                                       "Future changes:\nvoluntary additions"))
  tmp2$type <- type

  return(tmp2)
}


# run the computation
d_ch <- rbind(plot_ch(m_list$rich, newdf, "Species richness"),
              plot_ch(m_list$edge, newdf, "Edge-breeders"),
              plot_ch(m_list$field, newdf, "Field-breeders"),
              plot_ch(m_list$low, newdf, "Foraging visitors")) 

d_ch$type <- factor(d_ch$type, levels = c("Species richness", "Field-breeders",
                                          "Edge-breeders", "Foraging visitors"))
# make the plot
gg_2 <- ggplot(d_ch, aes(x=edge_std, y=Estimate, ymin=Q2.5, ymax=Q97.5)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(fill=fallow_cat), color = NA, alpha=0.25) +
  geom_line(aes(color=fallow_cat), size = 1.5) +
  facet_wrap(~type, scales = "free_y") +
  scale_x_continuous(breaks = (c(25, 50, 75, 100) - mean(bird$edge_mha)) / sd(bird$edge_mha),
                     labels = c(25, 50, 75, 100)) +
  scale_fill_brewer(name = "Change in share\nof fallow land:", palette = "Accent") +
  scale_color_brewer(name = "Change in share\nof fallow land:", palette = "Accent") +
  labs(x = "Landscape structural complexity (m edges / ha)",
       y = "Predicted change in bird\nrichness and scaled abundance")


ggsave("figures/03_change_landscape.png", gg_2,
       width = 8.34, height = 5.81, units = "in")


############ Table 2 ########
# comparison of the strength of fallow land effect
# between the different bird categories

# levels of edges to test
edge_std = matrix(c(rep(1, 3), c(quantile(bird$edge_std, probs = 0.025),
                                 0, quantile(bird$edge_std, probs = 0.975))),
                    ncol = 3, byrow = TRUE)

# grab the posterior samples and compute the effect
# conditional on the landscape complexity
pp_field <- as.matrix(posterior_samples(m_list$field, pars = c("fallow_std", "edge_std:fallow_std"))) %*% edge_std
pp_edge <- as.matrix(posterior_samples(m_list$edge, pars = c("fallow_std", "edge_std:fallow_std"))) %*% edge_std
pp_low <- as.matrix(posterior_samples(m_list$low, pars = c("fallow_std", "edge_std:fallow_std"))) %*% edge_std

# run the comparison
comp_groups <- data.frame(comparison = c("edge_field", "edge_low", "field_low"),
                          low = c(sum(pp_edge[,1] > pp_field[,1]) / 4000,
                                  sum(pp_edge[,1] > pp_low[,1]) / 4000,
                                  sum(pp_field[,1] > pp_low[,1]) / 4000),
                          middle = c(sum(pp_edge[,2] > pp_field[,2]) / 4000,
                                     sum(pp_edge[,2] > pp_low[,2]) / 4000,
                                     sum(pp_field[,2] > pp_low[,2]) / 4000),
                          high = c(sum(pp_edge[,3] > pp_field[,3]) / 4000,
                                   sum(pp_edge[,3] > pp_low[,3]) / 4000,
                                   sum(pp_field[,3] > pp_low[,3]) / 4000))


####### Figure 4 #########
# predicted changes at district level

# merge the district in one geoms
dshape %>%
  group_by(NUTS_CODE) %>%
  summarise() -> dsh

dsh <- st_cast(dsh, "MULTIPOLYGON")

names(dsh)[1] <- "NUTS_code"

# add scenario of 4% fallow on arable land
tmp <- district[district$year == 2016,]
tmp$NUTS_year <- gsub("_2016", "_2020", tmp$NUTS_year)
tmp$year <- 2020
tmp2 <- tmp # for 10%
tmp$SEFA <- with(tmp, ifelse((SEFA / ARAB) < 0.04, 0.04 * ARAB, SEFA))
tmp$fallow_per <- with(tmp, (SEFA / (ARAB + GRAS)) * 100)

# now 10%
tmp2$SEFA <- with(tmp2, ifelse((SEFA / ARAB) < 0.1, 0.1 * ARAB, SEFA))
tmp2$fallow_per <- with(tmp2, (SEFA / (ARAB + GRAS)) * 100)

# comine with other data
district_pred <- rbind(district[district$year %in% c(2007, 2016), c("NUTS_code", "NUTS_year", "fallow_per", "year")],
                       tmp[,c("NUTS_code", "NUTS_year", "fallow_per", "year")])
# for 10%
district_pred10 <- rbind(district[district$year %in% c(2007, 2016), c("NUTS_code", "NUTS_year", "fallow_per", "year")],
                       tmp2[,c("NUTS_code", "NUTS_year", "fallow_per", "year")])

# set the other covariates to their average values
district_pred$agri_std <- district_pred10$agri_std <- 0
district_pred$edge_std <- district_pred10$edge_std <- 0

# square root and std the fallow values
district_pred$fallow_sqrt <- sqrt(district_pred$fallow_per)
district_pred$fallow_std <- (district_pred$fallow_sqrt - mean(bird$fallow_sqrt)) / sd(bird$fallow_sqrt)

district_pred <- arrange(district_pred, NUTS_code, year)

district_pred10$fallow_sqrt <- sqrt(district_pred10$fallow_per)
district_pred10$fallow_std <- (district_pred10$fallow_sqrt - mean(bird$fallow_sqrt)) / sd(bird$fallow_sqrt)

district_pred10 <- arrange(district_pred10, NUTS_code, year)

# a helper function
plot_district <- function(model, pred, type = "Field-breeding bird"){
  tt <- posterior_epred(model, newdata = pred,
                        allow_new_levels = TRUE)

  dpred <- cbind(as.data.frame(district_pred), posterior_summary(tt))
  dpred %>%
    arrange(NUTS_code, year) %>%
    group_by(NUTS_code) %>%
    mutate(dval = Estimate - lag(Estimate),
           dval2 = Estimate - lag(Estimate, 2)) %>%
    pivot_longer(dval:dval2, names_to = "time",
                 values_to = "changes") %>%
    filter(!is.na(changes)) %>%
    unite("code", time, year) %>%
    full_join(dsh, by = "NUTS_code") %>%
    mutate(type = type) -> dd
  
  # remove Flensburg
  dd$changes[dd$NUTS_code == "DEF01"] <- NA
  #dd$dval2[dd$NUTS_code == "DEF01"] <- NA
  dd <- st_as_sf(dd)
  return(dd)
}

# run through the response variables for 4%,
# same can be done for the 10% using district_pred10
edge_change <- plot_district(m_list$edge, district_pred, "Edge-breeders")
edge_change$code <- factor(edge_change$code,
                           labels = c("Edge-breeders:\n2007 - 2016",
                                      "Edge-breeders:\n2016 - post-2020 (4%)",
                                      "Edge-breeders:\n2007 - post-2020 (4%)"))

r_edge <- ggplot() +
  geom_sf(data=edge_change, aes(fill=changes), size = 0.05) +
  facet_wrap(~code) +
  scale_fill_gradient2(na.value = "grey50",
                       name = "Predicted\nchanges") +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())


field_change <- plot_district(m_list$field, district_pred, "Field-breeders")
field_change$code <- factor(field_change$code, 
                            labels = c("Field-breeders:\n2007 - 2016",
                                       "Field-breeders:\n2016 - post-2020 (4%)",
                                       "Field-breeders:\n2007 - post-2020 (4%)"))

r_field <- ggplot() +
  geom_sf(data=field_change, aes(fill=changes), size = 0.05) +
  facet_wrap(~code) +
  scale_fill_gradient2(na.value = "grey50",
                       name = "Predicted\nchanges") +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())

low_change <- plot_district(m_list$low, district_pred, "FF")
low_change$code <- factor(low_change$code, 
                          labels = c("Foraging visitors:\n2007 - 2016",
                                     "Foraging visitors:\n2016 - post-2020 (4%)",
                                     "Foraging visitors:\n2007 - post-2020 (4%)"))

r_low <- ggplot() +
  geom_sf(data=low_change, aes(fill=changes), size = 0.05) +
  facet_wrap(~code) +
  scale_fill_gradient2(na.value = "grey50",
                       name = "Predicted\nchanges") +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())

# put the abundance together
gg_a <- gridExtra::grid.arrange(r_rich, r_field, r_edge, r_low, nrow = 4)
ggsave("figures/0X_abundance_changes.png", gg_a, width = 12, height = 16)

# species richness changes in one plot
rich_change <- plot_district(m_list$rich, district_pred, "")
rich_change$code <- factor(rich_change$code, 
                           labels = c("2007 - 2016",
                                      "2016 - post-2020 (4%)",
                                      "2007 - post-2020 (4%)"))

r_rich <- ggplot() +
  geom_sf(data=rich_change, aes(fill=changes), size = 0.05) +
  facet_wrap(~code) +
  scale_fill_gradient2(na.value = "grey50",
                       name = "Predicted\nchanges") +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())

ggsave("figures/04_spat_rich.png", r_rich)


# appendix figs, effect of landscape complexity
land_rich <- conditional_effects(m_list$rich, "edge_std:fallow_std")
pp_rich <- plot(land_rich)[[1]] +
  scale_x_continuous(breaks = (c(0.05, 0.15, 0.25) - mean(bird$edge_kmha)) / sd(bird$edge_kmha),
                     labels = c(50, 150, 250)) +
  labs(x="",
       y="Species richness") +
  scale_fill_brewer(name = "Share of fallow\nland (%):",
                    palette = "Accent",
                    labels = c("5.10", "2.50", "0.82")) +
  scale_color_brewer(name = "Share of fallow\nland (%):",
                    palette = "Accent",
                    labels = c("5.10", "2.50", "0.82"))

land_field <- conditional_effects(m_list$field, "edge_std:fallow_std")
pp_field <- plot(land_field)[[1]] +
  scale_x_continuous(breaks = (c(0.05, 0.15, 0.25) - mean(bird$edge_kmha)) / sd(bird$edge_kmha),
                     labels = c(50, 150, 250)) +
  labs(x="",
       y="Field-breeders") +
  scale_fill_brewer(name = "Share of fallow\nland (%):",
                    palette = "Accent",
                    labels = c("5.10", "2.50", "0.82")) +
  scale_color_brewer(name = "Share of fallow\nland (%):",
                     palette = "Accent",
                     labels = c("5.10", "2.50", "0.82"))

land_edge <- conditional_effects(m_list$edge, "edge_std:fallow_std")
pp_edge <- plot(land_edge)[[1]] +
  scale_x_continuous(breaks = (c(0.05, 0.15, 0.25) - mean(bird$edge_kmha)) / sd(bird$edge_kmha),
                     labels = c(50, 150, 250)) +
  labs(x="",
       y="Edge-breeders") +
  scale_fill_brewer(name = "Share of fallow\nland (%):",
                    palette = "Accent",
                    labels = c("5.10", "2.50", "0.82")) +
  scale_color_brewer(name = "Share of fallow\nland (%):",
                     palette = "Accent",
                     labels = c("5.10", "2.50", "0.82"))

land_low <- conditional_effects(m_list$low, "edge_std:fallow_std")
pp_low <- plot(land_low)[[1]] +
  scale_x_continuous(breaks = (c(0.05, 0.15, 0.25) - mean(bird$edge_kmha)) / sd(bird$edge_kmha),
                     labels = c(50, 150, 250)) +
  labs(x="",
       y="Foraging visitors") +
  scale_fill_brewer(name = "Share of fallow\nland (%):",
                    palette = "Accent",
                    labels = c("5.10", "2.50", "0.82")) +
  scale_color_brewer(name = "Share of fallow\nland (%):",
                     palette = "Accent",
                     labels = c("5.10", "2.50", "0.82"))

gg_land <- ggpubr::ggarrange(pp_rich, pp_field, pp_edge, pp_low,
                             ncol = 2, nrow = 2, common.legend = TRUE)

gg_land2 <- ggpubr::annotate_figure(gg_land, 
                                    bottom = "Landscape structural complexity (m edges per ha)")
ggsave("figures/0X_landscape_eff.png", gg_land2)