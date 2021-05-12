# plotting script for the fallow land data

setwd("~/PostDoc_Thunen/ASE-stuff/")

# load libraries
library(tidyverse)
library(sf)
library(brms)

# load data
bird <- read.csv("data/preprocessed/bird_fallow_v8.csv")
#bird$land <- substr(bird$routcode, 1, 2)

# some transformation
bird$year_cat <- factor(bird$year)
bird$fallow_sqrt <- sqrt(bird$fallow_ha)
bird$fallow_std <- scale(bird$fallow_sqrt)
bird$edge_std <- scale(bird$edge_m)
bird$urban_std <- scale(bird$urban)
bird$forest_std <- scale(bird$forest)
bird$bkr <- factor(bird$bkr)

# load models
m_list <- readRDS("model_output/fitted_model.rds")

# plot the fallow land coefficients values
# across different landscape complexity and
# for the different responses
plot_param <- function(model, type = "field"){
  tmp <- posterior_samples(model, pars = "fallow")
  # effect for non-zeroes
  tmp2 <- cbind(tmp[,1] - tmp[,2], tmp[,1], tmp[,1] + tmp[,2])
  df_out <- as.data.frame(posterior_summary(tmp2, probs = c(0.025, 0.25, 0.75, 0.975)))
  
  df_out$type <- type
  df_out$cpx <- factor(c("low", "medium", "high"), levels = c("low", "medium", "high"))
  df_out$param <- "avg"
 
  return(df_out)
}

df_param <- rbind(plot_param(m_list$field, "Field-breeding bird"),
                  plot_param(m_list$edge, "Edge-breeding bird"),
                  plot_param(m_list$low, "Non-fallow breeding bird"))


dd1 <- ggplot(df_param, aes(x=type, y=Estimate, color = cpx, group = cpx)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.05) +
  geom_linerange(aes(ymin = Q2.5, ymax = Q97.5), size = 0.5, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = Q25, ymax = Q75), size = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_color_brewer(palette = "Accent",
                     name = "Structural complexity\nof the landscape : ") +
  labs(x = "",
       y = "Effect of fallow land")
  
# richness
df_param2 <- plot_param(m_list$rich, "Bird richness")

dd2 <- ggplot(df_param2, aes(x=type, y=Estimate, color = cpx, group = cpx)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.05) +
  geom_linerange(aes(ymin = Q2.5, ymax = Q97.5), size = 0.5, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = Q25, ymax = Q75), size = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_color_brewer(palette = "Accent",
                     name = "Structural complexity\nof the landscape : ") +
  labs(x = "",
       y = "Effect of fallow land")

gg_dd <- ggpubr::ggarrange(dd1, dd2, nrow = 2, common.legend = TRUE, heights = c(2, 1))
ggsave("figures/slope_fallow_land.png", gg_dd)


# derive first model prediction along the gradient
# of fallow land values for different landscape complexities
# and controlling for all other factors (year, bkr and agricultural surface)

# predictions are derived using the function brms::conditional_effects

## a helper function
plot_1_mu <- function(model, type = "field"){
  df <- conditional_effects(model, effects = "fallow_std:edge_std", dpar = "mu")
  # put back the fallow area values in the orginal scale (ha)
  df[[1]]$effect1__ <- (df[[1]]$effect1__ * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)) ** 2
  # turn the edge values into categories
  df[[1]]$cpx <- factor(df[[1]]$effect2__, labels = c("low", "medium", "high"),
                        levels = c(-1, 0, 1))
  df[[1]]$type <- type
  # re-naming the response variable column
  names(df[[1]])[3] <- "resp"
  
  df_out <- arrange(df[[1]], cpx)
  
  return(df_out)
}
# for the zero-outcome
plot_1_hu <- function(model, type = "field"){
  df <- conditional_effects(model, effects = "fallow_std:edge_std", dpar = "hu")
  # put back the fallow area values in the orginal scale (ha)
  df[[1]]$effect1__ <- (df[[1]]$effect1__ * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)) ** 2
  # turn the edge values into categories
  df[[1]]$cpx <- factor(df[[1]]$effect2__, labels = c("low", "medium", "high"),
                        levels = c(-1, 0, 1))
  df[[1]]$type <- type
  # re-naming the response variable column
  names(df[[1]])[3] <- "resp"
  # turn into prob of presence
  df[[1]]$estimate__ <- 1 - df[[1]]$estimate__
  df[[1]]$lower__ <- 1 - df[[1]]$lower__
  df[[1]]$upper__ <- 1 - df[[1]]$upper__
  
  
  df_out <- arrange(df[[1]], cpx)
  
  return(df_out)
}

# combine the groups
pl1 <- rbind(plot_1_mu(m_list$field, "Field-breeding bird"),
             plot_1_mu(m_list$edge, "Edge-breeding bird"),
             plot_1_mu(m_list$low, "Non-fallow breeding bird"))
# the plot
gg_plot1 <- ggplot(pl1, aes(x=effect1__, y = estimate__, ymin = lower__, ymax = upper__)) +
  geom_ribbon(aes(fill = cpx), alpha = 0.2, color = NA) +
  geom_line(aes(color = cpx)) +
  facet_wrap(~ type) +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  guides(fill = guide_legend(title = "Structural complexity\nof the landscape : "),
         color = guide_legend(title = "Structural complexity\nof the landscape : ")) +
  labs(x = "Area of fallow in 1km buffer (ha)",
       y = "Predicted scaled abundance (with 95% CrI)") +
  theme(panel.background = element_blank(),
        axis.line = element_line())


ggsave("figures/fallow_effect_abund.png", gg_plot1)

# combine the groups
pl1_hu <- rbind(plot_1_hu(m_list$field, "Field-breeding bird"),
             plot_1_hu(m_list$edge, "Edge-breeding bird"),
             plot_1_hu(m_list$low, "Non-fallow breeding bird"))
# the plot
gg_plot1_hu <- ggplot(pl1_hu, aes(x=effect1__, y = estimate__, ymin = lower__, ymax = upper__)) +
  geom_ribbon(aes(fill = cpx), alpha = 0.2, color = NA) +
  geom_line(aes(color = cpx)) +
  facet_wrap(~ type) +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  guides(fill = guide_legend(title = "Structural complexity\nof the landscape : "),
         color = guide_legend(title = "Structural complexity\nof the landscape : ")) +
  labs(x = "Area of fallow in 1km buffer (ha)",
       y = "Probability of presence (with 95% CrI)") +
  theme(panel.background = element_blank(),
        axis.line = element_line())


ggsave("figures/fallow_effect_abund_hu2.png", gg_plot1_hu)

# for species richness
df_plot11 <- plot_1_mu(m_list$rich, type = "Species richness")

gg_plot11 <- ggplot(df_plot11, aes(x=effect1__, y = estimate__, ymin = lower__, ymax = upper__)) +
  geom_ribbon(aes(fill = cpx), alpha = 0.2, color = NA) +
  geom_line(aes(color = cpx)) +
  guides(fill = guide_legend(title = "Structural complexity\nof the landscape : "),
         color = guide_legend(title = "Structural complexity\nof the landscape : ")) +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  labs(x = "Area of fallow in 1km buffer (ha)",
       y = "Predicted bird richness (with 95% CrI)") +
  theme(panel.background = element_blank(),
        axis.line = element_line())



ggsave("figures/fallow_effect_rich_wgrass.png", gg_plot11)

## third set of figures predicting abundance and richness
## for various amount of fallow land in different years
## based on the observed data

# first compute lm to link amount in the past to
# amount in the present
bird %>%
  pivot_wider(id_cols = routcode, names_from = year_cat, 
              values_from = fallow_ha, names_prefix = "y") %>%
  mutate(d1 = y2010 - y2007,
         d2 = y2016 - y2010) -> ff

# model the first transition
mf <- lm(d1 ~ y2007, ff)
# predictions for changes and new values
f2007 <- seq(min(ff$y2007,na.rm=TRUE), 20, length.out = 20)
d2010 <- predict(mf, newdata = data.frame(y2007 = f2007))
f2010 <- d2010 + f2007

# for the plotting
dpred <- data.frame(y2007 = f2007, d1 = d2010)
gg_sc <- ggplot(ff, aes(x=y2007, y=d1)) +
  geom_point(alpha = 0.1) +
  stat_smooth(method = "lm") +
  geom_point(data=dpred, color = "red", size = 3) +
  labs(x = "Estimated area of fallow land in 2007 (in ha)",
       y = "Changes in estimated fallow land area between 2010 and 2007 (in ha)")

ggsave("figures/fig_prediction1.png", gg_sc)

# put together in a df
new_fallow <- data.frame(fallow_ha = c(f2007, f2010),
                         year = rep(c(2007, 2010), each = 20))

new_fallow$fallow_std <- (sqrt(new_fallow$fallow_ha) - mean(bird$fallow_sqrt)) / sd(bird$fallow_sqrt)

# expand with edges
edge_q <- c(-1, 0, 1)
new_df <- data.frame(year_cat = factor(2007),
                     agri_std = 0,
                     edge_std = rep(edge_q, each = 40),
                     fallow_std = rep(new_fallow$fallow_std, 3))

# the second transition
# model the second transition
mf2 <- lm(d2 ~ y2010, ff)

# predictions for changes and new values
f2010 <- seq(min(ff$y2010,na.rm=TRUE), 15, length.out = 20)
d2016 <- predict(mf2, newdata = data.frame(y2010 = f2010))
f2016 <- d2016 + f2010

# for the plotting
dpred <- data.frame(y2010 = f2010, d2 = d2016)
gg_sc <- ggplot(ff, aes(x=y2010, y=d2)) +
  geom_point(alpha = 0.1) +
  stat_smooth(method = "lm") +
  geom_point(data=dpred, color = "red", size = 3) +
  labs(x = "Estimated area of fallow land in 2010 (in ha)",
       y = "Changes in estimated fallow land area between 2016 and 2010 (in ha)")

ggsave("figures/fig_prediction2.png", gg_sc)

# put together in a df
new_fallow <- data.frame(fallow_ha = c(f2010, f2016),
                         year = rep(c(2007, 2010), each = 20))

new_fallow$fallow_std <- (sqrt(new_fallow$fallow_ha) - mean(bird$fallow_sqrt)) / sd(bird$fallow_sqrt)

# expand with edges
edge_q <- c(-1, 0, 1)
new_df2 <- data.frame(year_cat = factor(2007),
                     agri_std = 0,
                     edge_std = rep(edge_q, each = 40),
                     fallow_std = rep(new_fallow$fallow_std, 3))

# a helper function to fit this
# helper function
plot_2 <- function(model, newdata, delta_val, type = "field"){
  pp <- posterior_samples(model, pars = "year_cat")
  # model predictions
  predd <- posterior_epred(model, newdata = newdata, re.form = NA)
  # average out year-effect
  predd <- predd + apply(pp, 1, function(x) mean(c(0, x)))
  # make ratios
  pred_fr <- as.data.frame(posterior_summary(predd[,c(21:40, 61:80, 101:120)] / predd[,c(1:20, 41:60, 81:100)]))
  pred_fr$edge_std <- newdata$edge_std[c(21:40, 61:80, 101:120)]
  pred_fr$delta_fallow <- delta_val
  pred_fr$type <- type
  
  return(pred_fr)
}

# run this through the different groups
df_plot2.1 <- rbind(plot_2(m_list$edge, new_df, d2010, type = "Edge-breeding birds"),
                  plot_2(m_list$field, new_df, d2010, type = "Field-breeding birds"),
                  plot_2(m_list$low, new_df, d2010, type = "Non fallow-breeding birds"))
#df_plot2.1$transition <- "Transition 1:\nabolishment of set-aside"

df_plot2.2 <- rbind(plot_2(m_list$edge, new_df2, d2016, type = "Edge-breeding birds"),
                    plot_2(m_list$field, new_df2, d2016, type = "Field-breeding birds"),
                    plot_2(m_list$low, new_df2, d2016, type = "Non fallow-breeding birds"))
#df_plot2.2$transition <- "Transition 2:\nestablishment of EFA"

#df_plot2 <- rbind(df_plot2.1, df_plot2.2)
df_plot2.1$edge_cat <- factor(df_plot2.1$edge_std,
                              levels = c(1, 0, -1),
                            labels = c("high", "medium", "low"))
df_plot2.2$edge_cat <- factor(df_plot2.2$edge_std,
                              levels = c(1, 0, -1),
                              labels = c("high", "medium", "low"))

gg_1 <- ggplot(df_plot2.1, aes(x=delta_fallow, y=Estimate, ymin=Q2.5, ymax=Q97.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1.2) +
  geom_ribbon(aes(fill=edge_cat), alpha=0.1, color=NA) +
  geom_line(aes(color=edge_cat)) +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  facet_wrap(~type) +
  
  scale_y_continuous(breaks = c(0.90, 1, 1.1),
                     labels = c("-10%", "0%", "+10%")) +
  labs(x = "Changes in fallow land area between the 2010 and the 2007 agricultural census (in ha)",
       y = "Predicted changes\nin bird abundance (with 95% CrI)") +
  guides(color=guide_legend(title = "Structural complexity\nof the landscape"),
         fill = guide_legend(title = "Structural complexity\nof the landscape")) +
  theme(axis.line = element_line(),
        panel.background = element_blank())

ggsave("figures/bird_changes_transition1.png", gg_1)


gg_2 <- ggplot(df_plot2.2, aes(x=delta_fallow, y=Estimate, ymin=Q2.5, ymax=Q97.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1.2) +
  geom_ribbon(aes(fill=edge_cat), alpha=0.1, color=NA) +
  geom_line(aes(color=edge_cat)) +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  facet_wrap(~type) +
  scale_y_continuous(breaks = c(0.9, 1.0, 1.1, 1.2, 1.3),
                     labels = c("-10%", "0%", "+10%", "+20%", "+30%")) +
  labs(x = "Changes in fallow land area between the 2016 and the 2010 agricultural census (in ha)",
       y = "Predicted changes\nin bird abundance (with 95% CrI)") +
  guides(color=guide_legend(title = "Structural complexity\nof the landscape"),
         fill = guide_legend(title = "Structural complexity\nof the landscape")) +
  theme(axis.line = element_line(),
        panel.background = element_blank())

ggsave("figures/bird_changes_transition2.png", gg_2)

# other version with both together
gg_3 <- ggpubr::ggarrange(gg_1, gg_2, common.legend = TRUE, nrow = 2)
ggsave("figures/bird_changes_together_wgrass.png", gg_3)

# now for richness
df_plotr <- rbind(plot_2(m_list$rich, new_df, d2010, type = "Transition 1: abolishment\nof mandatory set-aside"),
      plot_2(m_list$rich, new_df2, d2016, type = "Transition 2: establishment\nof EFA"))
df_plotr$edge_cat <- factor(df_plotr$edge_std,
                            levels = c(1, 0, -1),
                            labels = c("high", "medium", "low"))

gg_r <- ggplot(df_plotr, aes(x=delta_fallow, y=Estimate, ymin=Q2.5, ymax=Q97.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1.2) +
  geom_ribbon(aes(fill=edge_cat), alpha=0.1, color=NA) +
  geom_line(aes(color=edge_cat)) +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  facet_wrap(~type) +
  scale_y_continuous(breaks = c(0.95, 1.0, 1.05),
                     labels = c("-5%", "0%", "+5%")) +
  labs(x = "Changes in fallow land area between the agricultural census (in ha)",
       y = "Predicted changes in bird richness (with 95% CrI)") +
  guides(color=guide_legend(title = "Structural complexity\nof the landscape"),
         fill = guide_legend(title = "Structural complexity\nof the landscape")) +
  theme(axis.line = element_line(),
        panel.background = element_blank())

ggsave("figures/rich_changes_transition_wgrass.png", gg_r)

##############
# other results: post. proba of effect being stronger in the different groups

# levels of edges to test
edge_std = matrix(c(rep(1, 3), c(-1, 0, 1)),
                    ncol = 3, byrow = TRUE)

pp_field <- as.matrix(posterior_samples(m_list$field, pars = c("fallow_std", "edge_std:fallow_std"))) %*% edge_std
pp_edge <- as.matrix(posterior_samples(m_list$edge, pars = c("fallow_std", "edge_std:fallow_std"))) %*% edge_std
pp_low <- as.matrix(posterior_samples(m_list$low, pars = c("fallow_std", "edge_std:fallow_std"))) %*% edge_std

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

