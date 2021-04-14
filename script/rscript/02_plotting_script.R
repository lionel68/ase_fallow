# plotting script for the fallow land data

setwd("~/PostDoc_Thunen/ASE-stuff/")

# load libraries
library(tidyverse)
library(sf)
library(brms)

# load data
bird <- read.csv("data/preprocessed/bird_fallow_v5.csv")
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
m_list <- readRDS("model_output/fitted_model_new.rds")


# new model predictions
## first derive predictions across the gradients of fallow land
## in different edge cat
pred_fallow <- expand.grid(fallow_std = seq(min(bird$fallow_std),
                                            max(bird$fallow_std),
                                            length.out = 10),
                           edge_std = quantile(bird$edge_std, 
                                               probs = c(0.25, 0.5, 0.75)),
                           urban_std = 0, forest_std = 0, year_cat = unique(bird$year_cat))
# helper function
plot_1 <- function(model, newdata, type = "field"){
  # expected model values
  predd <- posterior_epred(model, newdata, re.form = NA)
  # grab posterior samples to average over the years
  # and turn back to response scale
  predd <- (predd[,1:30] + predd[,31:60] + predd[,61:90]) / 3
  # compute quantiles
  pred_df <- cbind(newdata[1:30,], posterior_summary(predd))
  pred_df$edge_cat <- factor(rep(c("low", "middle", "high"), each = 10),
                             levels = c("low", "middle", "high"))
  pred_df$fallow_ha <- (pred_df$fallow_std * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)) ** 2
  pred_df$type <- type
  
  return(pred_df)
}

# run this through the groups
df_plot1 <- rbind(plot_1(m_list$edge, pred_fallow, type = "Edge-breeding birds"),
                  plot_1(m_list$field, pred_fallow, type = "Field-breeding birds"),
                  plot_1(m_list$low, pred_fallow, type = "Non fallow-breeding birds"))

# first plot
gg_plot1 <- ggplot(df_plot1, aes(x=fallow_ha, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_ribbon(aes(fill=edge_cat), color = NA, alpha = 0.15) +
  geom_line(aes(color=edge_cat)) +
  facet_grid(~type) +
  guides(fill = guide_legend(title = "Structural complexity\nof the landscape : "),
         color = guide_legend(title = "Structural complexity\nof the landscape : ")) +
  labs(x = "Area of fallow in 1km buffer (ha)",
       y = "Predicted scaled abundance (with 95% CrI)") +
  theme(panel.background = element_blank(),
        axis.line = element_line())

ggsave("figures/fallow_effect_scaled_new.png", gg_plot1)

# for species richness
df_plot11 <- plot_1(m_list$rich, pred_fallow, type = "Species richness")

gg_plot11 <- ggplot(df_plot11, aes(x=fallow_ha, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_ribbon(aes(fill=edge_cat), color = NA, alpha = 0.15) +
  geom_line(aes(color=edge_cat)) +
  guides(fill = guide_legend(title = "Structural complexity\nof the landscape : "),
         color = guide_legend(title = "Structural complexity\nof the landscape : ")) +
  labs(x = "Area of fallow in 1km buffer (ha)",
       y = "Predicted richness (with 95% CrI)") +
  theme(panel.background = element_blank(),
        axis.line = element_line())

ggsave("figures/fallow_effect_rich.png", gg_plot11)

### explore how fallow land area changed over time
scl <- function(x){
  (x-mean(x)) / sd(x)
}

# get the full lists of plots
mhb_fallow <- read.csv("data/landusedata/mhb_1kmbuffer_fallow2.csv")

# get an idea of how much fallow land changed in ha
mhb_fallow %>%
  group_by(routcode) %>%
  mutate(fallow_delta = fallow_ha - lag(fallow_ha)) %>%
  # filter(year != 2007) %>%
  group_by(year) %>%
  summarise(M = median(fallow_delta, na.rm=TRUE),
            M_orig = median(fallow_ha),
            LCI = quantile(fallow_delta, probs = 0.1, na.rm=TRUE),
            LCI_orig = quantile(fallow_ha, probs =0.1),
            UCI = quantile(fallow_delta, probs = 0.9, na.rm=TRUE),
            UCI_orig = quantile(fallow_ha, probs=0.9, na.rm=TRUE)) -> fallow_ref


# checking how fallow land area varied in regard to previous area
mhb_fallow %>%
  group_by(routcode) %>%
  mutate(fallow_delta = lead(fallow_ha) - fallow_ha) -> tmp

# plot distribution of fallow land changes
gg_fdelta <- ggpubr::ggdensity(tmp, x = "fallow_delta", fill = "year",
                  rug = TRUE, color = "year") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Changes in fallow land area between time points (ha)")


ggsave("figures/fallow_delta_distrib.png", gg_fdelta)

# re-format fallow_ref
fallow_ref %>%
  pivot_longer(cols = 2:7) %>%
  mutate(varn = rep(c("fallow_delta", "fallow_ha"), 9)) %>%
  mutate(quant = rep(rep(c("50%", "10%", "90%"), each = 2), 3)) %>%
  pivot_wider(id_cols = c(1, 5), names_from = "varn", values_from = "value") -> ref_dd

# some more wrangling around
ref_dd$fallow_delta <- c(-2.05, 0.013, -8.19, 0.0588, 3.55, -2.35, NA, NA, NA)
ref_dd$quant <- rep(paste("Scenario", c(2, 1, 3)), 3)

gg_f <- ggplot(subset(tmp, year != 2016), aes(x=fallow_ha, y=fallow_delta)) +
  geom_point(alpha=0.1) +
  geom_point(data=subset(ref_dd, year != 2016), color = "grey50", size = 5) +
  geom_point(data = subset(ref_dd, year != 2016), aes(color = quant), size = 2) +
  scale_color_discrete(name = "Scenarios : ") +
  labs(x="Fallow land area before (ha)",
       y = "Change in fallow land area (ha)") +
  stat_smooth(method="lm", se = FALSE) +
  facet_grid(~ year, labeller = as_labeller(c(`2007` = "Transition 1\n2007 - 2010",
                                              `2010` = "Transition 2\n2010 - 2016")))

ggsave("figures/fallow_temporal.png", gg_f)


## the idea of the derived prediction  from the data explorationis that
## the largest decline in fallow land area
## between 2010-2007 and 2016-2010 were observed around Mhb plots with large quantities
## of fallow land area. The smallest decline or the slight increases were observed 
## in the plots with low amount of fallow lands
## below are the construction of the relevant data frame to derive the predictions
## three categories of fallow land are used: 10, 50 and 90% of both the fallow land
## in 2007 and 2010 as well as the computed changes (at the route level) between
## 2010 and 2007 and 2016 and 2010.
# construct the prediction df

# first derive the values
# deriving the scaled square rooted fallow land values to derive the model predictions
mhb_fallow %>%
  mutate(fallow_std = scl(sqrt(fallow_ha))) %>%
  group_by(routcode) %>%
  mutate(fallow_delta = fallow_std - lag(fallow_std)) %>%
  group_by(year) %>%
  summarise(M_orig = median(fallow_std),
            M = median(fallow_delta, na.rm=TRUE),
            LCI_orig = quantile(fallow_std, probs = 0.1),
            LCI = quantile(fallow_delta, probs = 0.1, na.rm=TRUE),
            UCI_orig = quantile(fallow_std, probs = 0.9),
            UCI = quantile(fallow_delta, probs = 0.9, na.rm=TRUE)) -> fallow_dd

# put these in a dataframe
pred_f <- expand.grid(edge_std = quantile(bird$edge_std, 
                                          probs = c(0.1, 0.5, 0.9)),
                      fallow_std = c(c(1.80, 0.426, -1.12, 0.919, -0.415, -1.35),
                                     c(1.80-1.80, 0.426-0.612, -1.12+0.0191, 0.919-0.648, -0.415-0.0627, -1.35+1.17)),
                      urban_std = 0, forest_std = 0, year_cat = 2007)
# the final data frame
predd_f <- data.frame(edge_cat = rep(c("low", "middle", "high"), 6),
                      fallow_cat = rep(c("Scenario 3", "Scenario 2", "Scenario 1"), each = 3),
                      transition = rep(1:2, each = 9))

# helper function
plot_2 <- function(model, newdata, newdata2, type = "field"){
  pp <- posterior_samples(model, pars = "year_cat")
  # model predictions
  predd <- posterior_epred(model, newdata = newdata, re.form = NA)
  # average out year-effect
  predd <- predd + apply(pp, 1, function(x) mean(c(0, x)))
  # make ratios
  pred_fr <- predd[,19:36] / predd[,1:18]
  predd_f <- cbind(newdata2, posterior_summary(pred_fr))
  predd_f$type <- type
  
  return(predd_f)
}

# run this through the different groups
df_plot2 <- rbind(plot_2(m_list$edge, pred_f, predd_f, type = "Edge-breeding birds"),
                  plot_2(m_list$field, pred_f, predd_f, type = "Field-breeding birds"),
                  plot_2(m_list$low, pred_f, predd_f, type = "Non fallow-breeding birds"))

df_plot2$x <- rep(1:3, each = 18) + rep(c(-0.15, 0, 0.15), times = 54/3) 
df_plot2$edge_cat <- factor(df_plot2$edge_cat, levels = c("low", "middle", "high"))
df_plot2$fallow_cat <- factor(df_plot2$fallow_cat, labels = c("Scenario 1: low\namount before,\nlow/positice changes",
                                                              "Scenario 2: typical\namount before,\ntypical changes",
                                                              "Scenario 3: large\namount before,\nlarge negative changes"))
df_plot2$transition <- factor(df_plot2$transition, labels = c("Transition 1: abolishment\nof mandatory fallow",
                                                              "Transition 2: establishment\nof EFA"))


gg_plot2 <- ggplot(df_plot2, aes(x=x, y=Estimate, ymin = Q2.5, ymax = Q97.5, color = edge_cat)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1.05) +
  geom_linerange() +
  geom_point() +
  facet_grid(transition ~ fallow_cat) +
  scale_x_continuous(breaks = 1:3, labels = c("Edge", "Field", "Non-\nfallow")) +
  scale_y_continuous(breaks = seq(0.8, 1.2, length.out = 5),
                     labels = c("-20%", "-10%", "No change", "+10%", "+20%")) +
  scale_color_discrete(name = "Structural complexity\nof the landscape:") +
  labs(x = "Bird group",
       y = "Predicted changes in bird scaled abundance (95% CrI)")

ggsave("figures/bird_temporal_scaled.png", gg_plot2)

df_plot22 <- plot_2(m_rich, pred_f, predd_f, type = "Species richness")
df_plot22$edge_cat <- factor(df_plot22$edge_cat, levels = c("low", "middle", "high"))
df_plot22$fallow_cat <- factor(df_plot22$fallow_cat, labels = c("Scenario 1: low\namount before,\nlow/positice changes",
                                                              "Scenario 2: typical\namount before,\ntypical changes",
                                                              "Scenario 3: large\namount before,\nlarge negative changes"))
df_plot22$transition <- factor(df_plot22$transition, labels = c("Transition 1: abolishment\nof mandatory fallow",
                                                              "Transition 2: establishment\nof EFA"))

gg_rich <- ggplot(df_plot22, aes(x=edge_cat, y = Estimate, ymin=Q2.5, ymax=Q97.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1.05) +
  geom_linerange() +
  geom_point() +
  facet_grid(transition ~ fallow_cat) +
  scale_y_continuous(breaks = seq(0.8, 1.2, length.out = 5),
                                                           labels = c("-20%", "-10%", "No change", "+10%", "+20%")) +
  labs(x = "Structural complexity of the landscape",
       y = "Predicted changes in bird richness (95% CrI)")

ggsave("figures/bird_temporal_rich.png", gg_rich)

# new based on lm's
bird %>%
  pivot_wider(id_cols = routcode, names_from = year_cat, 
              values_from = fallow_ha, names_prefix = "y") %>%
  mutate(d1 = y2010 - y2007,
         d2 = y2016 - y2010) -> ff

# model the first transition
mf <- lm(d1 ~ y2007, ff)
# predictions for changes and new values
f2007 <- seq(min(ff$y2007,na.rm=TRUE), max(ff$y2007,na.rm=TRUE), length.out = 20)
d2010 <- predict(mf, newdata = data.frame(y2007 = f2007))
f2010 <- d2010 + f2007

# put together in a df
new_fallow <- data.frame(fallow_ha = c(f2007, f2010),
                         year = rep(c(2007, 2010), each = 20))

new_fallow$fallow_std <- (sqrt(new_fallow$fallow_ha) - mean(bird$fallow_sqrt)) / sd(bird$fallow_sqrt)

# expand with edges
edge_q <- quantile(bird$edge_std, probs = c(0.25, 0.5, 0.75))
new_df <- data.frame(year_cat = factor(2007),
                     intensity_std = 0,
                     edge_std = rep(edge_q, each = 40),
                     fallow_std = rep(new_fallow$fallow_std, 3))

# the second transition
# model the second transition
mf2 <- lm(d2 ~ y2010, ff)

# predictions for changes and new values
f2010 <- seq(min(ff$y2010,na.rm=TRUE), max(ff$y2010,na.rm=TRUE), length.out = 20)
d2016 <- predict(mf2, newdata = data.frame(y2010 = f2010))
f2016 <- d2016 + f2010

# put together in a df
new_fallow <- data.frame(fallow_ha = c(f2010, f2016),
                         year = rep(c(2007, 2010), each = 20))

new_fallow$fallow_std <- (sqrt(new_fallow$fallow_ha) - mean(bird$fallow_sqrt)) / sd(bird$fallow_sqrt)

# expand with edges
edge_q <- quantile(bird$edge_std, probs = c(0.25, 0.5, 0.75))
new_df2 <- data.frame(year_cat = factor(2007),
                     intensity_std = 0,
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
df_plot2.1$transition <- "Transition 1:\nabolishment of set-aside"

df_plot2.2 <- rbind(plot_2(m_list$edge, new_df2, d2016, type = "Edge-breeding birds"),
                    plot_2(m_list$field, new_df2, d2016, type = "Field-breeding birds"),
                    plot_2(m_list$low, new_df2, d2016, type = "Non fallow-breeding birds"))
df_plot2.2$transition <- "Transition 2:\nestablishment of EFA"

df_plot2 <- rbind(df_plot2.1, df_plot2.2)
df_plot2.1$edge_cat <- factor(df_plot2.1$edge_std,
                            labels = c("Low", "Middle", "High"))
df_plot2.2$edge_cat <- factor(df_plot2.2$edge_std,
                              labels = c("Low", "Middle", "High"))

gg_1 <- ggplot(df_plot2.1, aes(x=delta_fallow, y=Estimate, ymin=Q2.5, ymax=Q97.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1.2) +
  geom_ribbon(aes(fill=edge_cat), alpha=0.1, color=NA) +
  geom_line(aes(color=edge_cat)) +
  facet_wrap(~type) +
  scale_y_continuous(breaks = seq(0.4, 1.2, by = 0.2),
                     labels = c("-60%", "-40%", "-20%", "0%", "+20%")) +
  labs(x = "Changes in fallow land area between the 2010 and the 2007 agricultural census (in ha)",
       y = "Predicted changes in bird abundance (with 95% CrI)") +
  guides(color=guide_legend(title = "Structural complexity\nof the landscape"),
         fill = guide_legend(title = "Structural complexity\nof the landscape")) +
  theme(axis.line = element_line(),
        panel.background = element_blank())

ggsave("figures/bird_changes_transition1.png", gg_1)


gg_2 <- ggplot(df_plot2.2, aes(x=delta_fallow, y=Estimate, ymin=Q2.5, ymax=Q97.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1.2) +
  geom_ribbon(aes(fill=edge_cat), alpha=0.1, color=NA) +
  geom_line(aes(color=edge_cat)) +
  facet_wrap(~type) +
  scale_y_continuous(breaks = c(0.9, 1.0, 1.1, 1.2, 1.3),
                     labels = c("-10%", "0%", "+10%", "+20%", "+30%")) +
  labs(x = "Changes in fallow land area between the 2016 and the 2010 agricultural census (in ha)",
       y = "Predicted changes in bird abundance (with 95% CrI)") +
  guides(color=guide_legend(title = "Structural complexity\nof the landscape"),
         fill = guide_legend(title = "Structural complexity\nof the landscape")) +
  theme(axis.line = element_line(),
        panel.background = element_blank())

ggsave("figures/bird_changes_transition2.png", gg_2)

# now for richness
df_plotr <- rbind(plot_2(m_list$rich, new_df, d2010, type = "Transition 1"),
      plot_2(m_list$rich, new_df2, d2016, type = "Transition 2"))
df_plotr$edge_cat <- factor(df_plotr$edge_std,
                            labels = c("Low", "Middle", "High"))

gg_r <- ggplot(df_plotr, aes(x=delta_fallow, y=Estimate, ymin=Q2.5, ymax=Q97.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1.2) +
  geom_ribbon(aes(fill=edge_cat), alpha=0.1, color=NA) +
  geom_line(aes(color=edge_cat)) +
  facet_wrap(~type) +
  scale_y_continuous(breaks = c(0.9, 1.0, 1.1),
                     labels = c("-10%", "0%", "+10%")) +
  labs(x = "Changes in fallow land area between the agricultural census (in ha)",
       y = "Predicted changes in bird richness (with 95% CrI)") +
  guides(color=guide_legend(title = "Structural complexity\nof the landscape"),
         fill = guide_legend(title = "Structural complexity\nof the landscape")) +
  theme(axis.line = element_line(),
        panel.background = element_blank())

ggsave("figures/rich_changes_transition12.png", gg_r)

##############
# other results: post. proba of effect being stronger in the different groups

# levels of edges to test
edge_std = matrix(c(rep(1, 3), quantile(bird$edge_std, 
                    probs = c(0.25, 0.5, 0.75))),
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

