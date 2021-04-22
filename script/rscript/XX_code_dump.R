# cpde dump ASE analysis
setwd("PostDoc_Thunen/ASE-stuff/")
# bird$id <- 1:nrow(bird)

# calculating the end of the prediction gradient
bird %>%
  group_by(routcode) %>%
  mutate(delta_fallow = lead(fallow_std) - fallow_std) %>%
  group_by(year_cat) %>%
  summarise(delta_m = mean(delta_fallow, na.rm=TRUE),
            mi = min(fallow_std),
            ma = max(fallow_std)) -> fallow_end

# the data frame for the prediction
pred_fallow <- expand.grid(fallow_gradient = 1:10,
                           edge_std = c(-1.12, -0.18, 1.36),
                           year_cat = unique(bird$year_cat))

# new for edge_std
pred_fallow$fallow_std <- c(rep(seq(-1.44+0.767, 3.11-0.0536+0.767, length.out = 10), 3),
                            rep(seq(-1.44, 3.11-0.0536, length.out = 10), 3),
                            rep(seq(-1.44+0.0536, 3.11, length.out = 10), 3))

# old for swf cat
pred_fallow$fallow_std <- c(seq(-1.44+0.804, 3.56, length.out = 10),  seq(-1.44+0.754, 3.96, length.out = 10),
                            seq(-1.44, 3.56-0.804, length.out = 10), seq(-1.44, 3.96-0.754, length.out = 10),
                            seq(-1.44+0.048, 3.13-0.804+0.048, length.out = 10),  seq(-1.44+0.055, 3.96-0.754+0.055, length.out = 10))


pred_fallow$fallow_std <- c(seq(-1.44+0.804, 3.56, length.out = 10), seq(-1.44+0.193, 0.9+0.193, length.out = 10), seq(-1.44+0.806, 3.96, length.out = 10),
                            seq(-1.44, 3.56-0.804, length.out = 10), seq(-1.44, 0.9, length.out = 10), seq(-1.44, 3.96-0.806, length.out = 10),
                            seq(-1.44+0.048, 3.13-0.804+0.048, length.out = 10), seq(-1.44, 0.85, length.out = 10), seq(-1.44+0.066, 3.96-0.806+0.066, length.out = 10))
pred_fallow <- pred_fallow[, -1]
pred_fallow$urban_std <- 0
pred_fallow$forest_std <- 0

# pred_fallow$id <- 1:60
# 
# bird$routcode_2 <- factor(bird$routcode,
#                          levels = unique(bird$routcode))
# bird$group_2 <- factor(bird$group, labels = 1:6)
#   
# # set up adjancency matrix
# d_mhb <- as.matrix(dist(bird_dd[,c("X_COORD", "Y_COORD")]))
# W <- array(0, c(nrow(bird_dd), nrow(bird_dd)))
# W[d_mhb < 2500] <- 1
# rownames(W) <- bird_dd$routcode


# model checks
pp_check(m_strong)

dd <- createDHARMa(t(posterior_predict(m_strong)), bird$edge, 
                   integerResponse = TRUE)

# spatial autocorr
dd <- createDHARMa(t(posterior_predict(m_strong))[bird_dd$year==2010,],
                   bird_dd$strong_ground[bird_dd$year==2010], 
                   integerResponse = TRUE)
testSpatialAutocorrelation(dd, x=bird_dd$X_COORD[bird_dd$year==2010], y=bird_dd$Y_COORD[bird_dd$year==2010])

# plot the spatial pattern in the residuals
spp <- data.frame(x=bird$long_[bird$year==2010], y=bird$lat[bird$year==2010],
                  res=dd$scaledResiduals)
spp_sf <- st_as_sf(spp, coords = 1:2)
st_crs(spp_sf) <- 4326 


gg_sp <- ggplot() +
  geom_sf(data=de) +
  geom_sf(data = spp_sf, aes(color=res)) +
  scale_color_gradient2(low = "red", mid="white", 
                        high="blue", midpoint = 0.5) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank()) +
  labs(x="", y="")

ggsave("figures/spat_autocorr.png", gg_sp)

# for the prediction
# pred_dat <- pivot_longer(pred_fallow, cols=3:5, names_to = "fallow_id",
#                         values_to = "fallow_std")

# new model predictions
## first derive predictions across the gradients of fallow land
## in different edge cat
pred_fallow <- expand.grid(fallow_std = seq(min(bird$fallow_std),
                                            max(bird$fallow_std),
                                            length.out = 10),
                           edge_std = quantile(bird$edge_std, 
                                               probs = c(0.1, 0.5, 0.9)),
                           urban_std = 0, forest_std = 0, year_cat = 2007)
# model posterior draws
pp <- posterior_samples(m_strong, pars = "year_cat")
# expected model values
predd <- posterior_epred(m_strong, pred_fallow, re.form = NA)
# grab posterior samples to average over the years
predd <- predd + apply(pp, 1, function(x) mean(c(0, x)))
# compute quantiles
pred_df <- cbind(pred_fallow, posterior_summary(predd))
pred_df$edge_cat <- factor(rep(c("low", "middle", "high"), each = 10),
                           levels = c("low", "middle", "high"))
pred_df$fallow_ha <- (pred_df$fallow_std * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)) ** 2
# first plot
gg_edge1 <- ggplot(pred_df, aes(x=fallow_ha, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_ribbon(aes(fill=edge_cat), color = NA, alpha = 0.3) +
  geom_line(aes(color=edge_cat)) +
  labs(x = "Area of fallow in 1km buffer (ha)",
       y = "Predicted abundance (with 95% CrI)",
       title = "(a) Edge-breeding birds") +
  theme(panel.background = element_blank(),
        axis.line = element_line())

ggsave("figures/edge_fallow.png", gg_edge1)

## second compute mean and quantile changes of fallow lands
## and derive predictions from those in different edge cat
scl <- function(x){
  (x-mean(x)) / sd(x)
}

# get the full lists of plots
mhb_fallow <- read.csv("data/landusedata/mhb_1kmbuffer_fallow.csv")

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

# checking how fallow land area varied in regard to previous area
mhb_fallow %>%
  group_by(routcode) %>%
  mutate(fallow_delta = lead(fallow_ha) - fallow_ha) -> tmp

# re-format fallow_ref
fallow_ref %>%
  pivot_longer(cols = 2:7) %>%
  mutate(varn = rep(c("fallow_delta", "fallow_ha"), 9)) %>%
  mutate(quant = rep(rep(c("50%", "10%", "90%"), each = 2), 3)) %>%
  pivot_wider(id_cols = c(1, 5), names_from = "varn", values_from = "value") -> ref_dd

# some more wrangling around
ref_dd$fallow_delta <- c(-2.05, 0.013, -8.19, 0.0588, 3.55, -2.35, NA, NA, NA)

gg_f <- ggplot(subset(tmp, year != 2016), aes(x=fallow_ha, y=fallow_delta)) +
  geom_point(alpha=0.1) +
  geom_point(data=subset(ref_dd, year != 2016), color = "grey50", size = 5) +
  geom_point(data = subset(ref_dd, year != 2016), aes(color = quant), size = 2) +
  stat_smooth(method="lm", se = FALSE) +
  facet_grid(~ year)

ggsave("figures/fallow_temporal.png", gg_f)

# count how many plots stayed zero
mhb_fallow %>%
  group_by(routcode) %>%
  summarise(S = sum(fallow_ha)) -> tt

## the idea of the derived prediction is that the largest decline in fallow land area
## between 2010-2007 and 2016-2010 were observed around Mhb plots with large quantities
## of fallow land area. The smallest decline or the slight increases were observed 
## in the plots with low amount of fallow lands
## below are the construction of the relevant data frame to derive the predictions
## three categories of fallow land are used: 10, 50 and 90% of both the fallow land
## in 2007 and 2010 as well as the computed changes (at the route level) between
## 2010 and 2007 and 2016 and 2010.
# construct the prediction df
pred_f <- expand.grid(edge_std = quantile(bird$edge_std, 
                                          probs = c(0.1, 0.5, 0.9)),
                      fallow_std = c(c(1.80, 0.426, -1.12, 0.919, -0.415, -1.35),
                                     c(1.80-1.80, 0.426-0.612, -1.12+0.0191, 0.919-0.648, -0.415-0.0627, -1.35+1.17)),
                      urban_std = 0, forest_std = 0, year_cat = 2007)
# the final data frame
predd_f <- data.frame(edge_cat = rep(c("low", "middle", "high"), 6),
                      fallow_cat = rep(c("high", "middle", "low"), each = 3),
                      transition = rep(1:2, each = 9))

# model predictions
predd <- posterior_epred(m_strong, newdata = pred_f, re.form = NA)
# average out year-effect
predd <- predd + apply(pp, 1, function(x) mean(c(0, x)))
# make ratios
pred_fr <- predd[,19:36] / predd[,1:18]
predd_f <- cbind(predd_f, posterior_summary(pred_fr))
predd_fw <- pivot_wider(predd_f, id_cols = c(1, 3), names_from = fallow_cat, values_from = Estimate)
predd_fw$X <- rep(1:2, each = 3) + rep(c(-0.15, 0, 0.15), 2)
predd_fw$edge_cat <- factor(predd_fw$edge_cat, levels = c("low", "middle", "high"))

gg_edge2 <- ggplot(predd_fw, aes(x=X, y = middle, ymin = low, ymax=high, color=edge_cat)) +
  geom_hline(yintercept = 1, linetype="dashed", size=1.25)+
  geom_linerange() +
  geom_point() +
  scale_x_continuous(breaks = 1:2, labels = c("2007 -> 2010", "2010 -> 2016")) +
  scale_y_continuous(breaks = seq(0.8, 1.2, length.out = 5),
                     labels = c("-20%", "-10%", "No change", "+10%", "+20%")) +
  labs(x="",
       y = "Bird abundance changes",
       title = "(a) Edge-breeding birds") +
  theme(panel.background = element_blank(),
        axis.line = element_line())

ggsave("figures/edge_temporal.png", gg_edge2)

# predd_df <- expand.grid(fallow_id = 1:10,
#                         edge_std = c("low", "middle", "high"),
#                         transition = factor(1:2, 
#                                             labels = c("Transition 1:\nabolishment of mandatory set-aside",
#                                                        "Transition 2:\nestablishment of Ecological Focus Area")))
#   
#   
# predd_df <- cbind(predd_df, 
#                  t(apply((predd[,31:90] / predd[,1:60]), 2, quantile, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))))
# names(predd_df)[4:8] <- c("LCI", "LCI2", "Med", "UCI2", "UCI")
# # put the fallow area back in scale
# predd_df$fallow_sqrt <- pred_fallow$fallow_std[31:90] * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)
# # put in % arable land
# predd_df$fallow <- predd_df$fallow_sqrt ** 2
# 
# gg_edge <- ggplot(predd_df, aes(x=fallow, y=Med)) +
#   geom_hline(yintercept = 1, linetype= "dashed", size = 1.05) +
#   geom_ribbon(aes(fill=edge_std, ymin=LCI, ymax=UCI),alpha=0.1) +
#   geom_ribbon(aes(fill=edge_std, ymin=LCI2, ymax=UCI2), alpha = 0.3) +
#   geom_line(aes(color=edge_std)) +
#   guides(color = guide_legend("Landscape\ncontext"),
#          fill = guide_legend("Landscape\ncontext")) +
#   facet_grid(~transition) +
#   scale_y_continuous(breaks = c(0.5, 1, 1.5),
#                      labels = c("-50%", "No change", "+50%"),
#                      name = "Changes in bird abundance (with credible intervals)") +
#   labs(x="Fallow land area (ha in 1km buffer)")
#   
# ggsave("figures/edge.png", gg_edge)

# make hypothesis



# try out with regional var
# bird %>%
#   filter(!is.na(swf_cat)) %>%
#   group_by(land, year_cat, swf_cat) %>%
#   summarise(q25 = quantile(fallow_std, probs = 0.25),
#             q50 = quantile(fallow_std, probs = 0.50),
#             q75 = quantile(fallow_std, probs = 0.75)) -> pred_fallow
# 
# # for the prediction
# pred_dat <- pivot_longer(pred_fallow, cols=4:6, names_to = "fallow_id",
#                          values_to = "fallow_std")
# 
# # model predictions
# predd <- posterior_linpred(m_strong, newdata = pred_dat, re.form = ~land)
# predd_df <- cbind(as.data.frame(pred_dat), t(predd))
# predd_df %>%
#   pivot_longer(6:4005) %>%
#   group_by(land, swf_cat, fallow_id, name) %>%
#   mutate(lrr = log(value / lag(value))) %>%
#   mutate(transition = ifelse(year_cat == "2010", 1, ifelse(year_cat == "2016", 2, NA))) %>%
#   group_by(land, swf_cat, fallow_id, transition) %>%
#   summarise(med = median(lrr, na.rm = TRUE),
#             lci = quantile(lrr, probs=0.1, na.rm = TRUE),
#             uci = quantile(lrr, probs = 0.9, na.rm=TRUE)) %>%
#   filter(!is.na(transition)) -> predd_df2
# 
# ggplot(predd_df2, aes(x=land, y=med, color=fallow_id)) +
#   geom_jitter() +
#   facet_grid(transition~swf_cat)
# 
# predd_df <- expand.grid(fallow_id = c("q25", "q50", "q75"),
#                         swf_cat = c("cleared", "simple"),
#                         transition = factor(1:2, 
#                                             labels = c("Transition 1:\nabolishment of mandatory set-aside",
#                                                        "Transition 2:\nestablishment of Ecological Focus Area")),
#                         land = unique(bird$land))
# 
# 
# predd_df <- cbind(predd_df, 
#                   t(apply(log(predd[,7:18] / predd[,1:12]), 2, quantile, probs = c(0.1, 0.5, 0.9))))
# names(predd_df)[4:6] <- c("LCI", "Med", "UCI")
# 
# predd_df$X <- rep(1:2, each = 3) + rep(c(-0.2, 0, 0.2), times = 2)
# 
# gg_sg <- ggplot(predd_df, aes(x=X, y=Med,ymin=LCI,ymax=UCI, color = fallow_id)) +
#   geom_linerange() +
#   geom_point() +
#   scale_x_continuous(breaks = c(1, 2), 
#                      labels = c("Cleared landscape\nSWF < 2%", "Simple landscape\nSWF > 2%")) +
#   facet_grid(~transition) +
#   scale_color_discrete(name = "Amount of fallow\nland:", 
#                        labels = c("25% quantile", "Median", "75% quantile")) +
#   scale_y_continuous(breaks = log(c(0.7, 0.8, 0.9, 1, 1.1)),
#                      labels = c("-30%", "-20%", "-10%", "No change", "+10%"),
#                      name = "Changes in bird abundance") +
#   labs(title = "(a) Strong association to fallow land, ground-breeders")
# 

# next group
m_field <- brm(bf(field ~ year_cat + edge_std * fallow_std + (1 | bkr),
                  zi ~ urban_std + forest_std),
               data = bird, 
               family = "zero_inflated_negbinomial")

# model checks
pp_check(m_field)
dd <- createDHARMa(t(posterior_predict(m_field)), bird$field,
                   integerResponse = TRUE)


# new model predictions
## first derive predictions across the gradients of fallow land
## in different edge cat
pred_fallow <- expand.grid(fallow_std = seq(min(bird$fallow_std),
                                            max(bird$fallow_std),
                                            length.out = 10),
                           edge_std = quantile(bird$edge_std, 
                                               probs = c(0.1, 0.5, 0.9)),
                           urban_std = 0, forest_std = 0, year_cat = 2007)
# model posterior draws
pp <- posterior_samples(m_field, pars = "year_cat")
# expected model values
predd <- posterior_epred(m_field, pred_fallow, re.form = NA)
# grab posterior samples to average over the years
predd <- predd + apply(pp, 1, function(x) mean(c(0, x)))
# compute quantiles
pred_df <- cbind(pred_fallow, posterior_summary(predd))
pred_df$edge_cat <- factor(rep(c("low", "middle", "high"), each = 10),
                           levels = c("low", "middle", "high"))
pred_df$fallow_ha <- (pred_df$fallow_std * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)) ** 2
# first plot
gg_field1 <- ggplot(pred_df, aes(x=fallow_ha, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_ribbon(aes(fill=edge_cat), color = NA, alpha = 0.3) +
  geom_line(aes(color=edge_cat)) +
  labs(x = "Area of fallow in 1km buffer (ha)",
       y = "Predicted abundance (with 95% CrI)",
       title = "(b) Field-breeding birds") +
  theme(panel.background = element_blank(),
        axis.line = element_line())

ggsave("figures/field_fallow.png", gg_field1)

# construct the prediction df
pred_f <- expand.grid(edge_std = quantile(bird$edge_std, 
                                          probs = c(0.1, 0.5, 0.9)),
                      fallow_std = c(c(1.80, 0.426, -1.12, 0.919, -0.415, -1.35),
                                     c(1.80-1.80, 0.426-0.612, -1.12+0.0191, 0.919-0.648, -0.415-0.0627, -1.35+1.17)),
                      urban_std = 0, forest_std = 0, year_cat = 2007)
# the final data frame
predd_f <- data.frame(edge_cat = rep(c("low", "middle", "high"), 6),
                      fallow_cat = rep(c("high", "middle", "low"), each = 3),
                      transition = rep(1:2, each = 9))

# model predictions
predd <- posterior_epred(m_field, newdata = pred_f, re.form = NA)
# average out year-effect
predd <- predd + apply(pp, 1, function(x) mean(c(0, x)))
# make ratios
pred_fr <- predd[,19:36] / predd[,1:18]
predd_f <- cbind(predd_f, posterior_summary(pred_fr))
predd_fw <- pivot_wider(predd_f, id_cols = c(1, 3), names_from = fallow_cat, values_from = Estimate)
predd_fw$X <- rep(1:2, each = 3) + rep(c(-0.15, 0, 0.15), 2)
predd_fw$edge_cat <- factor(predd_fw$edge_cat, levels = c("low", "middle", "high"))

gg_field2 <- ggplot(predd_fw, aes(x=X, y = middle, ymin = low, ymax=high, color=edge_cat)) +
  geom_hline(yintercept = 1, linetype="dashed", size = 1.25) +
  geom_linerange() +
  geom_point() +
  scale_x_continuous(breaks = 1:2, labels = c("2007 -> 2010", "2010 -> 2016")) +
  scale_y_continuous(breaks = seq(0.8, 1.2, length.out = 5),
                     labels = c("-20%", "-10%", "No change", "+10%", "+20%")) +
  labs(x="",
       y = "Bird abundance changes",
       title = "(b) Field-breeding bird") +
  theme(panel.background = element_blank(),
        axis.line = element_line())

ggsave("figures/field_temporal.png", gg_field2)


# the plot
# model predictions
predd <- posterior_epred(m_field, newdata = pred_fallow, re.form = NA)

predd_df <- expand.grid(fallow_id = 1:10,
                        edge_std = c("low", "middle", "high"),
                        transition = factor(1:2, 
                                            labels = c("Transition 1:\nabolishment of mandatory set-aside",
                                                       "Transition 2:\nestablishment of Ecological Focus Area")))


predd_df <- cbind(predd_df, 
                  t(apply(predd[,31:90] / predd[,1:60], 2, quantile, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))))
names(predd_df)[4:8] <- c("LCI", "LCI2", "Med", "UCI2", "UCI")
# put the fallow area back in scale
predd_df$fallow_sqrt <- pred_fallow$fallow_std[31:90] * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)
# put in % arable land
predd_df$fallow_percent <- ((predd_df$fallow_sqrt ** 2) * 10000 / (pi * 1000 ** 2)) * 100
predd_df$fallow <- predd_df$fallow_sqrt ** 2

gg_field <- ggplot(predd_df, aes(x=fallow, y=Med)) +
  geom_hline(yintercept = 1, linetype= "dashed", size = 1.05) +
  geom_ribbon(aes(fill=edge_std, ymin=LCI, ymax=UCI),alpha=0.1) +
  geom_ribbon(aes(fill=edge_std, ymin=LCI2, ymax=UCI2), alpha = 0.3) +
  geom_line(aes(color=edge_std)) +
  guides(color = guide_legend("Landscape\ncontext"),
         fill = guide_legend("Landscape\ncontext")) +
  facet_grid(~transition) +
  scale_y_continuous(breaks = c(0.5, 1, 1.5, 2),
                     labels = c("-50%", "No change", "+50%", "+100%"),
                     name = "Changes in bird abundance (with credible intervals)") +
  labs(x="Fallow land area (ha in 1km buffer)")

ggsave("figures/field_nester.png", gg_field)

# next model
m_low <- brm(bf(low ~ year_cat + edge_std * fallow_std + (1 | bkr),
                zi ~ urban_std + forest_std),
             data = bird, 
             family = "zero_inflated_negbinomial")

# model checks
pp_check(m_low)
dd <- createDHARMa(t(posterior_predict(m_low)), bird$low,
                   integerResponse = TRUE)


# new model predictions
## first derive predictions across the gradients of fallow land
## in different edge cat
pred_fallow <- expand.grid(fallow_std = seq(min(bird$fallow_std),
                                            max(bird$fallow_std),
                                            length.out = 10),
                           edge_std = quantile(bird$edge_std, 
                                               probs = c(0.1, 0.5, 0.9)),
                           urban_std = 0, forest_std = 0, year_cat = 2007)
# model posterior draws
pp <- posterior_samples(m_low, pars = "year_cat")
# expected model values
predd <- posterior_epred(m_low, pred_fallow, re.form = NA)
# grab posterior samples to average over the years
predd <- predd + apply(pp, 1, function(x) mean(c(0, x)))
# compute quantiles
pred_df <- cbind(pred_fallow, posterior_summary(predd))
pred_df$edge_cat <- factor(rep(c("low", "middle", "high"), each = 10),
                           levels = c("low", "middle", "high"))
pred_df$fallow_ha <- (pred_df$fallow_std * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)) ** 2
# first plot
gg_low1 <- ggplot(pred_df, aes(x=fallow_ha, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_ribbon(aes(fill=edge_cat), color = NA, alpha = 0.3) +
  geom_line(aes(color=edge_cat)) +
  labs(x = "Area of fallow in 1km buffer (ha)",
       y = "Predicted abundance (with 95% CrI)",
       title = "(c) Low-association birds") +
  theme(panel.background = element_blank(),
        axis.line = element_line())

ggsave("figures/low_fallow.png", gg_low1)

# construct the prediction df
pred_f <- expand.grid(edge_std = quantile(bird$edge_std, 
                                          probs = c(0.1, 0.5, 0.9)),
                      fallow_std = c(c(1.80, 0.426, -1.12, 0.919, -0.415, -1.35),
                                     c(1.80-1.80, 0.426-0.612, -1.12+0.0191, 0.919-0.648, -0.415-0.0627, -1.35+1.17)),
                      urban_std = 0, forest_std = 0, year_cat = 2007)
# the final data frame
predd_f <- data.frame(edge_cat = rep(c("low", "middle", "high"), 6),
                      fallow_cat = rep(c("high", "middle", "low"), each = 3),
                      transition = rep(1:2, each = 9))

# model predictions
predd <- posterior_epred(m_low, newdata = pred_f, re.form = NA)
# average out year-effect
predd <- predd + apply(pp, 1, function(x) mean(c(0, x)))
# make ratios
pred_fr <- predd[,19:36] / predd[,1:18]
predd_f <- cbind(predd_f, posterior_summary(pred_fr))
predd_fw <- pivot_wider(predd_f, id_cols = c(1, 3), names_from = fallow_cat, values_from = Estimate)
predd_fw$X <- rep(1:2, each = 3) + rep(c(-0.15, 0, 0.15), 2)
predd_fw$edge_cat <- factor(predd_fw$edge_cat, levels = c("low", "middle", "high"))

gg_low2 <- ggplot(predd_fw, aes(x=X, y = middle, ymin = low, ymax=high, color=edge_cat)) +
  geom_hline(yintercept = 1, linetype="dashed", size = 1.25) +
  geom_linerange() +
  geom_point() +
  scale_x_continuous(breaks = 1:2, labels = c("2007 -> 2010", "2010 -> 2016")) +
  scale_y_continuous(breaks = seq(0.8, 1.2, length.out = 5),
                     labels = c("-20%", "-10%", "No change", "+10%", "+20%")) +
  labs(x="",
       y = "Bird abundance changes",
       title = "(c) Low-association birds") +
  theme(panel.background = element_blank(),
        axis.line = element_line())

ggsave("figures/low_temporal.png", gg_low2)

# the plot
# model predictions
predd <- posterior_epred(m_low, newdata = pred_fallow, re.form = NA)

predd_df <- expand.grid(fallow_id = 1:10,
                        swf_cat = c("low", "middle", "high"),
                        transition = factor(1:2, 
                                            labels = c("Transition 1:\nabolishment of mandatory set-aside",
                                                       "Transition 2:\nestablishment of Ecological Focus Area")))


predd_df <- cbind(predd_df, 
                  t(apply(predd[,31:90] / predd[,1:60], 2, quantile, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))))
names(predd_df)[4:8] <- c("LCI", "LCI2", "Med", "UCI2", "UCI")
# put the fallow area back in scale
predd_df$fallow_sqrt <- pred_fallow$fallow_std[31:90] * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)
# put in % arable land
#predd_df$fallow_percent <- ((predd_df$fallow_sqrt ** 2) * 10000 / (pi * 1000 ** 2)) * 100
predd_df$fallow <- predd_df$fallow_sqrt ** 2

gg_low <- ggplot(predd_df, aes(x=fallow, y=Med)) +
  geom_hline(yintercept = 1, linetype= "dashed", size = 1.05) +
  geom_ribbon(aes(fill=swf_cat, ymin=LCI, ymax=UCI),alpha=0.1) +
  geom_ribbon(aes(fill=swf_cat, ymin=LCI2, ymax=UCI2), alpha = 0.3) +
  geom_line(aes(color=swf_cat)) +
  guides(color = guide_legend("Landscape\ncontext"),
         fill = guide_legend("Landscape\ncontext")) +
  facet_grid(~transition) +
  scale_y_continuous(breaks = c(0.5, 1, 1.5, 2),
                     labels = c("-50%", "No change", "+50%", "+100%"),
                     name = "Changes in bird abundance (with credible intervals)") +
  labs(x="Fallow land area (ha in 1km buffer)")

ggsave("figures/non_breeders.png", gg_low)





# next diversity
m_shannon <- brm(shannon ~ year_cat * swf_cat * fallow_std + (1 | land), bird,
                 family = "skew_normal", control = list(adapt_delta=0.99))

# model checks
pp_check(m_shannon)
dd <- createDHARMa(t(posterior_predict(m_shannon)), bird$shannon)
plot(dd, quantreg=TRUE)

# the plot
# model predictions
predd <- posterior_epred(m_shannon, newdata = pred_fallow, re.form = NA)

predd_df <- expand.grid(fallow_id = 1:10,
                        swf_cat = c("cleared", "simple"),
                        transition = factor(1:2, 
                                            labels = c("Transition 1:\nabolishment of mandatory set-aside",
                                                       "Transition 2:\nestablishment of Ecological Focus Area")))


predd_df <- cbind(predd_df, 
                  t(apply(log(predd[,21:60] / predd[,1:40]), 2, quantile, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))))
names(predd_df)[4:8] <- c("LCI", "LCI2", "Med", "UCI2", "UCI")
# put the fallow area back in scale
predd_df$fallow_sqrt <- pred_fallow$fallow_std[21:60] * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)
# put in % arable land
predd_df$fallow_percent <- ((predd_df$fallow_sqrt ** 2) * 10000 / (pi * 1000 ** 2)) * 100
predd_df$fallow <- predd_df$fallow_sqrt ** 2

gg_div <- ggplot(predd_df, aes(x=fallow, y=Med)) +
  geom_hline(yintercept = log(1), linetype= "dashed", size = 1.05) +
  geom_ribbon(aes(fill=swf_cat, ymin=LCI, ymax=UCI),alpha=0.1) +
  geom_ribbon(aes(fill=swf_cat, ymin=LCI2, ymax=UCI2), alpha = 0.3) +
  geom_line(aes(color=swf_cat)) +
  guides(color = guide_legend("Landscape\ncontext"),
         fill = guide_legend("Landscape\ncontext")) +
  facet_grid(~transition) +
  scale_y_continuous(breaks = log(c(0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3)),
                     labels = c("-30%", "-20%", "-10%", "No change", "+10%", "+20%", "+30%"),
                     name = "Changes in bird diversity (with credible intervals)") +
  labs(x="Fallow land area (ha in 1km buffer)")

ggsave("figures/shannon_diversity.png", gg_div)


# next group
m_middle <- brm(middle_ground ~ year_cat * swf_cat * fallow_std + (1 | land),
                bird, family = "zero_inflated_negbinomial", control = list(adapt_delta=0.99))

# check model
pp_check(m_middle)
dd <- createDHARMa(t(posterior_predict(m_middle)), bird$middle_ground,
                   integerResponse = TRUE)
plot(dd, quantreg=TRUE)

# model predictions
predd <- posterior_epred(m_middle, newdata = pred_fallow, re.form = NA)

predd_df <- expand.grid(fallow_id = 1:10,
                        swf_cat = c("cleared", "simple"),
                        transition = factor(1:2, 
                                            labels = c("Transition 1:\nabolishment of mandatory set-aside",
                                                       "Transition 2:\nestablishment of Ecological Focus Area")))


predd_df <- cbind(predd_df, 
                  t(apply(log(predd[,21:60] / predd[,1:40]), 2, quantile, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))))
names(predd_df)[4:8] <- c("LCI", "LCI2", "Med", "UCI2", "UCI")
# put the fallow area back in scale
predd_df$fallow_sqrt <- pred_fallow$fallow_std[21:60] * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)
# put in % arable land
predd_df$fallow_percent <- ((predd_df$fallow_sqrt ** 2) * 10000 / (pi * 1000 ** 2)) * 100
predd_df$fallow <- predd_df$fallow_sqrt ** 2

gg_middle <- ggplot(predd_df, aes(x=fallow, y=Med)) +
  geom_hline(yintercept = log(1), linetype= "dashed", size = 1.05) +
  geom_ribbon(aes(fill=swf_cat, ymin=LCI, ymax=UCI),alpha=0.1) +
  geom_ribbon(aes(fill=swf_cat, ymin=LCI2, ymax=UCI2), alpha = 0.3) +
  geom_line(aes(color=swf_cat)) +
  guides(color = guide_legend("Landscape\ncontext"),
         fill = guide_legend("Landscape\ncontext")) +
  facet_grid(~transition) +
  scale_y_continuous(breaks = log(c(0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3)),
                     labels = c("-30%", "-20%", "-10%", "No change", "+10%", "+20%", "+30%"),
                     name = "Changes in bird diversity (with credible intervals)") +
  labs(x="Fallow land area (ha in 1km buffer)")

ggsave("figures/shannon_diversity.png", gg_div)


# for the prediction
pred_dat <- data.frame(year_cat = factor(rep(c(2007, 2010, 2016), times = 4)),
                       swf_cat = factor(rep(c("cleared", "simple"), each = 6)),
                       fallow_std = c(0.564, -0.215, -0.072, rep(0.564, 3), 0.158, -0.416, -0.258, rep(0.158, 3)),
                       fallow = rep(rep(c("changing", "constant"), each = 3), times = 2))
# model predictions
predd <- posterior_linpred(m_nstrong, newdata = pred_dat, re.form = NA)
predd_df <- data.frame(transition = rep(c(1, 2), 4),
                       swf = rep(c("simple", "cleared"), each = 4),
                       fallow = rep(c("changing", "constant"), each = 2))
predd_df <- cbind(predd_df, 
                  t(apply(log(predd[,c(2, 3, 5, 6, 8, 9, 11, 12)] / predd[,c(1, 2, 4, 5, 7, 8, 10, 11)]), 2, quantile, probs = c(0.1, 0.5, 0.9))))
names(predd_df)[4:6] <- c("LCI", "Med", "UCI")
predd_df$X <- c(0.9, 1.9, 0.95, 1.95, 1.1, 2.1, 1.05, 2.05)
ggplot(predd_df, aes(x=X, y=Med,ymin=LCI,ymax=UCI, shape=swf, color = fallow)) +
  geom_linerange() +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2), 
                     labels = c("2007 -> 2010", "2010 -> 2016"))

# data generation part

### Part 1. : combine the municipality level data
### geometries, fallow and production cost

## load municipality geometries
gem <- st_read("data/geodata/gemeinde_correct.shp")
names(gem)[4] <- "reg"
gem <- gem[,"reg"]
gem <- st_transform(gem, 3035) # project to 3035

## load ASE data
ase_1999 <- read.csv("data/landusedata/ase_landuse_1999_2003_2007.csv", dec = ",",
                     stringsAsFactors = FALSE)
ase_2010 <- read.csv("data/landusedata/ase_landuse_2010.csv", dec = ",",
                     stringsAsFactors = FALSE)
ase_2016 <- read.csv("data/landusedata/ase_landuse_2016.csv", dec = ",",
                     stringsAsFactors = FALSE)

## clean negative values and NAs
ase_2010 <- ase_2010[which(ase_2010$value >= 0),]
ase_2016 <- ase_2016[which(!is.na(ase_2016$value)),]

## some re-naming to be consistent
names(ase_1999) <- c("reg", "landuse", "year", "value")
names(ase_2010) <- c("reg", "landuse", "year", "value")
names(ase_2016)[c(3, 7)] <- c("year", "landuse")

## rbind
ase_all <- rbind(ase_1999, ase_2010, ase_2016[,c("reg", "landuse", "year", "value")])

# seom gemeinde do not have three record, set a dummy df to fill it up with 0s
df_reg <- expand.grid(reg = unique(ase_all$reg),
                      year = c(2007, 2010, 2016),
                      valueX = 0)

## compute changes
# put stilllegung and brache together
ase_all %>%
  mutate(value = value * 1000) %>% # put values in ha
  filter(landuse %in% c("SETA", "FALL")) %>%
  filter(year > 2004) %>%
  group_by(reg, year) %>%
  summarise(fallow = sum(value, na.rm = FALSE)) %>%
  right_join(df_reg, by = c("reg", "year")) %>%
  mutate(fallow = ifelse(is.na(fallow), valueX, fallow)) %>%
  arrange(reg, year) %>%
  ungroup() %>%
  select(-valueX) -> ddb


## load production cost data
intensity <- read.csv("data/landusedata/intensity_costs.csv", sep = ";", dec = ",")
names(intensity)[1] <- "reg"
# a helper function to replace outliers with the quantile
rem_99 <- function(x,p=0.99){
  q99 <- quantile(x, probs = p)
  x_out <- ifelse(x>q99, q99, x)
  return(x_out)
}
## data formatting
intensity %>%
  pivot_wider(id_cols = c(reg, ha),
              names_from = agg_verfahren,
              values_from = kosten_pro_ha) %>% # widen up the practices
  mutate(across(3:7, ~ifelse(is.na(.x), 0, .x))) %>% # fill na rows with 0
  mutate(across(5:7, ~rem_99(.x))) %>% # the 1% top values in DAIR, OGRL and GRAN are replaced
  pivot_longer(3:7, names_to="type",
               values_to = "cost_per_ha") %>% # put back together
  group_by(reg) %>%
  summarise(ges_kost_ha = sum(cost_per_ha), # compute total costs and total cost per ha
            ges_kost = sum(ha * cost_per_ha)) -> int_dd

# add fallow and intensity infos to the shapefile
gem %>%
  left_join(ddb, by = "reg") %>%
  left_join(int_dd, by = "reg") %>%
  mutate(area_municipality = st_area(.)) -> gem_fallow

# plot this
gem_fallow %>%
  group_by(reg) %>%
  mutate(delta_fallow = fallow - lag(fallow)) %>%
  ungroup() %>%
  mutate(delta_cat = Hmisc::cut2(delta_fallow, g = 6)) -> gem_fallow

gem_fallow$delta_cat <- factor(gem_fallow$delta_cat,
                               labels = c("-1800; -46", "-46; -15", "-15; -1", "-1; 0", "0; 15", "15; 900"))
gem_fallow$year_cat <- factor(gem_fallow$year, 
                              labels = c("", "Transition 1 :\nabolishment in mandatory fallow",
                                         "Transition 2 :\nestablishment of EFA"))

g_g <- ggplot() +
  geom_sf(data = subset(gem_fallow, year != 2007), aes(fill = delta_cat), color = NA) +
  facet_grid(~year_cat) +
  scale_fill_brewer(palette = "PRGn", name = "Change in fallow\nland (in ha)") +
  theme(axis.text = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) 

g_b <- ggplot() +
  geom_sf(data = bkr, aes(fill= BKR_Nummer), color = "white") +
  scale_fill_continuous(type = "viridis") +
  theme(axis.text = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) 


### Part 2. : compute the fallow and production cost
### per Mhb plot

## load Mhb geometries (with swf)
mhb_gem <- st_read("data/preprocessed/mhb_swf_patch.shp")

## turn the polygons into points
# find intersection using plot centroids
mhb_pts <- st_as_sf(st_drop_geometry(mhb_gem), coords = 3:2)
st_crs(mhb_pts) <- 4326
mhb_pts <- st_transform(mhb_pts, 3035)

mhb_fallow_int <- st_intersection(st_buffer(mhb_pts, 1000), gem_fallow) # 1000m around plot centroids

# compute area of intersection
mhb_fallow_int$area_intersection <- st_area(mhb_fallow_int)
mhb_fallow_int <- st_drop_geometry(mhb_fallow_int)
# compute weighted sum
mhb_fallow_int %>%
  group_by(routcode, year) %>%
  summarise(fallow = sum((area_intersection / area_municipality) * fallow, na.rm=TRUE),
            intensity = sum((area_intersection / area_municipality) * ges_kost, na.rm=TRUE),
            swf = mean(swf_vec),
            patch_size = mean(patch_size),
            mesh_size = mean(mesh_size)) -> mhb_fallow

## Part 3. : derive bird metrices
## and join with other data

## load bird data
mhb <- read.csv("data/birddata/mhb_bird.csv")
# species traits to get association with fallow and breeding
traits <- read.csv("data/birddata/traits_final.csv",
                   sep = ";", dec = ",", stringsAsFactors = FALSE)
names(traits)[1] <- "Artname"
traits$breeding[traits$breeding=="non-ground"] <- "nground"

## compute shannon diversity
mhb %>%
  rename(ROUTCODE = ROUTENCODE) %>%
  left_join(traits) %>%
  group_by(ROUTCODE, Jahr) %>%
  mutate(p = Abund / sum(Abund)) %>%
  summarise(shannon = exp(-sum(p*log(p),na.rm = TRUE)))  -> sha_dd

## compute abundances per groups
mhb %>%
  rename(ROUTCODE = ROUTENCODE) %>%
  left_join(traits) %>%
  group_by(ROUTCODE, Jahr, fallow, breeding) %>%
  summarise(Abund = sum(Abund)) %>%
  pivot_wider(1:2, names_from=c("fallow", "breeding"),
              values_from="Abund") %>%
  left_join(sha_dd, by=c("ROUTCODE", "Jahr")) -> bird_dd

## save this for first models of temporal trends
write.csv(bird_dd, "data/preprocessed/bird_temporal.csv", row.names = FALSE)

## add fallow, intensity, swf ... information
bird_dd %>%
  rename(year = Jahr, routcode = ROUTCODE) %>%
  inner_join(mhb_fallow, by = c("routcode", "year")) %>%
  ungroup() %>%
  mutate(across(10:11, ~as.numeric(.x))) -> bird_fallow

## save this
write.csv(bird_fallow, "data/preprocessed/bird_fallow.csv", row.names = FALSE)

# comparison SWF-edge length
edge <- read.csv("data/landusedata/mhb_1kmbuffer_edge.csv")
swf <- read.csv("data/landusedata/mhb_1kmbuffer_swf.csv")

swf %>%
  group_by(routcode) %>%
  summarise(swf = sum(swf_sqm)) %>%
  full_join(edge) -> swf_edge

psych::pairs.panels(swf_edge[,-1])

# identify a bunch of plots

# model predictions
predd <- posterior_epred(m_strong, newdata = pred_f, re.form = NA)
# average out year-effect
predd <- predd + apply(pp, 1, function(x) mean(c(0, x)))
# make ratios
pred_fr <- predd[,19:36] / predd[,1:18]
predd_f <- cbind(predd_f, posterior_summary(pred_fr))
predd_fw <- pivot_wider(predd_f, id_cols = c(1, 3), names_from = fallow_cat, values_from = Estimate)
predd_fw$X <- rep(1:2, each = 3) + rep(c(-0.15, 0, 0.15), 2)
predd_fw$edge_cat <- factor(predd_fw$edge_cat, levels = c("low", "middle", "high"))

gg_edge2 <- ggplot(predd_fw, aes(x=X, y = middle, ymin = low, ymax=high, color=edge_cat)) +
  geom_hline(yintercept = 1, linetype="dashed", size=1.25)+
  geom_linerange() +
  geom_point() +
  scale_x_continuous(breaks = 1:2, labels = c("2007 -> 2010", "2010 -> 2016")) +
  scale_y_continuous(breaks = seq(0.8, 1.2, length.out = 5),
                     labels = c("-20%", "-10%", "No change", "+10%", "+20%")) +
  labs(x="",
       y = "Bird abundance changes",
       title = "(a) Edge-breeding birds") +
  theme(panel.background = element_blank(),
        axis.line = element_line())

ggsave("figures/edge_temporal.png", gg_edge2)

# compare ASE UAAR (estimation) and ATKIS agri (code 4101, 4102, 4103, 4104 and 4109)
mhb_uaar <- read.csv("K:/Dokumente/Bird-stuff/monvia_voegel/ro2_ase_fallow/data/landusedata/mhb_1kmbuffer_arab.csv")
mhb_agri <- read.csv("K:/Dokumente/Bird-stuff/monvia_voegel/ro2_ase_fallow/data/landusedata/mhb_1kmbuffer_atkisagriarea.csv")

# merge the two
mhb_all <- inner_join(mhb_uaar[,-1], mhb_agri[,-1])
mhb_all$land <- substr(mhb_all$routcode, 1, 2)
mhb_all <- subset(mhb_all, !is.na(mhb_all$year))
mhb_all$uaar_ha[mhb_all$uaar_ha * 1000 > 314] <- NA
# plot
gg_comp <- ggplot(mhb_all, aes(x = uaar_ha * 1000, y = agriarea_sqm * 1e-4, color = year)) +
  geom_point() +
  facet_wrap(~ land, scales = "free") +
  labs(x = "Estimated arable land area in 1km buffer from ASE-ARAB (ha)",
       y = "Arable land area in 1km buffer from ATKIS (ha)")

ggsave("~/PostDoc_Thunen/ASE-stuff/figures/ase_atkis_uaar.png", gg_comp)

# now compare estimation of percent fallow land area
mhb_f <- read.csv("K:/Dokumente/Bird-stuff/monvia_voegel/ro2_ase_fallow/data/landusedata/mhb_1kmbuffer_fallow2.csv")
# add to the df
mhb_all2 <- inner_join(mhb_all, mhb_f[,-1], by = c("routcode", "year"))
# compute percent
mhb_all2$percent_fallow <- ifelse(mhb_all2$uaar_ha > 0, mhb_all2$fallow_ha / mhb_all2$uaar_ha, NA)
mhb_all2$percent_fallow[mhb_all2$percent_fallow > 1] <- NA

# a plot
gg_f <- ggplot(mhb_all2, aes(x = factor(year), y=percent_fallow * 100, color=factor(year))) +
  geom_boxplot() +
  facet_wrap(~land, scales = "free") +
  labs(y = "Estimated proportion of fallow land area\nover arable land area in 1km buffer of MhB transects",
       x = "ASE year")

ggsave("~/PostDoc_Thunen/ASE-stuff/figures/fallow_temporal_mhb.png", gg_f)

# make a plot of percent fallow land area per municipality
aa <- read.csv("data/landusedata/agraratlas_2003_2016.csv")
muni <- st_read("~/PostDoc_Thunen/ASE-stuff/data/geodata/gemeinde_correct.shp")
names(muni)[4] <- "reg"

aa %>%
  filter(!is.na(allYEAR)) %>%
  filter(variable %in% c("UAAR", "ARAB", "FALL", "SETA")) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  pivot_wider(id_cols = c(reg, allYEAR),
              names_from = variable,
              values_from = value) %>%
  mutate(percent_fallow = ifelse(UAAR > 0, ((SETA + FALL) / UAAR) * 100, NA)) %>%
  mutate(percent_fallow_cat = cut(percent_fallow, 
                                 breaks = c(0, 1e-4, 2, 5, 10, 25, 100), 
                                 labels = c("0%", ">0 - 2%", ">2 - 5%", ">5 - 10%", ">10 - 25%", ">25%"),
                                 include.lowest = TRUE)) %>%
  select(reg, allYEAR, percent_fallow, percent_fallow_cat) %>%
  inner_join(muni, by = "reg") -> ff

ff <- st_as_sf(ff)
ff$year_cat <- paste0("Census year:\n", ff$allYEAR)

gg_1 <- ggplot(subset(ff, allYEAR > 2004 & allYEAR < 2015)) +
  geom_sf(aes(fill=percent_fallow_cat),color = NA) +
  facet_wrap(~year_cat, ncol = 1) +
  scale_fill_brewer(palette = "Greens",
                    name = "Proportion fallow land\nover arable land (%)") +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "left") 

ggsave("~/PostDoc_Thunen/ASE-stuff/figures/fallow_percent_changes.png", gg_1)

# another with the changes in percent changes
ff %>%
  arrange(reg, allYEAR) %>%
  group_by(reg) %>%
  mutate(fallow_diff = (percent_fallow - lag(percent_fallow))) %>%
  mutate(diff_cat = cut(fallow_diff,
                        breaks = c(-100, -10, -2, 0, 2, 10, 100),
                        labels = c("-100; -10%", "-10; -2%",
                                   "-2; 0%", "0; 2%", "2; 10%",
                                   "10; 100%"),
                        include.lowest = TRUE)) %>%
  filter(allYEAR > 2008) -> ff2
ff2$year_cat <- factor(ff2$allYEAR,
                       labels = c("Abolishment\nof mandatory\nset-aside:\n2010 - 2007",
                                  "Establishment\nof EFA:\n2016 - 2010"))
  
gg_2 <- ggplot(ff2) +
  geom_sf(aes(fill=diff_cat),color = NA) +
  facet_wrap(~year_cat, ncol = 1) +
  scale_fill_brewer(palette = "PRGn",
                    name = "Difference in proportion\nof fallow land\nbetween two census\nyears") +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

gg_a <- ggpubr::ggarrange(gg_1, gg_2, ncol = 2)

ggsave("~/PostDoc_Thunen/ASE-stuff/figures/fallow_spatiotemp_variation.png", gg_a)

# compare municipality area with UAAR
muni <- st_read("~/PostDoc_Thunen/ASE-stuff/data/geodata/gemeinde_correct.shp")
ase <- read.csv("~/PostDoc_Thunen/ASE-stuff/data/landusedata/agraratlas_2003_2016.csv")
ll <- st_read("~/PostDoc_Thunen/Bird-stuff/data_gis/1000_NUTS1.shp")

ll <- st_transform(ll, st_crs(muni))

names(muni)[4] <- "reg"
muni$area_sqm <- st_area(muni)

ase %>%
  filter(!is.na(allYEAR)) %>%
  filter(variable == "ARAB") %>%
  inner_join(muni, by = "reg") %>%
  mutate(prop_uaar = (value * 1000) / (as.numeric(area_sqm) * 1e-4)) -> ase_muni

ase_muni$prop_over <- ifelse(ase_muni$prop_uaar > 1, "Above", "Below")
ase_muni <- st_as_sf(ase_muni)


gg_arab <- ggplot() +
  geom_sf(data = ase_muni, aes(fill= prop_over), color = NA) +
  geom_sf(data = ll, fill = NA) +
  facet_wrap(~ allYEAR) +
  scale_fill_manual(values = c("red", "grey90"), name = "Proportion of arable\narea is larger than 1?") +
  labs(title = "Municipalities with larger ARAB area reported than the municipality area") +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

# plots of fallow land vs ARAB and UAAR
aa %>%
  filter(!is.na(allYEAR)) %>%
  filter(variable %in% c("UAAR", "ARAB", "FALL", "SETA")) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  pivot_wider(id_cols = c(reg, allYEAR),
              names_from = variable,
              values_from = value) %>%
  mutate(fallow = SETA + FALL) -> tmp

ggplot(tmp, aes(x=UAAR, y = ARAB)) +
  geom_point() +
  facet_wrap(~allYEAR)

ggplot(tmp, aes(x=UAAR, y = fallow)) +
  geom_point() +
  facet_wrap(~allYEAR) +
  stat_smooth()

ggplot(tmp, aes(x=ARAB, y = fallow)) +
  geom_point() +
  facet_wrap(~allYEAR) +
  stat_smooth()

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
