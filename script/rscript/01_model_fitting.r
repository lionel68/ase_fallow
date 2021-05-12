# analysis script
setwd("~/PostDoc_Thunen/ASE-stuff/")

# load libraries
library(brms)
options(mc.cores = parallel::detectCores())
library(tidyverse)
# library(gridExtra)
library(sf)
library(DHARMa)

# load data
bird <- read.csv("data/preprocessed/bird_fallow_v7.csv")
#bird$land <- substr(bird$routcode, 1, 2)

# some transformation
bird$year_cat <- factor(bird$year)
bird$fallow_sqrt <- sqrt(bird$fallow_ha)
bird$fallow_std <- scale(bird$fallow_sqrt)
bird$edge_std <- scale(bird$edge_m)
bird$urban_std <- scale(bird$urban)
bird$forest_std <- scale(bird$forest)
bird$agri_std <- scale(bird$agriculture)
bird$bkr <- factor(bird$bkr)
bird$intensity_std <- scale(bird$intensity)
bird$x <- scales::rescale(bird$X_COORD, c(0, 1))
bird$y <- scales::rescale(bird$Y_COORD, c(0, 1))


# set the priors
bprior <- c(prior(normal(0, 2.5), class = Intercept),
            prior(normal(0, 2.5), class = b),
            prior(student_t(3, 0, 2.5), class = sd),
            prior(gamma(0.01, 0.01), class = shape))


bprior_rich <- c(prior(normal(0, 2.5), class = Intercept),
                prior(normal(0, 2.5), class = b),
                prior(student_t(3, 0, 2.5), class = sd),
                prior(student_t(3, 0, 2.5), class = sigma))

# the models
m_edge <- brm(bf(edge ~ year_cat + agri_std + edge_std * fallow_std + (1 | bkr),
                 hu ~ year_cat + agri_std + edge_std * fallow_std + (1 | bkr)),
                data = bird, prior = bprior,
                family = "hurdle_gamma")

m_field <- brm(bf(field ~ year_cat + edge_std * fallow_std + agri_std + (1 | bkr),
                  hu ~ year_cat + agri_std + edge_std * fallow_std + (1 | bkr)),
               data = bird, prior = bprior,
               family = "hurdle_gamma")

m_low <- brm(bf(low ~ year_cat + agri_std + edge_std * fallow_std + (1 | bkr),
                hu ~ year_cat + agri_std + edge_std * fallow_std + (1 | bkr)),
             data = bird, prior = bprior,
             family = "hurdle_gamma")

m_rich <- brm(bf(rich | trunc(lb = 0) ~ year_cat + agri_std + edge_std * fallow_std + (1 | bkr)),
              data = bird, prior = bprior_rich,
              family = "normal")

# save the models
m_list <- list(edge = m_edge, field = m_field, low = m_low, rich = m_rich)
saveRDS(m_list, "model_output/fitted_model_new.rds")
saveRDS(m_list, "model_output/fitted_model_new_wgrass.rds")
saveRDS(m_list, "model_output/fitted_model_new_hu.rds")



# re-load the models
m_list <- readRDS("model_output/fitted_model_new.rds")

## model checks

### first posterior predictive distributions
pp_field <- pp_check(m_list$field, nsamples = 50) +
  labs(title = "Field-breeding birds") +
  scale_x_continuous(trans = "sqrt")
pp_edge <- pp_check(m_list$edge, nsamples = 50) +
  labs(title = "Edge-breeding birds") +
  scale_x_continuous(trans = "sqrt")
pp_low <- pp_check(m_list$low, nsamples = 50) +
  labs(title = "Non fallow-breeding birds") +
  scale_x_continuous(trans = "sqrt")
pp_rich <- pp_check(m_list$rich, nsamples = 50) +
  labs(title = "Bird species richness") +
  scale_x_continuous(trans = "sqrt")

gg_pp <- ggpubr::ggarrange(pp_field, pp_edge, pp_low, pp_rich, nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
ggsave(filename = "figures/pp_checks.png", gg_pp)

### then DHARMa checks
dd_edge <- createDHARMa(t(posterior_predict(m_list$edge)),
                        bird$edge, integerResponse = TRUE)
dd_field <- createDHARMa(t(posterior_predict(m_list$field)),
                         bird$field, integerResponse = TRUE)
dd_low <- createDHARMa(t(posterior_predict(m_list$low)),
                       bird$low, integerResponse = TRUE)
dd_rich <- createDHARMa(t(posterior_predict(m_list$rich)),
                       bird$rich, integerResponse = TRUE)

png("figures/dharma_edge.png")
plot(dd_edge)
dev.off()

png("figures/dharma_field.png")
plot(dd_field)
dev.off()

png("figures/dharma_low.png")
plot(dd_low)
dev.off()

png("figures/dharma_rich.png")
plot(dd_rich)
dev.off()

### the spatial autocorr checks
## load the coords
mhb <- st_read("data/geodata/MhB_plots_S2637.shp")
mhb %>%
  rename(routcode = ROUTCODE) %>%
  st_drop_geometry() %>%
  select(routcode, X_COORD, Y_COORD) %>%
  right_join(bird) -> bird2

## an helper function
spat_autocorr <- function(model, data, colname = "field"){
  spp <- data.frame(routcode = character(0), type = character(0),
                    year = numeric(0), x = numeric(0),
                    y = numeric(0) , res = numeric(0))
  aa <- data.frame(type = character(0), year = numeric(0),
                   moran = numeric(0), p_val = numeric(0))
  for(i in 1:3){
    year <- c(2007, 2010, 2016)[i]
    tmp <- createDHARMa(t(posterior_predict(model))[data$year==year,],
                                 data[,colname][data$year==year], 
                                 integerResponse = TRUE)
    
    autocorr <- testSpatialAutocorrelation(tmp, x = data$X_COORD[data$year == year],
                                           y = data$Y_COORD[data$year == year],
                                           plot = FALSE)
    
    aa <- rbind(aa, data.frame(type = colname, year = year, 
                               moran = autocorr$statistic["observed"],
                               p_val = autocorr$p.value))
    
    spp <- rbind(spp, data.frame(routcode = data$routcode[data$year==year],
                                 type = colname,
                                 year = year,
                                 x = data$X_COORD[data$year==year],
                                 y = data$Y_COORD[data$year==year],
                                 res = tmp$scaledResiduals))
    
  }
  
  return(list(auto = aa, df = spp))
}

field_auto <- spat_autocorr(m_list$field, bird2, colname = "field")
edge_auto <- spat_autocorr(m_list$edge, bird2, colname = "edge")
low_auto <- spat_autocorr(m_list$low, bird2, colname = "low")
rich_auto <- spat_autocorr(m_list$rich, bird2, colname = "rich")

# put things together for plotting
dff <- rbind(field_auto$df, edge_auto$df, low_auto$df, rich_auto$df)
aaa <- rbind(field_auto$auto, edge_auto$auto, low_auto$auto, rich_auto$auto)
dff <- left_join(dff, aaa, by = c("type", "year"))
dff$code <- paste0("Taxa: ", dff$type, " , year: ", dff$year,
                   "\nMoran's I: ", round(dff$moran, 2),
                   " , p-value: ", round (dff$p_val, 0))

gg_sp <- ggplot(dff, aes(x=x, y=y, color = res)) +
  geom_point() +
  facet_wrap(~ code) +
  scale_color_gradient2(low = "red", mid="white", 
                        high="blue", midpoint = 0.5) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank()) +
  labs(x="", y="")

ggsave("figures/spat_autocorr.png", gg_sp)


# model explo
bird %>%
  select(routcode, year, edge, field, low, rich, edge_m, fallow_ha, bkr, intensity_std) %>%
  group_by(routcode) %>%
  mutate(delta_fallow = ifelse(lag(fallow_ha) != 0,
                               ((fallow_ha - lag(fallow_ha)) / lag(fallow_ha)) * 100,
                               NA)) %>%
  ungroup() %>%
  mutate(edge_std = scale(edge_m)) %>%
  pivot_longer(cols = 3:6) -> dd2

# build the datafrmae for the first transition for edge
dd2 %>%
  mutate(delta_fallow_std = scale(delta_fallow)) %>%
  filter(name == "rich") %>%
  pivot_wider(names_from = year, values_from = c(value, delta_fallow_std),
              id_cols = c(routcode, edge_std, bkr, intensity_std)) -> dd3

mm <- brm(value_2010 | trunc(lb = 0) ~ value_2007 + intensity_std + delta_fallow_std_2010 * edge_std + (1 | bkr), dd3)

# generate new prediction df
pred_df <- expand.grid(value_2007 = mean(dd3$value_2007, na.rm = TRUE),
                       intensity_std = mean(dd3$intensity_std, na.rm = TRUE),
                       delta_fallow_std_2010 = seq(-3.5, 0.77, length.out = 10),
                       edge_std = quantile(dd3$edge_std,
                                           probs = c(0.1, 0.5, 0.9), na.rm = TRUE))

ppred <- posterior_epred(mm, newdata = pred_df, re.form = NA) / pred_df$value_2007[1]
pred_df <- cbind(pred_df, posterior_summary(ppred))

ggplot(pred_df, aes(x=delta_fallow_std_2010, y=Estimate, ymin=Q2.5, ymax=Q97.5)) +
  geom_ribbon(aes(fill=factor(edge_std)), alpha = 0.1, color=NA) +
  geom_line(aes(color=factor(edge_std)))

mm2 <- brm(value_2016 | trunc(lb=0) ~ value_2010 + intensity_std+ delta_fallow_std_2016 * edge_std + (1 | bkr), dd3)

dd2 %>%
  mutate(delta_fallow_std = scale(delta_fallow)) %>%
  filter(name == "edge") %>%
  pivot_wider(names_from = year, values_from = c(value, delta_fallow_std),
              id_cols = c(routcode, edge_std, bkr, intensity_std)) -> dd3

mm <- brm(value_2010 ~ value_2007 + intensity_std + delta_fallow_std_2010 * edge_std + (1 | bkr), dd3,
          family = "poisson")

ress <- createDHARMa(t(posterior_predict(mm)),
                     observedResponse = dd3$value_2010[],
                     integerResponse = TRUE)


# some doodling with scenario dev
bird %>%
  pivot_wider(id_cols = routcode, names_from = year_cat, 
              values_from = fallow_ha, names_prefix = "y") %>%
  mutate(d1 = y2010 - y2007,
         d2 = y2016 - y2010) -> ff
