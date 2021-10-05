##############
# R script to fit the models
# and run model checks
#######################


# set wd
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load libraries
library(rstan)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7') # only on win machine
rstan_options(auto_write = TRUE)
library(tidyverse)
# library(gridExtra)
library(sf)
library(DHARMa)
library(brms)

# load data
bird <- read.csv("../data/bird_data.csv")
district <- read.csv("../data/district_fallow.csv")

# some transformation / scaling of covariates
bird$fallow_sqrt <- sqrt(bird$fallow_per)
bird$fallow_std <- scale(bird$fallow_sqrt)
bird$edge_std <- scale(bird$edge_mha)
bird$agri_std <- scale(bird$agriculture)



# priors for the species richness model
bprior_rich <- c(prior(normal(0, 2.5), class = Intercept),
                prior(normal(0, 2.5), class = b),
                prior(student_t(3, 0, 2.5), class = sd),
                prior(student_t(3, 0, 2.5), class = sigma))

# fitting the species richness model
# note that due to z ~ N(mu, sd) is the same as 
# z = mu + sd * d where d ~ N(0, 1)
# we can write the hierarchical model with group-level
# covariates in brms
m_rich <- brm(bf(rich | trunc(lb = 0) ~ agri_std + (edge_std + I(edge_std ^ 2)) * fallow_std + 
                   (1 | NUTS_code / NUTS_year) + (0 + edge_std | NUTS_code / NUTS_year) +
                   (0 + I(edge_std ^ 2) | NUTS_code /NUTS_year)),
              data = bird, family = "normal", prior = bprior_rich,
              control = list(adapt_delta = 0.9, max_treedepth = 25))

# priors for the aubundance models
bprior <- c(prior(normal(0, 2.5), class = Intercept),
            prior(normal(0, 2.5), class = b),
            prior(student_t(3, 0, 2.5), class = sd),
            prior(gamma(0.01, 0.01), class = shape))

# fitting the three abundance models
m_edge <- brm(bf(edge ~ agri_std + edge_std * fallow_std + (1 | NUTS_code / NUTS_year) +
                   (0 + edge_std | NUTS_code / NUTS_year),
                 zi ~ agri_std + edge_std),
              data = bird, family = "zero_inflated_negbinomial",
              prior = bprior,
              control = list(adapt_delta = 0.9, max_treedepth = 25))

m_field <- brm(bf(field ~ agri_std + edge_std * fallow_std + (1 | NUTS_code / NUTS_year) +
                    (0 + edge_std | NUTS_code / NUTS_year),
                  zi ~ agri_std + edge_std),
              data = bird, family = "zero_inflated_negbinomial", 
              control = list(max_treedepth = 25, adapt_delta = 0.9))

m_low <- brm(bf(low ~ agri_std + edge_std * fallow_std + (1 | NUTS_code / NUTS_year) +
                  (0 + edge_std | NUTS_code / NUTS_year),
                  zi ~ agri_std + edge_std),
               data = bird, family = "zero_inflated_negbinomial", 
               control = list(max_treedepth = 25, adapt_delta = 0.9))

# save the models
m_list <- list(edge = m_edge, field = m_field, low = m_low, rich = m_rich)
saveRDS(m_list, "model_output/fitted_model.rds")


## model checks

### first posterior predictive distributions
pp_field <- pp_check(m_list$field, nsamples = 50) +
  labs(title = "Field-breeders") +
  scale_x_continuous(trans = "sqrt")
pp_edge <- pp_check(m_list$edge, nsamples = 50) +
  labs(title = "Edge-breeders") +
  scale_x_continuous(trans = "sqrt")
pp_low <- pp_check(m_list$low, nsamples = 50) +
  labs(title = "Foraging visitors") +
  scale_x_continuous(trans = "sqrt")
pp_rich <- pp_check(m_list$rich, nsamples = 50) +
  labs(title = "Species richness") +
  scale_x_continuous(trans = "sqrt")

gg_pp <- ggpubr::ggarrange(pp_rich, pp_field, pp_edge, pp_low, nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")
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
# ! this code is not reproducible due to data policy on the spatial location of the plots !

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
dff$type <- factor(dff$type, 
                   labels = c("Field-breeders",
                              "Edge-breeders",
                              "Foraging visitors",
                              "Species richness"))

dff$code <- paste0(dff$type, "\nyear: ", dff$year,
                   "\nMoran's I: ", round(dff$moran, 2),
                   " , p-value: ", round (dff$p_val, 0))

gg_sp <- ggplot(dff, aes(x=x, y=y, color = res)) +
  geom_point() +
  facet_wrap(~ code, nrow = 4) +
  scale_color_gradient2(low = "red", mid="white", 
                        high="blue", midpoint = 0.5) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 15)) +
  labs(x="", y="")

ggsave("figures/spat_autocorr.png", gg_sp,
       width = 12, height = 16)