# analysis script
setwd("PostDoc_Thunen/ASE-stuff/")

# the path to save the outputs, must have sub-folder: model, dharma_fig and fallow_fig
output_path <- "/data/monvia_birds/ase_analysis"

# load libraries
library(brms)
options(mc.cores = parallel::detectCores())
library(tidyverse)
# library(gridExtra)
library(sf)
library(DHARMa)

# load data
bird <- read.csv("data/preprocessed/bird_fallow_species.csv")
#bird$land <- substr(bird$routcode, 1, 2)

# some transformation
bird$year_cat <- factor(bird$year)
bird$fallow_sqrt <- sqrt(bird$fallow)
bird$fallow_std <- scale(bird$fallow_sqrt)
bird$edge_std <- scale(bird$edge_m)
bird$urban_std <- scale(bird$urban)
bird$forest_std <- scale(bird$forest)
bird$bkr <- factor(bird$bkr)

# helper function to grab the warnings
tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}

# helper function and data for the fallow plot
pred_fallow <- expand.grid(fallow_std = seq(min(bird$fallow_std),
                                            max(bird$fallow_std),
                                            length.out = 10),
                           edge_std = quantile(bird$edge_std, 
                                               probs = c(0.1, 0.5, 0.9)),
                           urban_std = 0, forest_std = 0, year_cat = unique(bird$year_cat))
# helper function
plot_1 <- function(model, newdata, species){
  # expected model values
  predd <- posterior_epred(model, newdata, re.form = NA)
  # average over the years
  predd <- (predd[,1:30] + predd[,31:60] + predd[,61:90]) / 3
  # compute quantiles
  pred_df <- cbind(newdata, posterior_summary(predd))
  pred_df$edge_cat <- factor(rep(c("low", "middle", "high"), each = 10),
                             levels = c("low", "middle", "high"))
  pred_df$fallow_ha <- (pred_df$fallow_std * sd(bird$fallow_sqrt) + mean(bird$fallow_sqrt)) ** 2

  gg_plot11 <- ggplot(pred_df, aes(x=fallow_ha, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
    geom_ribbon(aes(fill=edge_cat), color = NA, alpha = 0.15) +
    geom_line(aes(color=edge_cat)) +
    guides(fill = guide_legend(title = "Structural complexity\nof the landscape : "),
           color = guide_legend(title = "Structural complexity\nof the landscape : ")) +
    labs(x = "Area of fallow in 1km buffer (ha)",
         y = "Predicted abundance (with 95% CrI)",
         title = species) +
    theme(panel.background = element_blank(),
          axis.line = element_line())
  
  ggsave(paste0(output_path, "/fallow_fig/", species, ".png"), gg_plot11)
}

# helper function to fit the model and generate various outputs
fit_brms <- function(data, species, output_path, ...){
  
  print(paste0("Starting fitting for : ", species))
  # the model
  m <- tryCatch.W.E(brm(bf(Abundance ~ year_cat + edge_std * fallow_std + (1 | bkr)),
           data = subset(data, species == sp), 
           family = "zero_inflated_negbinomial"))
  
  # save the model with warnings / errors
  saveRDS(m, file = paste0(output_path, "/model/", species, ".rds"))
  
  # make a DHARMa object
  if(!is.null(m$warning)){
    dd <- createDHARMa(t(posterior_predict(m$value)),
                       bird$Abundance[bird$Artname == species],
                       integerResponse = TRUE)
    
    # save the DHARMa plot
    png(paste0(output_path, "/dharma_fig/", species, ".png"))
    plot(dd)
    dev.off()
  }
  
  # make the fallow prediction plot
  if(!is.null(m$warning)){
    plot_1(m$value, pred_fallow, species)
  }
  
  return(paste0(species, " finished just fine!"))
}

# the actual fitting
for(sp in unique(bird$Artname)){
  fit_brms(bird, sp, output_path)
}
