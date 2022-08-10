## helper function ASE analysis
# function to fit a brm model to the abundance data
fit_stan_abund <- function(data, resp_name, output_path, file_id){
  # the formula
  form <- paste0(resp_name, " ~ agri_std + edge_std * fallow_std +
                         (1 | NUTS_code) + (0 + edge_std | NUTS_code)")
  
  # the prior
  bprior <- c(prior(student_t(3, 0, 2.5), class = Intercept),
              prior(normal(0, 2.5), class = b),
              prior("exponential(1)", class = sd),
              prior(gamma(0.01, 0.01), class = shape))
  
  # fit the model
  mod <- brm(bf(form,
                zi ~ agri_std + edge_std),
             data = data, family = "zero_inflated_negbinomial",
             prior = bprior,
             control = list(adapt_delta = 0.9, max_treedepth = 25))
  # save the model
  saveRDS(mod, file = paste0(output_path, file_id, ".rds"))
  gc()
  rm(mod)
}

# the same for species richness
fit_stan_rich <- function(data, resp_name, output_path, file_id){
  # the formula
  form <- paste0(resp_name, "| trunc(lb=0) ~ agri_std + fallow_std + edge_std + I(edge_std ** 2) +
                         fallow_std:edge_std + fallow_std:I(edge_std ** 2) +
                         (1 | NUTS_code) + (0 + edge_std | NUTS_code) +
                         (0 + I(edge_std ** 2) | NUTS_code)")
  
bprior <- c(prior(student_t(3, 0, 2.5), class = Intercept),
                        prior(normal(0, 2.5), class = b),
                        prior("exponential(1)", class = sd))
  
  # fit the model
  mod <- brms::brm(form, prior = bprior,
                       data = data)
  # save the model
  saveRDS(mod, file = paste0(output_path, file_id, ".rds"))
  gc()
  rm(mod)
}


  
# function to run stan checks:
# rhat, effective sample size, divergent iteration and dharma checks
stan_check <- function(model, data, colname, write_file = FALSE, file_id = "id"){
  require(posterior)
  # check Rhat and ESS
  ss <- summarise_draws(as_draws(model), default_convergence_measures())
  rhat_check <- sum(ss$rhat > 1.1)
  ess_check <- sum(ss$ess_bulk < 400)
  
  # check divergent iteration
  np <- bayesplot::nuts_params(model)
  div_check <- sum(np$Value[np$Parameter=="divergent__"])
  
  # get DHARMa
  dd <- DHARMa::createDHARMa(t(posterior_predict(model)),
                             data[, colname],
                             integerResponse = TRUE)
  # K-S test of uniformity
  ks_test <- DHARMa::testUniformity(dd, plot = FALSE)$statistic
  # check dispersion
  dis_test <- DHARMa::testDispersion(dd, plot=FALSE)$p.value
  # check spat autocorr
  sp_test <- DHARMa::testSpatialAutocorrelation(dd, x = data$X_COORD,
                                                y = data$Y_COORD, plot = FALSE)$statistic[1]
  
  # put together to output
  out <- data.frame(id = file_id,
                    rhat = rhat_check,
                    ess = ess_check,
                    divergence = div_check,
                    ks = ks_test,
                    disp = dis_test,
                    spat = sp_test)
  
  # save this to file
  if(write_file){
    if(file.exists("model_output/model_checks.csv")){
      write.table(out, file = "model_output/model_checks.csv",
                  append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE)
    } else {
      write.table(out, file = "model_output/model_checks.csv",
                  append = FALSE, sep = ",",
                  row.names = FALSE)
    }
  }
  
  return(out)
}

# extract model coefficients and standard error
extract_beta <- function(model, file_name, file_id){
  betas <- fixef(model)[c("fallow_std", "edge_std:fallow_std"), "Estimate"]
  betas_se <- fixef(model)[c("fallow_std", "edge_std:fallow_std"), "Est.Error"]
  
  out <- data.frame(id = file_id,
                    param = c("fallow", "fallow:edge"),
                    beta = betas,
                    beta_se = betas_se)
  
  # save this to file
  if(file.exists(paste0("model_output/", file_name, ".csv"))){
    write.table(out, file = paste0("model_output/", file_name, ".csv"),
                append = TRUE, sep = ",",
                row.names = FALSE, col.names = FALSE)
  } else {
    write.table(out, file = paste0("model_output/", file_name, ".csv"),
                append = FALSE, sep = ",",
                row.names = FALSE)
  }
  
  return(out)
  
}

# compute fallow land effect along the gradient of edge density
# from estimated model coefficients
extract_interaction <- function(model, edge_val, id){
  pp <- posterior::as_draws_df(model,
                    variable = c("b_fallow_std", "b_edge_std:fallow_std"))
  pp_slp <- sapply(edge_val, function(x) pp[,1] + x * pp[,2])
  pp_df <- plyr::ldply(pp_slp, quantile, probs = c(0.025, 0.5, 0.975))
  pp_df$edge <- edge_val
  pp_df$is <- id
  names(pp_df)[2:4] <- c("lci", "med", "hci")
  pp_df[,2:4] <- round(pp_df[,2:4], 2)
  
  return(pp_df)
}

# get spline correlograms of the residuals
# this function takes time to run better save the output
spline_correlog <- function(model, data, coln){
  
  dd <- DHARMa::createDHARMa(t(brms::posterior_predict(model)),
                             data[,coln])
  
  # run spline correlogram, take some time ...
  nn <- ncf::spline.correlog(data$X_COORD,
                             data$Y_COORD,
                             dd$scaledResiduals)
  
  return(nn)
}

# grab the xy coords of the plot together with the dharma residuals
stan_spatres <- function(model, data, colname, file_id = "id"){
  require(DHARMa)
  
  # get DHARMa
  dd <- DHARMa::createDHARMa(t(posterior_predict(model)),
                             data[, colname],
                             integerResponse = TRUE)
  # check spat autocorr
  sp_test <- data.frame(x = data$X_COORD,
                        y = data$Y_COORD, 
                        res = dd$scaledResiduals,
                        id = file_id)
  
  return(sp_test)
}