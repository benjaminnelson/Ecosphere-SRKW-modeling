## Optimized version of model data
setwd("C:/Users/nelso/OneDrive/BNelson/Consulting/Oceans initiative/Working/srkw_modeling")

# Load libraries
library(R2jags)
library(dplyr)
library(ggplot2)

## Clear everything ##
rm(list=ls()) 
gc() 

# Set randoms #
set.seed(99)

## Read in the data
removal_data<- read.csv("ipm_data_ver_3.csv", header = TRUE) # Removal capture data, by year and stage
survival_data<- read.csv("mortality_data.csv", header = TRUE) # Mortality/survival data
stage_data<- read.csv("stage_data_ver_2.csv", header = TRUE) # Annual observed stage composition data
rem_stage_start_stop<- read.csv("stage_start_stop_ver_3.csv", header=TRUE) # Age ranges for stage (multinomial indexing)
surv_age_to_stage<- read.csv("survival_age_to_stage.csv", header=TRUE) # Age-to-stage vectors for survival
salmon_data<- read.csv("chinook_indices_ver_2.csv", header=TRUE) # Data frame of PSC-WCVI, Raincoast, and FRAM abundance indices
psc_sims<- read.csv("psc_sims.csv", header = FALSE) # PSC simulated abundances
raincoast_sims<- read.csv("raincoast_sims.csv", header = FALSE)
fram_sims<- read.csv("fram_sims.csv", header = FALSE) # FRAM simulated abundances
fec_data<- read.csv("fec_data.csv", header=TRUE) # Data frame of mothers and their births, ages
nrkws<- read.csv("nrkws.csv", header=FALSE)[,2] # Time-series of NRKW abundance, adapted from Towers et al. 2020 (first ~20 are imputed linearly)

## Z-score the salmon indices
sd_std<- function(x, sd_group) {std=(x-mean(sd_group, na.rm=TRUE))/(sd(sd_group, na.rm=TRUE)); return(std)}
# Z-score simulated salmon indices with the mean and sd from observed time-series
raincoast_sims<- sd_std(raincoast_sims, salmon_data[,2])
psc_sims<- sd_std(psc_sims, salmon_data[,3])
fram_sims<- sd_std(fram_sims, salmon_data[,4])

# Z-score observed salmon abundance indices
salmon_data[,2]<- sd_std(salmon_data[,2], salmon_data[,2])
salmon_data[,3]<- sd_std(salmon_data[,3], salmon_data[,3])
salmon_data[,4]<- sd_std(salmon_data[,4], salmon_data[,4])

# Salmon indices
salmon_obs<- salmon_data$psc
salmon_sims<- psc_sims

## Pre-process the survival and fecundity data, define survival- and fecundity-related variables
source("ipm_survival_optim.r")

# Pars for Wishart dist.
diagp<- diag(2)
p<- 2

## Non-survival-related variables from the data
# Prep removal data matrix
removal_data<- removal_data[,-1]
removal_data[,"total_caps"]<- apply(removal_data, 1, sum)
removal_data<- as.matrix(removal_data)

# Prep stage-composition data matrix
stage_data<- stage_data[,-1]
stage_data[,"total_n"]<- apply(stage_data, 1, sum)
stage_data<- as.matrix(stage_data)
n_stage<- ncol(stage_data)-1
n_stage_obs<- nrow(stage_data)

# Age-to-stage indices
max_male<- 22 # Male max age (first age of old males)
max_female<- 43 # Female max age (first age of old females)

# Age-to-stage vectors for male and female survivals
male_surv_stage<- as.vector(na.omit(surv_age_to_stage$male_stage))
female_surv_stage<- as.vector(na.omit(surv_age_to_stage$female_stage))

# Stage age ranges for removal/stage indexing
n_rem_stage<- nrow(rem_stage_start_stop)
rem_stage_starts<- rem_stage_start_stop$start
rem_stage_stops<- rem_stage_start_stop$stop
rem_stage_seq<- seq(from=1, to=n_rem_stage, by=2) # Stage sequence for indexing the male and female stages in the removal likelihood
n_rem_years<- nrow(removal_data)

# Model time-steps
model_years<- seq(from=1940, to=2020, by=1)
n_years<- length(model_years)

# Dirichlet prior vectors
alpha_par<- 1 # Prior for the Dirichlet distributions; alpha values of 1 are a uniform distribution
alpha_seq<- (rem_stage_start_stop$stop - rem_stage_start_stop$start)+1
alpha_vect<- rep(alpha_par, length=max(alpha_seq))

## JAGS model function
model_func<- function(lag_surv_in = 0, lag_fec_in = 0, nrkw = 0, density_dep = 1, low_k_in = 95, up_k_in = 175){
  # Set lag years between salmon and survival/fecundity
  lag_surv<- lag_surv_in
  lag_fec<- lag_fec_in
  # Set code for whether NRKWs are included in density-dep process (0 = Off; 1 = On)
  nrkw_code<- nrkw
  # Set code for density dep process (0 = Off; 1 = On)
  density_dep_code<- density_dep
  # Set lower and upper bounds for K
  low_k<- low_k_in
  up_k<- up_k_in
  
  ##################################################################################################################################################################
  ## JAGS modeling block
  ##################################################################################################################################################################
  # Define model
  source("jags_ipm_optim_ver_2.r")
  
  # MCMC settings
  #iter = 2000; burn = iter/2; thins = 1; chains = 6     
  iter = 25000; burn = iter/2; thins = 10; chains = 8
  #iter = 40000; burn = iter/2; thins = 10; chains = 5 
  #iter = 20000; burn = iter/2; thins = 10; chains = 10
  
  # Pars to monitor during sampling
  jags_pars=c("n_tot", "b_stage", "b_fec", "lp_z", "lp_fec", "lp_n")
  #jags_pars=c("n_tot", "stage_v", "b_stage", "b_fec", "prop_k", "K", "theta", "n_zero", "year_effect", "annual_surv", "annual_fec")
  
  # Data
  jags_data = list( 
    # Survival data
    "n_animals",
    "model_years",
    "first_surv_year",
    "last_surv_year_id", 
    "stages_mat", 
    "z_surv", 
    "start_year_add_one", 
    "male_surv_stage",
    "female_surv_stage",
    "p",
    "diagp",
    # Demographic data
    "stage_data",
    "n_years",
    "max_male",
    "max_female",
    "n_stage",
    "n_stage_obs",
    "nrkws",
    # Removal data
    "n_rem_years",
    "n_rem_stage",
    "removal_data",
    "rem_stage_starts",
    "rem_stage_stops",
    "rem_stage_seq",
    "alpha_vect",
    "alpha_seq",
    # Fecundity data
    "n_moms", 
    "y_fec", 
    "mom_age", 
    "year_start_fec", 
    "year_stop_fec",
    # Salmon data
    "salmon_obs",
    "salmon_sims",
    "lag_surv",
    "lag_fec",
    # Model settings
    "nrkw_code",
    "density_dep_code",
    "low_k",
    "up_k")
  
  # Model file
  model_loc = paste("model.txt", sep="")
  
  # Fit the model w/ JAGS
  start_time<- Sys.time(); print(start_time)
  model = do.call(jags.parallel, list(jags_data, inits = NULL,
                                      parameters.to.save= jags_pars,
                                      model.file=model_loc,
                                      n.chains = chains,
                                      n.burnin = burn,
                                      n.thin = thins,
                                      n.iter = iter,
                                      jags.seed = 99,
                                      DIC = TRUE))
  end_time<- Sys.time()
  print(end_time-start_time) # Duration of model-fitting in JAGS
  
  attach.jags(model)
  # Print JAGS model summary
  model$BUGSoutput$summary
  # Print DIC
  model$BUGSoutput$DIC
  # Object containing posterior summaries
  posts<- model$BUGSoutput$summary
  
  OUTS<- NULL
  OUTS$model<- model
  OUTS$posts<- posts
  OUTS$dic<- model$BUGSoutput$DIC
  return(OUTS)
}

## Run candidate models
#mod_1<- model_func(lag_surv = 0, lag_fec = 0, nrkw = 0, density_dep = 0, low_k_in = 95, up_k_in = 150); saveRDS(mod_1, "mod_1.rds")
#mod_2<- model_func(lag_surv = 0, lag_fec = 0, nrkw = 0, density_dep = 1, low_k_in = 95, up_k_in = 150); saveRDS(mod_2, "mod_2.rds")
#mod_3<- model_func(lag_surv = 0, lag_fec = 0, nrkw = 1, density_dep = 1, low_k_in = 350, up_k_in = (150 + 400)); saveRDS(mod_3, "mod_3.rds")
#mod_4<- model_func(lag_surv = 0, lag_fec = 0, nrkw = 1, density_dep = 1, low_k_in = 350, up_k_in = (150 + 500)); saveRDS(mod_4, "mod_4.rds")
#mod_5<- model_func(lag_surv = 0, lag_fec = 1, nrkw = 0, density_dep = 0, low_k_in = 95, up_k_in = 150); saveRDS(mod_5, "mod_5.rds")
#mod_6<- model_func(lag_surv = 0, lag_fec = 1, nrkw = 0, density_dep = 1, low_k_in = 95, up_k_in = 150); saveRDS(mod_6, "mod_6.rds")
#mod_7<- model_func(lag_surv = 0, lag_fec = 1, nrkw = 1, density_dep = 1, low_k_in = 350, up_k_in = (150 + 400)); saveRDS(mod_7, "mod_7.rds")
#mod_8<- model_func(lag_surv = 0, lag_fec = 1, nrkw = 1, density_dep = 1, low_k_in = 350, up_k_in = (150 + 500)); saveRDS(mod_8, "mod_8.rds")
#mod_9<- model_func(lag_surv = 1, lag_fec = 1, nrkw = 0, density_dep = 0, low_k_in = 95, up_k_in = 150); saveRDS(mod_9, "mod_9.rds")
#mod_10<- model_func(lag_surv = 1, lag_fec = 1, nrkw = 0, density_dep = 1, low_k_in = 95, up_k_in = 150); saveRDS(mod_10, "mod_10.rds")
#mod_11<- model_func(lag_surv = 1, lag_fec = 1, nrkw = 1, density_dep = 1, low_k_in = 350, up_k_in = (150 + 400)); saveRDS(mod_11, "mod_11.rds")
#mod_12<- model_func(lag_surv = 1, lag_fec = 1, nrkw = 1, density_dep = 1, low_k_in = 350, up_k_in = (150 + 500)); saveRDS(mod_12, "mod_12.rds")
#mod_13<- model_func(lag_surv = 1, lag_fec = 0, nrkw = 0, density_dep = 0, low_k_in = 95, up_k_in = 150); saveRDS(mod_13, "mod_13.rds")
#mod_14<- model_func(lag_surv = 1, lag_fec = 0, nrkw = 0, density_dep = 1, low_k_in = 95, up_k_in = 150); saveRDS(mod_14, "mod_14.rds")
#mod_15<- model_func(lag_surv = 1, lag_fec = 0, nrkw = 1, density_dep = 1, low_k_in = 350, up_k_in = (150 + 400)); saveRDS(mod_15, "mod_15.rds")
#mod_16<- model_func(lag_surv = 1, lag_fec = 0, nrkw = 1, density_dep = 1, low_k_in = 350, up_k_in = (150 + 500)); saveRDS(mod_16, "mod_16.rds")