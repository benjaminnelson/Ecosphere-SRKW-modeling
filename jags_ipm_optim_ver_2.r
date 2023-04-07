# Same as Ver 23, but with update to harvest sub-model
model = cat("
model{
  #######################################################################################################################################################
  # Survival and fecundity sub-model; Ward et al. 2016
  #######################################################################################################################################################
  ## Priors
  p_detect <- 1.0; # Detection probability = 1 for SRKW; may be different seasonally (e.g., winter off the coast)
  
  # Priors for survival coefficients
  b_unif[1] ~ dunif(0,1); # Calves
  b_unif[2] ~ dunif(0,1); # Juveniles
  b_unif[3] ~ dunif(0,1); # Young males
  b_unif[4] ~ dunif(0,1); # Young females
  b_unif[5] ~ dunif(0,1); # Old males
  b_unif[6] ~ dunif(0,1); # Old females
  b_unif[7] ~ dunif(0,1); # Salmon effect

  # Transform survival coefficients; loop over the six stages and salmon effect
  for(i in 1:(n_stage+1)){b_stage[i] <-log(b_unif[i] / (1-b_unif[i]));}

  # Priors for fecundity coefficients
  for(i in 1:6){b_fec[i] ~ dnorm(0,1);} # Normal prior on regression coefs, incl. salmon effect b_fec[6]
  
  ## Pre-data year effects
  for(j in 1:(first_surv_year-model_years[1]-1)){temp_effect[j]<-0;}
  
  ## Year effects
  # Prior for rho par
  rho ~ dunif(0,1);
  # Prior for year effect variance/sd
  # Wishart
  tau_year[1:p,1:p] ~ dwish(diagp,p);
  sigma_year[1:p,1:p] <- inverse(tau_year[1:p,1:p]);
  rho_year <- sigma_year[1,2]/sqrt(sigma_year[1,1]*sigma_year[2,2]);

  # Year effects as correlated random walk
  # First year is zero
  year_effect[1,1]<- 0; 
  year_effect[1,2]<- 0;
  for(j in 2:last_surv_year_id){year_effect[j,1:2] ~ dmnorm(rho * year_effect[j-1,1:2], tau_year[1:2,1:2]);}

  ## Annual survival estimates
  # Pre survival/fec data
  for(j in 1:(first_surv_year-model_years[1]-1)){
    for(s in 1:n_stage){
      annual_surv[j,s]<- 1/(1+exp(-1*(temp_effect[j] + b_stage[s])));
    }
  }

  # During data period
  for(j in (first_surv_year-model_years[1]):(n_years-1)){
    for(s in 1:n_stage){
      annual_surv[j,s]<- 1/(1+exp(-1*(year_effect[j-(first_surv_year-model_years[1]-1), 1] + (salmon_surv[j-(first_surv_year-model_years[1]-1)]*b_stage[7]) + 
                          b_stage[s])));
    }
  }
  
  ## Annual fecundity estimates
  # Pre survival/fec data
  for(j in 1:(first_surv_year-model_years[1]-1)){
    for(a in 1:(max_female-1)){
      annual_fec[j,a]<- 1/(1+exp(-1*(temp_effect[j] + fec_age_vect[a])));
    }
  }
  # During data period
  for(j in (first_surv_year-model_years[1]):(n_years-1)){
    for(a in 1:(max_female-1)){
      annual_fec[j,a]<- 1/(1+exp(-1*(year_effect[j-(first_surv_year-model_years[1]-1), 2] + (salmon_fec[j-(first_surv_year-model_years[1]-1)]*b_fec[6]) + 
                        fec_age_vect[a])));   
    }
  }

  ## SURVIVAL DATA
  ## Calculate likelihood of all the individual animal encounter histories
  # Draw salmon dataset for survival
  salmon_rand ~ dunif(1,100);
  salmon_surv[1:(5+lag_surv)]<- salmon_sims[1:(5+lag_surv),round(salmon_rand)]
  salmon_surv[(6+lag_surv):last_surv_year_id]<- salmon_obs[1:((last_surv_year_id-5)-lag_surv)];
  # Draw salmon dataset for fecundity
  salmon_fec[1:(5+lag_fec)]<- salmon_sims[1:(5+lag_fec),round(salmon_rand)]
  salmon_fec[(6+lag_fec):last_surv_year_id]<- salmon_obs[1:((last_surv_year_id-5)-lag_fec)];

  # Loop over individual animals
  for(i in 1:n_animals){
    # Loop over years in each animal's encounter history
    for(j in start_year_add_one[i]:last_surv_year_id){
      # Each animal's sighting history starts at their birth
      mu_z[i,j]<- z_surv[i,j-1] * exp(year_effect[j, 1] + (salmon_surv[j]*b_stage[7]) + b_stage[stages_mat[i,j]]) / (1 + exp(year_effect[j, 1] + 
                  (salmon_surv[j]*b_stage[7]) + b_stage[stages_mat[i,j]]));
      z_surv[i,j] ~ dbern(mu_z[i,j]); # z = 0 if dead, 1 if alive
      y_surv[i,j] ~ dbern(z_surv[i,j] * p_detect);
      lp_z[i,j]<- logdensity.bern(z_surv[i,j], mu_z[i,j]);
    } # End j
  } # End i  

  ## FECUNDITY DATA
  ## Fit the model to calving data
  # Loop over individual moms
  for(i in 1:n_moms){
    # Loop over years in each mom's life
    for(j in year_start_fec[i]:year_stop_fec[i]){
      # probability of animal i in time t
      logit(prob[i, j])<- year_effect[j, 2] + (salmon_fec[j]*b_fec[6]) + (b_fec[1] + (b_fec[2] * mom_age[i, j]) + (b_fec[3] * pow(mom_age[i, j], 2)) + 
                          (b_fec[4]*pow(mom_age[i, j], 3)) + (b_fec[5] * pow(mom_age[i, j], 4)));
      y_fec[i,j] ~ dbern(prob[i, j]);
      lp_fec[i,j]<- logdensity.bern(y_fec[i,j], prob[i, j]);
    } # End t
  } # End i

  # Predict fec at age
  for(a in 1:(max_female-1)){
    fec_age_vect[a]<- b_fec[1] + (b_fec[2] * a) + (b_fec[3] * pow(a, 2)) + (b_fec[4]*pow(a, 3)) + (b_fec[5] * pow(a, 4));
  }

  #######################################################################################################################################################
  # Population process model; age- and sex- structured; leading pars: K; n_zero (initial abundance); theta (Pella-Tomlinson shape parameter)
  #######################################################################################################################################################
  ## Priors
  n_zero ~ dunif(80, 150); # Prior for abundance in 1940; must be <= to K, if NRKWs aren't included
  K ~ dunif(low_k, up_k); # Prior for K 
  theta ~ dnorm(5.0, 2.5) T(1,); # Truncated-normal prior for the Pella-Tomlinson parameter
  
  ## Initial age-distributions
  # Males
  sship_male[1]<- 1;
  for(a in 2:(max_male+1)){
    sship_male[a]<- sship_male[a-1] * b_unif[male_surv_stage[a-1]];
  }
  # Females
  sship_female[1]<- 1;
  for(a in 2:(max_female+1)){
    sship_female[a]<- sship_female[a-1] * b_unif[female_surv_stage[a-1]];
  }
  # Initial sex ratios
  male_rat<- sum(sship_male)/sum(sum(sship_male)+sum(sship_female)); # Initial proportion of males in the population
  female_rat<- 1-male_rat; # Initial proportion of females in the population
  
  ### Initialize numbers-at-age
  ## Fecundity
  # Fec females
  fec_females[1]<- sum(females[1,10:(max_female-1)] %*% annual_fec[1,10:(max_female-1)]); # Sumproduct of females and fecundity
  f[1]<- max(0.001, (1-((n_zero+nrkws[1]*nrkw_code)/K)^theta)) * density_dep_code + step(0-density_dep_code)*1; # Density-dep effects

  # Neonates
  neonates[1] <- fec_females[1] * f[1] * 0.50; # Fec females x density dep x sex ratio = births

  # Calves
  males[1,1] <- b_unif[1] * neonates[1];
  females[1,1] <- b_unif[1] * neonates[1];

  ## Initial juveniles
  # Ages 2-9
  for(a in 2:9){
    males[1, a] <- (male_rat*n_zero)*(b_unif[2]^(a-1))*(1-b_unif[2]);
    females[1, a] <- (female_rat*n_zero)*(b_unif[2]^(a-1))*(1-b_unif[2]);
  }
  
  ## Initial male adults
  # Age 10
  males[1, 10] <- (male_rat*n_zero)*(b_unif[2]^(10-1))*(1-b_unif[2]);
  # Ages 11-22
  for(a in 11:max_male) {
    males[1, a] <- (male_rat*n_zero)*(b_unif[3]^(a-1))*(1-b_unif[3]); # Young males
  } 
  # Age 23 (Plus group)
  males[1, (max_male+1)] <- (male_rat*n_zero)*(b_unif[5]^((max_male+1)-1)); # Old males > 22

  ## Initial female adults
  # Age 10
  females[1, 10] <- (female_rat*n_zero)*(b_unif[2]^(10-1))*(1-b_unif[2]);
  # Ages 11-43
  for(a in 11:max_female) {
    females[1, a] <- (female_rat*n_zero)*(b_unif[4]^(a-1))*(1-b_unif[4]); # Young females and first old female (age 43)
  } 
  # Age 44 (plus group)
  females[1, (max_female+1)] <- (female_rat*n_zero)*(b_unif[6]^((max_female+1)-1)); # Old females > 43
  
  ## Total initial abundances
  total_males[1]<- sum(males[1,]);
  total_females[1]<- sum(females[1,]);
  calves[1]<- males[1,1] + females[1,1];
  juveniles[1]<- sum(females[1, 2:9]) + sum(males[1, 2:9]);
  young_males[1]<- sum(males[1, 10:(max_male-1)]);
  young_females[1]<- sum(females[1, 10:(max_female-1)]);
  old_males[1]<- sum(males[1, max_male:(max_male+1)]);
  old_females[1]<- sum(females[1, max_female:(max_female+1)]);
  
  ## Stage abundance matrix
  stage_v[1,1]<- calves[1];
  stage_v[1,2]<- juveniles[1];
  stage_v[1,3]<- young_males[1];
  stage_v[1,4]<- young_females[1];
  stage_v[1,5]<- old_males[1];
  stage_v[1,6]<- old_females[1];
  
  ## Initial total population size
  n_tot[1] <- total_males[1] + total_females[1];
  prop_k[1]<- n_tot[1]/K;

  ### Project population forward in time; include removals 
  # No removals until 1962
  for(t in 1:22){
    for(a in 1:(max_male-1)) {male_removals[t,a]<-0;}
    for(aa in 1:(max_female-1)) {female_removals[t,aa]<-0;}
  }
  
  ## Loop over years (1940-2020)
  for(t in 1:(n_years-1)){
    ## Fecundity
    # Fec females
    fec_females[t+1]<- sum(females[t+1,10:(max_female-1)] %*% annual_fec[t,10:(max_female-1)]);
    f[t+1]<- max(0.001, (1-((n_tot[t]+nrkws[t]*nrkw_code)/K)^theta)) * density_dep_code + step(0-density_dep_code)*1;

    # Neonates
    neonates[t+1] <- fec_females[t+1] * f[t+1] * 0.50;
    
    # Calves
    males[t+1, 1] <- annual_surv[t,1] * neonates[t+1];
    females[t+1, 1] <- annual_surv[t,1] * neonates[t+1];
    
    # Male and female juveniles
    # Age 2
    males[t+1, 2] <- annual_surv[t,2] * males[t,1];
    females[t+1, 2] <- annual_surv[t,2] * females[t,1];
    # Ages 3-9
    for(a in 3:9){
      males[t+1, a] <- annual_surv[t,2] * max(0, males[t, a-1] - min(males[t, a-1], male_removals[t, a-1]));
      females[t+1, a] <- annual_surv[t,2] * max(0, females[t, a-1] - min(females[t, a-1], female_removals[t, a-1]));
    }
    
    ## Adult males (ages 10-22)
    males[t+1, 10] <- annual_surv[t,2] * max(0, males[t, 9] - min(males[t, 9], male_removals[t, 9]));
    for(a in 11:max_male){
      males[t+1, a] <- annual_surv[t,3] * max(0, males[t, a-1] - min(males[t, a-1], male_removals[t, a-1]));
    }
    # Age 23 (plus group)
    males[t+1, (max_male+1)] <- annual_surv[t,5] * (max(0, males[t, max_male]) + max(0, males[t, (max_male+1)])); 

    ## Adult females (ages 10-43)
    females[t+1, 10] <- annual_surv[t,2] * max(0, females[t, 9] - min(females[t, 9], female_removals[t, 9]));
    for(a in 11:max_female){
      females[t+1, a] <- annual_surv[t,4] * max(0, females[t, a-1] - min(females[t, a-1], female_removals[t, a-1]));
    }
    # Age 44 (plus group)
    females[t+1, (max_female+1)] <- annual_surv[t,6] * (max(0, females[t, max_female]) + max(0, females[t, (max_female+1)])); 
    
    ## Total annual abundances
    total_males[t+1]<- sum(males[t+1,]);
    total_females[t+1]<- sum(females[t+1,]);
    calves[t+1]<- males[t+1,1] + females[t+1,1];
    juveniles[t+1]<- sum(females[t+1, 2:9]) + sum(males[t+1, 2:9]);
    young_males[t+1]<- sum(males[t+1, 10:(max_male-1)]);
    young_females[t+1]<- sum(females[t+1, 10:(max_female-1)]);
    old_males[t+1]<- sum(males[t+1, max_male:(max_male+1)]);
    old_females[t+1]<- sum(females[t+1, max_female:(max_female+1)]);
    
    ## Stages matrix
    stage_v[t+1, 1]<- calves[t+1];
    stage_v[t+1, 2]<- juveniles[t+1];
    stage_v[t+1, 3]<- young_males[t+1];
    stage_v[t+1, 4]<- young_females[t+1];
    stage_v[t+1, 5]<- old_males[t+1];
    stage_v[t+1, 6]<- old_females[t+1];
    
    # Predicted total annual population sizes
    n_tot[t+1]<- total_males[t+1] + total_females[t+1]; 
    prop_k[t+1]<- n_tot[t+1]/K;

  } # End t (population process)
  
  # No removals after 1977
  for(t in 39:n_years){
    for(a in 1:(max_male-1)) {male_removals[t,a]<-0;}
    for(aa in 1:(max_female-1)) {female_removals[t,aa]<-0;}
  }
  
  ## Stage-comp/abundance likelihood calcs (obs start in 1977)
  # Loop over years with stage-comp data
  for(t in 1:n_stage_obs){
    # Poisson likelihood
    stage_data[t,(n_stage+1)] ~ dpois(n_tot[(t+((first_surv_year+2)-model_years[1]-1))]);
    lp_n[t]<-logdensity.pois(stage_data[t,(n_stage+1)], n_tot[(t+((first_surv_year+2)-model_years[1]-1))]);
  } # End t (likelihood years)

  #######################################################################################################################################################
  # Live removals/captures sub-model
  #######################################################################################################################################################
  ## Predicted removals/captures
  # Priors for obs equation
  # Loop over years with removal/capture data
  for(t in 1:n_rem_years){
  # No calf removals
    male_age_comp[t,1]<-0; 
    female_age_comp[t,1]<-0;
    male_removals[t+22,1]<-0;
    female_removals[t+22,1]<-0;
    for(s in rem_stage_seq){
      # Male removals, by age
      male_age_comp[t,(rem_stage_starts[s]:rem_stage_stops[s])] ~ ddirch(alpha_vect[1:alpha_seq[s]]);
      #male_removals[(t+22), (rem_stage_starts[s]:rem_stage_stops[s])] ~ dmulti(male_age_comp[t,(rem_stage_starts[s]:rem_stage_stops[s])], removal_data[t,s]);
      for(a in (rem_stage_starts[s]:rem_stage_stops[s])){
        male_removals[t+22,a]<- male_age_comp[t,a] * removal_data[t,s]; 
      } # End a
      #pred_removals[t,s]<- sum(male_removals[(t+22), (rem_stage_starts[s]:rem_stage_stops[s])]);
      # Female removals, by age
      female_age_comp[t,(rem_stage_starts[s+1]:rem_stage_stops[s+1])] ~ ddirch(alpha_vect[1:alpha_seq[s+1]]);
      #female_removals[(t+22), (rem_stage_starts[s+1]:rem_stage_stops[s+1])] ~ dmulti(female_age_comp[t,(rem_stage_starts[s+1]:rem_stage_stops[s+1])], removal_data[t,s+1]);
      for(a in (rem_stage_starts[s+1]:rem_stage_stops[s+1])){
        female_removals[t+22,a]<- female_age_comp[t,a] * removal_data[t,s+1]; 
      } # End a
      #pred_removals[t,(s+1)]<- sum(female_removals[(t+22), (rem_stage_starts[s+1]:rem_stage_stops[s+1])]);
    } # End s
    # Multinomial likelihood for removal/capture stage-composition data
    #removal_data[t, 1:n_rem_stage] ~ dmulti(pred_removals[t, 1:n_rem_stage], removal_data[t, (n_rem_stage+1)]); # Obs and predicted stage-comp vectors; col 9 in stage_data is [observed] total removals in year t
  } # End t
  
} # End JAGS model

", file = "model.txt")