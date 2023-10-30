## Functions for IPM model
# Survival data
## Clean/tidy the data
survival_data<- filter(survival_data, birth > 1973) # Filter for individuals born after a certain year
# Add stages
stages<- c("C", "J", "YM", "YF", "OM", "OF")
survival_data$stage<- NA
# Males
survival_data[survival_data$sex==2 & survival_data$age > 21, ]$stage<- "OM"
survival_data[survival_data$sex==2 & survival_data$age <= 21 & survival_data$age >= 10, ]$stage<- "YM"
# Females
survival_data[survival_data$sex==1 & survival_data$age > 42, ]$stage<- "OF"
survival_data[survival_data$sex==1 & survival_data$age <= 42 & survival_data$age >= 10, ]$stage<- "YF"
# Calves and juveniles
survival_data[survival_data$age==0, ]$stage<- "C"
survival_data[survival_data$age >=1 & survival_data$age < 10, ]$stage<- "J"
# Remove animals that can't be staged
survival_data<- survival_data[!is.na(survival_data$stage), ]

## Variables and indicies 
survival_last_year<- 2020
animals<- unique(survival_data$id)
n_animals<- length(animals)

# Create a df of individual animals, and their status each year
survival_df<- NULL
# Loop over individuals
for(i in 1:n_animals){
  # Df for animal i
  temp_animal<- filter(survival_data, id==animals[i]) %>%
    select(-age, -sex, -index)
  temp_years<- seq(from=temp_animal$birth[1], to=survival_last_year, by=1) # Vect of years
  temp_status<- rep(1, length=length(temp_years)) # Vector of status (1 = alive, 0 = dead)
  temp_stages<- temp_animal$stage
  last_stage<- temp_stages[length(temp_stages)]
  
  # If the animal died, create a new status vect with added zeros for years dead
  if(1 %in% temp_animal$death){ 
    temp_death<- temp_animal[temp_animal$death==1, ]$year
    temp_death_diff<- survival_last_year - temp_death
    temp_status[(length(temp_status) - temp_death_diff):length(temp_status)]<- 0
    temp_stage_diff<- rep(last_stage, length=temp_death_diff)
    temp_stages<- c(temp_stages, temp_stage_diff)
  }
  # Create the output for animal i to attach to the df
  temp_df<- data.frame(animal = rep(animals[i], length=length(temp_status)),
                       animal_id = rep(i, length=length(temp_status)),
                       year = temp_years,
                       status = temp_status,
                       stage = temp_stages)
  
  if(i==1){survival_df<- temp_df} else{
    survival_df<- rbind(survival_df, temp_df)
  }
}

# Replace stage labels with numeric ids
for(i in 1:6){survival_df[survival_df$stage==stages[i], ]$stage<- i}

## Variables and indicies
# Create matrix of animals and stages, and of the sighting data
first_surv_year<- min(survival_df$year)
last_surv_year<- max(survival_df$year)
surv_years<- seq(from=first_surv_year, to=last_surv_year, by=1)
last_surv_year_id<- length(surv_years)
# Add year IDs
survival_df$year_id<- NA
for(i in 1:length(surv_years)){survival_df[survival_df$year==surv_years[i], ]$year_id<- i}
start_year_add_one<- NULL
stages_mat<- matrix(NA, nrow = n_animals, ncol = length(surv_years))
z_surv<- matrix(NA, nrow = n_animals, ncol = length(surv_years))
for(i in 1:n_animals){
  temp_animal<- survival_df[survival_df$animal_id==i, ]
  temp_year_add_one<- temp_animal$year_id[1] + 1
  start_year_add_one<- c(start_year_add_one, temp_year_add_one)
  stages_mat[i,(temp_animal$year_id[1]:temp_animal$year_id[nrow(temp_animal)])]<- as.numeric(temp_animal$stage)
  z_surv[i,(temp_animal$year_id[1]:temp_animal$year_id[nrow(temp_animal)])]<- as.numeric(temp_animal$status)
}

## Fecundity data pre-processing
# Fecundity variables
n_moms<- length(unique(fec_data$id))
moms_indx<- fec_data$id

# Create data matrix
y_fec<- matrix(NA, nrow = n_moms, ncol = max(fec_data$year))
mom_age<- matrix(NA, nrow = n_moms, ncol = max(fec_data$year))
year_start_fec<- rep(NA, length=n_moms)
year_stop_fec<- rep(NA, length=n_moms)

# Fill birth and age matricies
for(i in 1:n_moms){
  temp_mom<- fec_data[fec_data$id==i, ] # Individual mom data
  temp_first_year<- min(temp_mom$year) # First year in her life
  temp_last_year<- max(temp_mom$year) # Last year in her life
  y_fec[i,(temp_first_year:temp_last_year)]<- temp_mom$birth
  mom_age[i,(temp_first_year:temp_last_year)]<- temp_mom$age
  year_start_fec[i]<- temp_first_year
  year_stop_fec[i]<- temp_last_year
}