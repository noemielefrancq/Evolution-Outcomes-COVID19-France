##################################################################
##                                                              ## 
##                                                              ##
##    Evolution of outcomes for patients hospitalised           ##
##      during the first 9 months of the SARS-CoV-2             ##
##      pandemic in France: A retrospective national            ## 
##                    surveillance data analysis,               ##
##         the first SARS-CoV-2 pandemic wave in France         ##
##                                                              ##
##         https://doi.org/10.1016/j.lanepe.2021.100087         ##
##                    Estimation framework                      ##
##                                                              ##
## Author: Noemie Lefrancq                                      ##
## Email: ncmjl2@cam.ac.uk                                      ##
## Last update: 21/12/2020                                      ##
##                                                              ##
##                                                              ##  
##################################################################

## Load libraries
library(questionr)
library(rstan)
library(lubridate)

#### Load stan model ####################################
model_MCMC <- stan_model(file = 'Model_indiv_data_changeprobsdelays_plostdischarge_final.stan')
#########################################################

## Load data ############################################
data=readRDS('Simulated_individual_trajectories.rds')

## Data format: Large dataframe, each row corresponds to the trajectory of one COVID-19 patient admitted at hospital
## The dataframe contains, for each patient:
###### Sexe 
###### Age
###### first_date (date of admittance in dd/mm/yy)
###### first_date_day (date of admittance in number of days since 01/03/2020 (01/03/2020 = day 1))
###### went_ICU (was the patient admitted to ICU at some point, yes or no)
###### hospi_ICU (time in days between hospital and ICU admissions, if no ICU admission, NA)
###### hospi_arrival_death (time in days between hospital admission and death, if no death, NA)
###### hospi_arrival_discharge (time in  days between hospital admission and discharge, if no discharge, NA)
###### ICU_death (time in days between ICU admission and death, if no death, NA)
###### ICU_discharge (time in days between ICU admission and discharge, if no death, NA)

## Set parameters #######################################

sexes = c('Female', 'Male')

## Age classes for delays distribution
min_age_classes_delay = c(0, 61, 71, 81)
max_age_classes_delay = c(60, 70, 80, 120)

## Age classes for outcome probabilities
min_age_classes_cases = c(0,(seq(40,80,10)+1))
max_age_classes_cases = c(seq(40,80,10), 120)

nb_age_classes_cases = length(min_age_classes_cases) ## number age classes for the cases
nb_age_classes_delays = length(min_age_classes_delay) ## number age classes for the delays
match_cases_delay_age = c(1,1,1,2,3,4) ## match between age classes of cases and delays

max_time = max(as.numeric(data$first_date_day), na.rm = T) # maximum outcome delay that could be reached in the dataset (usually = number of days since beginning of the epidemic)
wind = seq(1, max_time, 1) ## number of windows of observation, here: everyday
nb_obs = length(wind)-1  # number of observation points
time_obs_max = wind[-1]; # for each observation point, max date
time_obs_min = wind[-length(wind)];# for each observation point, min dat

## Delay truncation time
t_truncation = 120  ## truncation time

## Changes in p_lost (prob of losing a record in hospital or icu)
nb_changes_plost = 3
match_changes_windows = match_changes_windows
match_changes_windows_plost = c(1,1,1,1,2,2,2,3,3,3)

## Set periods of time for the changes in outcome probabilities
nb_changes = 10
# v = c(1, 32, 52, 72, 92, max_time)
v = c(1, 32, 52, 72, 92, 122, 153, 184, 214, 245, max_time)
match_changes_windows = c(rep(1, v[2]-v[1]),
                          rep(2, v[3]-v[2]),
                          rep(3, v[4]-v[3]),
                          rep(4, v[5]-v[4]),
                          rep(5, v[6]-v[5]),
                          rep(6, v[7]-v[6]),
                          rep(7, v[8]-v[7]),
                          rep(8, v[9]-v[8]),
                          rep(9, v[10]-v[9]),
                          rep(10, v[11]-v[10]))

## Check dates of switch
for(i in 1:(length(v)-1)){
  print(paste0('Time period ', i, ': ', ymd('2020-03-01')+v[i]-1, ' to ', ymd('2020-03-01')+(v[i+1]-1)))
}

## Set vectors and matrices to organize data #########################
## 1) Data for each observation point
# From hosp to ICU, Death or Discharge
delays_indiv_icu_F = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs))
delays_indiv_icu_M = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs))
delays_indiv_deaths_non_icu_F = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs))
delays_indiv_deaths_non_icu_M = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs))
delays_indiv_discharges_non_icu_F = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs))
delays_indiv_discharges_non_icu_M = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs))
delays_indiv_cases_unobs_distrib_hosp_F = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs)) ## number of patients F in hospital with unknown outcomes for each age, length of stay and observation point
delays_indiv_cases_unobs_distrib_hosp_M = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs)) ## number of patients M in hospital with unknown outcomes for each age, length of stay and observation point
delays_indiv_cases_unobs_lastday_hosp_F = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs)) ## number of patients F in hospital with unknown outcomes for each age, the last day of each observation point
delays_indiv_cases_unobs_lastday_hosp_M = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs)) ## number of patients M in hospital with unknown outcomes for each age, the last day of each observation point
delays_indiv_cases_hosp_F = matrix(NA, nb_age_classes_cases, nb_obs)
delays_indiv_cases_hosp_M = matrix(NA, nb_age_classes_cases, nb_obs)

# From ICU to Death or Discharge
delays_indiv_deaths_icu_F = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs))
delays_indiv_deaths_icu_M = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs))
delays_indiv_discharges_icu_F = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs))
delays_indiv_discharges_icu_M = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs))
delays_indiv_cases_unobs_distrib_icu_F = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs)) ## number of patients F in icu with unknown outcomes for each age, length of stay and observation point
delays_indiv_cases_unobs_distrib_icu_M = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs)) ## number of patients M in icu with unknown outcomes for each age, length of stay and observation point
delays_indiv_cases_unobs_lastday_icu_F = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs)) ## number of patients F in icu with unknown outcomes for each age, the last day of each observation point
delays_indiv_cases_unobs_lastday_icu_M = array(NA, dim = c(nb_age_classes_cases, max_time, nb_obs)) ## number of patients M in icu with unknown outcomes for each age, the last day of each observation point
delays_indiv_cases_icu_F = matrix(NA, nb_age_classes_cases, nb_obs)
delays_indiv_cases_icu_M = matrix(NA, nb_age_classes_cases, nb_obs)

## 2) Overall delay distribution (over the course of the pandemic)
####### This data is not used in the estimation, but useful for subsequent analysis
# From hosp to ICU, Death or Discharge
delay_icu_F = matrix(NA, nb_age_classes_delays, max_time)
delay_icu_M = matrix(NA, nb_age_classes_delays, max_time)
delay_deaths_non_icu_F = matrix(NA, nb_age_classes_delays, max_time)
delay_deaths_non_icu_M = matrix(NA, nb_age_classes_delays, max_time)
delay_discharges_non_icu_F = matrix(NA, nb_age_classes_delays, max_time)
delay_discharges_non_icu_M = matrix(NA, nb_age_classes_delays, max_time)

# From ICU to Death or Discharge
delay_deaths_icu_F = matrix(NA, nb_age_classes_delays, max_time)
delay_deaths_icu_M = matrix(NA, nb_age_classes_delays, max_time)
delay_discharges_icu_F = matrix(NA, nb_age_classes_delays, max_time)
delay_discharges_icu_M = matrix(NA, nb_age_classes_delays, max_time)

## Fill vectors and matrices #############################
################ WOMEN ###################################
data_tmp = data[which(data$sexe == 'Female'),]

## 1) Data for each observation point
for (kkk in 1:length(min_age_classes_cases)){
  for (jjj in 1:nb_obs){
    data_tmp_tmp = data_tmp[which(data_tmp$age>=min_age_classes_cases[kkk] & 
                                    data_tmp$age<= max_age_classes_cases[kkk] &
                                    data_tmp$first_date_day >= time_obs_min[jjj] &
                                    data_tmp$first_date_day < time_obs_max[jjj]),]
    
    ## ICU Death
    hosp_to_death_icu = data_tmp_tmp$ICU_death[which(data_tmp_tmp$went_to_ICU == 'yes')]
    dates_hosp_to_death_icu = as.numeric(levels(as.factor(hosp_to_death_icu)))
    hosp_to_death_tmp = freq(hosp_to_death_icu)$n
    if (length(which(rownames(freq(hosp_to_death_icu)) == 'NA')) >0) hosp_to_death_tmp = hosp_to_death_tmp[-length(hosp_to_death_tmp)] ## enlever le NA de la fin
    hosp_to_death_icu = hosp_to_death_tmp
    tmp = NULL
    for (i in 0:(max_time-1)){
      tmp2 = 0
      if (length(which(dates_hosp_to_death_icu == i))>0) tmp2 = hosp_to_death_icu[which(dates_hosp_to_death_icu == i)]
      tmp = c(tmp, tmp2)
    }
    hosp_to_death_icu = tmp
    
    ## ICU Discharge
    hosp_to_discharge_icu = data_tmp_tmp$ICU_discharge[which(data_tmp_tmp$went_to_ICU == 'yes')]
    dates_hosp_to_discharge_icu = as.numeric(levels(as.factor(hosp_to_discharge_icu)))
    hosp_to_dicharge_tmp = freq(hosp_to_discharge_icu)$n
    if (length(which(rownames(freq(hosp_to_discharge_icu)) == 'NA')) >0) hosp_to_dicharge_tmp = hosp_to_dicharge_tmp[-length(hosp_to_dicharge_tmp)] ## enlever le NA de la fin
    hosp_to_discharge_icu = hosp_to_dicharge_tmp
    tmp = NULL
    for (i in 0:(max_time-1)){
      tmp2 = 0
      if (length(which(dates_hosp_to_discharge_icu == i))>0) tmp2 = hosp_to_discharge_icu[which(dates_hosp_to_discharge_icu == i)]
      tmp = c(tmp, tmp2)
    }
    hosp_to_discharge_icu = tmp
    
    ## Hosp ICU
    hosp_to_icu = data_tmp_tmp$hospi_ICU
    dates_hosp_to_icu = as.numeric(levels(as.factor(hosp_to_icu)))
    hosp_to_icu_tmp = freq(hosp_to_icu)$n
    if (length(which(rownames(freq(hosp_to_icu)) == 'NA')) >0) hosp_to_icu_tmp = hosp_to_icu_tmp[-length(hosp_to_icu_tmp)] ## enlever le NA de la fin
    hosp_to_icu = hosp_to_icu_tmp
    tmp = NULL
    for (i in 0:(max_time-1)){
      tmp2 = 0
      if (length(which(dates_hosp_to_icu == i))>0) tmp2 = hosp_to_icu[which(dates_hosp_to_icu == i)]
      tmp = c(tmp, tmp2)
    }
    hosp_to_icu = tmp
    
    ## Hosp Death no icu
    hosp_to_death_non_icu = data_tmp_tmp$hospi_arrival_death[which(data_tmp_tmp$went_to_ICU == 'no')]
    dates_hosp_to_death_non_icu = as.numeric(levels(as.factor(hosp_to_death_non_icu)))
    hosp_to_death_tmp = freq(hosp_to_death_non_icu)$n
    if (length(which(rownames(freq(hosp_to_death_non_icu)) == 'NA')) >0) hosp_to_death_tmp = hosp_to_death_tmp[-length(hosp_to_death_tmp)] ## enlever le NA de la fin
    hosp_to_death_non_icu = hosp_to_death_tmp
    tmp = NULL
    for (i in 0:(max_time-1)){
      tmp2 = 0
      if (length(which(dates_hosp_to_death_non_icu == i))>0) tmp2 = hosp_to_death_non_icu[which(dates_hosp_to_death_non_icu == i)]
      tmp = c(tmp, tmp2)
    }
    hosp_to_death_non_icu = tmp
    
    ## Hosp discharge no icu
    hosp_to_discharge_non_icu = data_tmp_tmp$hospi_arrival_discharge[which(data_tmp_tmp$went_to_ICU == 'no')]
    dates_hosp_to_discharge_non_icu = as.numeric(levels(as.factor(hosp_to_discharge_non_icu)))
    hosp_to_discharge_tmp = freq(hosp_to_discharge_non_icu)$n
    if (length(which(rownames(freq(hosp_to_discharge_non_icu)) == 'NA')) >0) hosp_to_discharge_tmp = hosp_to_discharge_tmp[-length(hosp_to_discharge_tmp)] ## enlever le NA de la fin
    hosp_to_discharge_non_icu = hosp_to_discharge_tmp
    tmp = NULL
    for (i in 0:(max_time-1)){
      tmp2 = 0
      if (length(which(dates_hosp_to_discharge_non_icu == i))>0) tmp2 = hosp_to_discharge_non_icu[which(dates_hosp_to_discharge_non_icu == i)]
      tmp = c(tmp, tmp2)
    }
    hosp_to_discharge_non_icu = tmp
    
    
    ## Fill the matrices
    # From hosp to ICU, Death or Discharge
    delays_indiv_icu_F[kkk, ,jjj] = hosp_to_icu
    delays_indiv_deaths_non_icu_F[kkk, ,jjj] = hosp_to_death_non_icu
    delays_indiv_discharges_non_icu_F[kkk, ,jjj] = hosp_to_discharge_non_icu
    delays_indiv_cases_hosp_F[kkk, jjj] = nrow(data_tmp_tmp)
    # From ICU to Death or Discharge
    delays_indiv_deaths_icu_F[kkk, ,jjj] = hosp_to_death_icu
    delays_indiv_discharges_icu_F[kkk, ,jjj] = hosp_to_discharge_icu
    delays_indiv_cases_icu_F[kkk, jjj] = nrow(data_tmp_tmp[which(data_tmp_tmp$went_to_ICU == 'yes'),])    
    
    ## Fill the matrices for unknown outcomes
    # At hospital
    d = rep(nrow(data_tmp_tmp), max_time) - cumsum(hosp_to_icu) - cumsum(hosp_to_death_non_icu) - cumsum(hosp_to_discharge_non_icu) ## compute distribution of unknown outcome
    d = c(nrow(data_tmp_tmp), d)
    d = d[-length(d)]
    delays_indiv_cases_unobs_distrib_hosp_F[kkk, ,jjj] = d
    delays_indiv_cases_unobs_lastday_hosp_F[kkk, ,jjj] = rep(0, max_time)
    delays_indiv_cases_unobs_lastday_hosp_F[kkk, max_time - time_obs_min[jjj]+1,jjj] = delays_indiv_cases_unobs_distrib_hosp_F[kkk,max_time - time_obs_min[jjj]+1,jjj]
    if (time_obs_min[jjj] >1){
      delays_indiv_cases_unobs_distrib_hosp_F[kkk, (max_time - time_obs_min[jjj]+2): max_time,jjj] = NA ## After the observation time we don't know the outcomes
    }
    
    # In ICU
    d = rep(nrow(data_tmp_tmp[which(data_tmp_tmp$went_to_ICU == 'yes'),]), max_time)- cumsum(hosp_to_death_icu) - cumsum(hosp_to_discharge_icu) ## compute distribution of unknown outcome
    d = c(nrow(data_tmp_tmp[which(data_tmp_tmp$went_to_ICU == 'yes'),]), d) 
    d = d[-length(d)] 
    delays_indiv_cases_unobs_distrib_icu_F[kkk, ,jjj] = d
    delays_indiv_cases_unobs_lastday_icu_F[kkk, ,jjj] = rep(0, max_time)
    delays_indiv_cases_unobs_lastday_icu_F[kkk, max_time - time_obs_min[jjj] +1,jjj] = delays_indiv_cases_unobs_distrib_icu_F[kkk, max_time - time_obs_min[jjj]+1,jjj]
    if (time_obs_min[jjj] >1){
      delays_indiv_cases_unobs_distrib_icu_F[kkk, (max_time - time_obs_min[jjj]+2): max_time,jjj] = NA ## After the observation time we don't know the outcomes
    }
  }
}

## 2) Overall delay distribution (over the course of the pandemic)
for (kkk in 1:length(min_age_classes_delay)){
  data_tmp_tmp = data_tmp[which(data_tmp$age>=min_age_classes_delay[kkk] & 
                                  data_tmp$age<=max_age_classes_delay[kkk]),]
  
  ## ICU Death
  hosp_to_death_icu = data_tmp_tmp$ICU_death[which(data_tmp_tmp$went_to_ICU == 'yes')]
  dates_hosp_to_death_icu = as.numeric(levels(as.factor(hosp_to_death_icu)))
  hosp_to_death_tmp = freq(hosp_to_death_icu)$n
  if (length(which(rownames(freq(hosp_to_death_icu)) == 'NA')) >0) hosp_to_death_tmp = hosp_to_death_tmp[-length(hosp_to_death_tmp)] ## enlever le NA de la fin
  hosp_to_death_icu = hosp_to_death_tmp
  tmp = NULL
  for (i in 0:(max_time-1)){
    tmp2 = 0
    if (length(which(dates_hosp_to_death_icu == i))>0) tmp2 = hosp_to_death_icu[which(dates_hosp_to_death_icu == i)]
    tmp = c(tmp, tmp2)
  }
  hosp_to_death_icu = tmp
  
  ## ICU Discharge
  hosp_to_discharge_icu = data_tmp_tmp$ICU_discharge[which(data_tmp_tmp$went_to_ICU == 'yes')]
  dates_hosp_to_discharge_icu = as.numeric(levels(as.factor(hosp_to_discharge_icu)))
  hosp_to_dicharge_tmp = freq(hosp_to_discharge_icu)$n
  if (length(which(rownames(freq(hosp_to_discharge_icu)) == 'NA')) >0) hosp_to_dicharge_tmp = hosp_to_dicharge_tmp[-length(hosp_to_dicharge_tmp)] ## enlever le NA de la fin
  hosp_to_discharge_icu = hosp_to_dicharge_tmp
  tmp = NULL
  for (i in 0:(max_time-1)){
    tmp2 = 0
    if (length(which(dates_hosp_to_discharge_icu == i))>0) tmp2 = hosp_to_discharge_icu[which(dates_hosp_to_discharge_icu == i)]
    tmp = c(tmp, tmp2)
  }
  hosp_to_discharge_icu = tmp
  
  ## Hosp ICU
  hosp_to_icu = data_tmp_tmp$hospi_ICU
  dates_hosp_to_icu = as.numeric(levels(as.factor(hosp_to_icu)))
  hosp_to_icu_tmp = freq(hosp_to_icu)$n
  if (length(which(rownames(freq(hosp_to_icu)) == 'NA')) >0) hosp_to_icu_tmp = hosp_to_icu_tmp[-length(hosp_to_icu_tmp)] ## enlever le NA de la fin
  hosp_to_icu = hosp_to_icu_tmp
  tmp = NULL
  for (i in 0:(max_time-1)){
    tmp2 = 0
    if (length(which(dates_hosp_to_icu == i))>0) tmp2 = hosp_to_icu[which(dates_hosp_to_icu == i)]
    tmp = c(tmp, tmp2)
  }
  hosp_to_icu = tmp
  
  ## Hosp Death no icu
  hosp_to_death_non_icu = data_tmp_tmp$hospi_arrival_death[which(data_tmp_tmp$went_to_ICU == 'no')]
  dates_hosp_to_death_non_icu = as.numeric(levels(as.factor(hosp_to_death_non_icu)))
  hosp_to_death_tmp = freq(hosp_to_death_non_icu)$n
  if (length(which(rownames(freq(hosp_to_death_non_icu)) == 'NA')) >0) hosp_to_death_tmp = hosp_to_death_tmp[-length(hosp_to_death_tmp)] ## enlever le NA de la fin
  hosp_to_death_non_icu = hosp_to_death_tmp
  tmp = NULL
  for (i in 0:(max_time-1)){
    tmp2 = 0
    if (length(which(dates_hosp_to_death_non_icu == i))>0) tmp2 = hosp_to_death_non_icu[which(dates_hosp_to_death_non_icu == i)]
    tmp = c(tmp, tmp2)
  }
  hosp_to_death_non_icu = tmp
  
  ## Hosp discharge no icu
  hosp_to_discharge_non_icu = data_tmp_tmp$hospi_arrival_discharge[which(data_tmp_tmp$went_to_ICU == 'no')]
  dates_hosp_to_discharge_non_icu = as.numeric(levels(as.factor(hosp_to_discharge_non_icu)))
  hosp_to_discharge_tmp = freq(hosp_to_discharge_non_icu)$n
  if (length(which(rownames(freq(hosp_to_discharge_non_icu)) == 'NA')) >0) hosp_to_discharge_tmp = hosp_to_discharge_tmp[-length(hosp_to_discharge_tmp)] ## enlever le NA de la fin
  hosp_to_discharge_non_icu = hosp_to_discharge_tmp
  tmp = NULL
  for (i in 0:(max_time-1)){
    tmp2 = 0
    if (length(which(dates_hosp_to_discharge_non_icu == i))>0) tmp2 = hosp_to_discharge_non_icu[which(dates_hosp_to_discharge_non_icu == i)]
    tmp = c(tmp, tmp2)
  }
  hosp_to_discharge_non_icu = tmp
  
  ## Fill the matrices
  # From hosp to ICU, Death or Discharge
  delay_icu_F[kkk,] = hosp_to_icu
  delay_deaths_non_icu_F[kkk,] = hosp_to_death_non_icu
  delay_discharges_non_icu_F[kkk,] = hosp_to_discharge_non_icu
  # From ICU to Death or Discharge
  delay_deaths_icu_F[kkk,] = hosp_to_death_icu
  delay_discharges_icu_F[kkk,] = hosp_to_discharge_icu
}

################ MEN ####################################################################
data_tmp = data[which(data$sexe == 'Male'),]
## 1) Data for each observation point
for (kkk in 1:length(min_age_classes_cases)){
  for (jjj in 1:nb_obs){
    data_tmp_tmp = data_tmp[which(data_tmp$age>=min_age_classes_cases[kkk] & 
                                    data_tmp$age<= max_age_classes_cases[kkk] &
                                    data_tmp$first_date_day >= time_obs_min[jjj] &
                                    data_tmp$first_date_day < time_obs_max[jjj]),]
    
    ## ICU Death
    hosp_to_death_icu = data_tmp_tmp$ICU_death[which(data_tmp_tmp$went_to_ICU == 'yes')]
    dates_hosp_to_death_icu = as.numeric(levels(as.factor(hosp_to_death_icu)))
    hosp_to_death_tmp = freq(hosp_to_death_icu)$n
    if (length(which(rownames(freq(hosp_to_death_icu)) == 'NA')) >0) hosp_to_death_tmp = hosp_to_death_tmp[-length(hosp_to_death_tmp)] ## enlever le NA de la fin
    hosp_to_death_icu = hosp_to_death_tmp
    tmp = NULL
    for (i in 0:(max_time-1)){
      tmp2 = 0
      if (length(which(dates_hosp_to_death_icu == i))>0) tmp2 = hosp_to_death_icu[which(dates_hosp_to_death_icu == i)]
      tmp = c(tmp, tmp2)
    }
    hosp_to_death_icu = tmp
    
    ## ICU Discharge
    hosp_to_discharge_icu = data_tmp_tmp$ICU_discharge[which(data_tmp_tmp$went_to_ICU == 'yes')]
    dates_hosp_to_discharge_icu = as.numeric(levels(as.factor(hosp_to_discharge_icu)))
    hosp_to_dicharge_tmp = freq(hosp_to_discharge_icu)$n
    if (length(which(rownames(freq(hosp_to_discharge_icu)) == 'NA')) >0) hosp_to_dicharge_tmp = hosp_to_dicharge_tmp[-length(hosp_to_dicharge_tmp)] ## enlever le NA de la fin
    hosp_to_discharge_icu = hosp_to_dicharge_tmp
    tmp = NULL
    for (i in 0:(max_time-1)){
      tmp2 = 0
      if (length(which(dates_hosp_to_discharge_icu == i))>0) tmp2 = hosp_to_discharge_icu[which(dates_hosp_to_discharge_icu == i)]
      tmp = c(tmp, tmp2)
    }
    hosp_to_discharge_icu = tmp
    
    ## Hosp ICU
    hosp_to_icu = data_tmp_tmp$hospi_ICU
    dates_hosp_to_icu = as.numeric(levels(as.factor(hosp_to_icu)))
    hosp_to_icu_tmp = freq(hosp_to_icu)$n
    if (length(which(rownames(freq(hosp_to_icu)) == 'NA')) >0) hosp_to_icu_tmp = hosp_to_icu_tmp[-length(hosp_to_icu_tmp)] ## enlever le NA de la fin
    hosp_to_icu = hosp_to_icu_tmp
    tmp = NULL
    for (i in 0:(max_time-1)){
      tmp2 = 0
      if (length(which(dates_hosp_to_icu == i))>0) tmp2 = hosp_to_icu[which(dates_hosp_to_icu == i)]
      tmp = c(tmp, tmp2)
    }
    hosp_to_icu = tmp
    
    ## Hosp Death no icu
    hosp_to_death_non_icu = data_tmp_tmp$hospi_arrival_death[which(data_tmp_tmp$went_to_ICU == 'no')]
    dates_hosp_to_death_non_icu = as.numeric(levels(as.factor(hosp_to_death_non_icu)))
    hosp_to_death_tmp = freq(hosp_to_death_non_icu)$n
    if (length(which(rownames(freq(hosp_to_death_non_icu)) == 'NA')) >0) hosp_to_death_tmp = hosp_to_death_tmp[-length(hosp_to_death_tmp)] ## enlever le NA de la fin
    hosp_to_death_non_icu = hosp_to_death_tmp
    tmp = NULL
    for (i in 0:(max_time-1)){
      tmp2 = 0
      if (length(which(dates_hosp_to_death_non_icu == i))>0) tmp2 = hosp_to_death_non_icu[which(dates_hosp_to_death_non_icu == i)]
      tmp = c(tmp, tmp2)
    }
    hosp_to_death_non_icu = tmp
    
    ## Hosp discharge no icu
    hosp_to_discharge_non_icu = data_tmp_tmp$hospi_arrival_discharge[which(data_tmp_tmp$went_to_ICU == 'no')]
    dates_hosp_to_discharge_non_icu = as.numeric(levels(as.factor(hosp_to_discharge_non_icu)))
    hosp_to_discharge_tmp = freq(hosp_to_discharge_non_icu)$n
    if (length(which(rownames(freq(hosp_to_discharge_non_icu)) == 'NA')) >0) hosp_to_discharge_tmp = hosp_to_discharge_tmp[-length(hosp_to_discharge_tmp)] ## enlever le NA de la fin
    hosp_to_discharge_non_icu = hosp_to_discharge_tmp
    tmp = NULL
    for (i in 0:(max_time-1)){
      tmp2 = 0
      if (length(which(dates_hosp_to_discharge_non_icu == i))>0) tmp2 = hosp_to_discharge_non_icu[which(dates_hosp_to_discharge_non_icu == i)]
      tmp = c(tmp, tmp2)
    }
    hosp_to_discharge_non_icu = tmp
    
    ## Fill the matrices
    # From hosp to ICU, Death or Discharge
    delays_indiv_icu_M[kkk, ,jjj] = hosp_to_icu
    delays_indiv_deaths_non_icu_M[kkk, ,jjj] = hosp_to_death_non_icu
    delays_indiv_discharges_non_icu_M[kkk, ,jjj] = hosp_to_discharge_non_icu
    delays_indiv_cases_hosp_M[kkk, jjj] = nrow(data_tmp_tmp)
    # From ICU to Death or Discharge
    delays_indiv_deaths_icu_M[kkk, ,jjj] = hosp_to_death_icu
    delays_indiv_discharges_icu_M[kkk, ,jjj] = hosp_to_discharge_icu
    delays_indiv_cases_icu_M[kkk, jjj] = nrow(data_tmp_tmp[which(data_tmp_tmp$went_to_ICU == 'yes'),])    
    
    ## Fill the matrices for unknown outcomes
    # At hospital
    d = rep(nrow(data_tmp_tmp), max_time) - cumsum(hosp_to_icu) - cumsum(hosp_to_death_non_icu) - cumsum(hosp_to_discharge_non_icu) ## compute distribution of unknown outcome
    d = c(nrow(data_tmp_tmp), d)
    d = d[-length(d)]
    delays_indiv_cases_unobs_distrib_hosp_M[kkk, ,jjj] = d
    delays_indiv_cases_unobs_lastday_hosp_M[kkk, ,jjj] = rep(0, max_time)
    delays_indiv_cases_unobs_lastday_hosp_M[kkk, max_time - time_obs_min[jjj]+1,jjj] = delays_indiv_cases_unobs_distrib_hosp_M[kkk,max_time - time_obs_min[jjj]+1,jjj]
    if (time_obs_min[jjj] >1){
      delays_indiv_cases_unobs_distrib_hosp_M[kkk, (max_time - time_obs_min[jjj]+2): max_time,jjj] = NA ## After the observation time we don't know the outcomes
    }
    
    # In ICU
    d = rep(nrow(data_tmp_tmp[which(data_tmp_tmp$went_to_ICU == 'yes'),]), max_time)- cumsum(hosp_to_death_icu) - cumsum(hosp_to_discharge_icu) ## compute distribution of unknown outcome
    d = c(nrow(data_tmp_tmp[which(data_tmp_tmp$went_to_ICU == 'yes'),]), d)
    d = d[-length(d)]
    delays_indiv_cases_unobs_distrib_icu_M[kkk, ,jjj] = d
    delays_indiv_cases_unobs_lastday_icu_M[kkk, ,jjj] = rep(0, max_time)
    delays_indiv_cases_unobs_lastday_icu_M[kkk, max_time - time_obs_min[jjj] +1,jjj] = delays_indiv_cases_unobs_distrib_icu_M[kkk, max_time - time_obs_min[jjj]+1,jjj]
    if (time_obs_min[jjj] >1){
      delays_indiv_cases_unobs_distrib_icu_M[kkk, (max_time - time_obs_min[jjj]+2): max_time,jjj] = NA ## After the observation time we don't know the outcomes
    }
  }
}

## 2) Overall delay distribution (over the course of the pandemic)
for (kkk in 1:length(min_age_classes_delay)){
  data_tmp_tmp = data_tmp[which(data_tmp$age>=min_age_classes_delay[kkk] & 
                                  data_tmp$age<=max_age_classes_delay[kkk]),]
  
  ## ICU Death
  hosp_to_death_icu = data_tmp_tmp$ICU_death[which(data_tmp_tmp$went_to_ICU == 'yes')]
  dates_hosp_to_death_icu = as.numeric(levels(as.factor(hosp_to_death_icu)))
  hosp_to_death_tmp = freq(hosp_to_death_icu)$n
  if (length(which(rownames(freq(hosp_to_death_icu)) == 'NA')) >0) hosp_to_death_tmp = hosp_to_death_tmp[-length(hosp_to_death_tmp)] ## enlever le NA de la fin
  hosp_to_death_icu = hosp_to_death_tmp
  tmp = NULL
  for (i in 0:(max_time-1)){
    tmp2 = 0
    if (length(which(dates_hosp_to_death_icu == i))>0) tmp2 = hosp_to_death_icu[which(dates_hosp_to_death_icu == i)]
    tmp = c(tmp, tmp2)
  }
  hosp_to_death_icu = tmp
  
  ## ICU Discharge
  hosp_to_discharge_icu = data_tmp_tmp$ICU_discharge[which(data_tmp_tmp$went_to_ICU == 'yes')]
  dates_hosp_to_discharge_icu = as.numeric(levels(as.factor(hosp_to_discharge_icu)))
  hosp_to_dicharge_tmp = freq(hosp_to_discharge_icu)$n
  if (length(which(rownames(freq(hosp_to_discharge_icu)) == 'NA')) >0) hosp_to_dicharge_tmp = hosp_to_dicharge_tmp[-length(hosp_to_dicharge_tmp)] ## enlever le NA de la fin
  hosp_to_discharge_icu = hosp_to_dicharge_tmp
  tmp = NULL
  for (i in 0:(max_time-1)){
    tmp2 = 0
    if (length(which(dates_hosp_to_discharge_icu == i))>0) tmp2 = hosp_to_discharge_icu[which(dates_hosp_to_discharge_icu == i)]
    tmp = c(tmp, tmp2)
  }
  hosp_to_discharge_icu = tmp
  
  ## Hosp ICU
  hosp_to_icu = data_tmp_tmp$hospi_ICU
  dates_hosp_to_icu = as.numeric(levels(as.factor(hosp_to_icu)))
  hosp_to_icu_tmp = freq(hosp_to_icu)$n
  if (length(which(rownames(freq(hosp_to_icu)) == 'NA')) >0) hosp_to_icu_tmp = hosp_to_icu_tmp[-length(hosp_to_icu_tmp)] ## enlever le NA de la fin
  hosp_to_icu = hosp_to_icu_tmp
  tmp = NULL
  for (i in 0:(max_time-1)){
    tmp2 = 0
    if (length(which(dates_hosp_to_icu == i))>0) tmp2 = hosp_to_icu[which(dates_hosp_to_icu == i)]
    tmp = c(tmp, tmp2)
  }
  hosp_to_icu = tmp
  
  ## Hosp Death no icu
  hosp_to_death_non_icu = data_tmp_tmp$hospi_arrival_death[which(data_tmp_tmp$went_to_ICU == 'no')]
  dates_hosp_to_death_non_icu = as.numeric(levels(as.factor(hosp_to_death_non_icu)))
  hosp_to_death_tmp = freq(hosp_to_death_non_icu)$n
  if (length(which(rownames(freq(hosp_to_death_non_icu)) == 'NA')) >0) hosp_to_death_tmp = hosp_to_death_tmp[-length(hosp_to_death_tmp)] ## enlever le NA de la fin
  hosp_to_death_non_icu = hosp_to_death_tmp
  tmp = NULL
  for (i in 0:(max_time-1)){
    tmp2 = 0
    if (length(which(dates_hosp_to_death_non_icu == i))>0) tmp2 = hosp_to_death_non_icu[which(dates_hosp_to_death_non_icu == i)]
    tmp = c(tmp, tmp2)
  }
  hosp_to_death_non_icu = tmp
  
  ## Hosp discharge no icu
  hosp_to_discharge_non_icu = data_tmp_tmp$hospi_arrival_discharge[which(data_tmp_tmp$went_to_ICU == 'no')]
  dates_hosp_to_discharge_non_icu = as.numeric(levels(as.factor(hosp_to_discharge_non_icu)))
  hosp_to_discharge_tmp = freq(hosp_to_discharge_non_icu)$n
  if (length(which(rownames(freq(hosp_to_discharge_non_icu)) == 'NA')) >0) hosp_to_discharge_tmp = hosp_to_discharge_tmp[-length(hosp_to_discharge_tmp)] ## enlever le NA de la fin
  hosp_to_discharge_non_icu = hosp_to_discharge_tmp
  tmp = NULL
  for (i in 0:(max_time-1)){
    tmp2 = 0
    if (length(which(dates_hosp_to_discharge_non_icu == i))>0) tmp2 = hosp_to_discharge_non_icu[which(dates_hosp_to_discharge_non_icu == i)]
    tmp = c(tmp, tmp2)
  }
  hosp_to_discharge_non_icu = tmp
  
  ## Fill the matrices
  # From hosp to ICU, Death or Discharge
  delay_icu_M[kkk,] = hosp_to_icu
  delay_deaths_non_icu_M[kkk,] = hosp_to_death_non_icu
  delay_discharges_non_icu_M[kkk,] = hosp_to_discharge_non_icu
  # From ICU to Death or Discharge
  delay_deaths_icu_M[kkk,] = hosp_to_death_icu
  delay_discharges_icu_M[kkk,] = hosp_to_discharge_icu
}


############# Dataset for MCMC #########################################################
data.MCMC = list(
  ## Data format
  max_time = max_time,  ## maximum outcome delay that could be reached in the dataset (usually = number of days since beginning of the epidemic)
  nb_obs = nb_obs,  ## number of observation points
  nb_age_classes_cases = nb_age_classes_cases,  ## number age classes for the cases
  nb_age_classes_delays = nb_age_classes_delays,  ## number age classes for the delays
  match_cases_delay_age = match_cases_delay_age,  ## match between age classes of cases and delays
  
  ## Delay truncation time
  t_truncation = t_truncation,  ## truncation time
  
  ## Changes over the course of the pandemic
  nb_changes = nb_changes,  ## number of windows of changes
  match_changes_windows = match_changes_windows,  ## match between observation points and windows of time
  
  nb_changes_plost = nb_changes_plost, ## number of windows of changes for p lost
  match_changes_windows_plost = match_changes_windows_plost, ## match between observation points and windows of time
  
  ## Cases data
  ## From hosp admission to ICU, Death or Discharge
  delays_indiv_icu_F = delays_indiv_icu_F,  ## number of patients F admitted to icu for each age, delay and observation point
  delays_indiv_icu_M = delays_indiv_icu_M,  ## number of patients M admitted to icu for each age, delay and observation point
  delays_indiv_deaths_non_icu_F = delays_indiv_deaths_non_icu_F,  ## number of patients F deaths in hosp (no icu) for each age, delay and observation point
  delays_indiv_deaths_non_icu_M = delays_indiv_deaths_non_icu_M,  ## number of patients M deaths in hosp (no icu) for each age, delay and observation point
  delays_indiv_discharges_non_icu_F = delays_indiv_discharges_non_icu_F,  ## number of patients F discharges from hosp (no icu) for each age, delay and observation point
  delays_indiv_discharges_non_icu_M = delays_indiv_discharges_non_icu_M, ## number of patients M discharges from hosp (no icu) for each age, delay and observation point
  delays_indiv_cases_unobs_hosp_F = delays_indiv_cases_unobs_lastday_hosp_F,  ## number of patients F in hospital with unknown outcomes for each age, length of stay and observation point
  delays_indiv_cases_unobs_hosp_M = delays_indiv_cases_unobs_lastday_hosp_M,  ## number of patients M in hospital with unknown outcomes for each age, length of stay and observation point

  
  ## From ICU admisson to Death or Discharge
  delays_indiv_deaths_icu_F = delays_indiv_deaths_icu_F,  ## number of patients F deaths after icu for each age, delay and observation point
  delays_indiv_deaths_icu_M = delays_indiv_deaths_icu_M,  ##  number of patients M deaths after icu for each age, delay and observation point
  delays_indiv_discharges_icu_F = delays_indiv_discharges_icu_F,  ## number of patients F discharges after icu for each age, delay and observation point
  delays_indiv_discharges_icu_M = delays_indiv_discharges_icu_M,  ## number of patients M discharges after icu for each age, delay and observation point
  delays_indiv_cases_unobs_icu_F = delays_indiv_cases_unobs_lastday_icu_F,  ## number of patients F in icu with unknown outcomes for each age, length of stay and observation point
  delays_indiv_cases_unobs_icu_M = delays_indiv_cases_unobs_lastday_icu_M  ## number of patients M in icu with unknown outcomes for each age, length of stay and observation point
)

############# Supplementary dataset, for subsequent analyses #########################
data.complementary = list(
  ## Number of cases entering the hospital at each observation point
  delays_indiv_cases_hosp_F = delays_indiv_cases_hosp_F,
  delays_indiv_cases_hosp_M = delays_indiv_cases_hosp_M,
  delays_indiv_cases_icu_F = delays_indiv_cases_icu_F,
  delays_indiv_cases_icu_M = delays_indiv_cases_icu_M,
  
  ## Global distribution of the delays
  delay_icu_F = delay_icu_F,
  delay_icu_M = delay_icu_M,
  delay_deaths_non_icu_F = delay_deaths_non_icu_F,
  delay_deaths_non_icu_M = delay_deaths_non_icu_M,
  delay_discharges_non_icu_F = delay_discharges_non_icu_F,
  delay_discharges_non_icu_M = delay_discharges_non_icu_M,
  delay_deaths_icu_F = delay_deaths_icu_F,
  delay_deaths_icu_M = delay_deaths_icu_M,
  delay_discharges_icu_F = delay_discharges_icu_F,
  delay_discharges_icu_M = delay_discharges_icu_M,
  
  ## Distribution of unknown outcomes
  delays_indiv_cases_unobs_distrib_hosp_F = delays_indiv_cases_unobs_distrib_hosp_F,
  delays_indiv_cases_unobs_distrib_hosp_M = delays_indiv_cases_unobs_distrib_hosp_M,
  delays_indiv_cases_unobs_distrib_icu_F = delays_indiv_cases_unobs_distrib_icu_F,
  delays_indiv_cases_unobs_distrib_icu_M = delays_indiv_cases_unobs_distrib_icu_M
)
######################################################################################

#########################################################
## Run MCMC
#########################################################
fit_delay =  sampling(model_MCMC, data = data.MCMC, chains = 3, 
                      cores = 3,iter= 5000, control = list(adapt_delta = 0.9, max_treedepth = 20))
#########################################################
## Save results
#########################################################
filename=paste0('Results', today(), '.rds', sep = '_')
fit = list(fit=fit_delay,
           data= data.MCMC,
           data_complementary = data.complementary)    
saveRDS(fit, filename)
#########################################################
