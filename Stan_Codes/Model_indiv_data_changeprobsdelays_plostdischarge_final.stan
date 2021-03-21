// ##################################################################
// ##                                                              ## 
// ##                                                              ##
// ##    Evolution of outcomes for patients hospitalised           ##
// ##      during the first 9 months of the SARS-CoV-2             ##
// ##      pandemic in France: A retrospective national            ## 
// ##                    surveillance data analysis,               ##
// ##         the first SARS-CoV-2 pandemic wave in France         ##
// ##                                                              ##
// ##         https://doi.org/10.1016/j.lanepe.2021.100087         ##
// ##                    Estimation framework                      ##
// ##                         stan code                            ##
// ##                                                              ##
// ## Author: Noemie Lefrancq                                      ##
// ## Email: ncmjl2@cam.ac.uk                                      ##
// ## Last update: 21/12/2020                                      ##
// ##                                                              ##
// ##                                                              ##  
// ##################################################################

functions {
  // Accesory functions
	vector linspaced_vector(int min, int max){
		vector[max-min+1] a;
		for (i in 1:(max-min+1)){
			a[i] = min+i-1;
		}
	return (a);
  }

  // Trunctated mixture distribution (lognormal + exp)
  // cdf
  vector truncatedlognormal_exp_cdf(vector y1, int size, real mean, real sigma, int max_value, real p, real lambda){
  	vector[size] a;
  	for (i in 1:size){
  	 if (y1[i] == 0) {a[i] = 0;}
  	 if (y1[i] > 0 )  {a[i] = lognormal_cdf(y1[i], mean, sigma)/lognormal_cdf(max_value, mean, sigma);}
  	 if(i>=max_value) a[i] = 1;
  	}
  	a = (1-p)*(1-exp(-lambda*(y1))) +  p*a;
  	return(a);
  }
  // pmf
  vector truncatedlognormal_exp_pmf(vector y1, int size, real mean, real sigma, int max_value, real p, real lambda){
  	vector[size] a;
  	for (i in 1:size){
  	 if (y1[i] > 0) {a[i] = exp(lognormal_lpdf(y1[i]| mean, sigma))/lognormal_cdf(max_value, mean, sigma);}
  	 if (i>max_value) {a[i] = 0;}
  	}
  	a = (1-p)*(lambda*exp(-lambda*to_vector(y1))) +  p*a;
  	if (y1[1] == 0) a[1] = 0;
  	return a;
   }

  // Trunctated lognormal
  // cdf
  vector truncatedlognormal_cdf(vector y1, int size, real mean, real sigma, int max_value){
	vector[size] a;
	  for (i in 1:size){
	    if (y1[i] == 0) {a[i] = 0;}
	    if (y1[i] > 0) {a[i] = lognormal_cdf(y1[i], mean, sigma)/lognormal_cdf(max_value, mean, sigma);}
	    if(i>=max_value) a[i] = 1;
	  }
	return(a);
  }
  //pmf
  vector truncatedlognormal_pmf(vector y1, int size, real mean, real sigma, int max_value){
  	vector[size] a;
  	for (i in 1:size){
  	 if (y1[i] == 0) {a[i] = 0;}
  	 if (y1[i] > 0) {a[i] = exp(lognormal_lpdf(y1[i]| mean, sigma))/lognormal_cdf(max_value, mean, sigma);}
  	 if (y1[i]>=max_value) {a[i] = 0;}
  	}
  	return a;
  }

  // Trunctated inflated exponential distribution
  // cdf
  vector truncatedexp_inflated_cdf(vector y1, int size, real lambda, real p0, int max_value){
	  vector[size] a;
  	a = (1-p0)*(1-exp(-lambda*(y1)));
  	if (y1[1] == 0) {a[2:size] = a[2:size] + p0;}
  	if (y1[1] >= 1) {a = a + p0;}
  	return(a);
  }
  //pmf
  vector truncatedexp_inflated_pmf(vector y1, int size, real lambda, real p0, int max_value){
  	vector[size] a;
  	a = (1-p0)*(lambda*exp(-lambda*to_vector(y1)));
  	// if (y1[2] == 1) {a[2] = a[2] + p0;}
	  if (y1[1] == 0) {a[1] = a[1] + p0;}
  	return a;
  }

  // Likelyhood functions
  // lccdf for hosp data
  real likelyhood_lpmf(int[] y, int size, vector delay, real proba){
    real total= 0.0;
    for(i in 1:size){ 
      if(delay[i] == 0) {
        total = total + y[i]*log(0.0000000001);
      }
      if(delay[i]*proba>0) {
        total = total + y[i]*log(delay[i]*proba);      
      }
  	}
  	return  total;
  }
  // lccdf for hosp data  
  real likelyhood_unknown_hosp_lccdf(real[] y, int size, vector delay_icu, vector delay_death, vector delay_discharge, real proba_icu, real proba_death, real proba_discharge){
    real total= 0.0;
    total = total + y[1]*log(1);
    for(i in 1:(size-1)){ 
        if(1 - delay_icu[i]*proba_icu - delay_death[i]*proba_death - delay_discharge[i]*proba_discharge <= 0){
  	      total = total + y[i+1]*log(0.0000000001);
  	    }
  	    if(1 - delay_icu[i]*proba_icu - delay_death[i]*proba_death - delay_discharge[i]*proba_discharge >0) {
          total = total + y[i+1]*log(1 - delay_icu[i]*proba_icu - delay_death[i]*proba_death - delay_discharge[i]*proba_discharge);
        }
  	}
  	return  total;
  }
  // lccdf for icu data
  real likelyhood_unknown_icu_lccdf(real[] y, int size, vector delay_death, vector delay_discharge, real proba_death, real proba_discharge){
    real total= 0.0;
    total = total + y[1]*log(1);
    for(i in 1:(size-1)){
  	    if(1 - delay_death[i]*proba_death - delay_discharge[i]*proba_discharge<= 0){
  	      total = total + y[i+1]*log(0.0000000001);
  	    }
  	    if(1 - delay_death[i]*proba_death - delay_discharge[i]*proba_discharge > 0) {
          total = total + y[i+1]*log(1 - delay_death[i]*proba_death - delay_discharge[i]*proba_discharge);
        }
  	}
  	return  total;
  }
}

data {
  // Data format
  int max_time; // number of days since beginning of the epidemic
  int nb_obs; // number of windows
  int nb_age_classes_cases; // number age classes for the cases
  int nb_age_classes_delays; // number age classes for the global delays
  int match_cases_delay_age[nb_age_classes_cases]; // could be updated with different vectorss for each differetn delays 

  // Delay truncation time
  int t_truncation; // truncation time
  
  // Changes over the course of the pandemic
  int nb_changes; // number of windows of changes, for the probablities and delays
  int match_changes_windows[nb_obs]; // match between observation points and windows of time
  int nb_changes_plost; // number of windows of changes, for p_lost
  int match_changes_windows_plost[nb_changes]; // match between observation points and windows of time
  
  // Cases data
  // From hosp admisson to ICU, Death or Discharge
  int delays_indiv_icu_F[nb_age_classes_cases, max_time, nb_obs];
  int delays_indiv_icu_M[nb_age_classes_cases, max_time, nb_obs];
  int delays_indiv_deaths_non_icu_F[nb_age_classes_cases, max_time, nb_obs];
  int delays_indiv_deaths_non_icu_M[nb_age_classes_cases, max_time, nb_obs];
  int delays_indiv_discharges_non_icu_F[nb_age_classes_cases, max_time, nb_obs];
  int delays_indiv_discharges_non_icu_M[nb_age_classes_cases, max_time, nb_obs];
  real delays_indiv_cases_unobs_hosp_F[nb_age_classes_cases, max_time, nb_obs];
  real delays_indiv_cases_unobs_hosp_M[nb_age_classes_cases, max_time, nb_obs];
  // int delays_indiv_cases_hosp_F[nb_age_classes_cases, nb_obs];
  // int delays_indiv_cases_hosp_M[nb_age_classes_cases, nb_obs];
  
  // From ICU admisson to Death or Discharge
  int delays_indiv_deaths_icu_F[nb_age_classes_cases, max_time, nb_obs];
  int delays_indiv_deaths_icu_M[nb_age_classes_cases, max_time, nb_obs];
  int delays_indiv_discharges_icu_F[nb_age_classes_cases, max_time, nb_obs];
  int delays_indiv_discharges_icu_M[nb_age_classes_cases, max_time, nb_obs];
  real delays_indiv_cases_unobs_icu_F[nb_age_classes_cases, max_time, nb_obs];
  real delays_indiv_cases_unobs_icu_M[nb_age_classes_cases, max_time, nb_obs];
  // int delays_indiv_cases_icu_F[nb_age_classes_cases, nb_obs];
  // int delays_indiv_cases_icu_M[nb_age_classes_cases, nb_obs];

}

parameters {
  // Delays 
  
  // Hosp admision to ICU: inflated exponential distribution
  real<lower = 0.05, upper=5> lambda_icu[nb_changes];  // exponential rate (not varying across age&sex) 
  real<lower = 0.0001, upper=1> p0[nb_changes]; // 0 inflation (not varying across age&sex, varying over the epidemic) 
  
  // Hosp admision to death, non ICU: mixture lognormal+expo
  real<lower = 1, upper=3.5> mu_death_non_icu[nb_changes]; // log mean for the lognormal (not varying across age&sex) 
  real<lower = 0.0001, upper=3> sigma_death_non_icu; // log sd for the lognormal (not varying across age&sex) 
  real<lower = 0.0001, upper=1> rho_death_non_icu[nb_changes]; // proportion of slow deaths (not varying across age&sex, varying over the epidemic) 
  // real<lower = 0.2, upper=5> lambda_death_non_icu; // exponential rate (fixed)
   
  // Hosp admision to discharge, non ICU: mixture lognormal+expo
  real<lower = 1, upper=3.5> mu_discharge_non_icu_F[nb_age_classes_delays, nb_changes]; // log mean for the lognormal F (varying across ages groups and time) 
  real<lower = 1, upper=3.5> mu_discharge_non_icu_M[nb_age_classes_delays, nb_changes]; // log mean for the lognormal M (varying across ages groups and time) 
  real<lower = 0.0001, upper=3> sigma_discharge_non_icu_F[nb_age_classes_delays]; // log sd for the lognormal F(varying across ages groups) 
  real<lower = 0.0001, upper=3> sigma_discharge_non_icu_M[nb_age_classes_delays]; // log sd for the lognormal M(varying across ages groups) 
  real<lower = 0.0001, upper=1> rho_discharge_non_icu[nb_changes]; // proportion of slow discharges (not varying across age&sex, varying over the epidemic) 
  // real<lower = 0.2, upper=5> lambda_discharge_non_icu; // exponential rate (fixed)
  
  // ICU admision to Death: lognormal distribution
  real<lower = 1, upper=5> mu_death_icu_F[nb_age_classes_delays, nb_changes]; // log mean for the lognormal F (varying across ages groups and time) 
  real<lower = 1, upper=5> mu_death_icu_M[nb_age_classes_delays, nb_changes]; // log mean for the lognormal M (varying across ages groups and time)
  real<lower = 0.0001, upper=3> sigma_death_icu_F[nb_age_classes_delays]; // log sd for the lognormal F(varying across ages groups)
  real<lower = 0.0001, upper=3> sigma_death_icu_M[nb_age_classes_delays]; // log sd for the lognormal M(varying across ages groups)
  
  // ICU admision to Discharge: lognormal distribution
  real<lower = 1, upper=5> mu_discharge_icu_F[nb_age_classes_delays, nb_changes]; // log mean for the lognormal F (varying across ages groups and time)
  real<lower = 1, upper=5> mu_discharge_icu_M[nb_age_classes_delays, nb_changes]; // log mean for the lognormal M (varying across ages groups and time) 
  real<lower = 0.0001, upper=3> sigma_discharge_icu_F[nb_age_classes_delays]; // log sd for the lognormal F(varying across ages groups)
  real<lower = 0.0001, upper=3> sigma_discharge_icu_M[nb_age_classes_delays]; // log sd for the lognormal M(varying across ages groups) 
  
  // Probabilities
  // hosp
  simplex[3] p_hosp_F_unscaled[nb_age_classes_cases, nb_changes]; // probs from hosp admission to icu, death and discharge F
  simplex[3] p_hosp_M_unscaled[nb_age_classes_cases, nb_changes]; // probs from hosp admission to icu, death and discharge M
  
  // ICU
  simplex[2] p_icu_F_unscaled[nb_age_classes_cases, nb_changes]; // probs from icu admission to death and discharge F
  simplex[2] p_icu_M_unscaled[nb_age_classes_cases, nb_changes]; // probs from icu admission to death and discharge M
  
  // plost
  real<lower = 0.000001, upper=0.2> p_lost_hosp[nb_changes_plost]; // prob of losing a record in hospital
  real<lower = 0.000001, upper=0.2> p_lost_icu[nb_changes_plost]; // prob of losing a record in icu
}

transformed parameters{
  // Distributions pmf 
  // From hosp admision to ICU, Death or Discharge
  vector[t_truncation] distrib_icu_F_exp[nb_age_classes_delays, nb_changes];
  vector[t_truncation] distrib_icu_M_exp[nb_age_classes_delays, nb_changes];
  vector[t_truncation] distrib_deaths_non_icu_F_exp[nb_age_classes_delays, nb_changes];
  vector[t_truncation] distrib_deaths_non_icu_M_exp[nb_age_classes_delays, nb_changes];
  vector[t_truncation] distrib_discharges_non_icu_F_exp[nb_age_classes_delays, nb_changes];
  vector[t_truncation] distrib_discharges_non_icu_M_exp[nb_age_classes_delays, nb_changes];

  // From ICU admission to Death or Discharge
  vector[t_truncation]  distrib_deaths_icu_F_exp[nb_age_classes_delays, nb_changes];
  vector[t_truncation]  distrib_deaths_icu_M_exp[nb_age_classes_delays, nb_changes];
  vector[t_truncation]  distrib_discharges_icu_F_exp[nb_age_classes_delays, nb_changes];
  vector[t_truncation]  distrib_discharges_icu_M_exp[nb_age_classes_delays, nb_changes];
  
  // Distributions cdf
  // From hosp admission to ICU, Death or Discharge
  vector[t_truncation] distrib_icu_F_exp_cdf[nb_age_classes_delays, nb_changes];
  vector[t_truncation] distrib_icu_M_exp_cdf[nb_age_classes_delays, nb_changes];
  vector[t_truncation] distrib_deaths_non_icu_F_exp_cdf[nb_age_classes_delays, nb_changes];
  vector[t_truncation] distrib_deaths_non_icu_M_exp_cdf[nb_age_classes_delays, nb_changes];
  vector[t_truncation] distrib_discharges_non_icu_F_exp_cdf[nb_age_classes_delays, nb_changes];
  vector[t_truncation] distrib_discharges_non_icu_M_exp_cdf[nb_age_classes_delays, nb_changes];

  // From ICU admission to Death or Discharge
  vector[t_truncation]  distrib_deaths_icu_F_exp_cdf[nb_age_classes_delays, nb_changes];
  vector[t_truncation]  distrib_deaths_icu_M_exp_cdf[nb_age_classes_delays, nb_changes];
  vector[t_truncation]  distrib_discharges_icu_F_exp_cdf[nb_age_classes_delays, nb_changes];
  vector[t_truncation]  distrib_discharges_icu_M_exp_cdf[nb_age_classes_delays, nb_changes];
  
  // Probabilities
  // hosp
  vector[3] p_hosp_F_scaled[nb_age_classes_cases, nb_changes];
  vector[3] p_hosp_M_scaled[nb_age_classes_cases, nb_changes];
  
  // ICU
  vector[2] p_icu_F_scaled[nb_age_classes_cases, nb_changes];
  vector[2] p_icu_M_scaled[nb_age_classes_cases, nb_changes];

  //Tmp variable
  real s;
  
  // Constant parameters
  real lambda_discharge_non_icu;
  real lambda_death_non_icu;
  
  lambda_death_non_icu = 2;
  lambda_discharge_non_icu = 2;
  
  // Scaling the probabilities accounting for prob of losing a record in hospital/icu
  for (age_class in 1:(nb_age_classes_cases)){
    for (n in 1:nb_changes){
        p_hosp_F_scaled[age_class,n] = p_hosp_F_unscaled[age_class,n];
        p_hosp_M_scaled[age_class,n] = p_hosp_M_unscaled[age_class,n]; 
        p_icu_F_scaled[age_class,n] = p_icu_F_unscaled[age_class,n];
        p_icu_M_scaled[age_class,n] = p_icu_M_unscaled[age_class,n]; 
        
        p_hosp_F_scaled[age_class,n,3] = p_hosp_F_unscaled[age_class,n,3]*(1-p_lost_hosp[match_changes_windows_plost[n]]);
        p_hosp_M_scaled[age_class,n,3] = p_hosp_M_unscaled[age_class,n,3]*(1-p_lost_hosp[match_changes_windows_plost[n]]); 
        p_icu_F_scaled[age_class,n,2] = p_icu_F_unscaled[age_class,n,2]*(1-p_lost_icu[match_changes_windows_plost[n]]);
        p_icu_M_scaled[age_class,n,2] = p_icu_M_unscaled[age_class,n,2]*(1-p_lost_icu[match_changes_windows_plost[n]]); 
    }
  }
  
  // Computing the distributions of the delays for age class
  for (age_class in 1:(nb_age_classes_delays)){
    
    // Distributions that are varying through the pandemic
    for (wind in 1:nb_changes){
      // Hosp to ICU
      distrib_icu_F_exp[age_class, wind] = truncatedexp_inflated_cdf(linspaced_vector(1, t_truncation), t_truncation, lambda_icu[wind], p0[wind], t_truncation) - truncatedexp_inflated_cdf(linspaced_vector(0, t_truncation-1), t_truncation, lambda_icu[wind], p0[wind], t_truncation);
      distrib_icu_M_exp[age_class, wind] = truncatedexp_inflated_cdf(linspaced_vector(1, t_truncation), t_truncation, lambda_icu[wind], p0[wind], t_truncation) - truncatedexp_inflated_cdf(linspaced_vector(0, t_truncation-1), t_truncation, lambda_icu[wind], p0[wind], t_truncation);
      distrib_icu_F_exp_cdf[age_class, wind] = truncatedexp_inflated_cdf(linspaced_vector(0, t_truncation-1), t_truncation, lambda_icu[wind], p0[wind], t_truncation);
      distrib_icu_M_exp_cdf[age_class, wind] = truncatedexp_inflated_cdf(linspaced_vector(0, t_truncation-1), t_truncation, lambda_icu[wind], p0[ wind], t_truncation);

     // Hosp to death no ICU
      s = sigma_death_non_icu;
      distrib_deaths_non_icu_F_exp_cdf[age_class, wind] = truncatedlognormal_exp_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_death_non_icu[wind], s, t_truncation, rho_death_non_icu[wind], lambda_death_non_icu);
      distrib_deaths_non_icu_F_exp[age_class, wind] = truncatedlognormal_exp_cdf(linspaced_vector(1, t_truncation), t_truncation, mu_death_non_icu[wind], s, t_truncation, rho_death_non_icu[wind], lambda_death_non_icu) - truncatedlognormal_exp_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_death_non_icu[wind], s, t_truncation, rho_death_non_icu[wind], lambda_death_non_icu);
      s = sigma_death_non_icu;
      distrib_deaths_non_icu_M_exp_cdf[age_class, wind] = truncatedlognormal_exp_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_death_non_icu[wind], s, t_truncation, rho_death_non_icu[wind], lambda_death_non_icu);
      distrib_deaths_non_icu_M_exp[age_class, wind] = truncatedlognormal_exp_cdf(linspaced_vector(1, t_truncation), t_truncation, mu_death_non_icu[wind], s, t_truncation, rho_death_non_icu[wind], lambda_death_non_icu) - truncatedlognormal_exp_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_death_non_icu[wind], s, t_truncation, rho_death_non_icu[wind], lambda_death_non_icu);
    
      // Hosp to discharge no ICU
      s = sigma_discharge_non_icu_F[age_class];
      distrib_discharges_non_icu_F_exp_cdf[age_class, wind] = truncatedlognormal_exp_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_discharge_non_icu_F[age_class, wind], s, t_truncation, rho_discharge_non_icu[wind], lambda_discharge_non_icu);
      distrib_discharges_non_icu_F_exp[age_class, wind] = truncatedlognormal_exp_cdf(linspaced_vector(1, t_truncation), t_truncation, mu_discharge_non_icu_F[age_class, wind], s, t_truncation, rho_discharge_non_icu[wind], lambda_discharge_non_icu) - truncatedlognormal_exp_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_discharge_non_icu_F[age_class, wind], s, t_truncation, rho_discharge_non_icu[wind], lambda_discharge_non_icu);
      s = sigma_discharge_non_icu_M[age_class];
      distrib_discharges_non_icu_M_exp_cdf[age_class, wind] = truncatedlognormal_exp_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_discharge_non_icu_M[age_class, wind], s, t_truncation, rho_discharge_non_icu[wind], lambda_discharge_non_icu);
      distrib_discharges_non_icu_M_exp[age_class, wind] = truncatedlognormal_exp_cdf(linspaced_vector(1, t_truncation), t_truncation, mu_discharge_non_icu_M[age_class, wind], s, t_truncation, rho_discharge_non_icu[wind], lambda_discharge_non_icu) - truncatedlognormal_exp_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_discharge_non_icu_M[age_class, wind], s, t_truncation, rho_discharge_non_icu[wind], lambda_discharge_non_icu);

      // ICU to death
      s = sigma_death_icu_F[age_class];
      distrib_deaths_icu_F_exp_cdf[age_class, wind] = truncatedlognormal_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_death_icu_F[age_class, wind], s, t_truncation);
      distrib_deaths_icu_F_exp[age_class, wind] = truncatedlognormal_cdf(linspaced_vector(1, t_truncation), t_truncation, mu_death_icu_F[age_class, wind], s, t_truncation) - truncatedlognormal_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_death_icu_F[age_class, wind], s, t_truncation);
      s = sigma_death_icu_M[age_class];
      distrib_deaths_icu_M_exp_cdf[age_class, wind] = truncatedlognormal_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_death_icu_M[age_class, wind], s, t_truncation);
      distrib_deaths_icu_M_exp[age_class, wind] = truncatedlognormal_cdf(linspaced_vector(1, t_truncation), t_truncation, mu_death_icu_M[age_class, wind], s, t_truncation) - truncatedlognormal_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_death_icu_M[age_class, wind], s, t_truncation);

      // ICU to discharge
      s = sigma_discharge_icu_F[age_class];
      distrib_discharges_icu_F_exp_cdf[age_class, wind] = truncatedlognormal_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_discharge_icu_F[age_class, wind], s, t_truncation);
      distrib_discharges_icu_F_exp[age_class, wind] = truncatedlognormal_cdf(linspaced_vector(1, t_truncation), t_truncation, mu_discharge_icu_F[age_class, wind], s, t_truncation) - truncatedlognormal_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_discharge_icu_F[age_class, wind], s, t_truncation);
      s = sigma_discharge_icu_M[age_class];
      distrib_discharges_icu_M_exp_cdf[age_class, wind] = truncatedlognormal_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_discharge_icu_M[age_class, wind], s, t_truncation);
      distrib_discharges_icu_M_exp[age_class, wind] = truncatedlognormal_cdf(linspaced_vector(1, t_truncation), t_truncation, mu_discharge_icu_M[age_class, wind], s, t_truncation) - truncatedlognormal_cdf(linspaced_vector(0, t_truncation-1), t_truncation, mu_discharge_icu_M[age_class, wind], s, t_truncation);

      // Making sure there is no 0
      for (i in 1:t_truncation){
        if(distrib_icu_F_exp[age_class, wind, i] <0.0000001) distrib_icu_F_exp[age_class, wind, i] = 0.0000001;
        if(distrib_icu_M_exp[age_class, wind, i] <0.0000001) distrib_icu_M_exp[age_class, wind, i] = 0.0000001;
        if(distrib_deaths_non_icu_F_exp[age_class, wind, i] <0.0000001) distrib_deaths_non_icu_F_exp[age_class, wind, i] = 0.0000001;
        if(distrib_deaths_non_icu_M_exp[age_class, wind, i] <0.0000001) distrib_deaths_non_icu_M_exp[age_class, wind, i] = 0.0000001;
        if(distrib_discharges_non_icu_F_exp[age_class, wind, i] <0.0000001) distrib_discharges_non_icu_F_exp[age_class, wind, i] = 0.0000001;;
        if(distrib_discharges_non_icu_M_exp[age_class, wind, i] <0.0000001) distrib_discharges_non_icu_M_exp[age_class, wind, i] = 0.0000001;
        
        if(distrib_deaths_icu_F_exp[age_class, wind, i] <0.0000001) distrib_deaths_icu_F_exp[age_class, wind, i] = 0.0000001;
        if(distrib_deaths_icu_M_exp[age_class, wind, i] <0.0000001) distrib_deaths_icu_M_exp[age_class, wind, i] = 0.0000001;
        if(distrib_discharges_icu_F_exp[age_class, wind, i] <0.0000001) distrib_discharges_icu_F_exp[age_class, wind, i] = 0.0000001;
        if(distrib_discharges_icu_M_exp[age_class, wind, i] <0.0000001) distrib_discharges_icu_M_exp[age_class, wind, i] = 0.0000001;
      }
            
      distrib_icu_F_exp[age_class, wind] = distrib_icu_F_exp[age_class,wind]/sum(distrib_icu_F_exp[age_class,wind]);
      distrib_icu_M_exp[age_class, wind] = distrib_icu_M_exp[age_class,wind]/sum(distrib_icu_M_exp[age_class,wind]);
      distrib_deaths_non_icu_F_exp[age_class, wind] = distrib_deaths_non_icu_F_exp[age_class,wind]/sum(distrib_deaths_non_icu_F_exp[age_class,wind]);
      distrib_deaths_non_icu_M_exp[age_class, wind] = distrib_deaths_non_icu_M_exp[age_class,wind]/sum(distrib_deaths_non_icu_M_exp[age_class,wind]);
      distrib_discharges_non_icu_F_exp[age_class, wind] = distrib_discharges_non_icu_F_exp[age_class,wind]/sum(distrib_discharges_non_icu_F_exp[age_class,wind]);
      distrib_discharges_non_icu_M_exp[age_class, wind] = distrib_discharges_non_icu_M_exp[age_class,wind]/sum(distrib_discharges_non_icu_M_exp[age_class,wind]);
      
      distrib_deaths_icu_F_exp[age_class, wind] = distrib_deaths_icu_F_exp[age_class, wind]/sum(distrib_deaths_icu_F_exp[age_class, wind]);
      distrib_deaths_icu_M_exp[age_class, wind] = distrib_deaths_icu_M_exp[age_class, wind]/sum(distrib_deaths_icu_M_exp[age_class, wind]) ;
      distrib_discharges_icu_F_exp[age_class, wind] = distrib_discharges_icu_F_exp[age_class, wind]/sum(distrib_discharges_icu_F_exp[age_class, wind]);
      distrib_discharges_icu_M_exp[age_class, wind] = distrib_discharges_icu_M_exp[age_class, wind]/sum(distrib_discharges_icu_M_exp[age_class, wind]) ;
    }
  }
}

model {
    // Priors
  
  // Hosp admision to ICU: inflated exponential distribution
  lambda_icu ~ uniform(0.05,5);
  p0 ~ uniform(0.0001,1);
  
  // Hosp admision to death, non ICU: mixture lognormal+expo
  mu_death_non_icu ~ uniform(1,3.5);
  sigma_death_non_icu ~ uniform(0.0001,3);
  rho_death_non_icu ~ uniform(0.0001,1);
  
  // Hosp admision to discharge, non ICU: mixture lognormal+expo
  for (i in 1:nb_changes){
    mu_discharge_non_icu_F[,i] ~ cauchy(2.5,1); // for hosp to death WITHOUT icu
    mu_discharge_non_icu_M[,i] ~ cauchy(2.5,1); // for hosp to death WITHOUT icu
  }
  sigma_discharge_non_icu_F ~ uniform(0.0001,3);
  sigma_discharge_non_icu_M ~ uniform(0.0001,3);
  rho_discharge_non_icu ~ uniform(0.0001,1);
  
  // ICU admision to Death: lognormal distribution
  for (i in 1:nb_changes){
    mu_death_icu_F[,i] ~ cauchy(2.5,1); // for hosp to death WITH icu F
    mu_death_icu_M[,i] ~ cauchy(2.5,1); // for hosp to death WITH icu M
  }
  sigma_death_icu_F ~ uniform(0.0001,3); 
  sigma_death_icu_M ~ uniform(0.0001,3);
  
  // ICU admision to Discharge: lognormal distribution
  for (i in 1:nb_changes){
    mu_discharge_icu_F[,i] ~ cauchy(2.5,1); // for hosp to death WITH icu F
    mu_discharge_icu_M[,i] ~ cauchy(2.5,1); // for hosp to death WITH icu M
  }
  sigma_discharge_icu_F ~ uniform(0.0001,3); // for hosp to death WITH icu F
  sigma_discharge_icu_M ~ uniform(0.0001,3); // for hosp to death WITH icu M


  // Likelihood
  // Individual cases
  for (obs_i in 1:nb_obs){ 
    for(age_class in 1:(nb_age_classes_cases)){
      // Hosp admissions
      // to ICU
      target += likelyhood_lpmf(delays_indiv_icu_F[age_class, 1:t_truncation, obs_i]| t_truncation, distrib_icu_F_exp[match_cases_delay_age[age_class], match_changes_windows[obs_i],], p_hosp_F_scaled[age_class, match_changes_windows[obs_i],1]);
      target += likelyhood_lpmf(delays_indiv_icu_M[age_class, 1:t_truncation, obs_i]| t_truncation, distrib_icu_M_exp[match_cases_delay_age[age_class], match_changes_windows[obs_i],], p_hosp_M_scaled[age_class, match_changes_windows[obs_i],1]);
      // to death
      target += likelyhood_lpmf(delays_indiv_deaths_non_icu_F[age_class, 1:t_truncation, obs_i]| t_truncation, distrib_deaths_non_icu_F_exp[match_cases_delay_age[age_class], match_changes_windows[obs_i],], p_hosp_F_scaled[age_class, match_changes_windows[obs_i],2]);
      target += likelyhood_lpmf(delays_indiv_deaths_non_icu_M[age_class, 1:t_truncation, obs_i]| t_truncation, distrib_deaths_non_icu_M_exp[match_cases_delay_age[age_class], match_changes_windows[obs_i],], p_hosp_M_scaled[age_class, match_changes_windows[obs_i],2]);
      // to discharge
      target += likelyhood_lpmf(delays_indiv_discharges_non_icu_F[age_class, 1:t_truncation, obs_i]| t_truncation, distrib_discharges_non_icu_F_exp[match_cases_delay_age[age_class], match_changes_windows[obs_i],], p_hosp_F_scaled[age_class, match_changes_windows[obs_i],3]);
      target += likelyhood_lpmf(delays_indiv_discharges_non_icu_M[age_class, 1:t_truncation, obs_i]| t_truncation, distrib_discharges_non_icu_M_exp[match_cases_delay_age[age_class], match_changes_windows[obs_i],], p_hosp_M_scaled[age_class, match_changes_windows[obs_i],3]);
      // Unknown outcomes 
      target += likelyhood_unknown_hosp_lccdf(delays_indiv_cases_unobs_hosp_F[age_class, 1:t_truncation, obs_i] |  t_truncation , distrib_icu_F_exp_cdf[match_cases_delay_age[age_class], match_changes_windows[obs_i]], distrib_deaths_non_icu_F_exp_cdf[match_cases_delay_age[age_class], match_changes_windows[obs_i]], distrib_discharges_non_icu_F_exp_cdf[match_cases_delay_age[age_class], match_changes_windows[obs_i]], p_hosp_F_scaled[age_class, match_changes_windows[obs_i],1], p_hosp_F_scaled[age_class, match_changes_windows[obs_i],2], p_hosp_F_scaled[age_class, match_changes_windows[obs_i],3]);
      target += likelyhood_unknown_hosp_lccdf(delays_indiv_cases_unobs_hosp_M[age_class, 1:t_truncation, obs_i] |  t_truncation , distrib_icu_M_exp_cdf[match_cases_delay_age[age_class], match_changes_windows[obs_i]], distrib_deaths_non_icu_M_exp_cdf[match_cases_delay_age[age_class], match_changes_windows[obs_i]], distrib_discharges_non_icu_M_exp_cdf[match_cases_delay_age[age_class], match_changes_windows[obs_i]], p_hosp_M_scaled[age_class, match_changes_windows[obs_i],1], p_hosp_M_scaled[age_class, match_changes_windows[obs_i],2], p_hosp_M_scaled[age_class, match_changes_windows[obs_i],3]);


      // ICU admissions
      // to death
      target += likelyhood_lpmf(delays_indiv_deaths_icu_F[age_class, 1:t_truncation, obs_i]| t_truncation, distrib_deaths_icu_F_exp[match_cases_delay_age[age_class], match_changes_windows[obs_i],], p_icu_F_scaled[age_class, match_changes_windows[obs_i],1]);
      target += likelyhood_lpmf(delays_indiv_deaths_icu_M[age_class, 1:t_truncation, obs_i]| t_truncation, distrib_deaths_icu_M_exp[match_cases_delay_age[age_class], match_changes_windows[obs_i],], p_icu_M_scaled[age_class, match_changes_windows[obs_i],1]);
      // to discharge
      target += likelyhood_lpmf(delays_indiv_discharges_icu_F[age_class, 1:t_truncation, obs_i]| t_truncation, distrib_discharges_icu_F_exp[match_cases_delay_age[age_class], match_changes_windows[obs_i],], p_icu_F_scaled[age_class, match_changes_windows[obs_i],2]);
      target += likelyhood_lpmf(delays_indiv_discharges_icu_M[age_class, 1:t_truncation, obs_i]| t_truncation, distrib_discharges_icu_M_exp[match_cases_delay_age[age_class], match_changes_windows[obs_i],], p_icu_M_scaled[age_class, match_changes_windows[obs_i],2]);
      // Unknown outcomes
      target += likelyhood_unknown_icu_lccdf(delays_indiv_cases_unobs_icu_F[age_class, 1:t_truncation, obs_i] |  t_truncation, distrib_deaths_icu_F_exp_cdf[match_cases_delay_age[age_class], match_changes_windows[obs_i],], distrib_discharges_icu_F_exp_cdf[match_cases_delay_age[age_class], match_changes_windows[obs_i],], p_icu_F_scaled[age_class, match_changes_windows[obs_i],1], p_icu_F_scaled[age_class, match_changes_windows[obs_i],2]);
      target += likelyhood_unknown_icu_lccdf(delays_indiv_cases_unobs_icu_M[age_class, 1:t_truncation, obs_i] |  t_truncation, distrib_deaths_icu_M_exp_cdf[match_cases_delay_age[age_class], match_changes_windows[obs_i],], distrib_discharges_icu_M_exp_cdf[match_cases_delay_age[age_class], match_changes_windows[obs_i],], p_icu_M_scaled[age_class, match_changes_windows[obs_i],1], p_icu_M_scaled[age_class, match_changes_windows[obs_i],2]);
    }
  }
}





