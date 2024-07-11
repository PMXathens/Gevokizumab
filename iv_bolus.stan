data {
  int<lower=1> N;            // Number of virtual patients  
  int<lower = 1> nt;         // Number of time points
  real<lower = 0> time[nt];  // Time points
  real<lower=0> Dose;        // Dose (mg)
  real means[nt];            // Means observation values
  real sds[nt];              // Sds' observation values
  real rnorm_1[N];           // N values that follow the normal distribution (0,1) for plasma clearance
  real rnorm_2[N];           // N values that follow the normal distribution (0,1) for volume of distribution
}

////////////////////////////////////////////////////////////////////////////

transformed data {

}

//////////////////////////////////////////////////////////////////////////

parameters{
    real<lower=0>  sigma_1;     // Residual error for the means
    real<lower=0>  sigma_2;     // Residual error for the SDs
    real<lower=0>  CLp_;        // Mean plasma clearance
    real<lower=0>  Vd_;         // Mean volume of distribution
    real<lower=0>  Omega_CLp_;  // IIV of plasma clearance
    real<lower=0>  Omega_Vd_;   // IIV of volume of distribution
}

////////////////////////////////////////////////////////////////////////////

transformed parameters{
  
  
}
////////////////////////////////////////////////////////////////////////////
        
model{
  
  real y_hat1[nt]; real y_hat2[nt];      // Prediction values for means and sds respectively
  real CLearance_p[N]; real Volume_d[N]; // Values of plasma clearance and volume of distribution for the virtual patients
  matrix [N,nt] sol;                     // matrix to store plasma concentrations versus time for each one of the virtual patients
  
 //priors
 sigma_1 ~ normal(0,1);
 sigma_2 ~ normal(0,1);
 CLp_ ~ normal(0,100);
 Vd_ ~ normal(0,100);
 Omega_CLp_ ~ normal(0,1);
 Omega_Vd_ ~ normal(0,1);
 
 for (i in 1:N){
    CLearance_p[i] = CLp_*exp(Omega_CLp_*rnorm_1[i]);
    Volume_d[i] = Vd_*exp(Omega_Vd_*rnorm_2[i]);
    
    for(j in 1:nt){
      sol[i,j] = (Dose/Volume_d[i])*exp(-(CLearance_p[i]/Volume_d[i])*time[j]);  // Plasma concentration according to the IV bolus model
    }
  }
  
  //predicted mean values
  for (j in 1:nt){
    y_hat1[j] = mean(sol[:,j]);
  }
 
 //predicted SD values
 for (j in 1:nt){
    y_hat2[j] = sd(sol[:,j]);
  }
 
 //likelihood~  
  log(means) ~ normal(log(y_hat1),sigma_1);
  log(sds) ~ normal(log(y_hat2),sigma_2);
  
}

generated quantities{
    
}

