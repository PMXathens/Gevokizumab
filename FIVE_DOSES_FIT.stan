functions{
  //This whole function returns the solution to the ODEs
  
  
  //The matrix_exp() function is used to compute a matrix exponential solution for our linear ODE system:
  real[] effCptModel1(real t0, real t, real[] init, real amt, int cmt, int evid,
			  real CLp, real s1, real s2, real Vp, real Vt, real Vl, real VL) {
    
    real L; real sL;
    matrix[4, 4] K;
    real x[4];

    x = rep_array(0, 4);
    L = 2.9/24; sL = 0.2; 
    K = rep_matrix(0, 4, 4);
    
    K[1, 1] = (-(CLp+(0.33*L)*(1-s1)+ (0.67*L)*(1-s2)))/Vp;
    K[1, 4] =  L/VL;
    K[2, 1] = (0.33*L)*(1-s1)/Vp;
    K[2, 2] = (-(0.33*L)*(1-sL))/Vt;
    K[3, 1] = (0.67*L)*(1-s2)/Vp;
    K[3, 3] = (-(0.67*L)*(1-sL))/Vl;
    K[4, 2] = (0.33*L)*(1-sL)/Vt;
    K[4, 3] = (0.67*L)*(1-sL)/Vl;
    K[4, 4] = (-L)/VL;

    x = to_array_1d(matrix_exp((t - t0) * K) * to_vector(init));
    if (evid == 1) x[cmt] = x[cmt] + amt;
    return x;
  }
  
  //In this next piece of code, a function is used to handle the event schedule and calls the previous function effCptModel1:
  matrix effCptModel(real[] time, real[] amt, int[] cmt, int[] evid, 
		     real CLp, real s1, real s2, real Vp, real Vt, real Vl, real VL) {
		       
    real init[4];
    matrix[size(time), 4] result;
    int nt;
    real t0;

    nt = size(time);

    init = rep_array(0, 4);
    t0 = time[1];
    for (i in 1:nt){
      init = effCptModel1(time[max(1, i - 1)], time[i], init, amt[i], cmt[i], evid[i],
			       CLp, s1, s2, Vp, Vt, Vl, VL);
      for (j in 1:4) result[i, j] = init[j];
      t0 = time[i];
    }
    return result;
  }
  
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
data {
  int<lower = 1> nt;          //Number of time points
  real<lower = 0> amt_a[nt];  //Dose's amount (kg) at each time point fot the first dose of 0.01mg/kg
  real<lower = 0> amt_b[nt];  //Dose's amount (kg) at each time point for the second dose of 0.03mg/kg
  real<lower = 0> amt_c[nt];  //Dose's amount (kg) at each time point for the third dose of 0.1mg/kg
  real<lower = 0> amt_d[nt];  //Dose's amount (kg) at each time point for the fourth dose of 0.3mg/kg
  real<lower = 0> amt_e[nt];  //Dose's amount (kg) at each time point for the fifth dose of 1.0mg/kg
  int<lower = 1> cmt[nt];     //Compartment number
  int<lower = 0> evid[nt];    //Event ID
  real<lower = 0> time[nt];   //Time points
  real lhs_1[60]; real lhs_2[60]; //LHS Number
  real means_a[12];              //Observations - Mean plasma concentrtaions for the first dose
  real sds_a[9];                 //Observations - SDs for the first dose
  real means_b[14];              //Observations - Mean plasma concentrtaions for the second dose
  real sds_b[4];                 //Observations - SDs for the second dose
  real means_c[14];              //Observations - Mean plasma concentrtaions for the third dose
  real sds_c[13];                //Observations - SDs for the third dose
  real means_d[14];              //Observations - Mean plasma concentrtaions for the fourth dose
  real sds_d[12];                //Observations - SDs for the fourth dose
  real means_e[14];              //Observations - Mean plasma concentrtaions for the fifth dose
  real sds_e[8];                 //Observations - SDs for the fifth dose
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

transformed data {

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

parameters{
    real<lower=0>  sigma_1;              //Residual error for the means
    real<lower=0>  sigma_2;              //Residual error for the SDs
    real<lower=0>  theta1_mean;          //Group mean plasma Clearance
    real<lower=0, upper=1> theta2_mean;  //Group reflection coefficient ó1
    real<lower=0, upper=1> theta3_mean;  //Group reflection coefficient ó2
    
    real<lower=0> theta1_sd;             //inter-group variability of plasma Clearance
    real<lower=0> theta2_sd;             //inter-group variability of reflection coefficient ó1
    real<lower=0> theta3_sd;             //inter-group variability of reflection coefficient ó2
    
    vector<lower=0>[5] Omega_CLp;        //IIV of plasma Clearance for the first,second,third,fourth and fifth Dose 
    vector<lower=0>[5] Omega_V;          //IIV of volume of Human Body for the first,second,third,fourth and fifth Dose 
    
    vector[5] log_theta1;                //logarith of plasma Clearance for the first,second,third,fourth and fifth Dose
    vector<upper=0>[5] log_theta2;       //logarith of reflection coefficient ó1 for the first,second,third,fourth and fifth Dose
    vector<upper=0>[5] log_theta3;       //logarith of reflection coefficient ó2 for the first,second,third,fourth and fifth Dose
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

transformed parameters{
  
  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
model{
  
  //First Dose
  real y_hat1_a[12]; real y_hat2_a[9];    //prediction values for means and sds respectively
  real CLp_a[60]; real s1_a[60]; real s2_a[60]; real Vp_a[60]; real Vt_a[60]; real Vl_a[60]; real VL_a[60]; //parameters for 60 virtual patients
  matrix[nt, 4] x_a;     //matrix to store solution of the ODEs for 60 virtual patients
  vector[nt] sol_a[60];  //vector to store plasma concentrations at each time point for 60 virtual patients
  
  //Second Dose
  real y_hat1_b[14]; real y_hat2_b[4];    //prediction values for means and sds respectively
  real CLp_b[60]; real s1_b[60]; real s2_b[60]; real Vp_b[60]; real Vt_b[60]; real Vl_b[60]; real VL_b[60];  //parameters for 60 virtual patients
  vector[nt] sol_b[60]; //matrix to store solution of the ODEs for 60 virtual patients
  matrix[nt, 4] x_b;    //vector to store plasma concentrations at each time point for 60 virtual patients
  
  //Third Dose
  real y_hat1_c[14]; real y_hat2_c[13];    //prediction values for means and sds respectively
  real CLp_c[60]; real s1_c[60]; real s2_c[60]; real Vp_c[60]; real Vt_c[60]; real Vl_c[60]; real VL_c[60];  //parameters for 60 virtual patients
  vector[nt] sol_c[60];  //matrix to store solution of the ODEs for 60 virtual patients
  matrix[nt, 4] x_c;     //vector to store plasma concentrations at each time point for 60 virtual patients
  
  //Fourth Dose
  real y_hat1_d[14]; real y_hat2_d[12];    //prediction values for means and sds respectively
  real CLp_d[60]; real s1_d[60]; real s2_d[60]; real Vp_d[60]; real Vt_d[60]; real Vl_d[60]; real VL_d[60];  //parameters for 60 virtual patients
  vector[nt] sol_d[60];  //matrix to store solution of the ODEs for 60 virtual patients
  matrix[nt, 4] x_d;     //vector to store plasma concentrations at each time point for 60 virtual patients
  
  //Fifth Dose
  real y_hat1_e[14]; real y_hat2_e[8];    // prediction values for means and sds respectively
  real CLp_e[60]; real s1_e[60]; real s2_e[60]; real Vp_e[60]; real Vt_e[60]; real Vl_e[60]; real VL_e[60];  //parameters for 60 virtual patients
  vector[nt] sol_e[60];  //matrix to store solution of the ODEs for 60 virtual patients
  matrix[nt, 4] x_e;     //vector to store plasma concentrations at each time point for 60 virtual patients
  
 //priors
 sigma_1 ~ normal(0,1);
 sigma_2 ~ normal(0,1);
 theta1_mean ~ normal(0,1);
 theta2_mean ~ beta(1,1);
 theta3_mean ~ beta(1,1);
 theta1_sd ~ normal(0,1);
 theta2_sd ~ normal(0,1);
 theta3_sd ~ normal(0,1);
 Omega_CLp ~ normal(0,1);
 Omega_V ~ normal(0,1);
 
 //The logarithms of mean plasma Clearance, reflection coefficients ó1 and ó2 of the five doses follow a Student's-t distribution
 //with 5 degrees of freedom, location equal to log(theta_mean) and scale sigma equal to sqrt((((theta_sd)^2)*3)/5).
 for(i in 1:5){
   log_theta1[i] ~ student_t(5,log(theta1_mean),sqrt((((theta1_sd)^2)*3)/5));
   log_theta2[i] ~ student_t(5,log(theta2_mean),sqrt((((theta2_sd)^2)*3)/5));
   log_theta3[i] ~ student_t(5,log(theta3_mean),sqrt((((theta3_sd)^2)*3)/5));
  }
  
  for (i in 1:60){
    //First Dosage-group
    CLp_a[i] = exp(log_theta1[1])*exp(Omega_CLp[1]*lhs_1[i]);
    s1_a[i] = exp(log_theta2[1]);
    s2_a[i] = exp(log_theta3[1]);
    Vp_a[i] = 2.6*exp(Omega_V[1]*lhs_2[i]);
    Vt_a[i] = 8.112*exp(Omega_V[1]*lhs_2[i]);
    Vl_a[i] = 4.368*exp(Omega_V[1]*lhs_2[i]);
    VL_a[i] = 5.2*exp(Omega_V[1]*lhs_2[i]);
   
    x_a = effCptModel(time, amt_a, cmt, evid, CLp_a[i], s1_a[i], s2_a[i], Vp_a[i], Vt_a[i], Vl_a[i], VL_a[i]);
    sol_a[i] = x_a[:,1]/Vp_a[i];
    
    //Second Dosage-group
    CLp_b[i] = exp(log_theta1[2])*exp(Omega_CLp[2]*lhs_1[i]);
    s1_b[i] = exp(log_theta2[2]);
    s2_b[i] = exp(log_theta3[2]);
    Vp_b[i] = 2.6*exp(Omega_V[2]*lhs_2[i]);
    Vt_b[i] = 8.112*exp(Omega_V[2]*lhs_2[i]);
    Vl_b[i] = 4.368*exp(Omega_V[2]*lhs_2[i]);
    VL_b[i] = 5.2*exp(Omega_V[2]*lhs_2[i]);
   
    x_b = effCptModel(time, amt_b, cmt, evid, CLp_b[i], s1_b[i], s2_b[i], Vp_b[i], Vt_b[i], Vl_b[i], VL_b[i]);
    sol_b[i] = x_b[:,1]/Vp_b[i];
    
    //Third Dosage-group
    CLp_c[i] = exp(log_theta1[3])*exp(Omega_CLp[3]*lhs_1[i]);
    s1_c[i] = exp(log_theta2[3]);
    s2_c[i] = exp(log_theta3[3]);
    Vp_c[i] = 2.6*exp(Omega_V[3]*lhs_2[i]);
    Vt_c[i] = 8.112*exp(Omega_V[3]*lhs_2[i]);
    Vl_c[i] = 4.368*exp(Omega_V[3]*lhs_2[i]);
    VL_c[i] = 5.2*exp(Omega_V[3]*lhs_2[i]);
   
    x_c = effCptModel(time, amt_c, cmt, evid, CLp_c[i], s1_c[i], s2_c[i], Vp_c[i], Vt_c[i], Vl_c[i], VL_c[i]);
    sol_c[i] = x_c[:,1]/Vp_c[i];
    
    //Fourth Dosage-group
    CLp_d[i] = exp(log_theta1[4])*exp(Omega_CLp[4]*lhs_1[i]);
    s1_d[i] = exp(log_theta2[4]);
    s2_d[i] = exp(log_theta3[4]);
    Vp_d[i] = 2.6*exp(Omega_V[4]*lhs_2[i]);
    Vt_d[i] = 8.112*exp(Omega_V[4]*lhs_2[i]);
    Vl_d[i] = 4.368*exp(Omega_V[4]*lhs_2[i]);
    VL_d[i] = 5.2*exp(Omega_V[4]*lhs_2[i]);
   
    x_d = effCptModel(time, amt_d, cmt, evid, CLp_d[i], s1_d[i], s2_d[i], Vp_d[i], Vt_d[i], Vl_d[i], VL_d[i]);
    sol_d[i] = x_d[:,1]/Vp_d[i];
    
    //Fifth Dosage-group
    CLp_e[i] = exp(log_theta1[5])*exp(Omega_CLp[5]*lhs_1[i]);
    s1_e[i] = exp(log_theta2[5]);
    s2_e[i] = exp(log_theta3[5]);
    Vp_e[i] = 2.6*exp(Omega_V[5]*lhs_2[i]);
    Vt_e[i] = 8.112*exp(Omega_V[5]*lhs_2[i]);
    Vl_e[i] = 4.368*exp(Omega_V[5]*lhs_2[i]);
    VL_e[i] = 5.2*exp(Omega_V[5]*lhs_2[i]);
   
    x_e = effCptModel(time, amt_e, cmt, evid, CLp_e[i], s1_e[i], s2_e[i], Vp_e[i], Vt_e[i], Vl_e[i], VL_e[i]);
    sol_e[i] = x_e[:,1]/Vp_e[i];
  }
  
  //Prediction values for means and SDs for the first dose concerning the specific time points that were available
  for (j in 1:12){
    y_hat1_a[j] = mean(sol_a[:,j+1]);
  }
  for (j in 1:9){
    y_hat2_a[j] = sd(sol_a[:,j+2]);
  }
  
  //Prediction values for means and SDs for the second dose concerning the specific time points that were available
  for (j in 1:14){
    y_hat1_b[j] = mean(sol_b[:,j+1]);
  }
  for (j in 1:4){
    y_hat2_b[j] = sd(sol_b[:,j+10]);
  }
  
  //Prediction values for means and SDs for the third dose concerning the specific time points that were available
  for (j in 1:14){
    y_hat1_c[j] = mean(sol_c[:,j+1]);
  }
  for (j in 1:13){
    y_hat2_c[j] = sd(sol_c[:,j+2]);
  }
  
  //Prediction values for means and SDs for the fourth dose concerning the specific time points that were available
  for (j in 1:14){
    y_hat1_d[j] = mean(sol_d[:,j+1]);
  }
  for (j in 1:12){
    y_hat2_d[j] = sd(sol_d[:,j+3]);
  }
  
  //Prediction values for means and SDs for the fifth dose concerning the specific time points that were available
  for (j in 1:14){
    y_hat1_e[j] = mean(sol_e[:,j+1]);
  }
  for (j in 1:8){
    y_hat2_e[j] = sd(sol_e[:,j+7]);
  }
 
  //likelihood ~  
  log(means_a) ~ normal(log(y_hat1_a),sigma_1);
  log(sds_a) ~ normal(log(y_hat2_a),sigma_2);
  
  log(means_b) ~ normal(log(y_hat1_b),sigma_1);
  log(sds_b) ~ normal(log(y_hat2_b),sigma_2);
 
  log(means_c) ~ normal(log(y_hat1_c),sigma_1);
  log(sds_c) ~ normal(log(y_hat2_c),sigma_2);
  
  log(means_d) ~ normal(log(y_hat1_d),sigma_1);
  log(sds_d) ~ normal(log(y_hat2_d),sigma_2);
  
  log(means_e) ~ normal(log(y_hat1_e),sigma_1);
  log(sds_e) ~ normal(log(y_hat2_e),sigma_2);
}

