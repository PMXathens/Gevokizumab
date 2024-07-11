D = 100    # Dose (mg)
t = c(1,2,4,6,8,12,16,24)   # Time points when plasma concentration was measured (hours)

N1 = 24   # Number of patients in a theoretical clinical study
CLp = 5   # Plasma Clearance (L/h)
Vd = 20   # Volume of distribution
Omega_CLp = 0.2 # Standard deviation of the lognormal distribution of plasma clearance 
Omega_V = 0.2   # Standard deviation of the lognormal distribution of volume of distribution
rv = 0.05    # Residual variability
   
e = rnorm(N1,log(CLp),Omega_CLp)
CLear_p = exp(e)   # Values of plasma clearance for N1 = 24 patients
  
f = rnorm(N1,log(Vd),Omega_V)
Vol_d = exp(f)     # Values of volume of distribution for N1 = 24 patients
  
y1<-matrix(0,length(t),N1)     # Matrix to store concentrations without residual variability 
y2<-rnorm(N1*(length(t)),0,rv) # Residual variability for each patient's concentration
for (j in 1:N1){
  y1[,j]<- (D/Vol_d[j])*exp(-(CLear_p[j]/Vol_d[j])*t)    # IV bolus model (one compartment, single dose)
}
y<-matrix(0,length(t),N1)      # Matrix to store concentrations with residual variability
y<-y1*exp(y2)

y1_data_arx = rowMeans(y)      # Mean plasma concentrations
#print(y1_data_arx)
y2_data_arx = apply(y,1,sd)    # Standard deviations of plasma concentrations
#print(y2_data_arx)
  
library(rstan)
setwd("C:/Users/vangelisk/Desktop/")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
  
rnorm_1 = rnorm(1000,0,1)  # 1000 virtual values that follow the normal distribution (0,1) - Plasma Clearance
rnorm_2 = rnorm(1000,0,1)  # 1000 virtual values that follow the normal distribution (0,1) - Volume of Distribution

Gevokizumab_dat<-list(
    N = 1000,      # Number of virtual patients      
    nt = 8,        # Number of time points
    time = t,      # Time points
    Dose = D,      # Dose (mg)
    means = y1_data_arx, # Mean plasma concentrations
    sds = y2_data_arx,   # Standard deviations of plasma concentrations
    rnorm_1 = rnorm_1,
    rnorm_2 = rnorm_2
)
  
tic = proc.time()
  
fit1 <- stan(file = 'iv_bolus.stan',iter = 2000, data = Gevokizumab_dat, chains = 4, control=list(adapt_delta=0.99))
  
monitor(extract(fit1, permuted = FALSE, inc_warmup = TRUE))
  
fit_summary <- summary(fit1)
print(fit_summary$summary)
  
#Check diagnostic tools
check_hmc_diagnostics(fit1)

