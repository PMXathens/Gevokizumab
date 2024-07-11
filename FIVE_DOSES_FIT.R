library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

############################################
#the order of compartments is as follows:
#[1]= plasma
#[2]= tight tissue
#[3]= leaky tissue
#[4]= lymph
##################################################
#load data

#time vector in hours
time = c(0,4,8,24,48,72,96,168,216,264,336,504,672,1008,1344)

#mean values of plasma concentrations at different time points for each dose that was used

#Dose = 0.01 mg/kg for a 70kg person
y1_a_arx = c(0.230857533572143,0.1900112060476800,0.1616604472341880,0.1338851566029320,0.1166386460531590,0.1063097784111840,
             0.0844289897155243,0.0800018038916830,0.0758067655136998,0.0749076942726890,0.057574996220708,0.050321593592600)
#Dose = 0.03 mg/kg for a 70kg person
y1_b_arx = c(0.800823810803821,0.649928373475374,0.583555103226454,0.483293023857174,0.411182940243581,0.369191303121177,0.297635144163131,
             0.274534223383870,0.253226278169878,0.232980653051069,0.185696209704508,0.141854780146239,0.081504109619910,0.065873911240799)
#Dose = 0.1 mg/kg for a 70kg person
y1_c_arx = c(2.4040991835100100000,2.1077035344735200000,1.9574505122221200000,1.6054547740153000000,1.3995725850794500000,1.3167561658979200000,
             1.0799723719988500000,0.9855474670725830000,0.9131980245186050000,0.8083228036422610000,0.6629673117915560000,0.5274172953434570000,
             0.3419898168128350000,0.2395829193571490000)
#Dose = 0.3 mg/kg for a 70kg person
y1_d_arx = c(8.122269947080080,6.695861340634240,5.914660544591750,5.310631887314350,5.341474791658700,4.398198780581130,3.359818286283780,
             3.270543147432730,2.858514179684470,2.636650898730360,2.125618788191990,1.623776739188720,1.175372265130630,0.759231786940751)
#Dose = 1.0 mg/kg for a 70kg person
y1_e_arx = c(26.656704127384400,23.669743321043600,20.971883035581500,17.368650727593700,15.594895040582800,14.002282338513400,11.913061720030600,
             10.600301490059700,9.777560191660020,8.623280529014950,6.750068224441020,5.614489978846360,3.892494659848540,2.858514179684470)

# sds values of plasma concnentrations at different time points for each dose that was used
#Dose = 0.01 mg/kg for a 70kg person
y2_a_arx = c(0.025432262955508,0.021637623849055,0.017920008926793,0.012516320448329,0.012820838789121,
             0.009059949906539,0.008584875149325,0.008134711790879,0.009033783031890)
#Dose = 0.03 mg/kg for a 70kg person
y2_b_arx = c(0.034258806502265,0.029747259298680,0.022724162945933,0.011984830002154)
#Dose = 0.1 mg/kg for a 70kg person
y2_c_arx = c(0.2525626310110500000,0.0959838335051799000,0.0608735757883900000,0.0613945046364999000,0.0459441401037699000,0.0605929162498200000,
             0.0458831890644871000,0.0408974518313870000,0.0466903363937691000,0.0345217719836141000,0.0327262291842989000,0.0426320873407860000,
             0.0372239545846345000)
#Dose = 0.3 mg/kg for a 70kg person
y2_d_arx = c(0.310145856347200,0.278472528054360,0.573185752933044,0.230627458144800,0.176178110597270,0.171496808233300,
             0.149891328750810,0.138257528848900,0.111460641631290,0.043146051810590,0.061632757272770,0.091562493021995)
#Dose = 1.0 mg/kg for a 70kg person
y2_e_arx = c(1.002434930118200,0.996213277554320,0.918892287559185,0.980807683490410,0.992568602370250,0.972901145233585,0.875794378877120,0.789376098569680)


#Latin Hypercube Sampling (LHS) for N=60 virtual patients

N = 60
library("lhs")

X_a <- randomLHS(N, 2)
Y_a = matrix(0, nrow=N, ncol=2)
Y_a[,1] <- qnorm(X_a[,1], mean=0, sd=1)
Y_a[,2] <- qnorm(X_a[,2], mean=0, sd=1)
e_a<-Y_a[,1]
f_a<-Y_a[,2]
#LHS numbers for plasma Clearance - lhs_1
lhs_1 = e_a 
#LHS numbers for Volume of Human Body - lhs_2
lhs_2 = f_a

#Fixed values of the mPBPK model
L = 2.9/24
sL = 0.2
#########################################################################

#Generetaing the list for the Stan fit
Gevokizumab_dat<-list(
  lhs_1 = lhs_1,
  lhs_2 = lhs_2,
  nt = 15, # number of time points             
  evid = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),         #event ID at each time point concerning the plasma compartment. 1=Dose and 0=observation
  amt_a =  c(0.7,0,0,0,0,0,0,0,0,0,0,0,0,0,0),     #dose's amount (kg) at each time point adapted to a 70kg person concerning the first dose of 0.01mg/kg
  amt_b =  c(2.1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),     #dose's amount (kg) at each time point adapted to a 70kg person concerning the second dose of 0.03mg/kg
  amt_c =  c(7,0,0,0,0,0,0,0,0,0,0,0,0,0,0),       #dose's amount (kg) at each time point adapted to a 70kg person concerning the third dose of 0.1mg/kg
  amt_d =  c(21,0,0,0,0,0,0,0,0,0,0,0,0,0,0),      #dose's amount (kg) at each time point adapted to a 70kg person concerning the fourth dose of 0.3mg/kg
  amt_e =  c(70,0,0,0,0,0,0,0,0,0,0,0,0,0,0),      #dose's amount (kg) at each time point adapted to a 70kg person concerning the fifth dose of 1.0mg/kg
  time = c(0,4,8,24,48,72,96,168,216,264,336,504,672,1008,1344),  #time points
  means_a = y1_a_arx,               #mean plasma concentrations for the first dose
  sds_a = y2_a_arx,                 #SDs of plasma concentrations for the first dose
  means_b = y1_b_arx,               #mean plasma concentrations for the second dose
  sds_b = y2_b_arx,                 #SDs of plasma concentrations for the second dose
  means_c = y1_c_arx,               #mean plasma concentrations for the third dose
  sds_c = y2_c_arx,                 #SDs of plasma concentrations for the third dose
  means_d = y1_d_arx,               #mean plasma concentrations for the fourth dose
  sds_d = y2_d_arx,                 #SDs of plasma concentrations for the fourth dose
  means_e = y1_e_arx,               #mean plasma concentrations for the fifth dose
  sds_e = y2_e_arx,                 #SDs of plasma concentrations for the fifth dose
  cmt = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)    #compartment number. 1=plasma
)

tic = proc.time()

fit5 <- stan(file = 'FIVE_DOSES_FIT.stan',iter = 4000, warmup = 2000, data = Gevokizumab_dat, chains = 4)

#The stan-fit fit5 can be loaded from the Supplementary Material via the RDATA File FIVE_DOSES

#Summary Statistics
fit_summary <- summary(fit5)
print(fit_summary$summary)

#Check diagnostic tools
check_hmc_diagnostics(fit5)

monitor(extract(fit5, permuted = FALSE, inc_warmup = TRUE))

#Traceplots
parametersToPlot = c("sigma_1","sigma_2 ","theta1_mean","theta2_mean","theta3_mean","theta1_sd","theta2_sd","theta3_sd",
                     "Omega_CLp[1]","Omega_CLp[2]","Omega_CLp[3]","Omega_CLp[4]","Omega_CLp[5]","Omega_V[1]","Omega_V[2]",
                     "Omega_V[3]","Omega_V[4]","Omega_V[5]","log_theta1[1]","log_theta1[2]","log_theta1[3]","log_theta1[4]",
                     "log_theta1[5]","log_theta2[1]","log_theta2[2]","log_theta2[3]","log_theta2[4]",
                     "log_theta2[5]","log_theta3[1]","log_theta3[2]","log_theta3[3]","log_theta3[4]",
                     "log_theta3[5]")
stan_trace(fit5, parametersToPlot)

#Use shinystan for inference
library(shinystan)
launch_shinystan(fit5)
