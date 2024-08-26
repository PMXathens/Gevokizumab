library(rstan)
#load "ONE_DOSE.RData"

# Summary statistics of one dosage-group's fit
fit_summary <- summary(fit)
print(fit_summary$summary)

#Inference for one dosage-group's fit: Rhat, Bulk_ESS and Tail_ESS
monitor(extract(fit, permuted = FALSE, inc_warmup = TRUE))

#Use of shinystan for inference
# First installation of shinystan is required by typing install.packages("shinystan")
library(shinystan)
launch_shinystan(fit)

#Check diagnostic tools
check_hmc_diagnostics(fit)
