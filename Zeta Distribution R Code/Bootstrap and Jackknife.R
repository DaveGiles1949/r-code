# Monte Carlo Simulation Experiment to Evaluate the Small-Sample Properties of the
# Non-Parametric Bootstrap and Leave-one-out Jacknife Bias Correction of the MLE
# of the Zeta Distribution's Scale Parameter
#
# Reference for R package: Tibshirani, R., F. Leisch, and S. Kostyshack, 2025. 
# Package ‘bootstrap’: Functions for the book "An Introduction to the Bootstrap". 
# https://gitlab.com/scottkosty/bootstrap
# ---------------------------------------

# The following code was written by David Giles, October 2025
# -----------------------------------------------------------

# Run the code separately for the Bootstrap and the Jackknife (with same seed)
# to make the computation time manageable
#
library(bootstrap)
library(VGAM)
library (tolerance)
library(UnivRNG)

set.seed(123)
repMC<- 10000	# Number of Monte Carlo replications
nboot<- 1000      # Number of Bootstrap samples
s<- 1.25    	# True value of the Zeta-distribution parameter 
N<- 200    		# Sample size

res_tilde<-  c()
res_jack<-   c()
res_tilde2<- c()
res_jack2<-  c()
res_boot<-   c()
res_boot2<-  c()

# Define the statistic of interest (namely the MLE of 's')for the Jackknife & Bootstrap

my_statistic <- function(data) {

out.zeta<- zm.ll(data, N = Inf, dist = "Zeta")   # Using the 'tolerance' package
smle<- stats4::coef(out.zeta)[[1]]               # Using the 'VGAM' package
smle
}

# Start the MC loop:
start_time <- Sys.time()

for (ii in 1:repMC) {

zeta.data<- draw.zeta(N, s)$y                    # Using the 'UnivRNG' package. It is very fast. 
                                                 # The empirical moments of the simulated data match the 
                                                 # theoretical moments (where they exist) very closely.
out<- zm.ll(zeta.data, N = Inf, dist = "Zeta")   
sest<- stats4::coef(out)[[1]]   # MLE of 's'
res_tilde[ii]<- sest
res_tilde2[ii]<- sest^2

jackknife_results <- jackknife(zeta.data, my_statistic)  # Leave-one-out Jackknife 
#boot_results<- bootstrap(zeta.data, nboot, my_statistic)  # Non-parametric Bootstrap (with 'nboot' boot-samples)

# "Leave-one-out" Jackknife bias-corrected estimator:
res_jack[ii]<- sest - jackknife_results$jack.bias[[1]]
res_jack2[ii]<- res_jack[ii]^2

# (Non-parametric) Bootstrap bias-corrected estimator:
#res_boot[ii]<- 2*sest - mean(boot_results$thetastar)
#res_boot2[ii]<- res_boot[ii]^2

}     # End of the MC loop

end_time <- Sys.time()
time_taken <- end_time - start_time

perc_bias_tilde<- 100*(mean(res_tilde) - s)/s
perc_bias_jack<- 100*(mean(res_jack) - s)/s
#perc_bias_boot<- 100*(mean(res_boot) - s)/s
perc_mse_tilde<- 100*(mean(res_tilde2)-2*s*mean(res_tilde)+s^2)/(s^2)
perc_mse_jack<- 100*(mean(res_jack2)-2*s*mean(res_jack)+s^2)/(s^2)
#perc_mse_boot<- 100*(mean(res_boot2)-2*s*mean(res_boot)+s^2)/(s^2)

# Print %Biases and %MSEs:

c(s,N,repMC, nboot)

perc_bias_tilde
perc_bias_jack
#perc_bias_boot
perc_mse_tilde
perc_mse_jack
#perc_mse_boot

time_taken