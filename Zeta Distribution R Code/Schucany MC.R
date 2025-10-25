# Monte Carlo Simulation Experiment to Investigate the Effectiveness of the 
# Schucany et al. "peudo-estimator" of the Scale Parameter of the Zeta Distribution
#
# We call it a pseudo estimator here because in order to construct the "R" ratio we need to know the 
# true biases, and here we use the estimated biases - and to get these from the moan MC experiment
# we actually have to know the value of 's', which is the parameter weare trying to estimate!
#
# The 2 biased estimators that are combined are the MLE and the Cox-Snell (first-order)
# bias-adjusted etimator. The biases used in the combination are those estimated
# in the main Monte Carlo experiment - see Table 1 of the paper, for example.
# Note - the ratio of the biases equals the raio of the percenage biases, for any given case.
#
# Code written by David Giles, October 2025
# -------------------------

library(VGAM)
library (tolerance)
library(UnivRNG)

set.seed(123)
repMC<- 100000

s<- 1.25    # True value of the Zeta-distribution parameter. Clauset et al. find values between 1.7 and 3.7 for a range of data-sets.
n<- 200     # Sample size.  
perc_bias_tilde<- 9.76       # From Table 1
perc_bias_hat<- 9.61         # From Table 1

R<-  perc_bias_tilde/perc_bias_hat       # Schucany et al.

sum_sdbar<- sum_sdbar2<- 0   # for the Schucany et al. pseudo-estimator

# Gamma constants -  from Choudhury (1995), Table 5 (Note the shift in decimal places, as indicated in that table.)
gamma1<-  -0.072815845483676724861
gamma2<-  -0.009690036031902870231084845
gamma3<-   0.0020538344203033458662
gamma4<-   0.0023253700654673000575
gamma5<-   0.00079332381730106270175
gamma6<-  -0.00023876934543019960987
gamma7<-  -0.00052728956705775104607
gamma8<-  -0.00035212335380303950960
gamma9<-  -0.000034394774418088048178
gamma10<-  0.00020533281490906479468
gamma11<-  0.00027018443954390352667
gamma12<-  0.00016727291210514019335
gamma13<- -0.000027463806603760158860
gamma14<- -0.00020920926205929994584
gamma15<- -0.00028346865532024144664
gamma16<- -0.00019969685830896977471
gamma17<-  0.000026277037109918336699
gamma18<-  0.00030736840814925282659
gamma19<-  0.00050360545304735562606
gamma20<-  0.00046634356151155944940

# Start the MC Loop:

for (ii in 1:repMC) {

# Generate Zeta dist'n. sample data

zeta.data<- draw.zeta(n, s)$y        # Using the 'UnivRNG' package. Very fast. The empirical moments of the simulated data match the theoretical moments (where they exist) very well

# Obtain the MLE of "s":

out.zeta<- zm.ll(zeta.data, N = Inf, dist = "Zeta")   # Using the 'tolerance' package
smle<- stats4::coef(out.zeta)
sest<- smle[[1]]                   # Keep value only, not name ("s") as well. This speeds up the subsequent summing.

zetas<- zeta(sest,0,1)  # using the VGAM package to evaluate the zeta function & also in the next 2 lines
zeta1<- zeta(sest,1,1)	# first derivative of the zeta function
zeta2<- zeta(sest,2,1)	# second derivative of the zeta function

# Use Choudhury's equation 20 to get the third derivative of the zeta function:
zeta3<- -(6/(sest-1)^4+gamma3-gamma4*(sest-1)+(gamma5*(sest-1)^2)/2-(gamma6*(sest-1)^3)/6+(gamma7*(sest-1)^4)/24-(gamma8*(sest-1)^5)/120
+(gamma9*(sest-1)^6)/720-(gamma10*(sest-1)^7)/5040+(gamma11*(sest-1)^8)/40320-(gamma12*(sest-1)^9)/362880+(gamma13*(sest-1)^10)/3628800
- (gamma14*(sest-1)^11)/39916800 + (gamma15*(sest-1)^12)/479001600 - (gamma16*(sest-1)^13)/6227020800 + (gamma17*(sest-1)^14)/87178291200 
- (gamma18*(sest-1)^15)/1307674368000 )
# If "s"=2, for example, the remainder term (RT) in the expansion for zeta3 satisfies |RT| < |gamma19/15!| = |0.00050360545304735562606/15!| = 3.85*10^-16
bias<- (zetas/(2*n))*(3*zetas*zeta1*zeta2 - 2*zeta1^3 - zetas^2*zeta3)/((zeta1^2 - zetas*zeta2)^2)   # This is just the second-order bias!
shat<- sest-bias                     # Compute the bias-corrected estimate

sdbar<- (sest - R*shat)/(1-R)        # compute Schucany et al. pseudo-estimator
sum_sdbar<- sum_sdbar+sdbar
sum_sdbar2<- sum_sdbar2+sdbar^2

}    
# End of MC loop

bias_sdbar<- sum_sdbar/repMC - s       # This is the MC estimate of the FULL bias of the estimator
perc_bias_sdbar<- 100*bias_sdbar/s

mse_sdbar<- sum_sdbar2/repMC-2*s*sum_sdbar/repMC+s^2
perc_mse_sdbar<- 100*mse_sdbar/(s^2)

c(s,n,repMC)
c(perc_bias_sdbar, perc_mse_sdbar)

# END (not run)