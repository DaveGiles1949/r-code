# Monte Carlo simulation Experiment to Investigate the Effectiveness of the (First-Order) Cox-Snell & Firth Bias-Corrections of the MLE 
# of the Parameter of the Zeta Distribution
#
# David Giles, September 2024
# -------------------------

library(VGAM)
library (tolerance)
library(UnivRNG)

set.seed(123)
repMC<- 100000		# For stable results, the number of MC replications should be at least 50,000
				# Our results were based on 100,000 MC replications

s<- 1.25    # True value of the Zeta-distribution parameter. Clauset et al. find values between 1.7 and 3.7 for a range of data-sets.
n<- 10    # Sample size.  Clauset et al. have sample sizes above xmin. ranging from 39 to 1,000 and then very much larger.

sum_stilde<- sum_shat<- sum_scup<- sum_stilde2<- sum_shat2<- sum_scup2<- sum_shat_times_scup<- 0
              
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
start_time <- Sys.time()

for (ii in 1:repMC) {

# Generate Zeta dist'n. sample data

#zeta.data <- rzeta(n,s)             # Uses the VGAM package
#zeta.data <- rzipfman(n = n, s = s, N = Inf) # Uses the tolerance package
zeta.data<- draw.zeta(n, s)$y        # Using the UnivRNG package. Very fast. The empirical moments of the simulated data match the theoretical moments (where they exist) very well
sumlogdata<- sum(log(zeta.data))     # We need this for Firth's estimator

# Obtain the MLE of "s":

out.zeta<- zm.ll(zeta.data, N = Inf, dist = "Zeta")
smle<- stats4::coef(out.zeta)
sest<- smle[[1]]                   # Keep value only, not name ("s") as well. This speeds up the subsequent summing.
sum_stilde<- sum_stilde+sest       # Keep cumulative info. for bias and MSE calculations at end of MC simulations.
sum_stilde2<- sum_stilde2+sest^2

zetas<- zeta(sest,0,1)  # using the VGAM package to evaluate the zeta function & also in the next 2 lines
zeta1<- zeta(sest,1,1)	# first derivative of the zeta function
zeta2<- zeta(sest,2,1)	# second derivative of the zeta function

# Use Choudhury's equation 20 to get the third derivative of the zeta function:
zeta3<- -(6/(sest-1)^4+gamma3-gamma4*(sest-1)+(gamma5*(sest-1)^2)/2-(gamma6*(sest-1)^3)/6+(gamma7*(sest-1)^4)/24-(gamma8*(sest-1)^5)/120
+(gamma9*(sest-1)^6)/720-(gamma10*(sest-1)^7)/5040+(gamma11*(sest-1)^8)/40320-(gamma12*(sest-1)^9)/362880+(gamma13*(sest-1)^10)/3628800
- (gamma14*(sest-1)^11)/39916800 + (gamma15*(sest-1)^12)/479001600 - (gamma16*(sest-1)^13)/6227020800 + (gamma17*(sest-1)^14)/87178291200 
- (gamma18*(sest-1)^15)/1307674368000 )
# If "s"=2, for example, the remainder term (R) in the expansion for zeta3 satisfies |R| < |gamma19/15!| = |0.00050360545304735562606/15!| = 3.85*10^-16
bias<- (zetas/(2*n))*(3*zetas*zeta1*zeta2 - 2*zeta1^3 - zetas^2*zeta3)/((zeta1^2 - zetas*zeta2)^2)   # This is just the second-order bias!
shat<- sest-bias                     # Compute the bias-corrected estimate
sum_shat<- sum_shat+shat             # Keep cumulative info. for bias and MSE calculations at end of MC simulations
sum_shat2<- sum_shat2+shat^2

# Compute the Firth estimator:
# The general form -
#f1<- function(x)  -sumlogdata -n*zeta(x,1,1)/zeta(x,0,1) + (3*zeta(x,0,1)*zeta(x,1,1)*zeta(x,2,1) - 2*zeta(x,1,1)^3 - zeta(x,0,1)^2*

#( -(6/(x-1)^4+gamma3-gamma4*(x-1)+(gamma5*(x-1)^2)/2-(gamma6*(x-1)^3)/6+(gamma7*(x-1)^4)/24-(gamma8*(x-1)^5)/120
+(gamma9*(x-1)^6)/720-(gamma10*(x-1)^7)/5040+(gamma11*(x-1)^8)/40320-(gamma12*(x-1)^9)/362880+(gamma13*(x-1)^10)/3628800
#- (gamma14*(x-1)^11)/39916800 + (gamma15*(x-1)^12)/479001600 - (gamma16*(x-1)^13)/6227020800 + (gamma17*(x-1)^14)/87178291200 
#- (gamma18*(x-1)^15)/1307674368000 )  ))/ (2*zeta(x,0,1)*(zeta(x,1,1)^2-zeta(x,0,1)*zeta(x,2,1)))

#sfirth<- uniroot(f1 , c(1.0001,7), tol = 1e-6)$root

#Taking account of the fact that we have a canonical linear exponential model -

f<- function (x) -sumlogdata -(n+1)*zeta(x,1,1)/zeta(x,0,1) +
(zeta(x,0,1)*
( -(6/(x-1)^4+gamma3-gamma4*(x-1)+(gamma5*(x-1)^2)/2-(gamma6*(x-1)^3)/6+(gamma7*(x-1)^4)/24-(gamma8*(x-1)^5)/120
+(gamma9*(x-1)^6)/720-(gamma10*(x-1)^7)/5040+(gamma11*(x-1)^8)/40320-(gamma12*(x-1)^9)/362880+(gamma13*(x-1)^10)/3628800
- (gamma14*(x-1)^11)/39916800 + (gamma15*(x-1)^12)/479001600 - (gamma16*(x-1)^13)/6227020800 + (gamma17*(x-1)^14)/87178291200 
- (gamma18*(x-1)^15)/1307674368000 )) 

-zeta(x,1,1)*zeta(x,2,1)) /(2*(zeta(x,0,1)*zeta(x,2,1) - zeta(x,1,1)^2))

scup<- uniroot(f , c(1.00001,8.0), tol = 1e-5)$root

# scup and sfirth are equal in value !!
# -------------------------------------

sum_scup<- sum_scup+scup
sum_scup2<- sum_scup2+scup^2
sum_shat_times_scup<- sum_shat_times_scup + shat*scup

}    
# End of MC loop

end_time <- Sys.time()

bias_stilde<- sum_stilde/repMC - s       # This is the MC estimate of the FULL bias of the MLE, etc.
bias_shat<- sum_shat/repMC - s
bias_scup<- sum_scup/repMC - s
perc_bias_stilde<- 100*bias_stilde/s
perc_bias_shat<- 100*bias_shat/s
perc_bias_scup<- 100*bias_scup/s

mse_stilde<- sum_stilde2/repMC-2*s*sum_stilde/repMC+s^2
mse_shat<- sum_shat2/repMC-2*s*sum_shat/repMC+s^2
mse_scup<- sum_scup2/repMC-2*s*sum_scup/repMC+s^2
perc_mse_stilde<- 100*mse_stilde/(s^2)
perc_mse_shat<- 100*mse_shat/(s^2)
perc_mse_scup<- 100*mse_scup/(s^2)

# Form a weighted unbiased estimator. (See Shucany et al., JASA, 1971)

R<- bias_shat/bias_scup    # This is actually just an estimate of R. (Of course,the ratio of the % biases also equals R.)
bias_sstar<- (bias_shat - R*bias_scup)/(1-R)
perc_bias_sstar<- 100*bias_sstar/s    # This will always be exactly zero, by construction!
# Now construct the MSE for the weighted estimator. (This is simply the variance of the estimator, as it is unbiased.)

var_shat<- mse_shat - bias_shat^2
var_scup<- mse_scup - bias_scup^2 
cov_shat_scup<- (sum_shat_times_scup/repMC) - (sum_shat*sum_scup/repMC^2)
mse_sstar<- (var_shat+R^2*var_scup-2*R*cov_shat_scup)/(1-R)^2              # MSE = Variance in this case!
perc_mse_sstar<- 100*mse_sstar/s^2

c(s,n,repMC)
c(perc_bias_stilde, perc_bias_shat, perc_bias_scup, perc_bias_sstar)
c(perc_mse_stilde, perc_mse_shat, perc_mse_scup, perc_mse_sstar)

end_time - start_time

# END (not run)