# Monte Carlo Simulation Experiment to Investigate the Effectiveness of the (First-Order) Cox-Snell & Firth Bias-Corrections of the MLE 
# of the Scale Parameter of the Zeta Distribution
#
# Code written by David Giles, October 2025
# -------------------------

library(VGAM)
library (tolerance)
library(UnivRNG)

set.seed(123)
repMC<- 100000		

s<- 1.25   # True value of the Zeta-distribution parameter. Clauset et al. find values between 1.7 and 3.7 for a range of data-sets.
n<- 10     # Sample size.  Clauset et al. have sample sizes above xmin. ranging from 39 to 1,000 and then very much larger.


sum_stilde<- sum_shat<- sum_scup<- sum_stilde2<- sum_shat2<- sum_scup2<- 0
sum_sdbar<- sum_sdbar2<-0   # Schucany

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
sumlogdata<- sum(log(zeta.data))     # We need this for Firth's estimator

# Obtain the MLE of "s":

out.zeta<- zm.ll(zeta.data, N = Inf, dist = "Zeta")   # Using the 'tolerance' package
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

#Taking account of the fact that we have a canonical linear exponential model -

f<- function (x) -sumlogdata -(n+1)*zeta(x,1,1)/zeta(x,0,1) +
(zeta(x,0,1)*
( -(6/(x-1)^4+gamma3-gamma4*(x-1)+(gamma5*(x-1)^2)/2-(gamma6*(x-1)^3)/6+(gamma7*(x-1)^4)/24-(gamma8*(x-1)^5)/120
+(gamma9*(x-1)^6)/720-(gamma10*(x-1)^7)/5040+(gamma11*(x-1)^8)/40320-(gamma12*(x-1)^9)/362880+(gamma13*(x-1)^10)/3628800
- (gamma14*(x-1)^11)/39916800 + (gamma15*(x-1)^12)/479001600 - (gamma16*(x-1)^13)/6227020800 + (gamma17*(x-1)^14)/87178291200 
- (gamma18*(x-1)^15)/1307674368000 )) 

-zeta(x,1,1)*zeta(x,2,1)) /(2*(zeta(x,0,1)*zeta(x,2,1) - zeta(x,1,1)^2))

scup<- uniroot(f , c(1.00001,8.0), tol = 1e-5)$root

sum_scup<- sum_scup+scup
sum_scup2<- sum_scup2+scup^2

}    
# End of MC loop

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

c(s,n,repMC)
c(perc_bias_stilde, perc_bias_shat, perc_bias_scup)
c(perc_mse_stilde, perc_mse_shat, perc_mse_scup)

# END (not run)