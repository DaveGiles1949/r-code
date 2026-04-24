library(VGAM)

################### FIRST APPLICATION #####################

rh40<- read.table("C:/Users/OEM/Sync/Size-Biased sampling/Weibull bias-adjusted paper/Rod-hours.txt",header = TRUE)
rod_hours40<- rh40$Rod_Hours
anglers40<- rh40$Anglers
hours_per_angler40<- rod_hours40/anglers40            
x<- hours_per_angler40                                
n<- length(x)

# Now test H0 against H1:
# %%%%%%%%%%%%%%%%%%%%%%

# 1. Exponential distribution
#
crit40<- c(0.6632, 0.6877, 0.7317)  # n=40
xbar<- mean(x)
theta_tilde<- 1/xbar	                             # MLE
theta_tilde_1<- 2/xbar                               # SB MLE
ltilde<-  -n*(log(xbar)+1)                           # log-likelihood, MLE
ltilde_1<- -2*n*(-log(2)+1+log(xbar))+sum(log(x))     # log-likelihood, SB MLE
AIC_tilde<- 2-2*ltilde
AIC_tilde_1<- 2-2*ltilde_1
lambda1_tilde<- lambda1_hat<- ((prod(x))^(1/n))/xbar
AIC_tilde
AIC_tilde_1
lambda1_tilde
lambda1_hat
crit40

# 2. Half-Normal distribution
#
crit40_tilde<- c(0.7565, 0.7774, 0.8148)  # n=40
crit40_hat<- c(0.7935,0.8251,0.8825)
xbar<- mean(x)
sigma_tilde<- sqrt(mean(sum(x^2)))   # MLE
sigma_tilde_1<- sqrt(sum(x^2)/(2*n))  # SB MLE
ltilde<-   -n*log(sigma_tilde)+0.5*n*log(2)-0.5*n*log(pi)-0.5*sum(x^2)/sigma_tilde_1^2  # log-likelihood, MLE
ltilde_1<- -2*n*log(sigma_tilde_1)-0.5*n*log(2)+sum(log(x))-(0.5/sigma_tilde_1^2)*sum(x^2)     # log-likelihood, SB MLE
AIC_tilde<- 2-2*ltilde
AIC_tilde_1<- 2-2*ltilde_1
lambda1_tilde<- ((prod(x))^(1/n))/xbar
lambda1_hat<- ((prod(x))^(1/n))/(sqrt(mean(x^2)*2/pi)) 
AIC_tilde
AIC_tilde_1
lambda1_tilde
crit40_tilde
lambda1_hat
crit40_hat

# 3. Rayleigh distribution
#
crit40_tilde<- c(0.8921, 0.9018, 0.9187)  # n=40
crit40_hat<- c(0.9191,0.9357,0.9653)
xbar<- mean(x)
sigma_tilde<- sqrt(sum(x^2)/(2*n))	      # MLE
sigma_tilde_1<- sqrt(sum(x^2)/(3*n))      # SB MLE
ltilde<- -2*n*log(sigma_tilde)+sum(log(x))-0.5*sum(x^2)/sigma_tilde^2    # log-likelihood, MLE
ltilde_1<- -3*n*log(sigma_tilde_1)-0.5*n*log(2)+sum(log(x^2))-0.5*sum(x^2)/sigma_tilde_1^2-0.5*log(pi)-log(0.5)   # log-likelihood, SB MLE
AIC_tilde<- 2-2*ltilde
AIC_tilde_1<- 2-2*ltilde_1
lambda1_tilde<- ((prod(x))^(1/n))/xbar
lambda1_hat<- ((prod(x))^(1/n))/(sigma_tilde*sqrt(pi/2))   # CHECK THIS!
AIC_tilde
AIC_tilde_1
lambda1_tilde
crit40_tilde
lambda1_hat
crit40_hat

# THE AIC RESULTS FAVOUR THE RAYLEIGH DISTRIBUTION
# FOR THAT DISTRIBUTION, BOTH TEST STATISTICS REJECT H0 IN FAVOUR OF H1 AT THE 5% LEVEL - 
# THAT IS, IN FAVOUR OF LENGTH-BIASED SAMPLING

# create density for length-biased Rayleigh

f1<- function(x){
(1/sqrt(2))*(1/sigma_tilde_1)^3*x^2*exp(-x^2/(2*sigma_tilde_1^2))/(0.5*sqrt(pi)) 
}

# MLE and asymptotic std. error for MLE in each Exponential case:
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

ase_sigma_tilde<- sigma_tilde/(2*n)
ase_sigma_tilde_1<- sigma_tilde_1/(3*n)

sigma_tilde
ase_sigma_tilde
sigma_tilde_1
ase_sigma_tilde_1

par(mfrow=c(2,1))

# Plot the preferred results:

hist(hours_per_angler40,freq=FALSE, ylim=c(0,0.16), cex.main=1,xlab="Rod-hours per Angler", main="Figure 4(a): Algonquin Park fishing survey")
curve(drayleigh(x, scale=sigma_tilde), 
      col="red",lty=1, lwd=2, add=TRUE, yaxt="n")
curve(f1(x), 
      col="blue", lty=2,lwd=2, add=TRUE, yaxt="n")
legend(10, 0.14, legend=c("Rayleigh", "Length-biased Rayleigh"),
       col=c("red", "blue"), lty=1:2, cex=0.8, bty="n")

############### SECOND APPLICATION ###########

tsla_data<- read.table("C:/Users/OEM/Sync/Size-Biased sampling/Weibull bias-adjusted paper/TSLAshort.txt",header = TRUE)
tsla<- tsla_data$TSLA
n<- length(tsla)     # n=23                      
x<- tsla/4
   
# Now test H0 against H1:
# %%%%%%%%%%%%%%%%%%%%%%

# 1. Exponential distribution
#
crit23<- c(0.7004,0.7311,0.7838)  # n=23
xbar<- mean(x)
theta_tilde<- 1/xbar	                             # MLE
theta_tilde_1<- 2/xbar                               # SB MLE
ltilde<-  -n*(log(xbar)+1)                           # log-likelihood, MLE
ltilde_1<- -2*n*(-log(2)+1+log(xbar))+sum(log(x))     # log-likelihood, SB MLE
AIC_tilde<- 2-2*ltilde
AIC_tilde_1<- 2-2*ltilde_1
lambda1_tilde<- lambda1_hat<- ((prod(x))^(1/n))/xbar
AIC_tilde
AIC_tilde_1
lambda1_tilde
lambda1_hat
crit23

# 2. Half-Normal distribution
#
crit23_tilde<- c(0.7868,0.8121,0.8553)  # n=23
crit23_hat<- c(0.8391,0.8799,0.9518)
xbar<- mean(x)
sigma_tilde<- sqrt(mean(sum(x^2)))   # MLE
sigma_tilde_1<- sqrt(sum(x^2)/(2*n))  # SB MLE
ltilde<-   -n*log(sigma_tilde)+0.5*n*log(2)-0.5*n*log(pi)-0.5*sum(x^2)/sigma_tilde_1^2  # log-likelihood, MLE
ltilde_1<- -2*n*log(sigma_tilde_1)-0.5*n*log(2)+sum(log(x))-(0.5/sigma_tilde_1^2)*sum(x^2)     # log-likelihood, SB MLE
AIC_tilde<- 2-2*ltilde
AIC_tilde_1<- 2-2*ltilde_1
lambda1_tilde<- ((prod(x))^(1/n))/xbar
lambda1_hat<- ((prod(x))^(1/n))/(sqrt(mean(x^2)*2/pi)) 
AIC_tilde
AIC_tilde_1
lambda1_tilde
crit23_tilde
lambda1_hat
crit23_hat

# 3. Rayleigh distribution
#
crit23_tilde<- c(0.9068,0.9184,0.9370)  # n=23
crit23_hat<- c(0.9438,0.9645,0.9993)
xbar<- mean(x)
sigma_tilde<- sqrt(sum(x^2)/(2*n))	      # MLE
sigma_tilde_1<- sqrt(sum(x^2)/(3*n))      # SB MLE
ltilde<- -2*n*log(sigma_tilde)+sum(log(x))-0.5*sum(x^2)/sigma_tilde^2    # log-likelihood, MLE
ltilde_1<- -3*n*log(sigma_tilde_1)-0.5*n*log(2)+sum(log(x^2))-0.5*sum(x^2)/sigma_tilde_1^2-0.5*log(pi)-log(0.5)   # log-likelihood, SB MLE
AIC_tilde<- 2-2*ltilde
AIC_tilde_1<- 2-2*ltilde_1
lambda1_tilde<- ((prod(x))^(1/n))/xbar
lambda1_hat<- ((prod(x))^(1/n))/(sigma_tilde*sqrt(pi/2))   # CHECK THIS!
AIC_tilde
AIC_tilde_1
lambda1_tilde
crit23_tilde
lambda1_hat
crit23_hat

# THE AIC RESULTS FAVOUR THE EXPONENTIAL DISTRIBUTION
# FOR THAT DISTRIBUTION, BOTH TEST STATISTICS REJECT H0 IN FAVOUR OF H1 AT THE 10% LEVEL, BUT NOT QUITE AT THE 5% LEVEL - 
# THAT IS, IN FAVOUR OF LENGTH-BIASED SAMPLING
#

# create densities for regular and length-biased exponentials

f<- function(x){
theta_tilde*exp(-theta_tilde*x)
}
f1<- function(x){
theta_tilde_1^2*x*exp(-theta_tilde_1*x)
}

# MLE and asymptotic std. error for MLE in each Exponential case:
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ase_theta_tilde<- 1/(n*xbar)
ase_theta_tilde_1<- sqrt(2)/(n*xbar)
theta_tilde
ase_theta_tilde
theta_tilde_1
ase_theta_tilde_1

# Plot the preferred results:

hist(x,freq=FALSE, ylim=c(0,1.4), main="Figure 4(b): Tesla share price data", xlab="Weeks / 4", cex.main=1)  
curve(f(x), 
      col="red",lty=1, lwd=2, add=TRUE, yaxt="n")
curve(f1(x), 
      col="blue", lty=2,lwd=2, add=TRUE, yaxt="n")
legend(1.4, 1.2, legend=c("Exponential", "Length-biased Exponential"),
       col=c("red", "blue"), lty=1:2, cex=0.8,bty="n")

