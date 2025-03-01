# Area-Biased Sampling: Weibull distribution
############################################
set.seed(1234)
require(stats)
require(ggamma)

n<-150   # Sample size
nsamp<- 50000	# Number of MC replications for generating critical values
k<- 0.5   # k is the shape parameter for the Weibull distributions

theta<- 1   # Scale parameter - results invariant to value of this parameter
            # If k = 1 then we have the exponential distribution.

lambda_wei_1<- c() 
log_lambda_wei_1<- c()
lambda_wei_2<- c() 
log_lambda_wei_2<- c()


# The critical values are independent of scale (and rate) and hence independent of "beta"

##########################################################
# Create the critical values - these depend on the values of n, and k.

# Here is the likelihood equation for the concentrated log-likelihood, as a function of just "k"   
# 'x' is the unknown (parameter)
f<- function(x, y) {(1/x)+mean(log(y))-sum(y^x*log(y))/sum(y^x) }

# The log-transform of the lambda statistic is needed to avoid numerical issues in 'n' is large

for (i in 1: nsamp) {
y<- rggamma(n,theta,k,1)		 # Weibull random variates
ybar<- mean(y) 
k_hat<- uniroot(f,c(0.1, 20),y=y )$root
theta_hat<- (mean(y^k_hat))^(1/k_hat)
log_lambda_wei_1[i]<- (1/n)*sum(log(y))+0.5*log(n)-0.5*log(sum(y^2))
log_lambda_wei_2[i]<- (1/n)*sum(log(y))-log(theta_hat)-0.5*log(gamma(1+2/k_hat))

}

lambda_wei_1<- exp(log_lambda_wei_1)	   # Invert the log transformation
lambda_wei_2<- exp(log_lambda_wei_2)
crit_wei_1<- quantile(lambda_wei_1,  probs = c(0.90, 0.95, 0.99))
crit_wei_2<- quantile(lambda_wei_2,  probs = c(0.90, 0.95, 0.99))

k
n
crit_wei_1
crit_wei_2