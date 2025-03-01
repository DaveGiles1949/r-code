# Area-Biased Sampling
######################
set.seed(1234)
require(stats)
require(ggamma)

n<- 40          # Sample size
nsamp<- 100000	# Number of MC replications for generating critical values

log_lambda_exp_1<- c()
log_lambda_exp_2<- c()
log_lambda_ray_1<- c()
log_lambda_ray_2<- c()
log_lambda_hn_1<- c()
log_lambda_hn_2<- c()
lambda_exp_1<- c() # Use the sample means in the denominator of lambda
lambda_exp_2<- c()
lambda_hn_1<- c()
lambda_ray_1<- c()
lambda_exp_2<- c() # Use the paarameter MLEs in the denominator of lambda
lambda_hn_2<- c()
lambda_ray_2<- c()

# The critical values are independent of scale and hence independent of "theta" and "sigma"
theta<- 2.0   # theta is the rate parameter for Exponential
sigma<- 2.0   # sigma is the scale parameter for the half-normal and Rayleigh distributions
a<- sqrt(2*sigma^2)

# Create the critical values - these depend on the value of n, and the distribution.  

for (i in 1: nsamp) {

#####################################################################################################
# The log transformation is needed to avoid numerical issues when "n" is large #
#####################################################################################################

x1<- rexp(n,theta)				  # Exponential
log_lambda_exp_1[i]<- sum(log(x1))/n+0.5*log(n)-0.5*log(sum(x1^2))
log_lambda_exp_2[i]<- sum(log(x1))/n-log(mean(x1))-0.5*log(2)

x2<- rggamma(n,a,2,0.5)                     # Half-Normal
log_lambda_hn_1[i]<- sum(log(x2))/n+0.5*log(n)-0.5*log(sum(x2^2))
#log_lambda_hn_2[i]<- sum(log(x2))/n+0.5*log(n)-0.5*log(sum(x2^2))
# Same in both cases

x3<- rggamma(n,a,2,1)                       # Rayleigh
log_lambda_ray_1[i]<- sum(log(x3))/n+0.5*log(n)-0.5*log(sum(x3^2))
#log_lambda_ray_2[i]<- sum(log(x3))/n+0.5*log(n)-0.5*log(sum(x3^2))
# Same in both cases

}

lambda_exp_1<- exp(log_lambda_exp_1)	   # Invert the log transformations
lambda_exp_2<- exp(log_lambda_exp_2)
lambda_hn_1<-  exp(log_lambda_hn_1)
#lambda_hn_2<-  exp(log_lambda_hn_2)
lambda_ray_1<- exp(log_lambda_ray_1)
#lambda_ray_2<- exp(log_lambda_ray_2)

crit_exp_1<- quantile(lambda_exp_1,  probs = c(0.90, 0.95, 0.99))
crit_hn_1 <- quantile(lambda_hn_1,  probs = c(0.90, 0.95, 0.99))
crit_ray_1<- quantile(lambda_ray_1,  probs = c(0.90, 0.95, 0.99))

crit_exp_2<- quantile(lambda_exp_2,  probs = c(0.90, 0.95, 0.99))
#crit_hn_2<-  quantile(lambda_hn_2,  probs = c(0.90, 0.95, 0.99))
#crit_ray_2<- quantile(lambda_ray_2,  probs = c(0.90, 0.95, 0.99))

n
crit_exp_1
crit_exp_2
crit_hn_1
#crit_hn_2
crit_ray_1
#crit_ray_2

# _1 implies that sample mean used in the denominator of lambda
# _2 implies that the MLE of the population mean is used in the denominator of lambda
