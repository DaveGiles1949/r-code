# Area-Biased Sampling: Gamma distribution
############################################
set.seed(1234)
require(stats)
require(ggamma)

n<- 90  # Sample size
nsamp<- 50000	# Number of MC replications for generating critical values
b<- 3   # b is the shape parameter for the Gamma distributions

# If b = 1 then we have the exponential distribution.

lambda_gam_1<- c() 
log_lambda_gam_1<- c()
lambda_gam_2<- c() 
log_lambda_gam_2<- c()

# The critical values are independent of scale (or rate)
##########################################################
# Create the critical values - these depend on the value of n, and alpha.  

f<- function(x, n, y, ybar) {n*log(x/ybar)+sum(log(y))-n*digamma(x) }

# The log-transform of the lambda statistic is needed to avoid numerical issues in 'n' is large

for (i in 1: nsamp) {
y<- rgamma(n, shape=b, rate = 1)		 # Gamma variates
ybar<- mean(y)
b_hat<- uniroot(f, c(0.1, 200), n=n,y=y,ybar=ybar)$root
log_lambda_gam_1[i]<- (1/n)*sum(log(y))+0.5*log(n)-0.5*log(sum(y^2))
log_lambda_gam_2[i]<- (1/n)*sum(log(y))-log(ybar)-0.5*log(1+1/b_hat)
}

lambda_gam_1<- exp(log_lambda_gam_1)	   # Reverse the log transformation
lambda_gam_2<- exp(log_lambda_gam_2)	
crit_gam_1<- quantile(lambda_gam_1,  probs = c(0.90, 0.95, 0.99))
crit_gam_2<- quantile(lambda_gam_2,  probs = c(0.90, 0.95, 0.99))

b
n
crit_gam_1
crit_gam_2