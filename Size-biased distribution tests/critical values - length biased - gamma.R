# Length-Biased Sampling: Gamma distribution
############################################
set.seed(1234)
require(stats)
require(ggamma)

n<- 200    # Sample size
nsamp<- 100000	# Number of MC replications for generating critical values
alpha<- 3  # alpha is the shape parameter for the Gamma distributions

# If alpha = 1 then we have the exponential distribution.

lambda_gam<- c() 
log_lambda_gam<- c()

# The critical values are independent of scale (and rate) and hence independent of "beta"
# The two versions of the lambda test are the same, because:

# E[y]= (alpha/beta)
# Using MLE,beta_tilde = alpha_tilde/ybar
# So, the MLE of E[y] = (alpha_tilde/beta_tilde) = ybar
 
##########################################################
# Create the critical values - these depend on the value of n, and alpha.  
# The log-transform of the lambda statistic is needed to avoid numerical issues in 'n' is large

for (i in 1: nsamp) {
y<- rgamma(n, shape=alpha, rate = 1)		 # Gamma
log_lambda_gam[i]<- (1/n)*sum(log(y))-log(mean(y))

}

lambda_gam<- exp(log_lambda_gam)	   # Reverse the log transformation

crit_gam<- quantile(lambda_gam,  probs = c(0.90, 0.95, 0.99))

alpha
n
crit_gam