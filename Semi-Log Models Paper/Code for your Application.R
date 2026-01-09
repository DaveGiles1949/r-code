# Compute p-hat and a 95% confidence interval for "p" by bootstrapping, 
# rather than wrongly assuming that p-hat is normally distributed
#################################################################

# The following code can be adapted easily for any regression model
###################################################################

# Written by David Giles, <davegiles1949@gmail.com>; January 2026
#################################################################

library(tseries)   # for J-B test
library(stats)     # for quantiles

set.seed(1234)

X1<- c(1, 4, 3, 6, 2, 3, 3, 6, 1, 8)                     # Illustrative first regressor
X2<- c(0.5, 2, 1.5, 0.6, 1, 1.5, 2.5, 0.6, 3, 1.6)       # Illustrative second regressor
DUM<- c(1,0,0,0,1,1,0,0,1,0)                             # Illustrative dummy variable
y<- c(4.0, 1.6, 1.5, 3.3, 4.3, 4.5, 0.7, 3.5, 2.3, 3.8)  # Illustrative dependent variable

data<- cbind(y,X1,X2,DUM)
head(data)                          # Look at the illustrative data

OLS<- lm(log(y) ~ X1+X2+DUM)        # Basic OLS regression
summary(OLS)
jarque.bera.test(OLS$residuals)     # Test for normality of the errors (this DOES NOT imply normality of p-hat)

# Now bootstrap p-hat and get limits of a 95% Confidence Interval for "p":

nboot<- 5000
resid<- OLS$residuals
y_boot<- c()
resid_boot<- c()
p_hat_boot<- c()

for (ii in 1:nboot) {
resid_boot<- sample(resid)

y_boot<- OLS$coeff[[1]]+OLS$coeff[[2]]*X1+OLS$coeff[[3]]*X2+OLS$coeff[[4]]*DUM+resid_boot

ols_boot<- lm(y_boot ~ X1+X2+DUM)
p_hat_boot[ii]<- exp(ols_boot$coef[[4]]-0.5*vcov(ols_boot)[[4,4]])-1        # DUM is the 4th regressor (intercept is first)

}

# Compute the mean of the bootstrapped p-hat, and the 2.5% and 97.5% quantiles:

mean(p_hat_boot)
quantile(p_hat_boot,  probs = 0.025)
quantile(p_hat_boot,  probs = 0.975)  

# Compare these results with those based on just the original regression, and assuming normality 
# of the distribution of p-hat, namely:

Naive_p_hat<- exp(OLS$coef[[4]]-0.5*vcov(OLS)[[4,4]])-1
Naive_lower_limit<- p_hat-1.96*sqrt(vcov(OLS)[[4,4]])
Naive_upper_limit<- p_hat+1.96*sqrt(vcov(OLS)[[4,4]])

Naive_p_hat
Naive_lower_limit
Naive_upper_limit

# The naive estimator of "p" slightly understates the bootstrapped estimator.
# The naive confidence interval length = 0.6531; compared with the bootstrap confidence interval length of 1.8733
# Wrongly assuming normality of the distribution of p-hat results in a misleadingly short(& misleadingly "informative")
# confidence interval