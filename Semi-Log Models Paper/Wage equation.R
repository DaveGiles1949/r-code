# Code for Wage Equation Results in Tables 3 and 4
##################################################

library(lmtest)  # for ; RESET test
library(tseries)  # for J-B test
library(AFR)       # for BPG test
library(sandwich) # for N-W std. errors
library(stats)    # for quantiles

set.seed(1234)

# Results for Table 3
#####################
wage_data<- read.table("C:/Users/OEM/Sync/Semi-Log Models/Revision/CPS78.txt", header=TRUE)
head(wage_data)

# Get subset of data associated with HISP=1
subdat <- subset(wage_data, subset = HISP == 1)
OLS<- lm(LNWAGE ~ UNION+MANAG+PROF+SALES+SERV+FE+ED+EX+EXSQ, subdat)
NW<-coeftest(OLS, vcov.=NeweyWest(OLS, lag=0, prewhite=FALSE, adjust=TRUE, verbose=TRUE))
NW
jarque.bera.test(OLS$residuals)
bp(OLS)  # BPG "nR^2" test
resettest(OLS, power=c(2,3,4), "fitted")

# Results for Table 4
#####################

# First, the naive CI results

perc_impact<- c()

CI_naive<- matrix(nrow=7, ncol=2)
t<- qt(0.975,26) #   dof = (36-10) = 26; quantile for 2.5% in one tail (95% CI)

for (ii in 2:7) {
perc_impact[ii]<- 100*OLS$coef[[ii]]
CI_naive[ii,]<- 100*c(OLS$coef[[ii]]- t*NW[[ii,2]], OLS$coef[[ii]]+t*NW[[ii,2]])   # the intercept is first regressor in the NW results

}
perc_impact
CI_naive

# Second, the almost unbiased percentage impact, and bootstrapped CI's

nboot<- 5000


# Now bootstrap each mid-point and get limits of 95% CI:

resid<- OLS$residuals
y_boot<- c()
resid_boot<- c()
coeff_boot<- matrix(nrow=nboot, ncol=10)
phat<- c()
se_boot<- matrix(nrow=nboot, ncol=10)
perc_boot<- matrix(nrow=nboot, ncol=7)
lower<- c()
upper<- c()

UNION<- subdat$UNION
MANAG<- subdat$MANAG
PROF<- subdat$PROF
SALES<- subdat$SALES
SERV<- subdat$SERV
FE<- subdat$FE
ED<- subdat$ED
EX<- subdat$EX
EXSQ<- subdat$EXSQ

for (ii in 1:nboot) {
resid_boot<- sample(resid)

y_boot<- OLS$coeff[[1]]+OLS$coeff[[2]]*UNION+OLS$coeff[[3]]*MANAG+OLS$coeff[[4]]*PROF+
OLS$coeff[[5]]*SALES+OLS$coeff[[6]]*SERV+OLS$coeff[[7]]*FE+OLS$coeff[[8]]*ED+OLS$coeff[[9]]*EX+OLS$coeff[[10]]*EXSQ+resid_boot

ols_boot<- lm(y_boot ~ UNION+MANAG+PROF+SALES+SERV+FE+ED+EX+EXSQ)
coeff_boot[ii,]<- ols_boot$coeff
NW_boot<- coeftest(ols_boot, vcov.=NeweyWest(ols_boot, lag=0, prewhite=FALSE, adjust=TRUE, verbose=TRUE))
# Note: NW_boot comtains std. errors that need to be squared, below

for (jj in 2:7) {
perc_boot[ii,jj]<- 100*(exp(ols_boot$coef[[jj]]-0.5*NW_boot[[jj,2]]^2)-1)

}

}

# Compute the means of the bootstrapped phats, and the 2.5% and 97.5% percentiles:

for (jj in 2:7) {
phat[jj]<- mean(perc_boot[,jj])
lower[jj]<- quantile(perc_boot[,jj],  probs = 0.025)
upper[jj]<- quantile(perc_boot[,jj],  probs = 0.975)

}

phat
lower
upper
  