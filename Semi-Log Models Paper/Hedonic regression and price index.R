# Code for Hedonic Price Results in Tables 5 and 6
################################################## 

library(lmtest)  # for ; RESET test
library(tseries)  # for J-B test
library(AFR)       # for BPG test
library(sandwich) # for N-W std. errors
library(stats)    # for quantiles

set.seed(1234)

# Results for Table 5
#####################
data1<- read.table("C:/Users/OEM/Sync/Semi-Log Models/Revision/COLE1.txt", header=TRUE)
head(data1)
data2<- read.table("C:/Users/OEM/Sync/Semi-Log Models/Revision/COLE2.txt", header=TRUE)
head(data2)
my_data<- cbind(data1,data2)

# First use the full sample:


OLS<- lm(log(PRICE) ~ log(SPEED)+log(CAP)+D73+D74+D75+D76+D77+D78+D79+D80+D81+D82+D83+D84, my_data)
summary(OLS)
NW<-coeftest(OLS, vcov.=NeweyWest(OLS, lag=0, prewhite=FALSE, adjust=TRUE, verbose=FALSE))
NW
jarque.bera.test(OLS$residuals)
bp(OLS)  # BPG "nR^2" test
resettest(OLS, power=c(2,3,4), "fitted")

# Now get subset of data for 1973 - 1976:

subdat <- subset(my_data, subset = (D73 == 1 | D74==1 | D75==1 | D76==1))

OLS_sub<- lm(log(PRICE) ~ log(SPEED)+log(CAP)+D74+D75+D76, subdat)
summary(OLS_sub)
NW_sub<-coeftest(OLS_sub, vcov.=NeweyWest(OLS_sub, lag=0, prewhite=FALSE, adjust=TRUE, verbose=FALSE))
NW_sub
jarque.bera.test(OLS_sub$residuals)
bp(OLS_sub)  # BPG "nR^2" test
resettest(OLS_sub, power=c(2,3,4), "fitted")

# Results for Table 6
#####################

# 1. Normal Approximation:
##########################
# First, use the full sample, and then the sub-sample
# In what follows, the "-1" is not needed for "phat" because it is associated with the PROPORTIONAL (not ACTUAL) change
# We are multipyling by 100 becuase that is the base-period value of the price index

index<- c()
stdev<- c()
up<- c()
lo<- c()
index[3]<- 100

# Note, below, NW[,2] is a vector of standard errors, and needs to be squared

for (ii in 4:15) {
index[ii]<- 100*( exp(OLS$coef[ii] - 0.5*NW[ii,2]^2 ))
stdev[ii]<- sqrt( exp(2*OLS$coef[ii])*( exp(-NW[ii,2]^2)-exp(-2*NW[ii,2]^2) ) )
up[ii]<- 100*( index[ii]/100+1.96*stdev[ii])
lo[ii]<- 100*( index[ii]/100-1.96*stdev[ii])

}
lo
up
up-lo   # length of CI

# Now the sub-sample

index_sub<- c()
stdev_sub<- c()
up_sub<- c()
lo_sub<- c()
index_sub[3]<- 100
# Note, below, NW_sub[,2] is a vecor of standard errors, and needs to be squared

for (ii in 4:6) {
index_sub[ii]<- 100*( exp(OLS_sub$coef[ii] - 0.5*NW_sub[ii,2]^2 ))
stdev_sub[ii]<- sqrt( exp(2*OLS_sub$coef[ii])*( exp(-NW_sub[ii,2]^2)-exp(-2*NW_sub[ii,2]^2) ) )
up_sub[ii]<- 100*( index_sub[ii]/100+1.96*stdev_sub[ii])
lo_sub[ii]<- 100*( index_sub[ii]/100-1.96*stdev_sub[ii])
}

lo_sub
up_sub
up_sub-lo_sub  # length of CI

# 2. Bootstrap results:
#######################
# First, use the full sample, and then the sub-sample

nboot<- 5000

# Now bootstrap to get each mid-point and get limits of 95% CI:
# First, for the full sample:

y_boot<- c()
resid_boot<- c()
coeff_boot<- matrix(nrow=nboot, ncol=15)
index_boot<- matrix(nrow=nboot, ncol=15)
lower<- c()
upper<- c()
centre<- c()

SPEED<- my_data$SPEED
CAP<- my_data$CAP
D73<- my_data$D73
D74<- my_data$D74
D75<- my_data$D75
D76<- my_data$D76
D77<- my_data$D77
D78<- my_data$D78
D79<- my_data$D79
D80<- my_data$D80
D81<- my_data$D81
D82<- my_data$D82
D83<- my_data$D83
D84<- my_data$D84

# In what follows, the "-1" is not needed for "phat" because it is associated with the PROPORTIONAL (not ACTUAL) change
# We are multipyling by 100 becuase that is the base-period value of the price index

for (ii in 1:nboot) {
resid_boot<- sample(OLS$residuals)

y_boot<- OLS$coeff[[1]]+OLS$coeff[[2]]*log(SPEED)+OLS$coeff[[3]]*log(CAP)+OLS$coeff[[4]]*D73+
OLS$coeff[[5]]*D74+OLS$coeff[[6]]*D75+OLS$coeff[[7]]*D76+OLS$coeff[[8]]*D77+OLS$coeff[[9]]*D78+OLS$coeff[[10]]*D79+
OLS$coeff[[11]]*D80+OLS$coeff[[12]]*D81+OLS$coeff[[13]]*D82+OLS$coeff[[14]]*D83+OLS$coeff[[15]]*D84+resid_boot

ols_boot<- lm(y_boot ~ log(SPEED)+log(CAP)+D73+D74+D75+D76+D77+D78+D79+D80+D81+D82+D83+D84)
coeff_boot[ii,]<- ols_boot$coeff
NW_boot<- coeftest(ols_boot, vcov.=NeweyWest(ols_boot, lag=0, prewhite=FALSE, adjust=TRUE, verbose=FALSE))

# Note: NW_boot contains std. errors that need to be squared, below

index_boot[,3]<- 100
for (jj in 4:15) {
index_boot[ii,jj]<- 100*(exp(ols_boot$coef[[jj]]-0.5*NW_boot[[jj,2]]^2))

}

}

# Compute the 2.5% and 97.5% percentiles of the bootstrapped index:

for (jj in 4:15) {

centre[jj]<- mean(index_boot[,jj])
lower[jj]<- quantile(index_boot[,jj],  probs = 0.025)
upper[jj]<- quantile(index_boot[,jj],  probs = 0.975)

}

lower
centre
upper
upper-lower    # length of CI

# Now repeat with the sub-sample:
#################################

y_boot_sub<- c()
resid_boot_sub<- c()
coeff_boot_sub<- matrix(nrow=nboot, ncol=6)
se_boot_sub<- matrix(nrow=nboot, ncol=6)
index_boot_sub<- matrix(nrow=nboot, ncol=6)
lower_sub<- c()
upper_sub<- c()
centre_sub<- c()

SPEED<- subdat$SPEED
CAP<- subdat$CAP
D74<- subdat$D74
D75<- subdat$D75
D76<- subdat$D76

# In what follows, the "-1" is not needed for "phat" because it is associated with the PROPORTIONAL (not ACTUAL) change
# We are multipyling by 100 becuase that is the base-period value of the price index

for (ii in 1:nboot) {
resid_boot_sub<- sample(OLS_sub$residuals)

y_boot_sub<- OLS_sub$coeff[[1]]+OLS_sub$coeff[[2]]*log(SPEED)+OLS_sub$coeff[[3]]*log(CAP)+OLS_sub$coeff[[4]]*D74+
OLS_sub$coeff[[5]]*D75+OLS_sub$coeff[[6]]*D76+resid_boot_sub
ols_boot_sub<- lm(y_boot_sub ~ log(SPEED)+log(CAP)+D74+D75+D76)
coeff_boot_sub[ii,]<- ols_boot_sub$coeff
NW_boot_sub<- coeftest(ols_boot_sub, vcov.=NeweyWest(ols_boot_sub, lag=0, prewhite=FALSE, adjust=TRUE, verbose=FALSE))

# Note: NW_boot_sub contains std. errors that need to be squared, below

index_boot_sub[,3]<- 100
for (jj in 4:6) {
index_boot_sub[ii,jj]<- 100*(exp(ols_boot_sub$coef[[jj]]-0.5*NW_boot_sub[[jj,2]]^2))

}

}

# Compute the 2.5% and 97.5% percentiles of the bootstrapped index:

for (jj in 4:6) {

centre_sub[jj]<- mean(index_boot_sub[,jj])
lower_sub[jj]<- quantile(index_boot_sub[,jj],  probs = 0.025)
upper_sub[jj]<- quantile(index_boot_sub[,jj],  probs = 0.975)

}

lower_sub
centre_sub
upper_sub
upper_sub-lower_sub    # length of CI
