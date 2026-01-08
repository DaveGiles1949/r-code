# R Code for the Confidence Interval Results in Table 2
#######################################################

set.seed(321)

n<- 10
nrep<- 5000
nboot<- 1000
x1<- rchisq(n,df=1)
x2<- rnorm(n)
a<- 1.0
b1<- b2<- 0.1
c<- 0.25
sigsq<- 2
p<- exp(c)-1   # p = 0.2840
D1<- rep(1,n/5)
D2<- rep(0,0.8*n)
D<- c(D1,D2)
phat<- c()
vtilde<- c()
CI_naive<- matrix(nrow=nrep,ncol=2)
num<- c()
CL<- c()
CU<- c()
CI_approx<- c()
y_boot<- c()
resid_boot<- c()
p_boot<- c()

# Start of MC loop:
##################

for(ii in 1:nrep) {

logy<- a+b1*x1+b2*x2+c*D+rnorm(n)*sqrt(sigsq)
OLS<- lm(logy ~ x1+x2+D)
v<- (diag(vcov(OLS)))
phat[ii]<- (exp(OLS$coef[[4]])/exp(0.5*v[[4]]))-1       # the intercept is the 1st. regressor, D is the 4th. 

# First, create the naive CI results, assuming normality
########################################################

vtilde[ii]<- exp(2*OLS$coef[[4]])*(exp(-v[[4]])-exp(-2*v[[4]]))    
CI_naive[ii,1]<- phat[ii]- 1.96*sqrt(vtilde[ii])
CI_naive[ii,2]<- phat[ii]+ 1.96*sqrt(vtilde[ii])

# Does the interval that cover the true value of "p"?
num[ii]<- ( p>= CI_naive[ii,1] & p<= CI_naive[ii,2])


# Second, obtain bootstrapped CI's
######################################

resid<- OLS$residuals

# Now start the bootstrap simulation:


for (jj in 1:nboot) {

resid_boot<- sample(resid)
y_boot<- OLS$coef[[1]]+OLS$coef[[2]]*x1+OLS$coef[[3]]*x2+OLS$coef[[4]]*D+resid_boot
ols_boot<- lm(y_boot ~ x1+x2+D)
v_boot<- (diag(vcov(ols_boot)))
p_boot[jj]<- (exp(ols_boot$coef[[4]])/exp(0.5*v_boot[[4]]))-1

}   # end of the bootstrap loop

# Compute the 2.5% and 97.5% bootstrapped percentiles for p:

CL[ii]<- quantile(p_boot,  probs = 0.025)
CU[ii]<- quantile(p_boot,  probs = 0.975)

}             

# end of the MC loop

# Approximate rseults based on normality assumption:

n

CI_approx<- cbind(mean(CI_naive[,1]), mean(CI_naive[,2]))
CI_approx
# The coverage probability is the proportion of the simulated intervals that cover the truie value of "p"
CP<-sum(num)/nrep   
CP

# Results based on the bootstrap simulation:

CI_boot<- cbind(mean(CL), mean(CU))
CI_boot
