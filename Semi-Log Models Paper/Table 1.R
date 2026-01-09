# Code for the Monte Carlo Simulation Results in Table 1
####################################################

library(moments)
library(tseries)
set.seed(321)

n<- 10
nrep<- 5000
x1<- rchisq(n,df=1)
x2<- rnorm(n)
a<- 1.0
b1<- b2<- 0.1
c<- 0.25
sigsq<- 2
p<- exp(c)-1   # p = 0.2840
r_squared<- c()

# Case 1:

D1<- rep(1,n/5)
D2<- rep(0,0.8*n)
D<- c(D1,D2)



# Case 2:

#D<- rep(0,n)
#for(jj in 0:4) {
#D[n-jj]<- 1    # Just the last 5 obs. are "1"
# }
phat<- c()

for(ii in 1:nrep) {
logy<- (a+b1*x1+b2*x2+c*D+rnorm(n)*sqrt(sigsq))
OLS<- lm(logy~ x1+x2+D)
v<- (diag(vcov(OLS)))
phat[ii]<- (exp(OLS$coef[[4]])/exp(0.5*v[[4]]))-1
r_squared[ii]<- summary(OLS)$adj.r.squared
}
summary(phat)
100*(mean(phat)-p)/p
sd(phat)
skewness(phat)
kurtosis(phat)
hist(phat)

# Form inv(X'X) to get value of "d":
one<-rep(1,n)
X<- matrix(c(one, x1, x2, D), ncol = 4)
XPXINV<- solve(t(X)%*%X)
d<- XPXINV[4,4]
d
jarque.bera.test(phat)
mean(r_squared)