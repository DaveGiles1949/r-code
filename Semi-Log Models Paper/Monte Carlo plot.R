# Code for Figure 3
##########################

library(moments)
library(tseries)

nrep<- 5000
a<- 1.0
b1<- b2<- 0.1
c<- 0.25
sigsq<- 2
p<- exp(c)-1   # p = 0.2840

# Case 1:
#########
# Fig 3(a)

set.seed(321)
n<-100
x1<- rchisq(n,df=1)
x2<- rnorm(n)
D1<- rep(1,n/5)
D2<- rep(0,0.8*n)
D<- c(D1,D2)
phat_a<- c()

for(ii in 1:nrep) {
logy<- (a+b1*x1+b2*x2+c*D+rnorm(n)*sqrt(sigsq))
OLS<- lm(logy~ x1+x2+D)
v<- diag(vcov(OLS))
phat_a[ii]<- (exp(OLS$coef[[4]])/exp(0.5*v[[4]]))-1

}

# Fig 3(b)

set.seed(321)
n<- 1000
x1<- rchisq(n,df=1)
x2<- rnorm(n)
D1<- rep(1,n/5)
D2<- rep(0,0.8*n)
D<- c(D1,D2)
phat_b<- c()

for(ii in 1:nrep) {
logy<- (a+b1*x1+b2*x2+c*D+rnorm(n)*sqrt(sigsq))
OLS<- lm(logy~ x1+x2+D)
v<- diag(vcov(OLS))
phat_b[ii]<- (exp(OLS$coef[[4]])/exp(0.5*v[[4]]))-1

}

# Case 2:
# Fig 3(c)

set.seed(321)
n<- 100
x1<- rchisq(n,df=1)
x2<- rnorm(n)
D<- rep(0,n)

for(jj in 0:4) {
D[n-jj]<- 1    # Just the last 5 obs. are "1"
}
phat_c<- c()

for(ii in 1:nrep) {
logy<- (a+b1*x1+b2*x2+c*D+rnorm(n)*sqrt(sigsq))
OLS<- lm(logy~ x1+x2+D)
v<- diag(vcov(OLS))
phat_c[ii]<- (exp(OLS$coef[[4]])/exp(0.5*v[[4]]))-1

}

# Fig 3(d)

set.seed(321)
n<- 1000
x1<- rchisq(n,df=1)
x2<- rnorm(n)
D<- rep(0,n)
for(jj in 0:4) {
D[n-jj]<- 1    # Just the last 5 obs. are "1"
phat_d<- c()

}
for(ii in 1:nrep) {
logy<- (a+b1*x1+b2*x2+c*D+rnorm(n)*sqrt(sigsq))
OLS<- lm(logy~ x1+x2+D)
v<- diag(vcov(OLS))
phat_d[ii]<- (exp(OLS$coef[[4]])/exp(0.5*v[[4]]))-1

}

par(mfrow=c(2,2))
hist(phat_a,  prob=TRUE,col="red",main=expression(paste("Figure 3(a): Simulated Density of  ",  hat(p))),sub="Case 1: n = 100",cex.sub=0.8,
 ylab=expression(paste("f ", (hat(p)))), breaks=30,cex.main=0.8,cex.axis=0.6,cex.lab=0.7,  xlab=expression(hat(p)),xlim=c(-1,3))
hist(phat_b, prob=TRUE, col="red", main=expression(paste("Figure 3(b): Simulated Density of  ",  hat(p))),sub="Case 1: n = 1,000",cex.sub=0.8,
 ylab=expression(paste("f ", (hat(p)))), xlim=c(-0.2,0.8), breaks=30, cex.main=0.8,cex.axis=0.6,cex.lab=0.7, xlab=expression(hat(p)))
hist(phat_c,  prob=TRUE,col="blue",main=expression(paste("Figure 3(c): Simulated Density of  ",  hat(p))),sub="Case 2: n = 100",cex.sub=0.8,
 ylab=expression(paste("f ", (hat(p)))), breaks=60,cex.main=0.8,cex.axis=0.6,cex.lab=0.7,  xlab=expression(hat(p)),xlim=c(-1,3))
hist(phat_d, prob=TRUE, col="blue", main=expression(paste("Figure 3(d): Simulated Density of  ",  hat(p))),sub="Case 2: n = 1,000",cex.sub=0.8,
 ylab=expression(paste("f ", (hat(p)))), xlim=c(-1,3), breaks=60, cex.main=0.8,cex.axis=0.6,cex.lab=0.7, xlab=expression(hat(p)))

