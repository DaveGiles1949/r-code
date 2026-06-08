# To check Kossler et al. testing code - reproduce their results for for Fibonacci numbers
#
#############################################
#
# Code written by David Giles for the paper "Benford’s Law and Regression Results Published in
# Articles in New Zealand Economic Papers", last updated June 2026.
#
# Contact: David Giles; davegiles1949@gmail.com; davegiles.ca
# -------
############################################################
#  First digit tests
####################
# Kossler et al. invariance test (2024)
# ------------------------------------

library(benford.analysis)
library(BenfordTests)
order_of_magnitude <- function(x){
  if (x==0){
    return(0)
  }
  else if (x< 0){
    x = -1 * x
  }
  return(floor(log10(x)))
}

order<- c()

# Fibonacci Test:
x<- c()
n<- 1000
x[1]<- 1
x[2]<- 1
order[1]<- 0
order[2]<- 0
for (i in 3:n){
x[i]<- x[i-1]+x[i-2]
order[i]<- order_of_magnitude(x[i])
}
S1<- x/(10^order)

d1<- 1:9
E1<- 1/rep(log(10),9)        # Expected value of S1(S)
V1<- (2*d1+1)/(2*log(10))-(1/log(10))^2       # variance of S1(S)
cov1<-  -1/(log(10))^2		# covariance between any S1(S) terms
CORR1<- matrix(nrow=9,ncol=9)
fd<-signifd(x = x, digits = 1)
Sum1<- c()

for (j in 1:9){
Sum1[j]<- sum(S1[fd==j])  
for (i in 1:9){
CORR1[i,j]<- cov1/sqrt(V1[i]*V1[j])
} 
CORR1[j,j]<- 1                       # correlation matrix of Sum - matches Table 10
}

R1<- (Sum1 - n*E1)/sqrt(n*V1)       
IS_1E<- t(R1)%*%R1
IS_1E
IS_1M<- t(R1)%*%solve(CORR1)%*%R1
IS_1M

# These values for the test statistics match those in Table 7 of the 2024 paper
###################################

