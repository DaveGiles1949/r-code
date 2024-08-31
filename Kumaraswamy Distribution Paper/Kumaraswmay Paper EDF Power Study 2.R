require(GoFKernel)    	 # For inverting a cdf
require(univariateML)  	 # For the Kumaraswamy MLE      
require(VGAM)         	 # For the Kumaraswamy distribution 
require(stats)           # For the beta distribution
require(trapezoid)  	 # For trapezoidal
require(truncnorm)       # For truncated Normal
require(cascsim)         # For truncated Gamma
require(EnvStats)        # For truncated Log-Normal & Triangular
require(ReIns)           # For truncated Weibull

# Reference for A-D test in the case of the Beta distribution: 
# Mathias Raschke (2011), "Empirical behaviour of tests for the beta distribution and their application 
# in environmental research", Stochastic Environmental Research and Risk Assessment, 25:79–89
# (See the reference therein to his earlier paper)
#--------------------------------------------------------------

set.seed (1234)   
nrep<- 10000
n<-  100
alpha<- 3       # Values of the Kumaraswamy shape parameters used to generate critical values  
beta<- 8

y<- c()
AD<- c()
KS<- c()
CvM<- c()
Watson<- c()
Kuiper<- c()

Dplus<- function(x,N,px) {
max(seq(1:N)/N - px)
}

Dminus<- function(x,N,px) {
max(px-(seq(1:N)-1)/N)
}

s<- seq(1:n)

#  Start of simulation loop to get the percentiles (critical values) to ensure no size-distortion:
# (Thses values will depend on the values of n, a, and b)
 
for (i in 1:nrep) {

x<- rkumar(n, alpha, beta)                    # Null is true

n<- length(x) 
kumar<- mlkumar(x, a0=.1)
shape1<- kumar[[1]]
shape2<- kumar[[2]]
w<- pkumar(x, shape1, shape2)                  # Start of Raschke's testing procedure
f <- function(x) pnorm(x, mean=0, sd=1)         
f.inv <- inverse(f,lower=-25,upper=25)         

for (j in 1:n){
y[j]<- f.inv(w[j])                             

}

mu<- mean(y)
sigma<- sd(y)*sqrt((n-1)/n)
ys<- sort(y)
pys<- pnorm(ys, mean=mu, sd=sigma)
pysr<-  rev(pys)    # reverse the order of ths pdf

Dp<- Dplus(ys,n, pys)
Dm<- Dminus(ys,n, pys)
D<- max(Dp,Dm)
KS[i]<- Dstar<- D*(sqrt(n)-0.01+0.85/sqrt(n))

V<- Dp+Dm
Kuiper[i]<- Vstar<- V*(sqrt(n)+0.05+0.82/sqrt(n))

W2<- sum((pys -(2*s-1)/(2*n) )^2) + 1/(12*n)
CvM[i]<- W2star<- W2*(1+0.5/n)

U2<- W2-n*(sum(pys/n) - 0.5)^2   
Watson[i]<-U2star<- U2*(1+0.5/n)

A2<- -n-sum( (2*s-1)*(log(pys)+log(1-pysr) ) )/n
AD[i]<- A2star<- A2*(1+0.75/n+2.25/(n^2))

}

# End of first simulation loop
# Get the simulated critical values:

crit_KS<- quantile(KS, c(0.9,0.95))
crit_Watson<- quantile(Watson, c(0.9, 0.95))
crit_Kuiper<- quantile(Kuiper, c(.9, .95))
crit_CvM<- quantile(CvM, c(.9, .95))
crit_AD<- quantile(AD, c(.9, .95))

# Now start the simulation loop for the powers:
# (The results will depend on n, a. b, and the form of the alternative distribution)

for (i in 1:nrep) {

#x<- rtri(n, min = 0, max = 1, mode = 1/3)     # Null is false
x<- rbeta(n,20,20)                          # Null is false 
#x<- rtruncnorm(n, a=0, b=1, mean=0.2, sd=5)    # Null is false
#x<- rtrapezoid(n, min = 0, mode1 = 1/4, mode2 = 3/4, max = 1, n1 = 3, n3 = 3)  # Null is false
#x<- rtgamma(n, 2, 3, min = 0.0, max = 0.99)    # Null is false
#x<- rlnormTrunc(n, meanlog = 0.5, sdlog = 0.5, min = 0, max = 1)  # Null is false
#x<- rtweibull(n, shape = 2, scale=1, endpoint=1)               # Null is false  

n<- length(x) 
kumar<- mlkumar(x, a0=.1)
shape1<- kumar[[1]]
shape2<- kumar[[2]]
w<- pkumar(x, shape1, shape2)                  # Start of Raaschke's testing procedure
f <- function(x) pnorm(x, mean=0, sd=1)         
f.inv <- inverse(f,lower=-25,upper=25)         

for (j in 1:n){
y[j]<- f.inv(w[j])                             

}

mu<- mean(y)
sigma<- sd(y)*sqrt((n-1)/n)
ys<- sort(y)
pys<- pnorm(ys, mean=mu, sd=sigma)
pysr<-  rev(pys)    # reverse the order of ths pdf

Dp<- Dplus(ys,n, pys)
Dm<- Dminus(ys,n, pys)
D<- max(Dp,Dm)
KS[i]<- Dstar<- D*(sqrt(n)-0.01+0.85/sqrt(n))

V<- Dp+Dm
Kuiper[i]<- Vstar<- V*(sqrt(n)+0.05+0.82/sqrt(n))

W2<- sum((pys -(2*s-1)/(2*n) )^2) + 1/(12*n)
CvM[i]<- W2star<- W2*(1+0.5/n)

U2<- W2-n*(sum(pys/n) - 0.5)^2   
Watson[i]<-U2star<- U2*(1+0.5/n)

A2<- -n-sum( (2*s-1)*(log(pys)+log(1-pysr) ) )/n
AD[i]<- A2star<- A2*(1+0.75/n+2.25/(n^2))

}

# End of second simulation loop to enable calculation of powers

power_KS5<- sum(KS>crit_KS[2])/nrep
power_Kuiper5<- sum(Kuiper>crit_Kuiper[2])/nrep
power_CvM5<- sum(CvM>crit_CvM[2])/nrep
power_Watson5<- sum(Watson>crit_Watson[2])/nrep
power_AD5<- sum(AD>crit_AD[2])/nrep
power_KS10<- sum(KS>crit_KS[1])/nrep
power_Kuiper10<- sum(Kuiper>crit_Kuiper[1])/nrep
power_CvM10<- sum(CvM>crit_CvM[1])/nrep
power_Watson10<- sum(Watson>crit_Watson[1])/nrep
power_AD10<- sum(AD>crit_AD[1])/nrep

c(power_KS5, power_Kuiper5, power_CvM5, power_Watson5, power_AD5)
c(power_KS10, power_Kuiper10, power_CvM10, power_Watson10, power_AD10)
