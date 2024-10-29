require(GoFKernel)    	# For inverting a cdf
require(AcceptReject)	# For Akash distribution
require(VGAM)		# For Lindley distribution
require(stats)		# For Weibull and log-normal distributions
require (nakagami)	# For the Nakagami distribution
require(fdrtool)		# For the half-normal distribution

# Reference for A-D test in the case of the Beta distribution: 
# Mathias Raschke (2011), "Empirical behaviour of tests for the beta distribution and their application 
# in environmental research", Stochastic Environmental Research and Risk Assessment, 25:79–89
# (See the reference therein to his earlier paper)
#--------------------------------------------------------------

set.seed (1234)   

nrep<- 10000
n<- 10
#lambda<- 3.0      # Need this when null is true (size simulations)

y<- c()
AD<- c()
KS<- c()
CvM<- c()
Watson<- c()
Kuiper<- c()

crit_KS_10<- 0.819    # These are from Stephens (1986). See Raschke, 2011, p.81
crit_Kuiper_10<- 1.386
crit_AD_10<-  0.631   
crit_Watson_10<-  0.096 
crit_CvM_10<- 0.104  
crit_KS_5<- 0.895  
crit_Kuiper_5<- 1.489
crit_AD_5<-  0.752   
crit_Watson_5<-  0.117 
crit_CvM_5<- 0.126 

Dplus<- function(x,N,px) {
max(seq(1:N)/N - px)
}

Dminus<- function(x,N,px) {
max(px-(seq(1:N)-1)/N)
}

Akash<- function (lambda,x) {
pdf<- (lambda^3/(lambda^2+2))*(1+x^2)*exp(-lambda*x)
pdf
}

pakash<- function(lambda, x) {
pak<- 1-(1+lambda*x*(lambda*x+2)/(lambda^2+2))*exp(-lambda*x)
pak
}

leq<- function(lambda,xbar){
le<- lambda^3*xbar - lambda^2 + 2*xbar*lambda - 6
le
}

f <- function(x) pnorm(x, mean=0, sd=1) 
f.inv <- inverse(f,lower=-25,upper=25) 

#  Start of simulation loop:

for (i in 1:nrep) {

# Akash distribution: (Null is True)
#yy <- accept_reject(
#n = n,
#f = Akash,
#args_f=list(lambda=lambda),
#continuous = TRUE,
#xlim = c(0, 200)
#) 

#x<- as.vector(yy)
# End of Akash distribution generation

# Alternative Hypotheses:

#x<- rexp(n, rate=0.5)  				# Exponential distribution
#x<- rlind (n, theta=2.0)			# Lindley distribution
#x<- rweibull(n, shape=1.5, scale=1)	# Weibull distribution
#x<- rgamma(n, shape= 0.5,rate=1.0)  	# Gamma distribution
#x<- rnaka(n, shape = 1.0 , scale = 1.0 ) # Nakagami distribution
#x<- rhalfnorm(n, theta=3/2)			#Half-normal distribution	Mean = 1/theta
x<- rlnorm(n, meanlog = 1, sdlog = 2)	# Log-normal distribution

xbar<- mean(x)

lam_tilde<- uniroot(leq, c(0.0001,250), xbar=xbar)$root

w<- pakash(lam_tilde,x)                  # Start of Raschke's testing procedure
     
for (j in 1:n){
y[j]<- f.inv(w[j])                             

}

mu<- mean(y)
sigma<- sd(y)*sqrt((n-1)/n)
ys<- sort(y)
pys<- pnorm(ys, mean=mu, sd=sigma)
pysr<-  rev(pys)    # reverse the order of ths pdf
seqn<- seq(1:n)

Dp<- Dplus(ys,n, pys)
Dm<- Dminus(ys,n, pys)
D<- max(Dp,Dm)
KS[i]<- Dstar<- D*(sqrt(n)-0.01+0.85/sqrt(n))

V<- Dp+Dm
Kuiper[i]<- Vstar<- V*(sqrt(n)+0.05+0.82/sqrt(n))

W2<- sum((pys -(2*seqn-1)/(2*n) )^2) + 1/(12*n)
CvM[i]<- W2star<- W2*(1+0.5/n)

U2<- W2-n*(sum(pys/n) - 0.5)^2   
Watson[i]<-U2star<- U2*(1+0.5/n)

A2<- -n-sum( (2*seqn-1)*(log(pys)+log(1-pysr) ) )/n
AD[i]<- A2star<- A2*(1+0.75/n+2.25/(n^2))

}

# End of simulation loop

power_KS5<- sum(KS>crit_KS_5)/nrep
power_Kuiper5<- sum(Kuiper>crit_Kuiper_5)/nrep
power_CvM5<- sum(CvM>crit_CvM_5)/nrep
power_Watson5<- sum(Watson>crit_Watson_5)/nrep
power_AD5<- sum(AD>crit_AD_5)/nrep
power_KS10<- sum(KS>crit_KS_10)/nrep
power_Kuiper10<- sum(Kuiper>crit_Kuiper_10)/nrep
power_CvM10<- sum(CvM>crit_CvM_10)/nrep
power_Watson10<- sum(Watson>crit_Watson_10)/nrep
power_AD10<- sum(AD>crit_AD_10)/nrep

c(power_KS5, power_Kuiper5, power_CvM5, power_Watson5, power_AD5)
c(power_KS10, power_Kuiper10, power_CvM10, power_Watson10, power_AD10)
#------------------------------------------------------------------------------  