require(GoFKernel)    	 # For inverting a cdf
require(univariateML)  	 # For the Kumaraswamy MLE      
require(VGAM)         	 # For the Kumaraswamy distribution 
require(stats)           # For the beta distribution
require(trapezoid)  	 # For trapezoidal
require(truncnorm)       # For truncated Normal
require(cascsim)         # For truncated Gamma
require(EnvStats)        # For truncated Log-Normal & Triangular
require(ReIns)           # For truncated Weibull

# Reference for A-D test in case the of the Beta distribution: 
# Mathias Raschke (2011), "Empirical behaviour of tests for the beta distribution and their application 
# in environmental research", Stochastic Environmental Research and Risk Assessment, 25:79–89
# (See the reference therein to his earlier paper)
#--------------------------------------------------------------

set.seed (1234)   
nrep<- 10000
n<-  100

y<- c()
AD<- c()
KS<- c()
CvM<- c()
Watson<- c()
Kuiper<- c()

crit_KS_10<- 0.819  		# See Raschke, 2011, p.81
crit_Kuiper_10<- 1.386
crit_AD_10<-  0.631   
crit_Watson_10<-  0.096 
crit_CvM_10<- 0.104  
crit_KS_5<- 0.895  
crit_Kuiper_5<- 1.489
crit_AD_5<-  0.752   
crit_Watson_5<-  0.117 
crit_CvM_5<- 0.126 

alpha<- 0.5          # Only need values when null is true (size simulations)
beta<- 0.5

Dplus<- function(x,N,px) {
max(seq(1:N)/N - px)
}

Dminus<- function(x,N,px) {
max(px-(seq(1:N)-1)/N)
}

s<- seq(1:n)

#  Start of simulation loop:

for (i in 1:nrep) {

x<- rkumar(n, alpha, beta)                    # Null is true
#x<- rtri(n, min = 0, max = 1, mode = 7/8)     # Null is false
#x<- rbeta(n,0.5,0.5)                          # Null is false 
#x<- rtruncnorm(n, a=0, b=1, mean=0.2, sd=5)    # Null is false
#x<- rtrapezoid(n, min = 0, mode1 = 5/8, mode2 = 7/8, max = 1, n1 = 2, n3 = 2)  # Null is false
#x<- rtgamma(n, 2, 3, min = 0.0, max = 0.99)    # Null is false
#x<- rlnormTrunc(n, meanlog = 0, sdlog = 1, min = 0, max = 1)  # Null is false
x<- rtweibull(n, shape = 2, scale=1, endpoint=1)               # Null is false  

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
# Figure 1:
# --------

curve(dkumar(x, 2,20), col="blue", lwd=2,ylab="p.d.f.", ylim=c(0,4.5),main ="Figure 1: Kumaraswamy densities")
curve(dkumar(x, 1,3), col="black", lwd=2, lty=2, add=TRUE)
curve(dkumar(x,3,8), col="blue", lwd=2,lty=3,add=TRUE)
curve(dkumar(x,3,4), col="red", lty=3, lwd=2, add=TRUE)
curve(dkumar(x, 8,2), col="red", lwd=2, add=TRUE)
curve(dkumar(x, 2,2.5), col="black", lwd=2, add=TRUE)
curve(dkumar(x,0.5,0.5), col= "orange",lty=2, lwd=2, add=TRUE)
legend(0.4,4.1,legend = c(expression(paste(a, " = ", "2.0 ; ", b, " = ", "20.0")),
                           expression(paste(a, " = ", "1.0 ; ", b, " = ", "3.0")),
                           expression(paste(a, " = ", "3.0 ; ", b, " = ", "8.0")),
                           expression(paste(a, " = ", "3.0 ; ", b, " = ", "4.0")),
                           expression(paste(a, " = ", "8.0 ; ", b, " = ", "2.0")), 
                           expression(paste(a, " = ", "2.0 ; ", b, " = ", "2.5")),
				   expression(paste(a, " = ", "0.5 ; ", b, " = ", "0.5"))),
                           box.col="white", col=c("blue","black","blue", "red", "red","black", "orange"),
                           lty=c(1,2,3,3,1,1,2), ncol=1)
# -------------------------------------------------------------------------------
# Figure 2:
# --------

curve(dbeta(x, 3,3), col="blue", lwd=2,ylab="p.d.f.", ylim=c(0,6.5),main ="Figure 2: Beta densities")
curve(dbeta(x, 20,20), col="black", lwd=2, lty=2, add=TRUE)
curve(dbeta(x,4,2), col="blue", lwd=2,lty=3,add=TRUE)
curve(dbeta(x,2,4), col="red", lty=3, lwd=2, add=TRUE)
curve(dbeta(x, 3,20), col="red", lwd=2, add=TRUE)
curve(dbeta(x, 0.5,0.5), col="black", lwd=2, add=TRUE)
legend(0.65,6,legend = c(expression(paste(a, " = ", "3.0 ; ", b, " = ", "3.0")),
                           expression(paste(a, " = ", "20.0 ; ", b, " = ", "20.0")),
                           expression(paste(a, " = ", "4.0 ; ", b, " = ", "2.0")),
                           expression(paste(a, " = ", "2.0 ; ", b, " = ", "4.0")),
                           expression(paste(a, " = ", "3.0 ; ", b, " = ", "20.0")), 
                           expression(paste(a, " = ", "0.5 ; ", b, " = ", "0.5"))),
                           box.col="white", col=c("blue","black","blue", "red", "red","black"),
                           lty=c(1,2,3,3,1,1), ncol=1)
# -------------------------------------------------------------------------------


# Plots of pdf's of some alternative distributions

curve(dtri(x, min = 0, max = 1, mode = 1/4))   
curve(dtri(x, min = 0, max = 1, mode = 7/8))   

curve(dlnormTrunc(x, meanlog = 0, sdlog = 1, min = 0, max = 1))
curve(dlnormTrunc(x, meanlog = 0.5, sdlog = 0.5, min = 0, max = 1))

curve(dtrapezoid(x, min = 0, mode1 = 1/8, mode2 = 3/8, max = 1, n1 = 2, n3 = 2, alpha = 1))
curve(dtrapezoid(x, min = 0, mode1 = 5/8, mode2 = 7/8, max = 1, n1 = 2, n3 = 2, alpha = 1))

curve(dtruncnorm(x, a=0, b=1, mean=0.5 , sd=0.1))
curve(dtruncnorm(x, a=0, b=1, mean= 0.8 , sd=0.8))

curve(dbeta(x,shape1=3,shape2=3))
curve(dbeta(x,shape1=20,shape2=20))
curve(dbeta(x,shape1=2,shape2=4))
curve(dbeta(x,shape1=4,shape2=2))
curve(dbeta(x,shape1=3,shape2=20))
curve(dbeta(x,shape1=0.5,shape2=0.5))

curve(dtgamma(x, 2, 3, min = 0.0, max = 1))
curve(dtgamma(x, 2, 6, min = 0.0, max = 1))

curve(dtweibull(x, shape = 2, scale=1, endpoint=1))
