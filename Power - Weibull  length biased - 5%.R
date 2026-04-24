# Weibull: Length-biased sampling - Power
########################################

set.seed(1234)
require(stats)
require(ggamma)

N<- c(10,20,30,40,50,60,70,80,90,100)
m<- 20000		# Number of replications for computing the powers

laml_wei_1_hat<- c()		# 4 values of shape parameter: k = 0.5, 1.5, 2.0, 3.0
laml_wei_2_hat<- c()
laml_wei_3_hat<- c()
laml_wei_4_hat<- c()
laml_wei_1_tilde<- c()
laml_wei_2_tilde<- c()
laml_wei_3_tilde<- c()
laml_wei_4_tilde<- c()
power_wei_1_hat<- c()
power_wei_2_hat<- c()
power_wei_3_hat<- c()
power_wei_4_hat<- c()
power_wei_1_tilde<- c()
power_wei_2_tilde<- c()
power_wei_3_tilde<- c()
power_wei_4_tilde<- c()

sigma<- 1.0 # scale parameter. The results are invariant to the value of this parameter

# Read the critical values - these depend on the value of n, and the value of the shape parameter  
# and the distribution

# 5% critical values:  shape = k = 0.5, 1.5, 2.0, 3.0 (numbered 1, 2, 3, 4, in next few lines)

crit_wei_1_hat<-   c(0.5038,0.3657,0.3160,0.2894,0.2708,0.2593,0.2498,0.2421,0.2368,0.2317)
crit_wei_1_tilde<- c(0.5047,0.3662,0.3153,0.2886,0.2705,0.2591,0.2496,0.2420,0.2365,0.2315)
crit_wei_2_hat<-   c(0.9136,0.8694,0.8490,0.8369,0.8283,0.8219,0.8167,0.8127,0.8092,0.8064)
crit_wei_2_tilde<- c(0.9110,0.8670,0.8471,0.8350,0.8268,0.8205,0.8154,0.8114,0.8081,0.8053)
crit_wei_3_hat<-   c(0.9497,0.9224,0.9095,0.9018,0.8962,0.8919,0.8885,0.8858,0.8835,0.8815)
crit_wei_3_tilde<- c(0.9481,0.9208,0.9081,0.9003,0.8949,0.8907,0.8873,0.8847,0.8825,0.8806)
crit_wei_4_hat<-   c(0.9771,0.9640,0.9575,0.9537,0.9508,0.9487,0.9469,0.9455,0.9443,0.9432)
crit_wei_4_tilde<- c(0.9772,0.9640,0.9575,0.9535,0.9507,0.9485,0.9467,0.9454,0.9442,0.9431)

# Set up outer loop for different sample sizes:

for (j in 1:length(N)) {
n<- N[j]

# Here is the likelihood equation for the concentrated log-likelihood, as a function of just "k" ('x' denotes the parameter, 'k')   
f<- function(x, y) {(1/x)+mean(log(y))-sum(y^x*log(y))/sum(y^x) }

for ( i in 1:m) {

y1<- rggamma(n,sigma,0.5,3)  # length-biased Weibull; k = 0.5 ; last parameter = (k+1)/k
k_tilde<- uniroot(f,c(0.1, 40),y=y1 )$root
sigma_tilde<- (mean(y1^k_tilde))^(1/k_tilde)
laml_wei_1_hat[i]<- ((prod(y1))^(1/n))/(mean(y1))
laml_wei_1_tilde[i]<- (prod(y1))^(1/n)/(sigma_tilde*gamma(1+1/k_tilde))

y2<- rggamma(n,sigma,1.5,1.6667)  # length-biased weibull: k = 1.5
k_tilde<- uniroot(f,c(0.1, 40),y=y2 )$root
sigma_tilde<- (mean(y2^k_tilde))^(1/k_tilde)
laml_wei_2_hat[i]<- ((prod(y2))^(1/n))/(mean(y2))
laml_wei_2_tilde[i]<- (prod(y2))^(1/n)/(sigma_tilde*gamma(1+1/k_tilde))

y3<- rggamma(n,sigma,2,1.5)  # length-biased weibull; k = 2.0
k_tilde<- uniroot(f,c(0.1, 40),y=y3 )$root
sigma_tilde<- (mean(y3^k_tilde))^(1/k_tilde)
laml_wei_3_hat[i]<- ((prod(y3))^(1/n))/mean(y3)
laml_wei_3_tilde[i]<- (prod(y3))^(1/n)/(sigma_tilde*gamma(1+1/k_tilde))

y4<- rggamma(n,sigma,3,1.3333)  # length-biased weibull: k = 3.0
k_tilde<- uniroot(f,c(0.1, 40),y=y4 )$root
sigma_tilde<- (mean(y4^k_tilde))^(1/k_tilde)
laml_wei_4_hat[i]<- (prod(y4))^(1/n)/mean(y4)
laml_wei_4_tilde[i]<- ((prod(y4))^(1/n))/(sigma_tilde*gamma(1+1/k_tilde))

}   # end of loop for one value of n

power_wei_1_hat[j]<-   sum(laml_wei_1_hat>crit_wei_1_hat[j])/m     # 5% significance level
power_wei_1_tilde[j]<- sum(laml_wei_1_tilde>crit_wei_1_tilde[j])/m 
power_wei_2_hat[j]<-   sum(laml_wei_2_hat>crit_wei_2_hat[j])/m
power_wei_2_tilde[j]<- sum(laml_wei_2_tilde>crit_wei_2_tilde[j])/m
power_wei_3_hat[j]<-   sum(laml_wei_3_hat>crit_wei_3_hat[j])/m
power_wei_3_tilde[j]<- sum(laml_wei_3_tilde>crit_wei_3_tilde[j])/m
power_wei_4_hat[j]<-   sum(laml_wei_4_hat>crit_wei_4_hat[j])/m
power_wei_4_tilde[j]<- sum(laml_wei_4_tilde>crit_wei_4_tilde[j])/m


}   # end of loop for all n values

plot (N,power_wei_1_hat, type="l", lty=1, lwd=2, ylim=c(0,1),xlab="Sample Size (n)", ylab="Power", main=expression(paste("Figure 3a: Powers against H1: c = 1 (", alpha, "= 5%)")))
lines(N,power_wei_1_tilde, type="l", lty=2, lwd=2, col="black")
lines(N,power_wei_2_hat, type="l", lty=1, lwd=2, col="red")
lines(N,power_wei_2_tilde, type="l", lty=2, lwd=2, col="red")
lines(N,power_wei_3_hat, type="l", lty=1, lwd=2, col="blue")
lines(N,power_wei_3_tilde, type="l", lty=2, lwd=2, col="blue")
lines(N,power_wei_4_hat, type="l", lty=1, lwd=2, col="green")
lines(N,power_wei_4_tilde, type="l", lty=2, lwd=2, col="green")

legend(62, 0.45,legend = c(expression(paste("Weibull (k = 0.5): ",hat(lambda)[1])),
				expression(paste("Weibull (k = 0.5): ",tilde(lambda)[1])),
				expression(paste("Weibull (k = 1.5): ",hat(lambda)[1])),
                        expression(paste("Weibull (k = 1.5): ",tilde(lambda)[1])),
                        expression(paste("Weibull (k = 2.0): ",hat(lambda)[1])),
                        expression(paste("Weibull (k = 2.0): ",tilde(lambda)[1])),
				expression(paste("Weibull (k = 3.0): ",hat(lambda)[1])),
                        expression(paste("Weibull (k = 3.0): ",tilde(lambda)[1]))),
 				box.col="white", col=c("black","black","red","red","blue","blue","green","green"),lty=c(1,2,1,2,1,2,1,2),lwd=c(2,2,2,2,2,2,2,2),ncol=1)


