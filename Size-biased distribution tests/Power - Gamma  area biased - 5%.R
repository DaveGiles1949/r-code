# Gamma: Area-biased sampling - Powers
#######################################
set.seed(1234)
require(stats)
require(ggamma)

N<- c(10,20,30,40,50,60,70,80,90,100)
m<- 20000		# Number of replications for computing the powers

laml_gam_1_hat<- c()		# 4 values of shape parameter: b = 0.5, 1.5, 2.0, 3.0
laml_gam_2_hat<- c()
laml_gam_3_hat<- c()
laml_gam_4_hat<- c()
laml_gam_1_tilde<- c()		
laml_gam_2_tilde<- c()
laml_gam_3_tilde<- c()
laml_gam_4_tilde<- c()
power_gam_1_hat<- c()
power_gam_2_hat<- c()
power_gam_3_hat<- c()
power_gam_4_hat<- c()
power_gam_1_tilde<- c()
power_gam_2_tilde<- c()
power_gam_3_tilde<- c()
power_gam_4_tilde<- c()

theta<- 1
a<- 1/theta  # theta is the rate parameter for gamma distribution

# The results are invariant to the value of this parameter

# Read the critical values - these depend on the value of n, and the value of the shape parameter  
# and the null distribution

# 5% critical values:  numbered according to k = 0.5, k = 1.5, etc.

crit_gam_1_hat<-   c(0.4921,0.3607,0.3128,0.2887,0.2722,0.2603,0.2502,0.2440,0.2390,0.2340)
crit_gam_2_hat<-   c(0.7957,0.7166,0.6817,0.6595,0.6459,0.6356,0.6277,0.6214,0.6163,0.6119)
crit_gam_3_hat<-   c(0.8427,0.7785,0.7501,0.7310,0.7194,0.7105,0.7039,0.6983,0.6942,0.6904)
crit_gam_4_hat<-   c(0.8942,0.8467,0.8247,0.8107,0.8018,0.7957,0.7896,0.7859,0.7820,0.7791)
crit_gam_1_tilde<- c(0.4738,0.3484,0.3030,0.2796,0.2649,0.2536,0.2448,0.2380,0.2339,0.2290)
crit_gam_2_tilde<- c(0.7910,0.7110,0.6758,0.6545,0.6412,0.6314,0.6237,0.6175,0.6126,0.6084)
crit_gam_3_tilde<- c(0.8393,0.7744,0.7457,0.7274,0.7162,0.7072,0.7012,0.6953,0.6915,0.6877)
crit_gam_4_tilde<- c(0.8926,0.8441,0.8224,0.8089,0.7999,0.7938,0.7880,0.7842,0.7804,0.7775)

f<- function(x, n, y, ybar) {n*log(x/ybar)+sum(log(y))-n*digamma(x) }   # concentrated log-Likelihood equation for gamma distribution
                                                                        # remaining parameter (called 'x') is the shape parameter ("b")
# Set up outer loop for different sample sizes:

for (j in 1:length(N)) {
n<- N[j]

# The area-biased gamma is just the generalized Gamma density with a = 1/theta, b = 1,k = (b+2)
# We need to consider both lambda_hat and lambda_tilde

for ( i in 1:m) {

y1<- rggamma(n,a,1,2.5)    # length-biased gamma; b = 0.5
y1bar<- mean(y1)
b_tilde<- uniroot(f, c(0.1, 100), n=n,y=y1,ybar=y1bar)$root
laml_gam_1_hat[i]<- (prod(y1))^(1/n)/sqrt(mean(y1^2))
laml_gam_1_tilde[i]<- (prod(y1))^(1/n)/(sqrt(1+1/b_tilde)*y1bar)

y2<- rggamma(n,a,1,3.5)    # length-biased gamma; b = 1.5
y2bar<- mean(y2)
b_tilde<- uniroot(f, c(1, 100), n=n,y=y2,ybar=y2bar)$root
laml_gam_2_hat[i]<- (prod(y2))^(1/n)/sqrt(mean(y2^2))
laml_gam_2_tilde[i]<- (prod(y2))^(1/n)/(sqrt(1+1/b_tilde)*y2bar)

y3<- rggamma(n,a,1,4.0)    # length-biased gamma; b = 2.0
y3bar<- mean(y3)
b_tilde<- uniroot(f, c(1, 100), n=n,y=y3,ybar=y3bar)$root
laml_gam_3_hat[i]<- (prod(y3))^(1/n)/sqrt(mean(y3^2))
laml_gam_3_tilde[i]<- (prod(y3))^(1/n)/(sqrt(1+1/b_tilde)*y3bar)

y4<- rggamma(n,a,1,5.0)    # length-biased gamma; b = 3.0
y4bar<- mean(y4)
b_tilde<- uniroot(f, c(1, 100), n=n,y=y4,ybar=y4bar)$root
laml_gam_4_hat[i]<- (prod(y4))^(1/n)/sqrt(mean(y4^2))
laml_gam_4_tilde[i]<- (prod(y4))^(1/n)/(sqrt(1+1/b_tilde)*y4bar)


}   # end of loop for one value of n

power_gam_1_hat[j]<- sum(laml_gam_1_hat>crit_gam_1_hat[j])/m     # 5% significance level
power_gam_2_hat[j]<- sum(laml_gam_2_hat>crit_gam_2_hat[j])/m 
power_gam_3_hat[j]<- sum(laml_gam_3_hat>crit_gam_3_hat[j])/m
power_gam_4_hat[j]<- sum(laml_gam_4_hat>crit_gam_4_hat[j])/m
power_gam_1_tilde[j]<- sum(laml_gam_1_tilde>crit_gam_1_tilde[j])/m     # 5% significance level
power_gam_2_tilde[j]<- sum(laml_gam_2_tilde>crit_gam_2_tilde[j])/m 
power_gam_3_tilde[j]<- sum(laml_gam_3_tilde>crit_gam_3_tilde[j])/m
power_gam_4_tilde[j]<- sum(laml_gam_4_tilde>crit_gam_4_tilde[j])/m

}   # end of loop for all n values

plot (N,power_gam_1_hat, type="l", lty=1, lwd=2, ylim=c(0.2,1),xlab="Sample Size (n)", ylab="Power", main=expression(paste("Figure 2b: Powers against H1: c = 2 (", alpha, "= 5%)")))
lines(N,power_gam_2_hat, type="l", lty=1, lwd=2, col="red")
lines(N,power_gam_3_hat, type="l", lty=1, lwd=2, col="blue")
lines(N,power_gam_4_hat, type="l", lty=1, lwd=2, col="green")
lines(N,power_gam_1_tilde, type="l", lty=2, lwd=2, col="black")
lines(N,power_gam_2_tilde, type="l", lty=2, lwd=2, col="red")
lines(N,power_gam_3_tilde, type="l", lty=2, lwd=2, col="blue")
lines(N,power_gam_4_tilde, type="l", lty=2, lwd=2, col="green")

legend(60, 0.6,legend = c(expression(paste("Gamma (b = 0.5): ",hat(lambda)[2])),
				  expression(paste("Gamma (b = 0.5): ",tilde(lambda)[2])),
				  expression(paste("Gamma (b = 1.5): ",hat(lambda)[2])),
				  expression(paste("Gamma (b = 1.5): ",tilde(lambda)[2])),
				  expression(paste("Gamma (b = 2.0): ",hat(lambda)[2])),
                          expression(paste("Gamma (b = 2.0): ",tilde(lambda)[2])),
 				  expression(paste("Gamma (b = 3.0): ",hat(lambda)[2])),
				  expression(paste("Gamma (b = 3.0): ",tilde(lambda)[2]))),
 				  box.col="white", col=c("black","black", "red", "red", "blue","blue","green","green"),lty=c(1,2,1,2,1,2,1,2),lwd=c(2,2,2,2,2,2,2,2),ncol=1)


