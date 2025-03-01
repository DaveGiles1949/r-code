# Weibull Area-biased sampling - Powers
########################
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

sigma<- 1.0

# sigma is the scale parameter for Weibull distribution
# The results are invariant to the value of this parameter

# Read the critical values - these depend on the value of n, and the value of the shape parameter  
# and the null distribution

# 5% critical values:

crit_wei_1_hat<-   c(0.3546,0.2291,0.1875,0.1647,0.1502,0.1406,0.1332,0.1279,0.1235,0.1196)  # k = 0.5
crit_wei_1_tilde<- c(0.3386,0.2114,0.1707,0.1504,0.1372,0.1288,0.1223,0.1171,0.1134,0.1099)  # k = 0.5
crit_wei_2_hat<-   c(0.8462,0.7782,0.7487,0.7315,0.7194,0.7111,0.7041,0.6984,0.6942,0.6904)  # k = 1.5
crit_wei_2_tilde<- c(0.8458,0.7783,0.7491,0.7319,0.7200,0.7115,0.7045,0.6990,0.6945,0.6907)  # k = 1.5
crit_wei_3_hat<-   c(0.9067,0.8610,0.8411,0.8287,0.8206,0.8143,0.8092,0.8053,0.8021,0.7992)  # k = 2.0
crit_wei_3_tilde<- c(0.9053,0.8605,0.8404,0.8282,0.8201,0.8140,0.8090,0.8051,0.8018,0.7990)  # k = 2.0
crit_wei_4_hat<-   c(0.9558,0.9324,0.9214,0.9148,0.9101,0.9066,0.9037,0.9015,0.8996,0.8980)  # k = 3.0
crit_wei_4_tilde<- c(0.9545,0.9311,0.9204,0.9138,0.9093,0.9058,0.9030,0.9008,0.8989,0.8974)  # k = 3.0


# Set up outer loop for different sample sizes:

for (j in 1:length(N)) {
n<- N[j]

# The area-biased Weibull is just a generalized Gamma density
# We need to consider both lambda_hat and lambda_tilde

# Here is the likelihood equation for the concentrated log-likelihood, as a function of just "k" 
  
f<- function(x, y) {(1/x)+mean(log(y))-sum(y^x*log(y))/sum(y^x) } # 'x' is the unknown parameter, 'k'

for ( i in 1:m) {

y1<- rggamma(n,sigma,0.5,5)    # length-biased gamma; k = 0.5 ; last parameter = (k+2)/k
k_tilde<- uniroot(f,c(0.1, 20),y=y1 )$root
theta_tilde<- (mean(y1^k_tilde))^(1/k_tilde)
laml_wei_1_hat[i]<- (prod(y1))^(1/n)/sqrt(mean(y1^2))
laml_wei_1_tilde[i]<- (prod(y1))^(1/n)/(theta_tilde*sqrt(gamma(1+2/k_tilde)))

y2<- rggamma(n,theta,1.5,2.333)    # length-biased gamma; k = 1.5
k_tilde<- uniroot(f,c(0.1, 20),y=y2 )$root
theta_tilde<- (mean(y2^k_tilde))^(1/k_tilde)
laml_wei_2_hat[i]<- (prod(y2))^(1/n)/sqrt(mean(y2^2))
laml_wei_2_tilde[i]<- (prod(y2))^(1/n)/(theta_tilde*sqrt(gamma(1+2/k_tilde)))

y3<- rggamma(n,theta,2,2.0)    # length-biased gamma; k = 2.0
k_tilde<- uniroot(f,c(0.1, 20),y=y3 )$root
theta_tilde<- (mean(y3^k_tilde))^(1/k_tilde)
laml_wei_3_hat[i]<- (prod(y3))^(1/n)/sqrt(mean(y3^2))
laml_wei_3_tilde[i]<- (prod(y3))^(1/n)/(theta_tilde*sqrt(gamma(1+2/k_tilde)))

y4<- rggamma(n,theta,3,1.6667)    # length-biased gamma; k = 3.0
k_tilde<- uniroot(f,c(0.1, 50),y=y4 )$root
theta_tilde<- (mean(y4^k_tilde))^(1/k_tilde)
laml_wei_4_hat[i]<- (prod(y4))^(1/n)/sqrt(mean(y4^2))
laml_wei_4_tilde[i]<- (prod(y4))^(1/n)/(theta_tilde*sqrt(gamma(1+2/k_tilde)))

}   # end of loop for one value of n

power_wei_1_hat[j]<- sum(laml_wei_1_hat>crit_wei_1_hat[j])/m     # 5% significance level
power_wei_2_hat[j]<- sum(laml_wei_2_hat>crit_wei_2_hat[j])/m 
power_wei_3_hat[j]<- sum(laml_wei_3_hat>crit_wei_3_hat[j])/m
power_wei_4_hat[j]<- sum(laml_wei_4_hat>crit_wei_4_hat[j])/m
power_wei_1_tilde[j]<- sum(laml_wei_1_tilde>crit_wei_1_tilde[j])/m     # 5% significance level
power_wei_2_tilde[j]<- sum(laml_wei_2_tilde>crit_wei_2_tilde[j])/m 
power_wei_3_tilde[j]<- sum(laml_wei_3_tilde>crit_wei_3_tilde[j])/m
power_wei_4_tilde[j]<- sum(laml_wei_4_tilde>crit_wei_4_tilde[j])/m

}   # end of loop for all n values

plot (N,power_wei_1_hat, type="l", lty=1, lwd=2, ylim=c(0.2,1),xlab="Sample Size (n)", ylab="Power", main=expression(paste("Figure 3b: Powers against H1: c = 2 (", alpha, "= 5%)")))
#axis(1, at = c(25, 50, 75,100,125,150,175,200))
lines(N,power_wei_2_hat, type="l", lty=1, lwd=2, col="red")
lines(N,power_wei_3_hat, type="l", lty=1, lwd=2, col="blue")
lines(N,power_wei_4_hat, type="l", lty=1, lwd=2, col="green")
lines(N,power_wei_1_tilde, type="l", lty=2, lwd=2, col="black")
lines(N,power_wei_2_tilde, type="l", lty=2, lwd=2, col="red")
lines(N,power_wei_3_tilde, type="l", lty=2, lwd=2, col="blue")
lines(N,power_wei_4_tilde, type="l", lty=2, lwd=2, col="green")

legend(65, 0.5,legend = c(expression(paste("Weibull (k = 0.5): ",hat(lambda)[2])),
				expression(paste("Weibull (k = 0.5): ",tilde(lambda)[2])),
				expression(paste("Weibull (k = 1.5): ",hat(lambda)[2])),
				expression(paste("Weibull (k = 1.5): ",tilde(lambda)[2])),
				expression(paste("Weibull (k = 20): ",hat(lambda)[2])),
                        expression(paste("Weibull (k = 2.0): ",tilde(lambda)[2])),
 				expression(paste("Weibull (k = 3.0): ",hat(lambda)[2])),
				expression(paste("Weibull (k = 3.0): ",tilde(lambda)[2]))),
 				box.col="white", col=c("black","black", "red", "red", "blue","blue","green","green"),lty=c(1,2,1,2,1,2,1,2),lwd=c(2,2,2,2,2,2,2,2),ncol=1)


