# Area-biased sampling - Exponential, Half-normal and Rayleigh - Powers
#######################################################################
set.seed(1234)
require(stats)
require(ggamma)

N<- c(10,15,20,25,30,35,40,45,50)
m<- 20000		# Number of replications for computing the powers

laml_exp_1<- c()		# 1 = sample mean; 2 = MLE of pop mean
laml_exp_2<- c()
laml_hn_1<- c()
laml_ray_1<- c()
power_exp_1<- c()
power_exp_2<- c()
power_hn_1<- c()
power_ray_1<- c()

theta<- 1.0   # theta is the rate parameter for Exponential
sigma<- 1.0   # sigma is the scale parameter for the Half-normal and Rayleigh distributions
a<- sqrt(2*sigma^2)
# The results are invariant to the values of these parameters

# Read the critical values - these depend on the value of n, and the distribution.  
# and the distribution

# 5% critical values:

# There are 2 versions of the exponential, but only one of the H-H and Rayleigh
crit_exp_1<- c(0.7076,0.6426,0.6074,0.5814,0.5625,0.5504,0.5382,0.5297,0.5101)
crit_exp_2<- c(0.5799,0.5456,0.5260,0.5109,0.5007,0.4933,0.4863,0.4811,0.4692)
crit_hn_1<-  c(0.7975,0.7468,0.7152,0.6945,0.6800,0.6671,0.6584,0.6508,0.6336)
crit_ray_1<- c(0.9065,0.8786,0.8617,0.8506,0.8415,0.8346,0.8292,0.8242,0.8147 )


# Set up outer loop for different sample sizes:

for (j in 1:length(N)) {
n<-N[j]

# The area-biased exponential is just the generalized Gamma density with a = 1/theta, b=1,k=3
# Exponential: k  = 1 implies regular distribution; k = 3 implies area-biased


for ( i in 1:m) {

y1<- rggamma(n,1/theta,1,3)  # area-biased exponential
laml_exp_1[i]<- (prod(y1))^(1/n)/sqrt(mean(y1^2))
laml_exp_2[i]<- (prod(y1))^(1/n)/(sqrt(2)*mean(y1))
y2<- rggamma(n,a,2,1.5)   # area-biased H-N
laml_hn_1[i]<- (prod(y2))^(1/n)/sqrt(mean(y2^2))
y3<- rggamma(n,a,2,2)    # area-biased Rayleigh
laml_ray_1[i]<- (prod(y3))^(1/n)/sqrt(mean(y3^2))

}   # end of loop for one n

power_exp_1[j]<- sum(laml_exp_1>crit_exp_1[j])/m     # 5% significance level
power_exp_2[j]<- sum(laml_exp_2>crit_exp_2[j])/m 
power_hn_1[j]<- sum(laml_hn_1>crit_hn_1[j])/m
power_ray_1[j]<- sum(laml_ray_1>crit_ray_1[j])/m

}   # end of loop for all n

plot(N,power_exp_1,type="l",lty=1, lwd=2, ylim=c(0.35,1),xlab="Sample Size (n)", ylab="Power", main=expression( paste("Figure 1b: Powers against H1: c = 2 (", alpha, "= 5%)")))
lines(N,power_exp_2, type="l", lty=2, lwd=2, col="black")
lines(N,power_hn_1, type="l", lty=1, lwd=2, col="red")
lines(N,power_ray_1, type="l", lty=1, lwd=2, col="blue")
legend(35, 0.55,legend = c(expression(paste("Exponential: ",hat(lambda)[2])),
				expression(paste("Exponential: ",tilde(lambda)[2])),
                        expression(paste("Half-Normal: ",hat(lambda)[2]," = ",tilde(lambda)[2])),
 				expression(paste("Rayleigh: ",hat(lambda)[2]," = ",tilde(lambda)[2]))),
           			  box.col="white", col=c("black","black", "red", "blue"),lty=c(1,2,1,1),lwd=c(2,2,2,2),ncol=1)


