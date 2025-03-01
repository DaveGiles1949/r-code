# Length-biased sampling - Exponential, Half-Normal & Rayleigh - Powers
#######################################################################
set.seed(1234)
require(stats)
require(ggamma)

N<- c(10,20,30,40,50, 60, 70 ,80, 90, 100, 110, 120, 130)
m<- 20000		# Number of replications for computing the powers

laml_exp_1<- c()		# 1 = sample mean; 2 = MLE of pop mean
laml_exp_3<- c()
laml_hn_1<- c()
laml_hn_2<- c()
laml_hn_3<- c()
laml_ray_1<- c()
laml_ray_2<- c()
power_exp_1<- c()
power_exp_3<- c()
power_hn_1<- c()
power_hn_2<- c()
power_hn_3<- c()
power_ray_1<- c()
power_ray_2<- c()

theta<- 1.0   # theta is the rate parameter for Exponential
sigma<- 1.0   # sigma is the scale parameter for the Half-normal and Rayleigh distributions
a<- sqrt(2*sigma^2)
# The results are invariant to the values of these parameters

# Read the critical values - these depend on the value of n, and the distribution.  
# and the distribution

# 5% critical values: 

crit_exp_1<- c(0.8202,0.7439,0.7082,0.6877,0.6740,0.6636,0.6552,0.6493,0.6432,0.6393,0.6359,0.6326,0.6294) 
crit_hn_1<- c(0.8815,0.8224,0.7951,0.7774,0.7660,0.7570,0.7500,0.7444,0.7395,0.7356,0.7326,0.7294,0.7271)
crit_hn_2<- c(0.9996,0.8963,0.8523,0.8251,0.8076,0.7941,0.7838,0.7760,0.7686,0.7631,0.7582,0.7539,0.7507)
crit_ray_1<- c(0.9497,0.9228,0.9098,0.9018,0.8961,0.8921,0.8888,0.8859,0.8837,0.8819,0.8801,0.8789,0.8776)
crit_ray_2<- c(1.0229,0.9723,0.9495,0.9357,0.9258,0.9193,0.9129,0.9089,0.9054,0.9024, 0.8994,0.8975,0.8954)


# Set up outer loop for different sample sizes:

for (j in 1:length(N)) {
n<-N[j]

# The length-biased exponential is just the generalized Gamma density with a = 1/theta, b=1,k=2

#k<-  2   # Exponential: k  = 1 implies regular distribution; k = 2 implies length-biased


for ( i in 1:m) {

y1<- rggamma(n,1/theta,1,2)  # length-biased exponential
laml_exp_1[i]<- (prod(y1))^(1/n)/(mean(y1))
y2<- rggamma(n,a,2,1)   # length-biased H-N = Rayleigh
laml_hn_1[i]<- (prod(y2))^(1/n)/(mean(y2))
laml_hn_2[i]<- (prod(y2))^(1/n)/sqrt(2*sum(y2^2)/(n*pi))
y3<- rggamma(n,a,2,1.5)    # Length-biased Rayleigh
laml_ray_1[i]<- (prod(y3))^(1/n)/(mean(y3))
laml_ray_2[i]<- 2*(prod(y3))^(1/n)/sqrt(pi*sum(y3^2)/n)

}   # end of loop for one n

power_exp_1[j]<- sum(laml_exp_1>crit_exp_1[j])/m     # 5% significance level
power_hn_1[j]<- sum(laml_hn_1>crit_hn_1[j])/m
power_hn_2[j]<- sum(laml_hn_2>crit_hn_2[j])/m
power_ray_1[j]<- sum(laml_ray_1>crit_ray_1[j])/m
power_ray_2[j]<- sum(laml_ray_2>crit_ray_2[j])/m

}   # end of loop for all n

plot(N,power_exp_1,type="l",lty=1, lwd=2, ylim=c(0.2,1),xlab="Sample Size (n)", ylab="Power", main=expression( paste("Figure 1a: Powers against H1: c = 1 (", alpha, "= 5%)")))
lines(N,power_hn_1, type="l", lty=1, lwd=2, col="red")
lines(N,power_hn_2, type="l", lty=2, lwd=2, col="red")
lines(N,power_ray_1, type="l", lty=1, lwd=2, col="blue")
lines(N,power_ray_2, type="l", lty=2, lwd=2, col="blue")

legend(85, 0.5,legend = c(expression(paste("Exponential: ",hat(lambda)[1]," = ",tilde(lambda)[1])),
				expression(paste("Half-Normal: ",hat(lambda)[1])),
                         expression(paste("Half-Normal: ",tilde(lambda)[1])),
				expression(paste("Rayleigh: ",hat(lambda)[1])),
 				expression(paste("Rayleigh: ",tilde(lambda)[1]))),
				  box.col="white", col=c("black","red", "red", "blue","blue"),lty=c(1,1,2,1,2),lwd=c(2,2,2,2,2),ncol=1)


