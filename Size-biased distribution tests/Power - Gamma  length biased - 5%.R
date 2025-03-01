# Gamma distribution - Length-biased sampling - Power
#####################################################
set.seed(1234)
require(stats)
require(ggamma)

N<- c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200)
m<- 20000		# Number of replications for computing the powers

laml_gam_1<- c()		# 5 values of shape parameter: b = 0.5, 1.0, 1.5, 2.0, 3.0
laml_gam_2<- c()
laml_gam_3<- c()
laml_gam_4<- c()
laml_gam_5<- c()
power_gam_1<- c()
power_gam_2<- c()
power_gam_3<- c()
power_gam_4<- c()
power_gam_5<- c()

theta<- 1.0
a<- 1/theta  # theta is the scale parameter for gamma distribution

# The results are invariant to the value of this parameter

# Read the critical values - these depend on the value of n, and the value of the shape parameter  
# and the distribution

# 5% critical values: numbering is fo k = 0.5, k = 1.0, etc. 

crit_gam_1<- c(0.6359,0.5113,0.4611,0.4327,0.4143,0.4003,0.3901,0.3817,0.3760,0.3700,0.3657,0.3614,0.3582,0.3548,0.3522,0.3499,0.3474,0.3456,0.3440,0.3419)
crit_gam_2<- c(0.8292,0.7439,0.7082,0.6877,0.6740,0.6636,0.6552,0.6493,0.6432,0.6393,0.6359,0.6326,0.6294,0.6267,0.6246,0.6219,0.6203,0.6186,0.6172,0.6159)			                                                                                                                      
crit_gam_3<- c(0.8822,0.8279,0.8025,0.7873,0.7772,0.7697,0.7639,0.7591,0.7552,0.7520,0.7491,0.7467,0.7444,0.7422,0.7404,0.7390,0.7373,0.7361,0.7349,0.7337)
crit_gam_4<- c(0.9124,0.8708,0.8521,0.8398,0.8319,0.8257,0.8211,0.8172,0.8144,0.8118,0.8094,0.8076,0.8057,0.8041,0.8026,0.8012,0.8001,0.7991,0.7980,0.7971)
crit_gam_5<- c(0.9427,0.9145,0.9013,0.8929,0.8871,0.8832,0.8797,0.8772,0.8749,0.8731,0.8716,0.8701,0.8688,0.8677,0.8667,0.8659,0.8651,0.8644,0.8638,0.8630)

# Set up outer loop for different sample sizes:

for (j in 1:length(N)) {
n<- N[j]

# The length-biased gamma is just the generalized Gamma density with a = 1/theta, b = 1,k = (b+1)

for ( i in 1:m) {

y1<- rggamma(n,a,1,1.5)  # length-biased gamma; b = 0.5
laml_gam_1[i]<- (prod(y1))^(1/n)/(mean(y1))
y2<- rggamma(n,a,1,2.0)  # length-biased gamma; b = 1.0 ; this is equal to exponential
laml_gam_2[i]<- (prod(y2))^(1/n)/mean(y2)
y3<- rggamma(n,a,1,2.5)  # length-biased gamma; b = 1.5
laml_gam_3[i]<- (prod(y3))^(1/n)/mean(y3)
y4<- rggamma(n,a,1,3.0)  # length-biased gamma; b = 2.0
laml_gam_4[i]<- (prod(y4))^(1/n)/mean(y4)
y5<- rggamma(n,a,1,4.0)  # length-biased gamma; b = 3.0
laml_gam_5[i]<- (prod(y5))^(1/n)/mean(y5)

}   # end of loop for one value of n

power_gam_1[j]<- sum(laml_gam_1>crit_gam_1[j])/m     # 5% significance level
power_gam_2[j]<- sum(laml_gam_2>crit_gam_2[j])/m 
power_gam_3[j]<- sum(laml_gam_3>crit_gam_3[j])/m
power_gam_4[j]<- sum(laml_gam_4>crit_gam_4[j])/m
power_gam_5[j]<- sum(laml_gam_5>crit_gam_5[j])/m

}   # end of loop for all n values

# Don't include exponential case (b = 1) in the plots

plot (N,power_gam_1, type="l", lty=1, lwd=2, ylim=c(0.1,1),xlab="Sample Size (n)", ylab="Power", main=expression(paste("Figure 2a: Powers against H1: c = 1 (", alpha, "= 5%)")))
axis(1, at = c(25, 50, 75,100,125,150,175,200))

lines(N,power_gam_3, type="l", lty=1, lwd=2, col="red")
lines(N,power_gam_4, type="l", lty=1, lwd=2, col="blue")
lines(N,power_gam_5, type="l", lty=1, lwd=2, col="green")

legend(100, 0.4,legend = c(expression(paste("Gamma (b = 0.5): ",hat(lambda)[1]," = ",tilde(lambda)[1])),
				expression(paste("Gamma (b = 1.5): ",hat(lambda)[1]," = ",tilde(lambda)[1])),
				expression(paste("Gamma (b = 2.0): ",hat(lambda)[1]," = ",tilde(lambda)[1])),
				expression(paste("Gamma (b = 3.0): ",hat(lambda)[1]," = ",tilde(lambda)[1]))),
 				box.col="white", col=c("black", "red", "blue","green"),lty=c(1,1,1,1),lwd=c(2,2,2,2),ncol=1)



 			


