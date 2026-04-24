# Figure 4: Gamma distribution - Length-biased sampling
#####################################################
set.seed(1234)
require(moments)
require(ggamma)

m<- 50000		# Number of replications

laml0_gam_1<- c()     # 5 values of shape parameter: b = 0.5, 1.0, 1.5, 2.0, 3.0
laml0_gam_2<- c()
laml0_gam_3<- c()
laml0_gam_4<- c()
laml0_gam_5<- c()
laml1_gam_1<- c()    
laml1_gam_2<- c()
laml1_gam_3<- c()
laml1_gam_4<- c()
laml1_gam_5<- c()

power_gam_1<- c()
power_gam_2<- c()
power_gam_3<- c()
power_gam_4<- c()
power_gam_5<- c()

theta<- 1.0
a<- 1/theta  # theta is the scale parameter for gamma distribution

# The results are invariant to the value of this parameter

# The length-biased gamma is just the generalized Gamma density with a = 1/theta, b = 1,k = (b+1)

# generate the data, and simulate the distributions of the test statistic, under both the null hypothesis
# and the alternative hypothesis of length-biased sampling:

n<- 10

for ( i in 1:m) {

y10<- rgamma(n,0.5,rate=a) # Null hypothesis: gamma; b = 0.5
y11<- rggamma(n,a,1,1.5)  # Alternative hypothesis: length-biased gamma; b = 0.5; etc.
laml0_gam_1[i]<- (prod(y10))^(1/n)/(mean(y10))
laml1_gam_1[i]<- (prod(y11))^(1/n)/(mean(y11))

y20<- rgamma(n,1.0,rate=a)
y21<- rggamma(n,a,1,2.0)  # length-biased gamma; b = 1.0 
laml0_gam_2[i]<- (prod(y20))^(1/n)/mean(y20)
laml1_gam_2[i]<- (prod(y21))^(1/n)/mean(y21)

y30<- rgamma(n,1.5,rate=a)
y31<- rggamma(n,a,1,2.5)  # length-biased gamma; b = 1.5
laml0_gam_3[i]<- (prod(y30))^(1/n)/mean(y30)
laml1_gam_3[i]<- (prod(y31))^(1/n)/mean(y31)

y40<- rgamma(n,2.0,rate=a)
y41<- rggamma(n,a,1,3.0)  # length-biased gamma; b = 2.0
laml0_gam_4[i]<- (prod(y40))^(1/n)/mean(y40)
laml1_gam_4[i]<- (prod(y41))^(1/n)/mean(y41)

y50<- rgamma(n,3.0,rate=a)
y51<- rggamma(n,a,1,4.0)  # length-biased gamma; b = 3.0
laml0_gam_5[i]<- (prod(y50))^(1/n)/mean(y50)
laml1_gam_5[i]<- (prod(y51))^(1/n)/mean(y51)

}   # end of loop for one value of n


par(mfrow=c(5,2))
# the ablines mark the position of the 5% critical value (for n=10, and each value of the shape parameter)

# First,when the null is true:

hist(laml0_gam_1, prob=TRUE)
abline(v=0.6359, lwd=2, col="red")
hist(laml0_gam_2, prob=TRUE)
abline(v=0.8292, lwd=2, col="red")
hist(laml0_gam_3, prob=TRUE)
abline(v=0.8822, lwd=2, col="red")
hist(laml0_gam_4, prob=TRUE)
abline(v=0.9124, lwd=2,col="red")
hist(laml0_gam_5, prob=TRUE)
abline(v=0.9427, lwd=2, col="red")

mean(laml0_gam_1)
mean(laml0_gam_2)
mean(laml0_gam_3)
mean(laml0_gam_4)
mean(laml0_gam_5)

moments::skewness(laml0_gam_1)
moments::skewness(laml0_gam_2)
moments::skewness(laml0_gam_3)
moments::skewness(laml0_gam_4)
moments::skewness(laml0_gam_5)

moments::kurtosis(laml0_gam_1)
moments::kurtosis(laml0_gam_2)
moments::kurtosis(laml0_gam_3)
moments::kurtosis(laml0_gam_4)
moments::kurtosis(laml0_gam_5)

# next when the null is false:

hist(laml1_gam_1, prob=TRUE)
abline(v=0.6359, lwd=2, col="blue")
hist(laml1_gam_2, prob=TRUE)
abline(v=0.8292, lwd=2, col="blue")
hist(laml1_gam_3, prob=TRUE)
abline(v=0.8822, lwd=2, col="blue")
hist(laml1_gam_4, prob=TRUE)
abline(v=0.9124, lwd=2,col="blue")
hist(laml1_gam_5, prob=TRUE)
abline(v=0.9427, lwd=2, col="blue")

mean(laml1_gam_1)
mean(laml1_gam_2)
mean(laml1_gam_3)
mean(laml1_gam_4)
mean(laml1_gam_5)

moments::skewness(laml1_gam_1)
moments::skewness(laml1_gam_2)
moments::skewness(laml1_gam_3)
moments::skewness(laml1_gam_4)
moments::skewness(laml1_gam_5)

moments::kurtosis(laml1_gam_1)
moments::kurtosis(laml1_gam_2)
moments::kurtosis(laml1_gam_3)
moments::kurtosis(laml1_gam_4)
moments::kurtosis(laml1_gam_5)

par(mfrow=c(2,2))
hist(laml0_gam_1, prob=TRUE, main="Gamma distribution
(Shape parameter = 0.5)", xlab="Test statistic", cex.main=0.8,cex.lab=0.8,cex.axis=0.8)
abline(v=0.6359, lwd=3, col="red")
mtext("Figure 4: Test statistic densities (n = 10)", side = 3, line = -1, outer = TRUE)
hist(laml0_gam_5, prob=TRUE,main="Gamma Distribution
(Shape parameter = 3.0)", xlab="Test statistic", cex.main=0.8,cex.lab=0.8,cex.axis=0.8)
abline(v=0.9427, lwd=3, col="red")

hist(laml1_gam_1, prob=TRUE, main="Length-biased Gamma distribution
(Shape parameter = 0.5)", xlab="Test statistic", cex.main=0.8,cex.lab=0.8,cex.axis=0.8)
abline(v=0.6359, lwd=3, col="blue")
hist(laml1_gam_5, prob=TRUE,main="Length-biased Gamma distribution
(Shape parameter = 3.0)", xlab="Test statistic", cex.main=0.8,cex.lab=0.8,cex.axis=0.8)
abline(v=0.9427, lwd=3, col="blue")

################################################

# Now repeat with n = 50:

n<- 50

for ( i in 1:m) {

y10<- rgamma(n,0.5,rate=a) # Null hypothesis: gamma; b = 0.5
y11<- rggamma(n,a,1,1.5)  # Alternative hypothesis: length-biased gamma; b = 0.5; etc.
laml0_gam_1[i]<- (prod(y10))^(1/n)/(mean(y10))
laml1_gam_1[i]<- (prod(y11))^(1/n)/(mean(y11))

y20<- rgamma(n,1.0,rate=a)
y21<- rggamma(n,a,1,2.0)  # length-biased gamma; b = 1.0 
laml0_gam_2[i]<- (prod(y20))^(1/n)/mean(y20)
laml1_gam_2[i]<- (prod(y21))^(1/n)/mean(y21)

y30<- rgamma(n,1.5,rate=a)
y31<- rggamma(n,a,1,2.5)  # length-biased gamma; b = 1.5
laml0_gam_3[i]<- (prod(y30))^(1/n)/mean(y30)
laml1_gam_3[i]<- (prod(y31))^(1/n)/mean(y31)

y40<- rgamma(n,2.0,rate=a)
y41<- rggamma(n,a,1,3.0)  # length-biased gamma; b = 2.0
laml0_gam_4[i]<- (prod(y40))^(1/n)/mean(y40)
laml1_gam_4[i]<- (prod(y41))^(1/n)/mean(y41)

y50<- rgamma(n,3.0,rate=a)
y51<- rggamma(n,a,1,4.0)  # length-biased gamma; b = 3.0
laml0_gam_5[i]<- (prod(y50))^(1/n)/mean(y50)
laml1_gam_5[i]<- (prod(y51))^(1/n)/mean(y51)

}   # end of loop for one value of n


par(mfrow=c(5,2))
# the ablines mark the position of the 5% critical value (for n=10, and each value of the shape parameter)

# First,when the null is true:

hist(laml0_gam_1, prob=TRUE)
abline(v=0.4143, lwd=2, col="red")   # "ab-line" is to left of lower limit of x-axis
hist(laml0_gam_2, prob=TRUE)
abline(v=0.6740, lwd=2, col="red")
hist(laml0_gam_3, prob=TRUE)
abline(v=0.7772, lwd=2, col="red")
hist(laml0_gam_4, prob=TRUE)
abline(v=0.8319, lwd=2,col="red")
hist(laml0_gam_5, prob=TRUE)
abline(v=0.8871, lwd=2, col="red")

mean(laml0_gam_1)
mean(laml0_gam_2)
mean(laml0_gam_3)
mean(laml0_gam_4)
mean(laml0_gam_5)

moments::skewness(laml0_gam_1)
moments::skewness(laml0_gam_2)
moments::skewness(laml0_gam_3)
moments::skewness(laml0_gam_4)
moments::skewness(laml0_gam_5)

moments::kurtosis(laml0_gam_1)
moments::kurtosis(laml0_gam_2)
moments::kurtosis(laml0_gam_3)
moments::kurtosis(laml0_gam_4)
moments::kurtosis(laml0_gam_5)

# next when the null is false:

hist(laml1_gam_1, prob=TRUE)
abline(v=0.4143, lwd=2, col="blue")
hist(laml1_gam_2, prob=TRUE)
abline(v=0.6740, lwd=2, col="blue")
hist(laml1_gam_3, prob=TRUE)
abline(v=0.7772, lwd=2, col="blue")
hist(laml1_gam_4, prob=TRUE)
abline(v=0.8319, lwd=2,col="blue")
hist(laml1_gam_5, prob=TRUE)
abline(v=0.8871, lwd=2, col="blue")

mean(laml1_gam_1)
mean(laml1_gam_2)
mean(laml1_gam_3)
mean(laml1_gam_4)
mean(laml1_gam_5)

moments::skewness(laml1_gam_1)
moments::skewness(laml1_gam_2)
moments::skewness(laml1_gam_3)
moments::skewness(laml1_gam_4)
moments::skewness(laml1_gam_5)

moments::kurtosis(laml1_gam_1)
moments::kurtosis(laml1_gam_2)
moments::kurtosis(laml1_gam_3)
moments::kurtosis(laml1_gam_4)
moments::kurtosis(laml1_gam_5)

################################################

