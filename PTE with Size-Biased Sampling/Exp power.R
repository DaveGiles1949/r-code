set.seed(1234)
require(stats)
require(ggamma)

N<- c(10,20,30,40,50,60,70,80,90,100,110, 120,130,140,150,160,170,180,190,200)
nsamp<- 100000	# Number of replications for generating critical values
m<- 20000		# Number of replications for simulating the biases and risks
lambda_exp<- c()
lambda_hn<- c()
lambda_ray<- c()
laml_exp<- c()
laml_hn<- c()
laml_ray<- c()
power_exp<- c()
power_hn<- c()
power_ray<- c()
theta<- 1.0   # theta is the rate parameter for Exponential
sigma<- 1.0   # sigma is the scale parameter for the half-normal and Rayleigh distributions

# Set up outer loop for different sample sizes:

for (j in 1:length(N)) {
n<-N[j]
power_exp[j]<-0
power_hn[j]<- 0
power_ray[j]<- 0

# Create the critical values - these depend on the value of n, and the distribution.  
# and the distribution

for (i in 1: nsamp) {

x1<- rexp(n,theta)				  # Exponential
lambda_exp[i]<- (prod(x1))^(1/n)/(mean(x1))   
a<- sqrt(2*sigma^2)
x2<- rggamma(n,a,2,0.5)                     # Half-normal
lambda_hn[i]<- (prod(x2))^(1/n)/(mean(x2))
x3<- rggamma(n,a,2,1)                       # Rayleigh
lambda_ray[i]<- (prod(x3))^(1/n)/(mean(x3))

}

crit_exp<- quantile(lambda_exp,  probs = c(0.90, 0.95, 0.99))
crit_hn<- quantile(lambda_hn,  probs = c(0.90, 0.95, 0.99))
crit_ray<-quantile(lambda_ray,  probs = c(0.90, 0.95, 0.99))

# The length-biased exponential is just the generalized Gamma density with a = 1/theta, b=1,k=2

k<-  2   # Exponential: k  = 1 implies regular distribution; k = 2 implies length-biased

for ( i in 1:m) {

y1<- rggamma(n,1/theta,1,k)
laml_exp[i]<- (prod(y1))^(1/n)/(mean(y1))
y2<- rggamma(n,a,2,1)   # length-biased H-N = Rayleigh
laml_hn[i]<- (prod(y2))^(1/n)/(mean(y2))
y3<- rggamma(n,a,2,1.5)    # Length-biased Rayleigh
laml_ray[i]<- (prod(y3))^(1/n)/(mean(y3))

}   # end of loop for one n

power_exp[j]<- sum(laml_exp>crit_exp[2])/m     # 5% significance level
power_hn[j]<- sum(laml_hn>crit_hn[2])/m
power_ray[j]<- sum(laml_ray>crit_ray[2])/m

}   # end of loop for all n

plot(N,power_exp,type="l",lty=2, lwd=2, ylim=c(0,1),xlab="Sample Size (n)", ylab="Power", main=" " )
lines(N,power_hn, type="l", lty=3, lwd=3, col="red")
lines(N,power_ray, type="l", lty=1, lwd=2, col="blue")
legend(150, 0.2,legend = c(expression("Exponential"),
                         expression("Half-normal"),
				 expression("Rayleigh")), 
                         box.col="white", col=c("black","red","blue"),lty=c(2,3,1),lwd=c(2,3,2),ncol=1)


