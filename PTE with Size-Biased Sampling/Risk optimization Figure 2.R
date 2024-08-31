set.seed(1234)
require(stats)
require(ggamma)
n<- 50  # sample size
nsamp<- 100000	# Number of replications for generating critical values
m<- 20000		# Number of replications for simulating the biases and risks
lambda<- c()
lam1<- c()

theta_tilde<- c()
theta_hat<- c()
theta_pte10<- c()
theta_pte1<- c()

perc_mse_tilde<- c()
perc_mse_hat<- c()
perc_mse_pte10<- c()
perc_mse_pte1<- c()
true_pmse_tilde<- c()
true_pmse_hat<- c()
perc_bias_tilde<- c()
perc_bias_hat<- c()
perc_bias_pte10<- c()
perc_bias_pte1<- c()

K<- c()
theta<- 1.0   # theta is the rate parameter for Exponential
a<- theta
b<- 1
k<- 1

# Create the critical values - these depend only on n 

for (i in 1: nsamp) {
x<- rggamma(n,a,b,k)				# Exponential with rate parameter theta
		
lambda[i]<- (prod(x))^(1/n)/mean(x)    # This is for the case of length-biased

}

#############################################
crit<- quantile(lambda,  probs = c(0.999, 0.9, 0.987, 0.9965))
crit1<- crit[1]
crit10<- crit[2]

# Set up outer loop for different bias sizes:


for (j in 1:21) {
k<- 1+0.05*(j-1)
K[j]<- k-1
  # Exponential: k  = 1 implies regular distribution; k = 2 implies length-biased; k = 3 implies area-biased

for ( i in 1:m) {

y<- rggamma(n,a,b,k)		
lam1[i]<- (prod(y))^(1/n)/mean(y)
theta_tilde[i]<- 1/mean(y)
theta_hat[i]<- 2/mean(y)

if (lam1[i]<=crit10) {		
	theta_pte10[i]<- 1/mean(y)
} 	else {
	theta_pte10[i]<- 2/mean(y)
} 
if (lam1[i]<=crit1) {		
	theta_pte1[i]<- 1/mean(y)
} 	else {
	theta_pte1[i]<- 2/mean(y)
} 

}      # end of loop for a particular k

perc_bias_tilde[j]<- 100*(mean(theta_tilde)-theta)/theta
perc_bias_hat[j]<- 100*(mean(theta_hat)-theta)/theta
perc_bias_pte10[j]<- 100*(mean(theta_pte10)-theta)/theta
perc_bias_pte1[j]<- 100*(mean(theta_pte1)-theta)/theta


perc_mse_tilde[j]<- 100*var(theta_tilde)/theta^2 + (perc_bias_tilde[j])^2/100
perc_mse_hat[j]<- 100*var(theta_hat)/theta^2 + (perc_bias_hat[j])^2/100
perc_mse_pte10[j]<- 100*var(theta_pte10)/theta^2 + (perc_bias_pte10[j])^2/100
perc_mse_pte1[j]<- 100*var(theta_pte1)/theta^2 + (perc_bias_pte1[j])^2/100

}      # end of loop for all k

# The PTE risk is lower on the left, and higher on the right, as alpha decreases.
# This is the same as for the regression problem.
# As alpha decreases, Pr[Reject Ho when Ho is True] decreases. So, more likely to 
# incorporate tilde estimator into the PTE. So the PTE risk pivots to become closer to the "unweighted" estimator's
# at each of the extremes

par(mfrow=c(1,1))
plot(K,perc_mse_tilde,type="l",lty=4, lwd=2, ylim=c(0,25),xlab="c", ylab="%MSE", main=" ")
lines(K,perc_mse_hat, type="l", lty=5, lwd=2,col="red")
lines(K,perc_mse_pte1, type="l", lty=1, lwd=2,col="blue")
lines(K, perc_mse_pte10, type="l", lty=3, lwd=3, col="blue")
legend(0.5,25,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("PTE (0.1%)"),
                         expression ("PTE (10%)")), 
                         box.col="white", col=c("black","red","blue", "blue"),lwd=c(2,2,2,3),lty=c(4,5,1,3),ncol=1)
text(x=0.55, y=-0.5, "c*")