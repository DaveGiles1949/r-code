set.seed(1234)
require(stats)
require(ggamma)
n<- 50   # sample size
nsamp<- 5000	# Number of replications for generating critical values
m<- 5000		# Number of replications for simulating the biases and risks
lambda<- c()
lam1<- c()
theta_tilde<- c()
theta_hat<- c()
theta_pte<- c()
perc_bias_tilde<- c()
perc_bias_hat<- c()
perc_bias_pte<- c()
perc_mse_tilde<- c()
perc_mse_hat<- c()
perc_mse_pte<- c()
true_pbias_tilde<- c()
true_pbias_hat<- c()
true_pmse_tilde<- c()
true_pmse_hat<- c()
C<- c()

theta<- 1.0   # theta is the rate parameterfor Exponential

# Create the critical values - these depend on n and the particular parameter values

# NEED TO SET UP MLE'S FOR THE GAMMA PARAMETERS, AND SIMILARLY FOR THE OTHER DISTRIBUTIONS
# ****************************************************************************************

for (i in 1: nsamp) {
x<- rexp(n,theta)				# Exponential with rate parameter theta
		
lambda[i]<- prod(x)/mean(x)^n    # This is for the case of length-biased

}

crit<- quantile(lambda,  probs = c(0.90, 0.95, 0.99))
crit		# The critical values are specific to the Exponential distribution; the value of the parameter; and the sample size

# Set up outer loop for different bias sizes:


for (j in 1:21) {
k<- 1+0.05*(j-1)
C[j]<- k
  # Exponential: k  = 1 implies regular distribution; k = 2 implies length-biased; k = 3 implies area-biased

for ( i in 1:m) {

y<- rggamma(n,1/theta,1,k)		# Exponential
lam1[i]<- prod(y)/mean(y)^n
theta_tilde[i]<- 1/mean(y)
theta_hat[i]<- 2/mean(y)

if (lam1[i]<=crit[3]) {		# Alter element from [1], [2], [3] for different sig. levels (10%, 5%, 1%)											
	theta_pte[i]<- 1/mean(y)
} 	else {
	theta_pte[i]<- 2/mean(y)
} 

}      # end of loop for a particular k

#pp<-c(sum(lam1>crit[1]), sum(lam1>crit[2]), sum(lam1>crit[3]))/m
#pp
perc_bias_tilde[j]<- 100*(mean(theta_tilde)-theta)/theta
perc_bias_hat[j]<- 100*(mean(theta_hat)-theta)/theta
perc_bias_pte[j]<- 100*(mean(theta_pte)-theta)/theta
perc_mse_tilde[j]<- 100*var(theta_tilde)/theta^2 + (perc_bias_tilde[j])^2/100
perc_mse_hat[j]<- 100*var(theta_hat)/theta^2 + (perc_bias_hat[j])^2/100
perc_mse_pte[j]<- 100*var(theta_pte)/theta^2 + (perc_bias_pte[j])^2/100

}      # end of loop for all k


perc_bias_tilde
perc_bias_hat
perc_bias_pte

perc_mse_tilde
perc_mse_hat
perc_mse_pte

# The PTE risk is lower on the left, and higher on the right, as alpha decreases.
# This is the same aas for the regression problem.
# As alpha decreases, Pr[Reject Ho when Ho is True] decreases. So, more likely to 
# incorporate tilde estimator into the PTE. SO the PTE risk pivots to become closer to the "unweighted" estimator's
# at each of the extremes

plot(C,perc_mse_tilde,type="l",lty=2, ylim=c(0,130),xlab="c", ylab="% MSE", main="Exponential; % MSE; alpha = 1%")
lines(C,perc_mse_hat, type="l", lty=1, col="red")
lines(C,perc_mse_pte, type="l", lty=1, col="blue")
legend(1.6, 100,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lty=c(2,1,1),ncol=1)

