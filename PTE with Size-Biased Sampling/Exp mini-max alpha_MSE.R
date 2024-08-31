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
theta_pte<- c()
perc_mse_tilde<- c()
perc_mse_hat<- c()
perc_mse_pte<- c()
true_pmse_tilde<- c()
true_pmse_hat<- c()
perc_bias_tilde<- c()
perc_bias_hat<- c()
perc_bias_pte<- c()

K<- c()
diff<- c()

theta<- 1.0   # theta is the rate parameter for Exponential
a<- theta
b<- 1
k<- 1
# Create the critical values - these depend only on n 

for (i in 1: nsamp) {
x<- rggamma(n,a,b,k)				# Exponential with rate parameter theta
		
lambda[i]<- (prod(x))^(1/n)/mean(x)    # This is for the case of length-biased

}

crit<- quantile(lambda,  probs = c(0.95))
		# The critical values are specific to the Exponential distribution; and the sample size
# vary the first element of crit to minimize max(diff) at bottom of file
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#for (jj in 1 :3) {   
# this loop needs to be built to search over different alpha values
#}
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


if (lam1[i]<=crit[1]) {		
	theta_pte[i]<- 1/mean(y)
} 	else {
	theta_pte[i]<- 2/mean(y)
} 


}      # end of loop for a particular k

#pp<-c(sum(lam1>crit[1]), sum(lam1>crit[2]), sum(lam1>crit[3]))/m
#pp       # check that the significance levels are correct

perc_bias_tilde[j]<- 100*(mean(theta_tilde)-theta)/theta
perc_bias_hat[j]<- 100*(mean(theta_hat)-theta)/theta
perc_bias_pte[j]<- 100*(mean(theta_pte)-theta)/theta
perc_mse_tilde[j]<- 100*var(theta_tilde)/theta^2 + (perc_bias_tilde[j])^2/100
perc_mse_hat[j]<- 100*var(theta_hat)/theta^2 + (perc_bias_hat[j])^2/100
perc_mse_pte[j]<- 100*var(theta_pte)/theta^2 + (perc_bias_pte[j])^2/100

diff[j]<- (perc_mse_pte[j] - min(perc_mse_tilde[j], perc_mse_hat[j]))

}      # end of loop for all k

perc_bias_tilde
perc_bias_hat
perc_bias_pte

perc_mse_tilde
perc_mse_hat
perc_mse_pte

# The PTE risk is lower on the left, and higher on the right, as alpha decreases.
# This is the same as for the regression problem.
# As alpha decreases, Pr[Reject Ho when Ho is True] decreases. So, more likely to 
# incorporate tilde estimator into the PTE. So the PTE risk pivots to become closer to the "unweighted" estimator's
# at each of the extremes

plot(K,perc_mse_tilde,type="l",lty=2, ylim=c(0,25),xlab="c", ylab="% MSE", main="Exponential; %MSE 
n = 50")
lines(K,perc_mse_hat, type="l", lty=3, col="red")
lines(K,perc_mse_pte, type="l", lty=1, col="blue")
legend(0.65, 15,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lty=c(2,3,1),ncol=1)
max(diff)
crit
