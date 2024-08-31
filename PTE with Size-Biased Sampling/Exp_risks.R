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
theta_pte5<- c()
theta_pte_opt_mse<- c()
theta_pte_opt_mae<- c()
theta_pte_opt_LINEX<- c()

perc_mse_tilde<- c()
perc_mse_hat<- c()
perc_mse_pte5<- c()
perc_mse_pte_opt<- c()
true_pmse_tilde<- c()
true_pmse_hat<- c()
perc_mae_tilde<- c()
perc_mae_hat<- c()
perc_mae_pte5<- c()
perc_mae_pte_opt<- c()
perc_LINEX_risk_tilde<- c()
perc_LINEX_risk_hat<- c()
perc_LINEX_risk_pte5<- c()
perc_LINEX_risk_pte_opt<- c()
perc_bias_tilde<- c()
perc_bias_hat<- c()
perc_bias_pte5<- c()
perc_bias_pte_opt_mse<- c()
perc_bias_pte_opt_mae<- c()
perc_bias_pte_opt_LINEX<- c()

K<- c()
theta<- 1.0   # theta is the rate parameter for Exponential
a<- theta
b<- 1
k<- 1

w<- 2   # for the LINEX loss
# Create the critical values - these depend only on n 

for (i in 1: nsamp) {
x<- rggamma(n,a,b,k)				# Exponential with rate parameter theta
		
lambda[i]<- (prod(x))^(1/n)/mean(x)    # This is for the case of length-biased

}
#The following critical values are for n = 50
#############################################
crit<- quantile(lambda,  probs = c(0.95, 0.991, 0.987, 0.9965))
crit5<- crit[1]
crit_opt_mse<-    crit[2]
crit_opt_mae<-    crit[3]
crit_opt_LINEX<-   crit[4] # for w=2   



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

if (lam1[i]<=crit5) {		
	theta_pte5[i]<- 1/mean(y)
} 	else {
	theta_pte5[i]<- 2/mean(y)
} 
if (lam1[i]<=crit_opt_mse) {		
	theta_pte_opt_mse[i]<- 1/mean(y)
} 	else {
	theta_pte_opt_mse[i]<- 2/mean(y)
} 
if (lam1[i]<=crit_opt_mae) {		
	theta_pte_opt_mae[i]<- 1/mean(y)
} 	else {
	theta_pte_opt_mae[i]<- 2/mean(y)
} 
if (lam1[i]<=crit_opt_LINEX) {		
	theta_pte_opt_LINEX[i]<- 1/mean(y)
} 	else {
	theta_pte_opt_LINEX[i]<- 2/mean(y)
} 

}      # end of loop for a particular k

#pp<-c(sum(lam1>crit[1]), sum(lam1>crit[2]), sum(lam1>crit[3]))/m
#pp       # check that the significance levels are correct

perc_bias_tilde[j]<- 100*(mean(theta_tilde)-theta)/theta
perc_bias_hat[j]<- 100*(mean(theta_hat)-theta)/theta
perc_bias_pte5[j]<- 100*(mean(theta_pte5)-theta)/theta
perc_bias_pte_opt_mse[j]<- 100*(mean(theta_pte_opt_mse)-theta)/theta
perc_bias_pte_opt_mae[j]<- 100*(mean(theta_pte_opt_mae)-theta)/theta
perc_bias_pte_opt_LINEX[j]<- 100*(mean(theta_pte_opt_LINEX)-theta)/theta

perc_mse_tilde[j]<- 100*var(theta_tilde)/theta^2 + (perc_bias_tilde[j])^2/100
perc_mse_hat[j]<- 100*var(theta_hat)/theta^2 + (perc_bias_hat[j])^2/100
perc_mse_pte5[j]<- 100*var(theta_pte5)/theta^2 + (perc_bias_pte5[j])^2/100
perc_mse_pte_opt[j]<- 100*var(theta_pte_opt_mse)/theta^2 + (perc_bias_pte_opt_mse[j])^2/100

perc_mae_tilde[j]<- 100*(mean(abs(theta_tilde-theta)))/theta
perc_mae_hat[j]<- 100*(mean(abs(theta_hat-theta)))/theta
perc_mae_pte5[j]<- 100*(mean(abs(theta_pte5-theta)))/theta
perc_mae_pte_opt[j]<- 100*(mean(abs(theta_pte_opt_mae - theta)))/theta

perc_LINEX_risk_tilde[j]<- 100*mean(exp(w*(theta_tilde-theta))-w*(theta_tilde-theta)-1)/theta
perc_LINEX_risk_hat[j]<- 100*mean(exp(w*(theta_hat-theta))-w*(theta_hat-theta)-1)/theta
perc_LINEX_risk_pte5[j]<- 100*mean(exp(w*(theta_pte5-theta))-w*(theta_pte5-theta)-1)/theta
perc_LINEX_risk_pte_opt[j]<- 100*mean(exp(w*(theta_pte_opt_LINEX-theta))-w*(theta_pte_opt_LINEX-theta)-1)/theta
}      # end of loop for all k

perc_bias_tilde
#perc_bias_hat
#perc_bias_pte

#perc_mse_tilde
#perc_mse_hat
#perc_mse_pte

# The PTE risk is lower on the left, and higher on the right, as alpha decreases.
# This is the same as for the regression problem.
# As alpha decreases, Pr[Reject Ho when Ho is True] decreases. So, more likely to 
# incorporate tilde estimator into the PTE. So the PTE risk pivots to become closer to the "unweighted" estimator's
# at each of the extremes

par(mfrow=c(3,1))
plot(K,perc_mse_tilde,type="l",lwd=2,lty=2, ylim=c(0,25),xlab="c", ylab="%MSE", main=" ")
lines(K,perc_mse_hat, type="l", lwd=2,lty=5, col="red")
lines(K,perc_mse_pte_opt, type="l", lwd=2,lty=1, col="blue")
lines(K, perc_mse_pte5, type="l", lwd=3, lty=3, col="blue")
legend(0.78, 19.5,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("PTE opt)"),
                         expression ("PTE (5%)")), 
                         box.col="white", col=c("black","red","blue", "blue"),lwd=c(2,2,2,3),lty=c(2,5,1,3),ncol=1)
plot(K,perc_mae_tilde,type="l",lwd=2,lty=2, ylim=c(0,50),xlab="c", ylab="%MAE", main=" ")
lines(K,perc_mae_hat, type="l", lwd=2,lty=5, col="red")
lines(K,perc_mae_pte_opt, type="l", lwd=2, lty=1, col="blue")
lines(K,perc_mae_pte5, type="l", lwd=3, lty=3, col="blue")				
legend(0.74,40, legend = c(expression("Unweighted"),
				 expression("Length-Weighted"),
                         expression("PTE (opt)"), 							
                         expression ("PTE (5%)")), 
                         box.col="white", col=c("black","red","blue", "blue"),lwd=c(2,2,2,3),lty=c(2,5,1,3),ncol=1)
plot(K,perc_LINEX_risk_tilde,type="l",lwd=2,lty=2, ylim=c(0,75),xlab="c", ylab="%LINEX Risk", main=" ")
lines(K,perc_LINEX_risk_hat, type="l", lwd=2,lty=5, col="red")
lines(K,perc_LINEX_risk_pte_opt, type="l", lwd=2,lty=1, col="blue")
lines(K, perc_LINEX_risk_pte5, type="l", lwd=3,lty=3, col="blue")
legend(0.74, 75,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("PTE (opt)"),
                         expression ("PTE (5%)")), 
                         box.col="white", col=c("black","red","blue", "blue"),lwd=c(2,2,2,3),lty=c(2,5,1,3),ncol=1)