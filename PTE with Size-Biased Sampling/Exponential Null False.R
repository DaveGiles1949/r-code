set.seed(1234)
require(stats)
require(ggamma)

N<- c(10,15,20,25,30,35, 40,45,50,55,60, 65,70,75,80,85,90,95,100)
nsamp<- 100000	# Number of replications for generating critical values
m<- 20000		# Number of replications for simulating the biases and risks
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
perc_abs_risk_tilde<- c()
perc_abs_risk_hat<- c()
perc_abs_risk_pte<- c()
perc_LINEX_risk_tilde<- c()
perc_LINEX_risk_hat<- c()
perc_LINEX_risk_pte<- c()

true_pbias_tilde<- c()
true_pbias_hat<- c()
true_pmse_tilde<- c()
true_pmse_hat<- c()

theta<- 1.0   # theta is the rate parameterfor Exponential

# Set up outer loop for different sample sizes:
for (j in 1:19) {
n<-N[j]

# Create the critical values - these depend on the value of n 

# NEED TO SET UP MLE'S FOR THE GAMMA PARAMETERS, AND SIMILARLY FOR THE OTHER DISTRIBUTIONS
# ****************************************************************************************

for (i in 1: nsamp) {
x<- rexp(n,theta)				# Exponential with rate parameter theta
		
lambda[i]<- (prod(x))^(1/n)/(mean(x))    # This is for the case of length-biased

}

crit<- quantile(lambda,  probs = c(0.90, 0.95, 0.99))
crit		# The critical values are specific to the Exponential distribution; the value of the parameter; and the sample size

# The length-biased exponential is just the generalized Gamma density with a = 1/theta, b=1,k=2
# For area-biased, set a = 1/theta, b=1 and k=3
# (using the notation in the ggamma package)


k<-  2   # Exponential: k  = 1 implies regular distribution; k = 2 implies length-biased; k = 3 implies area-biased
for ( i in 1:m) {

y<- rggamma(n,1/theta,1,k)		# Exponential
lam1[i]<- (prod(y))^(1/n)/(mean(y))
theta_tilde[i]<- 1/mean(y)
theta_hat[i]<- 2/mean(y)

if (lam1[i]<=crit[2]) {		# Alter element from [1], [2], [3] for different sig. levels (10%, 5%, 1%)											
	theta_pte[i]<- 1/mean(y)
} 	else {
	theta_pte[i]<- 2/mean(y)
} 

}      # end of loop for a particular n

#pp<-c(sum(lam1>crit[1]), sum(lam1>crit[2]), sum(lam1>crit[3]))/m
#pp
perc_bias_tilde[j]<- 100*(mean(theta_tilde)-theta)/theta
true_pbias_tilde[j]<- -100*(N[j]-1)/(2*N[j]-1)
perc_bias_hat[j]<- 100*(mean(theta_hat)-theta)/theta
true_pbias_hat[j]<- 100/(2*N[j]-1)
perc_bias_pte[j]<- 100*(mean(theta_pte)-theta)/theta
perc_mse_tilde[j]<- 100*var(theta_tilde)/theta^2 + (perc_bias_tilde[j])^2/100
true_pmse_tilde[j]<- 100*(2*N[j]^3-5*N[j]^2+6*N[j]-2)/(2*(N[j]-1)*(2*N[j]-1)^2)
perc_mse_hat[j]<- 100*var(theta_hat)/theta^2 + (perc_bias_hat[j])^2/100
true_pmse_hat[j]<- 100*(N[j]+1)/((N[j]-1)*(2*N[j]-2))
perc_mse_pte[j]<- 100*var(theta_pte)/theta^2 + (perc_bias_pte[j])^2/100

# Consider other loss functions, apart from quadratic

# Absolute error loss:
# -------------------

perc_abs_risk_tilde[j]<- 100*mean(abs(theta_tilde-theta))/theta
perc_abs_risk_hat[j]<- 100*mean(abs(theta_hat-theta))/theta
perc_abs_risk_pte[j]<- 100*mean(abs(theta_pte-theta))/theta

# LINEX loss:
# -----------

w<- 2.0	#  Shape parameter for LINEX loss function (+ve or -ve)
# w>0 implies over-estimation more costly than under-estimation
# w<0 implies under-estimation more costly than over-estimation

perc_LINEX_risk_tilde[j]<- 100*mean(exp(w*(theta_tilde-theta))-w*(theta_tilde-theta)-1)/theta
perc_LINEX_risk_hat[j]<- 100*mean(exp(w*(theta_hat-theta))-w*(theta_hat-theta)-1)/theta
perc_LINEX_risk_pte[j]<- 100*mean(exp(w*(theta_pte-theta))-w*(theta_pte-theta)-1)/theta

}      # end of loop for all n


perc_bias_tilde
true_pbias_tilde
perc_bias_hat
true_pbias_hat
perc_bias_pte

perc_mse_tilde
true_pmse_tilde 
perc_mse_hat
true_pmse_hat
perc_mse_pte

perc_abs_risk_tilde
perc_abs_risk_hat
perc_abs_risk_pte

perc_LINEX_risk_tilde
perc_LINEX_risk_hat
perc_LINEX_risk_pte


par(mfrow=c(2,2))
plot(N,true_pbias_tilde,type="l",lwd=2,lty=2, ylim=c(-50,10),xlab="Sample Size (n)", ylab="% Bias", main=" ")
lines(N,true_pbias_hat, type="l", lwd=3,lty=3, col="red")
lines(N,perc_bias_pte, type="l", lwd=2,lty=1, col="blue")
legend(35, -20,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lwd=c(2,3,2),lty=c(2,3,1),ncol=1)

plot(N,true_pmse_tilde,type="l",lwd=2,lty=2, ylim=c(0,30),xlab="Sample Size (n)", ylab="%MSE", main=" ")
lines(N,true_pmse_hat, type="l", lwd=3,lty=3, col="red")
lines(N,perc_mse_pte, type="l", lwd=2,lty=1, col="blue")
legend(35, 15,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lwd=c(2,3,2),lty=c(2,3,1),ncol=1)

plot(N,perc_abs_risk_tilde,type="l",lwd=2,lty=2, ylim=c(0,50),xlab="Sample Size (n)", ylab="%MAE", main=" ")
lines(N,perc_abs_risk_hat, type="l", lwd=3,lty=3, col="red")
lines(N,perc_abs_risk_pte, type="l", lwd=2,lty=1, col="blue")
legend(35, 31,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lwd=c(2,3,2),lty=c(2,3,1),ncol=1)

plot(N,perc_LINEX_risk_tilde,type="l",lwd=2,lty=2, ylim=c(0,40),xlab="Sample Size (n)", ylab="%LINEX Risk", main=" ")
lines(N,perc_LINEX_risk_hat, type="l", lwd=3, lty=3, col="red")
lines(N,perc_LINEX_risk_pte, type="l", lwd=2,lty=1, col="blue")
legend(35, 22,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lwd=c(2,3,2),lty=c(2,3,1),ncol=1)


