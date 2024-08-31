set.seed(1234)
require(stats)
require(ggamma)

N<- c(10,15,20,25,30,35, 40,45,50,55,60, 65,70,75,80,85,90,95,100)
nsamp<- 100000	# Number of replications for generating critical values
m<- 20000		# Number of replications for simulating the biases and risks
lambda<- c()
lam1<- c()
sigma_tilde<- c()
sigma_hat<- c()
sigma_pte<- c()
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
perc_abs_risk_tilde<- c()
perc_abs_risk_hat<- c()
perc_abs_risk_pte<- c()
perc_LINEX_risk_tilde<- c()
perc_LINEX_risk_hat<- c()
perc_LINEX_risk_pte<- c()

sigma<- 1.0   # sigma is the scale parameter for Half-Normal
a<- sqrt(2*sigma^2)
b<- 2
k<- 1.0      # Length-Biased Half-Normal = Rayleigh

# Set up outer loop for different sample sizes:
for (j in 1:19) {

n<-N[j]

# Create the critical values - these depend on n 

# ****************************************************************************************
for (i in 1: nsamp) {
x<- rggamma(n,a,b,0.5)			   # Half-Normal with scale parameter sigma; k = 0.5 if Ho True
lambda[i]<- (prod(x))^(1/n)/sqrt(2*sum(x^2)/(n*pi))           # This is for the case of length-biased

}

crit<- quantile(lambda,  probs = c(0.90, 0.95, 0.99))
#crit		 # The critical values are specific to the Rayleigh distribution; the value of the scale parameter; and the sample size

for (i in 1:m) {
y<- rggamma(n,a,b,k)		# Half-Normal
lam1[i]<-  (prod(y))^(1/n)/sqrt(2*sum(y^2)/(n*pi))
sigma_tilde[i]<- sqrt(sum(y^2)/n)
sigma_hat[i]<- sqrt(sum(y^2)/(2*n))

if (lam1[i]<=crit[2]) {		# Alter element from [1], [2], [3] for different sig. levels (10%, 5%, 1%)											
	sigma_pte[i]<- sigma_tilde[i]
} 	else {
	sigma_pte[i]<- sigma_hat[i]
} 

}      # end of loop for a particular n

#pp<-c(sum(lam1>crit[1]), sum(lam1>crit[2]), sum(lam1>crit[3]))/m
#pp
perc_bias_tilde[j]<- 100*(mean(sigma_tilde)-sigma)/sigma
true_pbias_tilde[j]<- 100*(sqrt(2/N[j])*gamma(N[j]+1/2)/gamma(N[j]) - 1)
perc_bias_hat[j]<- 100*(mean(sigma_hat)-sigma)/sigma
true_pbias_hat[j]<- 100*((1/sqrt(N[j]))*gamma(N[j]+1/2)/gamma(N[j]) -1)
perc_bias_pte[j]<- 100*(mean(sigma_pte)-sigma)/sigma
perc_mse_tilde[j]<- 100*var(sigma_tilde)/sigma^2 + (perc_bias_tilde[j])^2/100
true_pmse_tilde[j]<- 100*(3-2*sqrt(2/n)*gamma(N[j]+1/2)/gamma(N[j]) )
perc_mse_hat[j]<- 100*var(sigma_hat)/sigma^2 + (perc_bias_hat[j])^2/100
true_pmse_hat[j]<- 100*2*(1 - (1/sqrt(N[j]))*gamma(N[j]+1/2)/gamma(N[j])  )
perc_mse_pte[j]<- 100*var(sigma_pte)/sigma^2 + (perc_bias_pte[j])^2/100


# Absolute error loss:
# -------------------

perc_abs_risk_tilde[j]<- 100*mean(abs(sigma_tilde-sigma))/sigma
perc_abs_risk_hat[j]<- 100*mean(abs(sigma_hat-sigma))/sigma
perc_abs_risk_pte[j]<- 100*mean(abs(sigma_pte-sigma))/sigma

# LINEX loss:
# -----------

w<- 2.0	#  Shape parameter for LINEX loss function (+ve or -ve)
# w>0 implies over-estimation more costly than under-estimation
# w<0 implies under-estimation more costly than over-estimation

perc_LINEX_risk_tilde[j]<- 100*mean(exp(w*(sigma_tilde-sigma))-w*(sigma_tilde-sigma)-1)/sigma
perc_LINEX_risk_hat[j]<- 100*mean(exp(w*(sigma_hat-sigma))-w*(sigma_hat-sigma)-1)/sigma
perc_LINEX_risk_pte[j]<- 100*mean(exp(w*(sigma_pte-sigma))-w*(sigma_pte-sigma)-1)/sigma


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
plot(N,true_pbias_tilde,type="l",lwd=2,lty=2, ylim=c(-5,50),xlab="Sample Size (n)", ylab="%Bias", main=" ")
lines(N,true_pbias_hat, type="l", lwd=3,lty=3, col="red")
lines(N,perc_bias_pte, type="l", lwd=2,lty=1, col="blue")
legend(35, 30,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lwd=c(2,3,2),lty=c(2,3,1),ncol=1)
plot(N,true_pmse_tilde,type="l",lwd=2,lty=2, ylim=c(0,25),xlab="Sample Size (n)", ylab="%MSE", main=" ")
lines(N,true_pmse_hat, type="l", lwd=3,lty=3, col="red")
lines(N,perc_mse_pte, type="l", lwd=2,lty=1, col="blue")
legend(35, 14,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lwd=c(2,3,2),lty=c(2,3,1),ncol=1)

plot(N,perc_abs_risk_tilde,type="l",lwd=2,lty=2, ylim=c(0,50),xlab="Sample Size (n)", ylab="%MAE", main=" ")
lines(N,perc_abs_risk_hat, type="l", lwd=3,lty=3, col="red")
lines(N,perc_abs_risk_pte, type="l", lwd=2,lty=1, col="blue")
legend(35, 30,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lwd=c(2,3,2),lty=c(2,3,1),ncol=1)

plot(N,perc_LINEX_risk_tilde,type="l",lwd=2,lty=2, ylim=c(0,70),xlab="Sample Size (n)", ylab="%LINEX Risk", main=" ")
lines(N,perc_LINEX_risk_hat, type="l", lwd=3,lty=3, col="red")
lines(N,perc_LINEX_risk_pte, type="l", lwd=2,lty=1, col="blue")
legend(35, 30,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lwd=c(2,3,2),lty=c(2,3,1),ncol=1)


