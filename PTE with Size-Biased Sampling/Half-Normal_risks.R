set.seed(1234)
require(stats)
require(ggamma)

n<-  25
nsamp<- 100000	# Number of replications for generating critical values
m<- 20000		# Number of replications for simulating the biases and risks
lambda<- c()
lam1<- c()
sigma_tilde<- c()
sigma_hat<- c()
sigma_pte5<- c()
sigma_pte_opt_mse<- c()
sigma_pte_opt_mae<- c()
sigma_pte_opt_LINEX<- c()

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

sigma<- 1.0   # sigma is the scale parameter for Half-Normal
a<- sqrt(2*sigma^2)
b<- 2
k<- 0.5      # Half-Normal
w<- 2    # for LINEX
# Create the critical values - these depend on n 

# ****************************************************************************************
for (i in 1: nsamp) {
x<- rggamma(n,a,b,k)			   # Half-Normal with scale parameter sigma
lambda[i]<- (prod(x))^(1/n)/sqrt(2*sum(x^2)/(n*pi))          # This is for the case of length-biased

}

#The following critical values are for n = 25
#############################################
crit<- quantile(lambda,  probs = c(0.95, 0.867, 0.882, 0.807))
crit5<- crit[1]
crit_opt_mse<-    crit[2]
crit_opt_mae<-    crit[3]
crit_opt_LINEX<-   crit[4] # for w=2   


# Set up outer loop for different bias sizes:

#   k=1/2 is half-nornal. Size-biased gives Rayleigh, or k =1

for (j in 1:21) {
k<- 0.5+0.025*(j-1)
K[j]<- 2*(k-0.5)

for (i in 1:m) {
y<- rggamma(n,a,b,k)		# Half-Normal
lam1[i]<- (prod(y))^(1/n)/sqrt(2*sum(y^2)/(n*pi))
sigma_tilde[i]<- sqrt(sum(y^2)/n)
sigma_hat[i]<- sqrt(sum(y^2)/(2*n))

if (lam1[i]<=crit5) {		
	sigma_pte5[i]<- sigma_tilde[i]
} 	else {
	sigma_pte5[i]<- sigma_hat[i]
} 
if (lam1[i]<=crit_opt_mse) {		
	sigma_pte_opt_mse[i]<- sigma_tilde[i]
} 	else {
	sigma_pte_opt_mse[i]<- sigma_hat[i]
} 
if (lam1[i]<=crit_opt_mae) {		
	sigma_pte_opt_mae[i]<- sigma_tilde[i]
} 	else {
	sigma_pte_opt_mae[i]<- sigma_hat[i]
} 
if (lam1[i]<=crit_opt_LINEX) {		
	sigma_pte_opt_LINEX[i]<- sigma_tilde[i]
} 	else {
	sigma_pte_opt_LINEX[i]<- sigma_hat[i]
} 


}      # end of loop for a particular k

#pp<-c(sum(lam1>crit[1]), sum(lam1>crit[2]), sum(lam1>crit[3]))/m
#pp
perc_bias_tilde[j]<- 100*(mean(sigma_tilde)-sigma)/sigma
perc_bias_hat[j]<- 100*(mean(sigma_hat)-sigma)/sigma
perc_bias_pte5[j]<- 100*(mean(sigma_pte5)-sigma)/sigma
perc_bias_pte_opt_mse[j]<- 100*(mean(sigma_pte_opt_mse)-sigma)/sigma
perc_bias_pte_opt_mae[j]<- 100*(mean(sigma_pte_opt_mae)-sigma)/sigma
perc_bias_pte_opt_LINEX[j]<- 100*(mean(sigma_pte_opt_LINEX)-sigma)/sigma

perc_mse_tilde[j]<- 100*var(sigma_tilde)/sigma^2 + (perc_bias_tilde[j])^2/100
perc_mse_hat[j]<- 100*var(sigma_hat)/sigma^2 + (perc_bias_hat[j])^2/100
perc_mse_pte5[j]<- 100*var(sigma_pte5)/sigma^2 + (perc_bias_pte5[j])^2/100
perc_mse_pte_opt[j]<- 100*var(sigma_pte_opt_mse)/sigma^2 + (perc_bias_pte_opt_mse[j])^2/100

perc_mae_tilde[j]<- 100*(mean(abs(sigma_tilde-sigma)))/sigma
perc_mae_hat[j]<- 100*(mean(abs(sigma_hat-sigma)))/sigma
perc_mae_pte5[j]<- 100*(mean(abs(sigma_pte5-sigma)))/sigma
perc_mae_pte_opt[j]<- 100*(mean(abs(sigma_pte_opt_mae - sigma)))/sigma

perc_LINEX_risk_tilde[j]<- 100*mean(exp(w*(sigma_tilde-sigma))-w*(sigma_tilde-sigma)-1)/sigma
perc_LINEX_risk_hat[j]<- 100*mean(exp(w*(sigma_hat-sigma))-w*(sigma_hat-sigma)-1)/sigma
perc_LINEX_risk_pte5[j]<- 100*mean(exp(w*(sigma_pte5-sigma))-w*(sigma_pte5-sigma)-1)/sigma
perc_LINEX_risk_pte_opt[j]<- 100*mean(exp(w*(sigma_pte_opt_LINEX-sigma))-w*(sigma_pte_opt_LINEX-sigma)-1)/sigma

}      # end of loop for all k


par(mfrow=c(3,1))
plot(K,perc_mse_tilde,type="l",lwd=2,lty=2, ylim=c(1,9),xlab="c", ylab="%MSE", main=" ")
lines(K,perc_mse_hat, type="l", lwd=2,lty=5, col="red")
lines(K,perc_mse_pte_opt, type="l", lwd=2,lty=1, col="blue")
lines(K, perc_mse_pte5, type="l", lwd=3,lty=3, col="blue")
legend("topright",inset=c(0.01,0.01),legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("PTE (opt)"),
                         expression ("PTE (5%)")), 
                         box.col="white", col=c("black","red","blue", "blue"),lwd=c(2,2,2,3),lty=c(2,5,1,3),ncol=1)
plot(K,perc_mae_tilde,type="l",lwd=2,lty=2, ylim=c(5,30),xlab="c", ylab="%MAE", main=" ")
lines(K,perc_mae_hat, type="l", lwd=2,lty=5, col="red")
lines(K,perc_mae_pte_opt, type="l", lwd=2,lty=1, col="blue")
lines(K, perc_mae_pte5, type="l", lwd=3,lty=3, col="blue")
legend("topright", inset=c(0.01,0.01),legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("PTE (opt)"),
                         expression ("PTE (5%)")), 
                         box.col="white", col=c("black","red","blue", "blue"),lwd=c(2,2,2,3),lty=c(2,5,1,3),ncol=1)
plot(K,perc_LINEX_risk_tilde,type="l",lwd=2,lty=2, ylim=c(0,30),xlab="c", ylab="%LINEX Risk", main=" ")
lines(K,perc_LINEX_risk_hat, type="l", lwd=2,lty=5, col="red")
lines(K,perc_LINEX_risk_pte_opt, type="l", lwd=2,lty=1, col="blue")
lines(K, perc_LINEX_risk_pte5, type="l", lwd=3,lty=3, col="blue")
legend("topright", inset=c(0.01,0.01),legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("PTE (opt)"),
                         expression ("PTE (5%)")), 
                         box.col="white", col=c("black","red","blue", "blue"),lwd=c(2,2,2,3),lty=c(2,5,1,3),ncol=1)