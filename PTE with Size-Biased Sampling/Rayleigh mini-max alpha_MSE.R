set.seed(1234)
require(stats)
require(ggamma)

n<-  15
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
K<- c()
diff<- c()
maxdiff<- c()

sigma<- 1.0   # sigma is the scale parameter for Rayleigh
a<- sqrt(2*sigma^2)
b<- 2
k<- 1.0      # Rayleigh
# Create the critical values - these depend on n and the particular parameter values

# ****************************************************************************************
for (i in 1: nsamp) {
x<- rggamma(n,a,b,k)			   # Rayleigh with scale parameter sigma
lambda[i]<- (prod(x))^(1/n)/sqrt(pi*sum(x^2)/(4*n))          # This is for the case of length-biased

}
crit<- quantile(lambda,  probs = c(0.656, 0.95, 0.99))
# The critical values are specific to the Rayleigh distribution; and the sample size

# Set up outer loop for different bias sizes:

#   k=1 is Rayleigh. Size-biased gives k = 3/2

for (j in 1:21) {
k<- 1.0+0.025*(j-1)
K[j]<- 2*(k-1.0)

for (i in 1:m) {
y<- rggamma(n,a,b,k)		# Rayleigh
lam1[i]<- (prod(y))^(1/n)/sqrt(pi*sum(y^2)/(4*n))
sigma_tilde[i]<- sqrt(sum(y^2)/(2*n))
sigma_hat[i]<-   sqrt(sum(y^2)/(3*n))

if (lam1[i]<=crit[1]) {		# Alter element from [1], [2], [3] for different sig. levels (10%, 5%, 1%)											
	sigma_pte[i]<- sigma_tilde[i]
} 	else {
	sigma_pte[i]<- sigma_hat[i]
} 

}      # end of loop for a particular k

#pp<-c(sum(lam1>crit[1]), sum(lam1>crit[2]), sum(lam1>crit[3]))/m
#pp
perc_bias_tilde[j]<- 100*(mean(sigma_tilde)-sigma)/sigma
perc_bias_hat[j]<- 100*(mean(sigma_hat)-sigma)/sigma
perc_bias_pte[j]<- 100*(mean(sigma_pte)-sigma)/sigma
perc_mse_tilde[j]<- 100*var(sigma_tilde)/sigma^2 + (perc_bias_tilde[j])^2/100
perc_mse_hat[j]<- 100*var(sigma_hat)/sigma^2 + (perc_bias_hat[j])^2/100
perc_mse_pte[j]<- 100*var(sigma_pte)/sigma^2 + (perc_bias_pte[j])^2/100

diff[j]<- (perc_mse_pte[j] - min(perc_mse_tilde[j], perc_mse_hat[j]))

}      # end of loop for all k

max(diff)
#maxdiff[jj]<- max(diff)

perc_bias_tilde
perc_bias_hat
perc_bias_pte

perc_mse_tilde
perc_mse_hat
perc_mse_pte


plot(K,perc_mse_tilde,type="l",lty=2, ylim=c(0,5),xlab="c", ylab="% MSE", main="Rayleigh; %MSE; n=25;alpha = 5%")
lines(K,perc_mse_hat, type="l", lty=3, col="red")
lines(K,perc_mse_pte, type="l", lty=1, col="blue")
legend(0.6, 3,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lty=c(2,3,1),ncol=1)
max(diff)
crit

