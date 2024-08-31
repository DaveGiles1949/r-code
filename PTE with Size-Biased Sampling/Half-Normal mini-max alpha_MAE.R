set.seed(1234)
require(stats)
require(ggamma)

n<- 250
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
perc_mae_tilde<- c()
perc_mae_hat<- c()
perc_mae_pte<- c()
K<- c()
diff<- c()
maxdiff<- c()

sigma<- 1.0   # sigma is the scale parameter for Half-Normal
a<- sqrt(2*sigma^2)
b<- 2
k<- 0.5      # Half-Normal
# Create the critical values - these depend on n 

# ****************************************************************************************
for (i in 1: nsamp) {
x<- rggamma(n,a,b,k)			   # Half-Normal with scale parameter sigma
lambda[i]<- (prod(x))^(1/n)/sqrt(2*sum(x^2)/(n*pi))          # This is for the case of length-biased

}

crit<- quantile(lambda,  probs = c(0.99993, 0.95, 0.99))
#crit		 # The critical values are specific to the Rayleigh distribution; the value of the scale parameter; and the sample size

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
perc_mae_tilde[j]<- 100*(mean(abs(sigma_tilde-sigma)))/sigma
perc_mae_hat[j]<- 100*(mean(abs(sigma_hat-sigma)))/sigma
perc_mae_pte[j]<- 100*(mean(abs(sigma_pte-sigma)))/sigma

diff[j]<- (perc_mae_pte[j] - min(perc_mae_tilde[j], perc_mae_hat[j]))

}      # end of loop for all k

max(diff)
#maxdiff[jj]<- max(diff)

perc_bias_tilde
perc_bias_hat
perc_bias_pte

perc_mae_tilde
perc_mae_hat
perc_mae_pte


plot(K,perc_mae_tilde,type="l",lty=2, ylim=c(5,40),xlab="c", ylab="% MAE", main="Half-Normal; %MAE; n=25;alpha = 5%")
lines(K,perc_mae_hat, type="l", lty=1, col="red")
lines(K,perc_mae_pte, type="l", lty=1, col="blue")
legend(0.6, 20,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lty=c(2,1,1),ncol=1)


