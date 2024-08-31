set.seed(1234)
require(stats)
require(ggamma)

n<- 250
nsamp<- 1000000	# Number of replications for generating critical values
m<- 20000		# Number of replications for simulating the biases and risks
lambda<- c()
lam1<- c()
sigma_tilde<- c()
sigma_hat<- c()
sigma_pte<- c()
perc_LINEX_tilde<- c()
perc_LINEX_hat<- c()
perc_LINEX_pte<- c()
K<- c()
diff<- c()

sigma<- 1.0   # sigma is the scale parameter for Half-Normal
a<- sqrt(2*sigma^2)
b<- 2
k<- 0.5    # Half-Normal
w<- -2     # For the LINEX risk function

# Create the critical values - these depend on n 

# ****************************************************************************************
for (i in 1: nsamp) {
x<- rggamma(n,a,b,k)			   # Half-Normal with scale parameter sigma
lambda[i]<- (prod(x))^(1/n)/sqrt(2*sum(x^2)/(n*pi))          # This is for the case of length-biased

}

crit<- quantile(lambda,  probs = c(0.999985, 0.95, 0.99))
crit		 # The critical values are specific to the Half-Normal distribution; the value of the scale parameter; and the sample size

# Set up outer loop for different bias sizes:

#   k = 0.5 is half-normal. Size-biased gives Rayleigh, or k = 1

for (j in 1:21) {
k<- 0.5+0.025*(j-1)
K[j]<- 2*(k-0.5)

for (i in 1:m) {
y<- rggamma(n,a,b,k)		
lam1[i]<- (prod(y))^(1/n)/sqrt(2*sum(y^2)/(n*pi))
sigma_tilde[i]<- sqrt(sum(y^2)/n)
sigma_hat[i]<- sqrt(sum(y^2)/(2*n))

if (lam1[i]<=crit[1]) {		# Alter element from [1], [2], [3] for different sig. levels (10%, 5%, 1%)											
	sigma_pte[i]<- sigma_tilde[i]
} 	else {
	sigma_pte[i]<- sigma_hat[i]
} 

}      # end of loop for a particular k

perc_LINEX_tilde[j]<- 100*mean(exp(w*(sigma_tilde-sigma))-w*(sigma_tilde-sigma)-1)/sigma
perc_LINEX_hat[j]<- 100*mean(exp(w*(sigma_hat-sigma))-w*(sigma_hat-sigma)-1)/sigma
perc_LINEX_pte[j]<- 100*mean(exp(w*(sigma_pte-sigma))-w*(sigma_pte-sigma)-1)/sigma

diff[j]<- (perc_LINEX_pte[j] - min(perc_LINEX_tilde[j], perc_LINEX_hat[j]))

}      # end of loop for all k


perc_LINEX_tilde
perc_LINEX_hat
perc_LINEX_pte


plot(K,perc_LINEX_tilde,type="l",lty=2, ylim=c(0,9),xlab="c", ylab="% LINEX risk", main="Half-Normal: %LINEX Risk 
(n = 25)")
lines(K,perc_LINEX_hat, type="l", lty=3, col="red")
lines(K,perc_LINEX_pte, type="l", lty=1, col="blue")
legend(0.58, 16,legend = c(expression("Unweighted"),
                         expression("Length-Weighted"),
				 expression("Pre-Test")), 
                         box.col="white", col=c("black","red","blue"),lty=c(2,3,1),ncol=1)
max(diff)
crit


