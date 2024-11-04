require(AcceptReject)
require(stats)
set.seed(1234)

lam_hat<- c()
lam_tilde<- c()
lam_cup<- c()
lam_dot<- c()

nrep<-50000
n<- 100
lambda<- 5
mu<- (lambda^2+6)/(lambda*(lambda^2+2))
sd<- sqrt((lambda^4+16*lambda^2+12)/(lambda^2*(lambda^2+2)^2))

Akash<- function (lambda,x) {
pdf<- (lambda^3/(lambda^2+2))*(1+x^2)*exp(-lambda*x)
pdf
}

leq<- function(lambda,xbar){
le<- lambda^3*xbar - lambda^2 + 2*xbar*lambda - 6
le
}

ckbias<- function(L) {
 L^4*(L^2+2)^4/(3*n*(L^2+2)^2 -4*n*L^4 +2*n*L^2*(L^2+2))^2*((3*n/L^3)-8*n*L^3/(L^2+2)^3 + 6*n*L/(L^2+2)^2 )

}

Godwin<- function (ldot, ltilde) {
gg<- ldot- ltilde+ckbias(ldot)
}

Firth<- function(L, xbar) {
grad<- 3/L - 2*L/(L^2+2)-xbar
I<- (3*(L^2+2)^2 -4*L^4 +2*L^2*(L^2+2))/(L^2*(L^2+2)^2)
B<- L^4*(L^2+2)^4/(3*n*(L^2+2)^2 -4*n*L^4 +2*n*L^2*(L^2+2))^2*((3*n/L^3)-8*n*L^3/(L^2+2)^3 + 6*n*L/(L^2+2)^2 )
ff<- grad-I*B
ff 
}

# Start MC Loop:

for(i in 1:nrep) {

y <- accept_reject(
n = n,
f = Akash,
args_f=list(lambda=lambda),
continuous = TRUE,
xlim = c(0, 50)
) 

x<- as.vector(y)
xbar<- mean(x)

lam_tilde[i]<- uniroot(leq, c(0.0001,50), xbar=xbar)$root
L<- lam_tilde[i]
bias_CS<- L^4*(L^2+2)^4/(3*n*(L^2+2)^2 -4*n*L^4 +2*n*L^2*(L^2+2))^2*((3*n/L^3)-8*n*L^3/(L^2+2)^3 + 6*n*L/(L^2+2)^2 )
lam_hat[i]<- L- bias_CS 
lam_cup[i]<- uniroot(Firth, c(0.0001,50), xbar=xbar)$root
lam_dot[i]<- uniroot(Godwin, c(0.0001,50), ltilde=lam_tilde[i])$root

}

perc_bias_tilde<- 100*(mean(lam_tilde)-lambda)/lambda
perc_mse_tilde<- (perc_bias_tilde^2)/100 + 100*(var(lam_tilde))/lambda^2
c(perc_bias_tilde, perc_mse_tilde)
perc_bias_hat<- 100*(mean(lam_hat)-lambda)/lambda
perc_mse_hat<- (perc_bias_hat^2)/100 + 100*(var(lam_hat))/lambda^2
c(perc_bias_hat, perc_mse_hat)
perc_bias_cup<- 100*(mean(lam_cup)-lambda)/lambda
perc_mse_cup<- (perc_bias_cup^2)/100 + 100*(var(lam_cup))/lambda^2
c(perc_bias_cup, perc_mse_cup)
perc_bias_dot<- 100*(mean(lam_dot)-lambda)/lambda
perc_mse_dot<- (perc_bias_dot^2)/100 + 100*(var(lam_dot))/lambda^2
c(perc_bias_dot, perc_mse_dot)
