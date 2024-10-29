require(stats)
require(GoFKernel)    	# For inverting a cdf
require(survival)		# For data-sets
require(AcceptReject)	# For Akash distribution r.v.'s


y<- c()
x<- read.table(file = 'C:/Users/David Giles/Dropbox/Akash distribution/Applications/application_4.txt',header = TRUE)$Cycles

#data(cancer, package="survival")   # second application (2 lines)
#x<- ovarian$futime
n<- length(x)
xbar<- mean(x)

logl<- function(L, x, xbar,n) {
ll<- 3*n*log(L)-n*log(L^2+2)-L*n*xbar+sum(log(1+x^2))
ll
}

Akash<- function (lambda,x) {
pdf<- (lambda^3/(lambda^2+2))*(1+x^2)*exp(-lambda*x)
pdf
}

hazard<- function(lambda,x) {
haz<- lambda^3*(1+x^2)/(lambda*x*(lambda*x+2)+(lambda^2+2))
}

pakash<- function(lambda, x) {
pak<- 1-(1+lambda*x*(lambda*x+2)/(lambda^2+2))*exp(-lambda*x)
pak
}

f <- function(x) pnorm(x, mean=0, sd=1) 
f.inv <- inverse(f,lower=-25,upper=25)

Dplus<- function(x,N,px) {
max(seq(1:N)/N - px)
}

Dminus<- function(x,N,px) {
max(px-(seq(1:N)-1)/N)
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

lam_tilde<- uniroot(leq, c(0.0001,20), xbar=xbar)$root
L<- lam_tilde
bias_CS<- L^4*(L^2+2)^4/(3*n*(L^2+2)^2 -4*n*L^4 +2*n*L^2*(L^2+2))^2*((3*n/L^3)-8*n*L^3/(L^2+2)^3 + 6*n*L/(L^2+2)^2 )
lam_hat<- L- bias_CS 
lam_cup<- uniroot(Firth, c(0.0001,20), xbar=xbar)$root
lam_dot<- uniroot(Godwin, c(0.0001,20), ltilde=lam_tilde)$root

c(lam_tilde, lam_hat, lam_cup, lam_dot)

AIC_tilde<- 2*(1-logl(lam_tilde,x, xbar,n))
AIC_hat<- 2*(1-logl(lam_hat,x, xbar,n))
AIC_cup<- 2*(1-logl(lam_cup,x, xbar,n))
AIC_dot<- 2*(1-logl(lam_dot,x, xbar,n))
BIC_tilde<- log(n)-2*logl(lam_tilde, x, xbar,n)
BIC_hat<- log(n)-2*logl(lam_hat, x, xbar,n)
BIC_cup<- log(n)-2*logl(lam_cup, x, xbar,n)
BIC_dot<- log(n)-2*logl(lam_dot, x, xbar,n)

c(AIC_tilde, AIC_hat, AIC_cup, AIC_dot)
c(BIC_tilde, BIC_hat, BIC_cup, BIC_dot)
#######################################
 # Start of Raschke's testing procedure
#######################################
 
w<- pakash(lam_tilde,x)                 
# CHANGE lam_dot TO lam_hat, lam_cup, AND lam_tilde
    
for (j in 1:n){
y[j]<- f.inv(w[j])                             

}

mu<- mean(y)
sigma<- sd(y)*sqrt((n-1)/n)
ys<- sort(y)
pys<- pnorm(ys, mean=mu, sd=sigma)
pysr<-  rev(pys)    # reverse the order of ths pdf
seqn<- seq(1:n)

Dp<- Dplus(ys,n, pys)
Dm<- Dminus(ys,n, pys)
D<- max(Dp,Dm)
KS<- Dstar<- D*(sqrt(n)-0.01+0.85/sqrt(n))

V<- Dp+Dm
Kuiper<- Vstar<- V*(sqrt(n)+0.05+0.82/sqrt(n))

W2<- sum((pys -(2*seqn-1)/(2*n) )^2) + 1/(12*n)
CvM<- W2star<- W2*(1+0.5/n)

U2<- W2-n*(sum(pys/n) - 0.5)^2   
Watson<-U2star<- U2*(1+0.5/n)

A2<- -n-sum( (2*seqn-1)*(log(pys)+log(1-pysr) ) )/n
AD<- A2star<- A2*(1+0.75/n+2.25/(n^2))

KS_10<- 0.819    # These are from Stephens (1986). See Raaschke, 2011, p.81
Kuiper_10<- 1.386
AD_10<-  0.631   
Watson_10<-  0.096 
CvM_10<- 0.104  
KS_5<- 0.895  
Kuiper_5<- 1.489
AD_5<-  0.752   
Watson_5<-  0.117 
CvM_5<- 0.126 

c(KS_5, Kuiper_5, Watson_5, CvM_5, AD_5)
c(KS_10, Kuiper_10, Watson_10, CvM_10, AD_10)

c(KS, Kuiper, Watson, CvM, AD)

par(mfrow=c(1,1))
p<- ecdf(x)  
plot(p, verticals=TRUE, do.points=FALSE, col="black", main="Figure 1 (a): Distribution Functions
(Yarn Failure Data)", xlab="Cycles", ylab="CDF")
xxfit<- seq(50,1250, length=26)
lam<- lam_tilde
yyfit<- 1-(1+(lam*xxfit*(lam*xxfit+2))/(lam^2+2))*exp(-lam*xxfit)

lines(xxfit, yyfit, col="red", lwd=2,lty=2)
legend(600, 0.2,legend = c(expression("Akash CDF"), 
                         expression("Empirical CDF")),
                         box.col="white", col=c("red","black"),lwd=c(2,1),lty=c(2,1),ncol=1, cex=0.8)

