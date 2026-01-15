require(stats)
require(GoFKernel)    	# For inverting a cdf
require(survival)		# For data-sets
require(AcceptReject)	# For Akash distribution r.v.'s

density<- function (lambda,x) {
pdf<- (lambda^3/(lambda^2+2))*(1+x^2)*exp(-lambda*x)
pdf
}

hazard<- function(lambda,x) {
haz<- lambda^3*(1+x^2)/(lambda*x*(lambda*x+2)+(lambda^2+2))
haz
}

par(mfrow=c(1,1))
#p<- ecdf(x)  
#plot(p, verticals=TRUE, do.points=FALSE, col="black", main="Figure 2 (b): Distribution Functions
#(Ovarian Cancer Data)", xlab="Days", ylab="CDF")
#xxfit<- seq(50,1250, length=26)
#lam<- lam_tilde
#yyfit<- 1-(1+(lam*xxfit*(lam*xxfit+2))/(lam^2+2))*exp(-lam*xxfit)

#lines(xxfit, yyfit, col="red", lwd=2,lty=2)
#legend(800, 0.2,legend = c(expression("Akash CDF"), 
#                         expression("Empirical CDF")),
#                         box.col="white", col=c("red","black"),lwd=c(2,1),lty=c(2,1),ncol=1, cex=0.8)

# p.d.f. and hazard function
#---------------------------

xx<- seq(0,20 , length=100)
pdf1<- density(0.05, xx)
pdf2<- density(0.1,xx)
pdf3<- density(0.5,xx)
pdf4<- density(1.0,xx)

plot(xx, pdf1, main = "Figure 1 (a): Akash Density Functions", type="l",lty=1, ylim=c(0,0.4))
lines(xx, pdf2, lty=1, col="red")
lines(xx, pdf3, lty=2, col="black")
lines(xx, pdf4, lty=2, col="red")

legend(7,0.3, legend = c(expression("lambda = 0.05"),
				expression("lambda = 0.1"),
				expression("lambda = 0.5"),
				expression("lambda = 1.0")),
				box.col="white", col=c("black", "red","black","red"), lty=c(1,1,2,2), ncol=1, cex=0.8)
haz1<- hazard(1,xx)
haz2<- hazard(2,xx)

plot(xx, haz1, main = "Figure 1 (b): Akash Hazard Functions", type="l",lty=1, ylim=c(0.2,2))
lines(xx, haz2, lty=2, col="red")
legend(7,0.7, legend = c(expression("lambda = 1"),
				expression("lambda = 2")),
				box.col="white", col=c("black", "red"), lty=c(1,2), ncol=1, cex=0.8)



