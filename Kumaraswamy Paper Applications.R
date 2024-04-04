require(GoFKernel)    # For inverting a cdf
require(nortest)      # For A-D test
require(fitdistrplus) # For fitting the Beta-distribution
require(univariateML)  # For the Kumaraswamy MLE      
require(VGAM)          # For the Kumaraswamy distribution 
require(stats)

# Note: The 2 data series that are used can be accessed at   https://github.com/DaveGiles1949/Data

HE_dat<- read.table(file = 'C:/Users/David Giles/Dropbox/GOF Testing/Kumaraswamy/Kumaraswamy Test_HE Data.txt',header = TRUE)
HE17<- HE_dat[,27]/100

gini<- read.table(file = 'C:/Users/David Giles/Dropbox/GOF Testing/Kumaraswamy/Kumaraswamy Test_Gini Data.txt',header = TRUE)$GINI/100

x<- c()
wk<- c()
wb<- c()
AD<- c()
KS<- c()
CvM<- c()
Watson<- c()
Kuiper<- c()
shape1<- c()    #[1] for Kumar; [2] for Beta
shape2<- c()    #[1] for Kumar; [2] for Beta

# Use the next 2 lines to select the data for each application
x<- gini     # Gini Index data
#x<- HE17      # Hidden Economy data 
n<- length(x)
y<- matrix(nrow=n, ncol=2)

#Kumar:
kumar<- mlkumar(x)
shape1[1]<- kumar[[1]]
shape2[1]<- kumar[[2]]

#Beta:
fbeta <- fitdist(x, "beta")
shape1[2]<- fbeta$estimate[[1]]
shape2[2]<- fbeta$estimate[[2]]


wk<- pkumar(x, shape1[1], shape2[1])  
wb<- pbeta(x, shape1[2], shape2[2])           # Start of Raaschke's testing procedure
f <- function(x) pnorm(x, mean=0, sd=1)
f.inv <- inverse(f,lower=-20,upper=20)

for (i in 1:n){
y[i,1]<- f.inv(wk[i])
y[i,2]<- f.inv(wb[i])

}

# Start the EDF tests

Dplus<- function(x,N,px) {
max(seq(1:N)/N - px)
}

Dminus<- function(x,N,px) {
max(px-(seq(1:N)-1)/N)
}

s<- seq(1:n)

for (j in 1:2){

mu<- mean(y[,j])
sigma<- sd(y[,j])*sqrt((n-1)/n)
ys<- sort(y[,j])
pys<- pnorm(ys, mean=mu, sd=sigma)
pysr<-  rev(pys)    # reverse the order of the pdf

Dp<- Dplus(ys,n, pys)
Dm<- Dminus(ys,n, pys)
D<- max(Dp,Dm)
KS[j]<- Dstar<- D*(sqrt(n)-0.01+0.85/sqrt(n))      # [1] for Kumar; [2] for Beta

V<- Dp+Dm
Kuiper[j]<- Vstar<- V*(sqrt(n)+0.05+0.82/sqrt(n))

W2<- sum((pys -(2*s-1)/(2*n) )^2) + 1/(12*n)
CvM[j]<- W2star<- W2*(1+0.5/n)

U2<- W2-n*(sum(pys/n) - 0.5)^2   
Watson[j]<-U2star<- U2*(1+0.5/n)

A2<- -n-sum( (2*s-1)*(log(pys)+log(1-pysr) ) )/n
AD[j]<- A2star<- A2*(1+0.75/n+2.25/(n^2))

}

shape1
shape2
crit_KS_10<- 0.819  		# See Raschke, 2011, p.81
crit_Kuiper_10<- 1.386
crit_AD_10<-  0.631   
crit_Watson_10<-  0.096 
crit_CvM_10<- 0.104  
crit_KS_5<- 0.895  
crit_Kuiper_5<- 1.489
crit_AD_5<-  0.752   
crit_Watson_5<-  0.117 
crit_CvM_5<- 0.126 

AD
c(crit_AD_5, crit_AD_10)
KS
c(crit_KS_5, crit_KS_10)
Kuiper
c(crit_Kuiper_5, crit_Kuiper_10)
Watson
c(crit_Watson_5, crit_Watson_10)
CvM
c(crit_CvM_5, crit_CvM_10)

# Plot the results:
#----------------
# First for the HE data:
# Note that the shape coefficients have been computed with x=HE data

par(mfrow=c(1,1))
#y<- rkumar(n,shape1[1], shape2[1])
#y<- rbeta(n, shape1[2], shape2[2])
h<- hist(x, prob=TRUE,main="Figure 3 (a): Hidden Economy Densities", xlab="HE (share of GDP)", col="seashell2",breaks=12)
xfit<- seq(0,1,length=100)
yfitk<- dkumar(xfit,shape1[1], shape2[1])
yfit<- dbeta(xfit, shape1[2], shape2[2])
lines(xfit, yfit, col="tomato3", lwd=1)
lines(xfit, yfitk, col="blue", lwd=1, lty=2)
legend(0.41, 3.36,legend = c(expression("Beta PDF"),
                         expression("Kuma. PDF")), 
                         box.col="white", col=c("red","blue"),lty=c(1,2),ncol=1, cex=0.8)
p<- ecdf(x)  
plot(p, verticals=TRUE, do.points=FALSE, col="black", main="Figure 3 (b): Distribution Functions", xlab="HE (share of GDP)", ylab="CDF")
xxfit<- seq(0,1, length=100)
yyfitk<- pkumar(xxfit, shape1[1], shape2[1])
yyfitb<- pbeta(xxfit, shape1[2], shape2[2])
lines(xxfit, yyfitb, col="red", lwd=1)
lines(xxfit, yyfitk, col="blue", lwd=1,lty=2)
legend(0.26, 0.4,legend = c(expression("Beta CDF"), expression("Kuma. CDF"),
                         expression("Empirical CDF")),
                         box.col="white", col=c("red","blue", "black"),lty=c(1,2,1),ncol=1, cex=0.8)

# Now for the Gini coefficient data:
# Note that the shape coefficients have been computed with x=GINI data

par(mfrow=c(1,1))
#y<- rkumar(n,shape1[1], shape2[1])
#y<- rbeta(n, shape1[2], shape2[2])
h<- hist(x, prob=TRUE,main="Figure 4 (a): Gini Index Densities", xlab="Gini Index", col="seashell2",breaks=12)
xfit<- seq(0,1,length=100)
yfitk<- dkumar(xfit,shape1[1], shape2[1])
yfit<- dbeta(xfit, shape1[2], shape2[2])
lines(xfit, yfit, col="tomato3", lwd=1)
lines(xfit, yfitk, col="blue", lwd=1, lty=2)
legend(0.43, 5,legend = c(expression("Beta PDF"),
                         expression("Kuma. PDF")), 
                         box.col="white", col=c("red","blue"),lty=c(1,2),ncol=1, cex=0.8)
p<- ecdf(x)  
plot(p, verticals=TRUE, do.points=FALSE, col="black", main="Figure 4 (b): Distribution Functions", xlab="Gini Index", ylab="CDF")
xxfit<- seq(0,1, length=100)
yyfitk<- pkumar(xxfit, shape1[1], shape2[1])
yyfitb<- pbeta(xxfit, shape1[2], shape2[2])
lines(xxfit, yyfitb, col="red", lwd=1)
lines(xxfit, yyfitk, col="blue", lwd=1,lty=2)
legend(0.31, 0.25,legend = c(expression("Beta CDF"), expression("Kuma. CDF"),
                         expression("Empirical CDF")),
                         box.col="white", col=c("red","blue", "black"),lty=c(1,2,1),ncol=1, cex=0.8)
