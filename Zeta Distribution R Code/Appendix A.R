# David Giles, October 2025
# -------------------------

library(VGAM)

# Gamma constants -  from Choudhury (1995), Table 5 (Note the shift in decimal places, as indicated in that table.)
gamma1<-  -0.072815845483676724861
gamma2<-  -0.009690036031902870231084845
gamma3<-   0.0020538344203033458662
gamma4<-   0.0023253700654673000575
gamma5<-   0.00079332381730106270175
gamma6<-  -0.00023876934543019960987
gamma7<-  -0.00052728956705775104607
gamma8<-  -0.00035212335380303950960
gamma9<-  -0.000034394774418088048178
gamma10<-  0.00020533281490906479468
gamma11<-  0.00027018443954390352667
gamma12<-  0.00016727291210514019335
gamma13<- -0.000027463806603760158860
gamma14<- -0.00020920926205929994584
gamma15<- -0.00028346865532024144664
gamma16<- -0.00019969685830896977471
gamma17<-  0.000026277037109918336699
gamma18<-  0.00030736840814925282659
gamma19<-  0.00050360545304735562606
gamma20<-  0.00046634356151155944940

#Taking account of the fact that we have a canonical linear exponential model -

s<- seq(from = 1, to = 5, by = 0.005)
x<- s
y<- c()
y<- -zeta(x,1,1)/zeta(x,0,1) +
(zeta(x,0,1)*
( -(6/(x-1)^4+gamma3-gamma4*(x-1)+(gamma5*(x-1)^2)/2-(gamma6*(x-1)^3)/6+(gamma7*(x-1)^4)/24-(gamma8*(x-1)^5)/120
+(gamma9*(x-1)^6)/720-(gamma10*(x-1)^7)/5040+(gamma11*(x-1)^8)/40320-(gamma12*(x-1)^9)/362880+(gamma13*(x-1)^10)/3628800
- (gamma14*(x-1)^11)/39916800 + (gamma15*(x-1)^12)/479001600 - (gamma16*(x-1)^13)/6227020800 + (gamma17*(x-1)^14)/87178291200 
- (gamma18*(x-1)^15)/1307674368000 )) 

-zeta(x,1,1)*zeta(x,2,1)) /(2*(zeta(x,0,1)*zeta(x,2,1) - zeta(x,1,1)^2))

par(mfrow=c(3,1))

y_sub1<- y[s>1 & s<=4.25]
s_sub1<- s[s>1 & s<=4.25]

plot(s_sub1,y_sub1, xaxt="n",col="red", type="l", xlab="s",ylab="f(s)",main= "Figure A.1(a): f(s)")
axis(1,xaxp=c(1,4.25,13), las=2)
abline(v = 1.25, col = "blue", lty = 2, lwd = 2)
text(x=1.4, y=-150,labels="s = 1.25", col = "blue")

y_sub2<- y[s>=1.25 & s<=5] 
s_sub2<- s[s>=1.25 & s<=5]

plot(s_sub2,y_sub2, xaxt="n",col="red", type="l",lty=1, xlab="s",ylab="f(s)",main= "Figure A.1(b): f(s)")
axis(1,xaxp=c(1.25,5,15), las=2)

y_sub3<- y[s>=3.5 & s<=5.0]
s_sub3<- s[s>=3.5 & s<= 5.0]

ols<- lm(y_sub3 ~ s_sub3)

plot(s_sub3,y_sub3, xaxt="n",col="red", type="l",lty=1, xlab="s",ylab="f(s)",main= "Figure A.1(c): f(s)")
axis(1,xaxp=c(3.5, 5.0,15),las=2)
lines(s_sub3,ols$fit, col ="black", type="l",lty=2)
legend("bottomright",legend=c("f(s)", "OLS"),col=c("red", "black"), lty=c(1,2))

summary(ols)