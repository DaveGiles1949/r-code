# Code for evaluating the density function of "p-hat"
####################################################

library(BAS)
# Needed for the confluent hypergeometric function

phat<- c(-0.99,-0.98, -0.96,-0.94,-0.92,-0.9,-0.8,-.7,-.6,-.5,-.4,-.3,-.2,-.1,0,0.1,.2,.3,.4,.5,
.6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,
4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,7.2,7.4,7.6,7.8,8,8.2,8.4,8.6,8.8,9,9.5,10,10.5,11,11.5,12,13,14,15,16,17,18,19,20)

# Case 1: sigsq=100

c<- 2
sigsq<- 100
d<- 0.022
v<-  5
alpha<- log(phat+1)-c
beta<- 1/(sigsq*d)  
bav22<- 0.5*beta*(alpha+v)^2
av2b<- (alpha+v)*sqrt(2*beta)
v4<- v/4
v2<- v/2
v24<- (2+v)/4
v12<- (v+1)/2
kpp<- c()
Dmv2<- c()
LM1<- c()
LM2<- c()
f1<- c()

# The "hypergeometric1F1" function requires a scalar third element, hence the loop

for (ii in 1:length(phat)){

kpp[ii]<- v^v2*beta^v12/(sqrt(2*pi)*gamma(v2))

LM1[ii]<- hypergeometric1F1(v4, 0.5, bav22[ii])   # these are LOGS of the Kummer function
LM2[ii]<- hypergeometric1F1(v24,1.5,bav22[ii])

Dmv2[ii]<- 1/(2^v4)*sqrt(pi)*exp(-(bav22[ii]/2))*(exp(LM1[ii])/gamma(v24) - (av2b[ii]/gamma(v4))*exp(LM2[ii]))

f1[ii]<- kpp[ii]*Dmv2[ii]*(phat[ii]+1)^(v*beta-1)*beta^(-v4)*gamma(v2)*exp(0.5*bav22[ii])*exp(-beta*(v*log(phat[ii]+1)+alpha[ii]^2/2))
# Simplify the line above
}

# Case 2: sigmasq=50

c<- 2
sigsq<- 50
d<- 0.022
v<-  5
alpha<- log(phat+1)-c
beta<- 1/(sigsq*d)  
bav22<- 0.5*beta*(alpha+v)^2
av2b<- (alpha+v)*sqrt(2*beta)
v4<- v/4
v2<- v/2
v24<- (2+v)/4
v12<- (v+1)/2
kpp<- c()
Dmv2<- c()
LM1<- c()
LM2<- c()
f2<- c()

# The "hypergeometric1F1" function requires a scalar third element, hence the loop

for (ii in 1:length(phat)){

kpp[ii]<- v^v2*beta^v12/(sqrt(2*pi)*gamma(v2))

LM1[ii]<- hypergeometric1F1(v4, 0.5, bav22[ii])   # these are LOGS of the Kummer function
LM2[ii]<- hypergeometric1F1(v24,1.5,bav22[ii])

Dmv2[ii]<- 1/(2^v4)*sqrt(pi)*exp(-(bav22[ii]/2))*(exp(LM1[ii])/gamma(v24) - (av2b[ii]/gamma(v4))*exp(LM2[ii]))

f2[ii]<- kpp[ii]*Dmv2[ii]*(phat[ii]+1)^(v*beta-1)*beta^(-v4)*gamma(v2)*exp(0.5*bav22[ii])*exp(-beta*(v*log(phat[ii]+1)+alpha[ii]^2/2))
# Simplify the line above
}

plot(phat,f1, xlab="p-hat", cex.main=0.8,cex.lab=0.7, cex.axis=0.6,type="l",col="red", xlim=c(-1,20), main="Figure 1: Density of p-hat
(v = 5, c = 2, d=0.022", ylab="f (p-hat)",ylim=c(0,0.4))
lines(phat,f2, type="l",col="blue")


