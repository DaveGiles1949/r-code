library(BAS)
library(VGAM)
# Check:
hypergeometric1F1(100, 102,10)  # log is TRUE answer should be 9.8127
################

z<- c(-3.5, -3,-2.5,-2,-1.5,-1,-0.5, 0,0.5,1,1.5,2,2.5,3,3.5)
pi<- 3.14159
Phi<- c()
# the "hypergeometric1F1" function needs a SCALR third element, and the default return is the LOG of the function
for(ii in 1:15){
Phi[ii]<- 0.5*(1+ (z[ii]*sqrt(2)/sqrt(pi))*(hypergeometric1F1(0.5, 1.5,-z[ii]^2/2, log=FALSE)))
}

phi<- 0.5*(1+erf(z/sqrt(2))) # Alternative formulation of std. normal CDF
plot(z,Phi, main="Std. Normal CDF: 2 Equivalent Calculations", ylab="CDF(z)")
lines(z, phi, col="red")

# Everything checks out
