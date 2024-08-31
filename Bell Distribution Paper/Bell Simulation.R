# NOTE: MLE = MOM

library(bellreg)           # Used to generate random "Bell" data
library(lamW)              # Used for the Lambert function to get MLE

set.seed(123)
nrep<- 50000
n<- 250
theta<- 1.0

theta_hat<- vector()            # MLE
theta_tilde<- vector()          # Cox-Snell
theta_cup<- vector()            # Firth
theta_star<- vector()           # Godwin-Giles

fun_firth<- function(x,dat,n) {
n*(dat/x-exp(x))+(x+2)/(2*(x+1))
}                                     # Used to get Firth estimator

fun_gg<- function(x,n,mle) {
x-x*(x+2)/(2*n*exp(x)*(1+x)^2)-mle
}                                     # Used to get Godwin-Giles estimator

for (i in 1:nrep) {
y<- rbell(n,theta)
ybar<-mean(y)
theta_tilde[i]<- lambertW0(ybar)    # Data are positive, so only the principal branch of the Lambert function exits.
est_bias<- -theta_tilde[i]*(theta_tilde[i]+2)/(2*n*exp(theta_tilde[i])*(1+theta_tilde[i])^2)
theta_hat[i]<- theta_tilde[1]-est_bias
theta_cup[i]<-uniroot(fun_firth, dat=ybar,n=n, lower=0.0001, upper=theta+1)$root     # Upper limit set acccording to true theta value
theta_star[i]<- uniroot(fun_gg, n=n, mle=theta_tilde[i], lower=0.0001, upper=theta+1)$root
}

summary(theta_tilde)
summary(theta_hat)
summary(theta_cup)
summary(theta_star)

perc_bias_tilde<- 100*(mean(theta_tilde)-theta)/theta
perc_bias_hat<- 100*(mean(theta_hat)-theta)/theta
perc_bias_cup<- 100*(mean(theta_cup)-theta)/theta
perc_bias_star<- 100*(mean(theta_star)-theta)/theta

perc_mse_tilde<- 100*((mean(theta_tilde)-theta)^2 + var(theta_tilde))/theta^2
perc_mse_hat<- 100*((mean(theta_hat)-theta)^2 + var(theta_hat))/theta^2
perc_mse_cup<- 100*((mean(theta_cup)-theta)^2 + var(theta_cup))/theta^2
perc_mse_star<- 100*((mean(theta_star)-theta)^2 + var(theta_star))/theta^2

c(perc_bias_tilde,perc_bias_hat,perc_bias_cup, perc_bias_star)
c(perc_mse_tilde,perc_mse_hat,perc_mse_cup,perc_mse_star)
c(theta, n, nrep)

par(mfrow=c(4,1))
hist(theta_tilde)
hist(theta_hat)
hist(theta_cup)
hist(theta_star)

# END OF THE BASIC SIMULATION EXPERIMENT
#############################################

# plot the bias function as a function of theta:
dev.off()
t<- seq(1, 600)/100
bias_func<- -t*(t+2)/(2*n*exp(t)*(1+t)^2)
plot(t,bias_func, type="l", xlab=expression(paste(theta)), ylab="Bias", main=c(expression(paste("Figure 1: First-Order Bias of")~ tilde(theta))))

# Plot the derivative of the bias function:
d<- seq(0,100)/100
dbias<- -2*(1+d)^2+d*(2+d)*(3+d)
plot (d, dbias)

# Where is the bias function's turning point?
f<- function(x) {
 db<-    x^3+3*x^2+2*x-2  
}
uniroot(f,lower=0.4, upper=0.6)   # The turning point is at theta = 0.5213

#############################################
# Generate the Bell numbers:

B0=1
B<- c(1,2,5,15,52,203,877,4140,21147,115975,678570,4213597,27644437, 190899322,1266867570)  # Note B[14] and B[15] calculated as below     

B3<- B0+choose(2,1)*B[1]+choose(2,2)*B[2]
BB15<- B0+B[1]*choose(14,1)+choose(14,2)*B[2]+choose(14,3)*B[3]+choose(14,4)*B[4]+choose(14,5)*B[5]+choose(14,6)*B[6]+choose(14,7)*B[7]+choose(14,8)*B[8]+choose(14,9)*B[9]+choose(1,10)*B[10]+choose(14,11)*B[11]+choose(14,12)*B[12]+choose(14,13)*B[13]+B[14]

########################################
# Plot the p.m.f.'s for various values of theta:

prob<- matrix(nrow=12, ncol=15)
Beta<- c(1,1,2,5,15,52,203,877,4140,21147,115975,678570,4213597,27644437,190899322,1266867570)    # starts at B0
# Theta will take values 0.75 to 3

for (i in 3:12) {
theta<- i/4
for (j in 0:14) {
prob[i,j+1]<- theta^j*exp(-exp(theta)+1)*Beta[j+1]/factorial(j)    # The "Beta" vector starts at B0, not B1

} 
} 

# Each column of p gives us the p.m.f. for a particular value of theta, starting at 0.25 (so the third is for 0.75)
dev.off()

plot(prob[3,], type="o", col="blue", pch="o", lty=1, ylim=c(0,0.35),xlim=c(0,15), xaxt="n",xlab="y", ylab="Pr.(Y=y)", main = "Bell Distribution p.m.f." )
axis(1, at = c(-2,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
     labels = c("0", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9","10","11","12","13","14"))

points(prob[4,], col="red", pch="*")       # theta=1.0
 lines(prob[4,], col="red",lty=2)

points(prob[5,], col="dark red",pch="+")   # theta=1.25
 lines(prob[5,], col="dark red", lty=3)

points(prob[6,], col="green",pch="o")      # theta=1.5
 lines(prob[6,], col="green", lty=3)

points(prob[8,], col="black",pch="+")      # theta=2.0
 lines(prob[8,], col="black", lty=3)

legend(10,0.325,legend=c("theta=0.75","theta=1.0","theta=1.25", "theta=1.5", "theta=2.0"), box.col="white", col=c("blue","red","dark red", "green", "black"),
   pch=c("o","*","+", "o","+"),lty=c(1,2,3,3), ncol=1)


############################################ 
# NOW LOOK AT PROPERTIES OF MLE'S OF THE MEAN, MU = THETA*EXP(THETA)
#################################################################

# FIRST, JUST TRANSFORM TO MU AFTER THETA IS ESTIMATED
# LATER - CONSIDER RE-PARAMETERIZING IN TERMS OF MU ITSELF

library(bellreg)           # Used to generate random "Bell" data
library(lamW)              # Used for the Lambert function to get MLE

set.seed(123)
nrep<- 50000
n<- 10
theta<- 1.0
mu<- theta*exp(theta)

theta_tilde<- vector()            # MLE
mu_tilde<- vector()
theta_hat<- vector()          # Cox-Snell
mu_hat<- vector()
theta_cup<- vector()            # Firth
mu_cup<- vector()
theta_star<- vector()           # Godwin-Giles
mu_star<- vector()

fun_firth<- function(x,dat,n) {
n*(dat-x*exp(x))/x -exp(x)*x*(x+2)/(2*(dat+x^2*exp(x)))
}                                     # Used to get Firth estimator

fun_gg<- function(x,n,mle) {
x-x*(x+2)/(2*n*exp(x)*(1+x)^2)-mle
}                                     # Used to get Godwin-Giles estimator

for (i in 1:nrep) {
y<- rbell(n,theta)
ybar<-mean(y)
theta_tilde[i]<- lambertW0(ybar)    # Data are positive, so only the principal branch of the Lambert function exits.
est_bias<- -theta_tilde[i]*(theta_tilde[i]+2)/(2*n*exp(theta_tilde[i])*(1+theta_tilde[i])^2)
theta_hat[i]<- theta_tilde[1]-est_bias
theta_cup[i]<-uniroot(fun_firth, dat=ybar,n=n, lower=0.0001, upper=theta+1)$root     # Upper limit set acccording to true theta value
theta_star[i]<- uniroot(fun_gg, n=n, mle=theta_tilde[i], lower=0.0001, upper=theta+1)$root
mu_tilde[i]<- theta_tilde[i]*exp(theta_tilde[i])
mu_hat[i]<- theta_hat[i]*exp(theta_hat[i])
mu_cup[i]<- theta_cup[i]*exp(theta_cup[i])
mu_star[i]<- theta_star[i]*exp(theta_star[i])

}

summary(mu_tilde)
summary(mu_hat)
summary(mu_cup)
summary(mu_star)

mu
perc_bias_tilde<- 100*(mean(mu_tilde)-mu)/mu
perc_bias_hat<- 100*(mean(mu_hat)-mu)/mu
perc_bias_cup<- 100*(mean(mu_cup)-mu)/mu
perc_bias_star<- 100*(mean(mu_star)-mu)/mu

perc_mse_tilde<- 100*((mean(mu_tilde)-mu)^2 + var(mu_tilde))/mu^2
perc_mse_hat<- 100*((mean(mu_hat)-mu)^2 + var(mu_hat))/mu^2
perc_mse_cup<- 100*((mean(mu_cup)-mu)^2 + var(mu_cup))/mu^2
perc_mse_star<- 100*((mean(mu_star)-mu)^2 + var(mu_star))/mu^2

c(perc_bias_tilde,perc_bias_hat,perc_bias_cup, perc_bias_star)
c(perc_mse_tilde,perc_mse_hat,perc_mse_cup,perc_mse_star)
c(theta, n, nrep)

par(mfrow=c(4,1))
hist(mu_tilde)
hist(mu_hat)
hist(mu_cup)
hist(mu_star)

# END OF THE BASIC SIMULATION EXPERIMENT FOR THE MEAN
####################################################
#
# NOW RE-PARAMETERIZE IN TERMS OF MU

# The MLE for Mu is ybar

library(bellreg)           # Used to generate random "Bell" data
library(lamW)              # Used for the Lambert function to get MLE

set.seed(123)
nrep<- 5000
n<- 10
theta<- 1.0
mu<- theta*exp(theta)

mu_tilde<- vector()
mu_hat<- vector()
mu_cup<- vector()
mu_star<- vector()

#     the 2 functions need to be fixed up for mu
#fun_firth<- function(x,dat,n) {
#n*(dat-x*exp(x))/x -exp(x)*x*(x+2)/(2*(dat+x^2*exp(x)))
#}                                     # Used to get Firth estimator

#fun_gg<- function(x,n,mle) {
#x-x*(x+2)/(2*n*exp(x)*(1+x)^2)-mle
#}                                     # Used to get Godwin-Giles estimator

for (i in 1:nrep) {

y<- rbell(n,theta)
ybar<-mean(y)
mu_tilde[i]<- ybar
W<- lambertW0(mu_tilde)
        # Note that the kappa terms have been evaluated at the MLE
k11<- -n*(mu_tilde[i]*(W^2+3*W+1)-W^2*exp(W)) / (mu_tilde[i]^2*(1+W)^3)
T1<-  ( mu_tilde*(1+W)*(W^2+3*W+1)+mu_tilde*(2*W^2+3*W-2*W^2*exp(W)-W^3*exp(W)))/(mu_tilde^3*(1+W)^4)
T2<- -(mu_tilde*(W^2+3*W+1) -W^2*exp(W))*(2*W^2+W+2) / (mu_tilde^3*(1+W)^5)
k11_1<- -n*(T1+T2)
T3<-   (mu_tilde*W*(2*W+3) - W^2*exp(W)*(W+2) ) /(mu_tilde^3*(1+W)^4)
T4<-  - (mu_tilde*(W^2+3*W+1) - W^2*exp(W) )   / (mu_tilde^3*(1+W)^4 )
k111<- -n*(T3+T4)
a11<- k11_1-0.5*k111
bias_est<- a11/k11^2
mu_hat[i]<- mu_tilde-bias_est
#mu_cup[i]<- 
mu_star[i]<-  uniroot(fun_gg, n=n, mle=mu_tilde[i], lower=0.0001, upper=theta+1)$root

}

summary(mu_tilde)
summary(mu_hat)

mu
perc_bias_tilde<- 100*(mean(mu_tilde)-mu)/mu
perc_bias_hat<- 100*(mean(mu_hat)-mu)/mu
perc_mse_tilde<- 100*((mean(mu_tilde)-mu)^2 + var(mu_tilde))/mu^2
perc_mse_hat<- 100*((mean(mu_hat)-mu)^2 + var(mu_hat))/mu^2

c(perc_bias_tilde, perc_bias_hat)
c(perc_mse_tilde, perc_mse_hat)

hist(mu_tilde)
hist(mu_hat)

########################################
#######
# plot the bias function as a function of theta:

t<- seq(1, 600)/100
bias_func<- -t*(t+2)/(2*n*exp(t)*(1+t)^2)
plot(t,bias_func, type="l", xlab=expression(paste(theta)), ylab="Bias", main=c(expression(paste("Figure 1: First-Order Bias of")~ tilde(theta))))

# Plot the derivative of the bias function:
d<- seq(0,100)/100
dbias<- -2*(1+d)^2+d*(2+d)*(3+d)
plot (d, dbias)

# Where is the bias function's turning point?
f<- function(x) {
 db<-    x^3+3*x^2+2*x-2  
}
uniroot(f,lower=0.4, upper=0.6)   # The turning point is at theta = 0.5213


