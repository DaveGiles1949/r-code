# David Giles, April 2023
# -------------------------

library(VGAM)            # For computing Riemann zeta function & its first 2 derivatives
library (tolerance)      # For computing the MLE of the zeta distribution parameter using the 'zm.ll' function
library(formattable)      # Used to display output          

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

#Predictions based on estimated zeta distributions
predtilde<- c()
predhat<- c()
predcup<- c()

groups<- 16
ones<- c(1,1,1,1,1,1,1,1,1,1)
twos<- c(2,2,2,2,2,2,2,2,2,2)
threes<- c(3,3,3,3,3,3,3,3,3,3)
fours<- c(4,4,4,4,4,4,4,4,4,4)

#Seal's data, 1952, Table 1, 17.5 year age-group
#zeta.data<- c(ones,ones,ones,1,1,1,1,1,1,2,2,2)
#actual<-c(36,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#Seal's data, 1952, Table 1, 22.5 year age-group
#zeta.data<- c(ones,ones,ones,ones,ones,1,1,1,1,1,1,1,1,2,2,2)
#actual<-c(58,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#Seal's data, 1952, Table 1, 27.5 year age-group
#zeta.data<- c(ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,1,2,2,2,2,2,2,2,2,3,3)
#actual<-c(101,8,2,0,0,0,0,0,0,0,0,0,0,0,0,0)

#Seal's data, 1952, Table 1, 32.5 year age-group
#zeta.data<- c(ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,1,twos,twos,2,2,2,2,2,2,3,3,3,4,4,4,6,6)
#actual<-c(241,26,3,3,0,2,0,0,0,0,0,0,0,0,0,0)

#Seal's data, 1952, Table 1, 37.5 year age-group
#zeta.data<- c(ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,1,1,1,twos,twos,twos,2,2,2,2,2,3,3,3,3,3,4,4,6,6,8,9)
#actual<-c(283,35,4,2,0,2,0,1,1,0,0,0,0,0,0,0)

#Seal's data, 1952, Table 1, 42.5 year age-group
#zeta.data<- c(ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,1,1,1,1,1,1,1,twos,twos,twos,2,threes,3,3,4,4,4,4,4,4,6,7,7,11,11)
#actual<-c(307,31,12,6,0,1,2,0,0,0,2,0,0,0,0,0)

#Seal's data, 1952, Table 1, 47.5 year age-group
#zeta.data<- c(ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,1,1,1,twos,twos,twos,2,2,2,2,2,threes,3,3,3,3,4,4,4,4,4,5,5,5,5,5,7,13)
#actual<-c(233,35,4,5,3,0,1,0,0,0,0,0,1,0,0,0)

#Seal's data, 1952, Table 1, 52.5 year age-group
#zeta.data<- c(ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,1,1,1,1,1,1,1,1,1,twos,twos,2,2,2,2,3,3,3,3,3,3,3,4,4,6,6,6,7,7,10)
#actual<-c(209,24,7,2,0,3,1,1,0,1,0,0,0,0,0,0)

#Seal's data, 1952, Table 1, 57.5 year age-group
#zeta.data<- c(ones,ones,ones,ones,ones,ones,ones,ones,ones,ones,1,1,1,1,1,1,1,1,twos,twos,3,3,3,3,3,4,4,5,5,5,5,8)
#actual<-c(108,20,5,2,4,0,0,1,0,0,0,0,0,0,0,0)

#Seal's data, 1952, Table 1, 62.5 year age-group
#zeta.data<- c(ones,ones,ones,ones,ones,ones,1,1,1,1,1,1,1,1,1,twos,3)
#actual<-c(69,10,1,0,0,0,0,0,0,0,0,0,0,0,0,0)

#Seal's data, 1952, Table 1, 67.5 year age-group
#zeta.data<- c(ones,ones,ones,1,1,1,2,2,2,2,2,2,2,3,3,3,3,4)
#actual<- c(33,7,4,1,0,0,0,0,0,0,0,0,0,0,0,0)

#Seal's data, 1952, Table 1, 72.5 year age-group
zeta.data<- c(ones,ones,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,5,5)
actual<- c(26,5,4,1,2,0,0,0,0,0,0,0,0,0,0,0)

# All data for all years
#groups<- 16     # Again, j=18 is omitted because we can't have zero values in zeta.data
#zeta.data<- c(rep(ones,169),1,1,1,1,1, rep(twos,20),2,2,2,2,2,2,2,threes, threes, threes, threes,3,3,3,3,3,3,fours,fours,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,7,7,7,7,8,8,8,9,10,11,11,13)
#actual<- c(1695,207,46,22,9,8,4,3,1,1,2,0,1,0,0,0)

N<- sum(actual)
r<- actual

sumlogdata<- sum(log(zeta.data))     # We need this for Firth's estimator

# Obtain the MLE of "s":

out.zeta<- zm.ll(zeta.data, N = Inf, dist = "Zeta")
stilde<- stats4::coef(out.zeta)[[1]]           # Keep value only, not name ("s") as well. This speeds up the subsequent summing.
se_tilde<- sqrt(stats4::vcov(out.zeta)[[1]])   # Asymptotic std. error of the MLE

zetas<- zeta(stilde,0,1)  # using the VGAM package in this & the next 2 lines
zeta1<- zeta(stilde,1,1)
zeta2<- zeta(stilde,2,1)

# Use Choudhury's equation 20 to get the third derivative:
zeta3<- -(6/(stilde-1)^4+gamma3-gamma4*(stilde-1)+(gamma5*(stilde-1)^2)/2-(gamma6*(stilde-1)^3)/6+(gamma7*(stilde-1)^4)/24-(gamma8*(stilde-1)^5)/120
+(gamma9*(stilde-1)^6)/720-(gamma10*(stilde-1)^7)/5040+(gamma11*(stilde-1)^8)/40320-(gamma12*(stilde-1)^9)/362880+(gamma13*(stilde-1)^10)/3628800
- (gamma14*(stilde-1)^11)/39916800 + (gamma15*(stilde-1)^12)/479001600 - (gamma16*(stilde-1)^13)/6227020800 + (gamma17*(stilde-1)^14)/87178291200 
- (gamma18*(stilde-1)^15)/1307674368000 )
bias<- (zetas/(2*N))*(3*zetas*zeta1*zeta2 - 2*zeta1^3 - zetas^2*zeta3)/((zeta1^2 - zetas*zeta2)^2)   # This is just the second-order bias!
shat<- stilde-bias                     # Compute the Cox-Snell bias-corrected estimate

#Taking account of the fact that we have a canonical linear exponential model -

f<- function (x) -sumlogdata -(N+1)*zeta(x,1,1)/zeta(x,0,1) +
(zeta(x,0,1)*
( -(6/(x-1)^4+gamma3-gamma4*(x-1)+(gamma5*(x-1)^2)/2-(gamma6*(x-1)^3)/6+(gamma7*(x-1)^4)/24-(gamma8*(x-1)^5)/120
+(gamma9*(x-1)^6)/720-(gamma10*(x-1)^7)/5040+(gamma11*(x-1)^8)/40320-(gamma12*(x-1)^9)/362880+(gamma13*(x-1)^10)/3628800
- (gamma14*(x-1)^11)/39916800 + (gamma15*(x-1)^12)/479001600 - (gamma16*(x-1)^13)/6227020800 + (gamma17*(x-1)^14)/87178291200 
- (gamma18*(x-1)^15)/1307674368000 )) 

-zeta(x,1,1)*zeta(x,2,1)) /(2*(zeta(x,0,1)*zeta(x,2,1) - zeta(x,1,1)^2))

scup<- uniroot(f , c(1.1,7), tol = 1e-6)$root

for(ii in 1:16) {
predtilde[ii]<- (ii^-stilde/zeta(stilde,0,1))*N
predhat[ii]<- (ii^-shat/zeta(shat,0,1))*N
predcup[ii]<- (ii^-scup/zeta(scup,0,1))*N
}

c(stilde, shat, scup, se_tilde)
df1<- data.frame(actual, predtilde,predhat,predcup)			# Print results n a tidy format
df1$actual<-formattable(df1$actual,format="f",digits=0)
df1$predtilde<-formattable(df1$predtilde,format="f",digits=2)
df1$predhat<-formattable(df1$predhat,format="f",digits=2)
df1$predcup<-formattable(df1$predcup,format="f",digits=2)

print(df1, row.names = FALSE)

##########################################################

# See my paper about GOF tests with discrete data
S<- cumsum(r/N-predtilde/N)
UN2_tilde<- (N/groups)*( (sum(S^2)-S[groups]^2 )- (sum(S)-S[groups])^2/groups )
S<- cumsum(r/N-predhat/N)
UN2_hat<- (N/groups)*( (sum(S^2)-S[groups]^2 )- (sum(S)-S[groups])^2/groups )
S<- cumsum(r/N-predcup/N)
UN2_cup<- (N/groups)*( (sum(S^2)-S[groups]^2 )- (sum(S)-S[groups])^2/groups )

c(UN2_tilde, UN2_hat, UN2_cup)

# See R file for p-value calculation
# 
# END (not run)
	