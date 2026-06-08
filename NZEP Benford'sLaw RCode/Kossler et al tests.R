# Kossler et al. Tests
#---------------------
#############################################
#
# Code written by David Giles for the paper "Benford’s Law and Regression Results Published in
# Articles in New Zealand Economic Papers", last updated June 2026.
#
# Contact: David Giles; davegiles1949@gmail.com; davegiles.ca
# -------
############################################################

library(benford.analysis)
library(stringr)
library(readxl)

set.seed(123)

coef<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/Digits_All.xlsx",sheet="coef")
coef<- coef[!is.na(coef)]
summary(coef)
se<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/Digits_All.xlsx",sheet="se")
se<- se[!is.na(se)]
summary(se)
tstat<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/Digits_All.xlsx",sheet="tstat")
tstat<- tstat[!is.na(tstat)]
summary(tstat)


# First digit Analysis
#====================
xc<- coef
xs<- se
xt<- tstat
nc<- length(xc)
ns<- length(xs)
nt<- length(xt)
bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p1<- bfdc$bfd$benford.dist     
rN_c1<- bfdc$bfd$data.dist   # these values are "r/N", or "phat"
rN_s1<- bfds$bfd$data.dist
rN_t1<- bfdt$bfd$data.dist

order_of_magnitude <- function(x){
  if (x==0){
    return(0)
  }
  else if (x< 0){
    x = -1 * x
  }
  return(floor(log10(x)))
}

order_c<- c()
order_s<- c()
order_t<- c()

Mcrit10<- qchisq(0.9,9,ncp=0, log.p = FALSE)
Mcrit5<- qchisq(0.95,9,ncp=0, log.p = FALSE)
Mcrit1<- qchisq(0.99,9,ncp=0, log.p = FALSE)

c(Mcrit10, Mcrit5, Mcrit1)

df1<- 7.84619
beta10<- 0.0780258
beta11<- 1.13711

Ecrit10<- beta10+beta11*qchisq(0.9,df1,ncp=0, log.p = FALSE)
Ecrit5<- beta10+beta11*qchisq(0.95,df1,ncp=0, log.p = FALSE)
Ecrit1<- beta10+beta11*qchisq(0.99,df1,ncp=0, log.p = FALSE)
c(Ecrit10, Ecrit5, Ecrit1)

d1<- 1:9
E1<- 1/rep(log(10),9)        # Expected value of S1(S)
V1<- (2*d1+1)/(2*log(10))-(1/log(10))^2       # variance of S1(S)
cov1<-  -1/(log(10))^2		# covariance between any S1(S) terms
CORR1<- matrix(nrow=9,ncol=9)
  for (j in 1:9){
for (i in 1:9){
CORR1[i,j]<- cov1/sqrt(V1[i]*V1[j])
} 
CORR1[j,j]<- 1                       # correlation matrix of Sum - matches Table 10
}

# Coefficients
# -------------

MAD_c<- sqrt(nc)*sum(abs(p1-rN_c1))   # Kossler et al
MAD_c

# Asy. critical values: 2.869 3.084 3.485   (90, 95, 99) # Table 4 K et al. (2024)

# Recall that we have both positive and negative coeffs.

x<- abs(xc)
n<- length(x)
for ( i in 1:n) {
order_c[i]<- order_of_magnitude(x[i])
}
S1<- x/(10^order_c)
fd<-signifd(x = x, digits = 1)
Sum1<- c()

for (j in 1:9){
Sum1[j]<- sum(S1[fd==j])  
} 

R1<- (Sum1 - n*E1)/sqrt(n*V1)       # Sum is scaled down, the expected value and variance of S1(S) are used
IS_1E<- t(R1)%*%R1

IS_1M<- t(R1)%*%solve(CORR1)%*%R1
IS_1E
c(Ecrit10, Ecrit5, Ecrit1)
IS_1M
c(Mcrit10, Mcrit5, Mcrit1)
pchisq(IS_1M,9, ncp = 0, lower.tail = FALSE, log.p = FALSE)


# Std. Errors
# -----------

MAD_s<- sqrt(ns)*sum(abs(p1-rN_s1))   # Kossler et al
MAD_s
# Asy. critical values: 2.869 3.084 3.485   (90, 95, 99) # Table 4 K et al. (2024)

x<- xs
n<- length(x)

for ( i in 1:n) {
order_s[i]<- order_of_magnitude(x[i])
}
S1<- x/(10^order_s)
fd<-signifd(x = x, digits = 1)
Sum1<- c()

for (j in 1:9){
Sum1[j]<- sum(S1[fd==j])  
                  
}

R1<- (Sum1 - n*E1)/sqrt(n*V1)       # Sum is scaled down, the expected value and variance of S1(S) are used
IS_1E<- t(R1)%*%R1

IS_1M<- t(R1)%*%solve(CORR1)%*%R1
IS_1E
c(Ecrit10, Ecrit5, Ecrit1)
IS_1M
c(Mcrit10, Mcrit5, Mcrit1)
pchisq(IS_1M,9, ncp = 0, lower.tail = FALSE, log.p = FALSE)

# t-Statistics
# ------------

MAD_t<- sqrt(nt)*sum(abs(p1-rN_t1))   # Kossler et al
MAD_t
# Asy. critical values: 2.869 3.084 3.485   (90, 95, 99)

x<- abs(xt)
n<- length(x)

for ( i in 1:n) {
order_t[i]<- order_of_magnitude(x[i])
}
S1<- x/(10^order_t)

fd<-signifd(x = x, digits = 1)
Sum1<- c()

for (j in 1:9){
Sum1[j]<- sum(S1[fd==j])  
} 

R1<- (Sum1 - n*E1)/sqrt(n*V1)       # Sum is scaled down, the expected value and variance of S1(S) are used
IS_1E<- t(R1)%*%R1

IS_1M<- t(R1)%*%solve(CORR1)%*%R1
IS_1E
c(Ecrit10, Ecrit5, Ecrit1)
IS_1M

c(Mcrit10, Mcrit5, Mcrit1)
pchisq(IS_1M,9, ncp = 0, lower.tail = FALSE, log.p = FALSE)

#############################################################################
#############################################################################
# Second digit tests
###################

digits<- c(0,1,2,3,4,5,6,7,8,9)
coef_2_ch <- str_sub(as.character(abs(coef)), 2, 2)   # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
# Convert back to numeric
coef_2 <- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

xc<- coef_2
xs<- se_2
xt<- tstat_2
nc<- length(xc)
ns<- length(xs)
nt<-length(xt)

# Benford 2nd-digit probs. (for i = 0 to 9)
p2<- c(0.1197,0.1139,0.1088,0.1043,0.1,0.0967, 0.0934,0.0904,0.0876,0.0850)

rN_c2<- c()
rN_s2<- c()
rN_t2<- c()
phat_c2<- c()
phat_s2<- c()
phat_t2<- c()

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
rN_c2[ii] <- sum(coef_2 == i1)/nc
rN_s2[ii] <- sum(se_2 == i1)/ns
rN_t2[ii] <- sum(tstat_2 == i1)/nt
}

Mcrit10<- qchisq(0.9,10,ncp=0, log.p = FALSE)
Mcrit5<- qchisq(0.95,10,ncp=0, log.p = FALSE)
Mcrit1<- qchisq(0.99,10,ncp=0, log.p = FALSE)

c(Mcrit10, Mcrit5, Mcrit1)

df2<- 9
beta20<- 0
beta21<- 1.77527      # Kossler et al. (2024), p.3795

Ecrit10<- beta20+beta21*qchisq(0.9,df2,ncp=0, log.p = FALSE)
Ecrit5<- beta20+beta21*qchisq(0.95,df2,ncp=0, log.p = FALSE)
Ecrit1<- beta20+beta21*qchisq(0.99,df2,ncp=0, log.p = FALSE)

c(Ecrit10, Ecrit5, Ecrit1)

d2<- 0:9
E2<- 9/rep(10*log(10),10)        # Expected value of S2(S)
V2<- 9*(2*d2+101)/(200*log(10))-E2^2       # variance of S2(S)
cov2<-  -81/(10*log(10))^2		# covariance between any S2(S) terms
CORR2<- matrix(nrow=10,ncol=10)

for (j in 1:10){
for (i in 1:10){
CORR2[i,j]<- cov2/sqrt(V2[i]*V2[j])
} 
CORR2[j,j]<- 1                       # correlation matrix of Sum 
}


# Coefficients
# ------------
MAD_c<- sqrt(nc)*sum(abs(p2-rN_c2))   # Kossler et al
MAD_c

# MAD Asy. critical values: 3.18 3.42 3.92   (90, 95, 99) # Table 3 K et al. (2024)
x<- c()
x<- abs(coef)       # need absolute value
#x<- c(2,5,1,1,8,7,9,0,0,4,3,5,6,6)
n<- length(x)
order_c<- c()
for ( i in 1:n) {
order_c[i]<- order_of_magnitude(x[i])
}

S2<- x/(10^order_c)
Sum2<- c()

for (j in 1:10){
jj<- j-1
Sum2[j]<- sum(S2[coef_2==jj])  
                  
}

R2<- (Sum2 - n*E2)/sqrt(n*V2)       # Sum is scaled down, the expected value and variance of S1(S) are used
IS_2E<- t(R2)%*%R2

IS_2M<- t(R2)%*%solve(CORR2)%*%R2
IS_2E
c(Ecrit10, Ecrit5, Ecrit1)
IS_2M
c(Mcrit10, Mcrit5, Mcrit1)
pchisq(IS_2M,10, ncp = 0, lower.tail = FALSE, log.p = FALSE)


# Std. Errors
# -----------
MAD_s<- sqrt(ns)*sum(abs(p2-rN_s2))   # Kossler et al
MAD_s

# MAD Asy. critical values: 3.18 3.42 3.92   (90, 95, 99) # Table 3 of Kossler et al. (2024)
x<- c()
x<- se  
n<- length(x)
order_c<- c()
for ( i in 1:n) {
order_c[i]<- order_of_magnitude(x[i])
}

S2<- x/(10^order_c)
Sum2<- c()

for (j in 1:10){
jj<- j-1
Sum2[j]<- sum(S2[se_2==jj])  
                  
}

R2<- (Sum2 - n*E2)/sqrt(n*V2)       # Sum is scaled down, the expected value and variance of S1(S) are used
IS_2E<- t(R2)%*%R2

IS_2M<- t(R2)%*%solve(CORR2)%*%R2
IS_2E
c(Ecrit10, Ecrit5, Ecrit1)
IS_2M
c(Mcrit10, Mcrit5, Mcrit1)
pchisq(IS_2M,10, ncp = 0, lower.tail = FALSE, log.p = FALSE)

# t-Statistics
# ------------
MAD_t<- sqrt(nt)*sum(abs(p2-rN_t2))   # Kossler et al
MAD_t

# MAD Asy. critical values: 3.18 3.42 3.92   (90, 95, 99) # Table 3 K et al. (2024)
x<- c()
x<- abs(tstat)   # abs value needed  
n<- length(x)
order_c<- c()
for ( i in 1:n) {
order_c[i]<- order_of_magnitude(x[i])
}

S2<- x/(10^order_c)
Sum2<- c()

for (j in 1:10){
jj<- j-1
Sum2[j]<- sum(S2[tstat_2==jj])  
                  
}

R2<- (Sum2 - n*E2)/sqrt(n*V2)       # Sum is scaled down, the expected value and variance of S1(S) are used
IS_2E<- t(R2)%*%R2

IS_2M<- t(R2)%*%solve(CORR2)%*%R2
IS_2E
c(Ecrit10, Ecrit5, Ecrit1)
IS_2M
c(Mcrit10, Mcrit5, Mcrit1)
pchisq(IS_2M,10, ncp = 0, lower.tail = FALSE, log.p = FALSE)


