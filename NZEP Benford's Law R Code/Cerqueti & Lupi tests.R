# Cerqueti and Lupi (2022)
# ----------------
library(benford.analysis)
library(stringr)
library(readxl)
library(Matrix)

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

# First digit analysis
# --------------------

xc<- coef   #sample(coef,1000)
xs<- se   #sample(se,1000)
xt<- tstat   #sample(tstat,1000)
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

k<- 9
MAD_c<- sum(abs(p1-rN_c1))/k
MAD_c
MAD_s<- sum(abs(p1-rN_s1))/k
MAD_s
MAD_t<- sum(abs(p1-rN_t1))/k
MAD_t

# These MAD values are below the conformity criterion of 0.015 for first digit case 
# suggested by Nigrini (2102, p.160)

# See authors, p.11 - delta_star = 0.01 for first-digit case
# The MAD values
# Obtain asy. distribution

i<- rep(1,k)
D<- matrix(0,nrow=k,ncol=k)
R<- matrix(nrow=k,ncol=k)
rho<- matrix(nrow=k,ncol=k)
for (ii in 1:k) {
D[ii,ii]<- sqrt(p1[ii]*(1-p1[ii]))
for (jj in 1:k){
rho[ii,jj]<- -sqrt(p1[ii]*p1[jj]/((1-p1[ii])*(1-p1[jj])))
rho[ii,ii]<- 1
R[ii,jj]<- (2/pi)*(rho[ii,jj]*asin(rho[ii,jj])+sqrt(1-rho[ii,jj]^2))-(2/pi)
}
}
asy_mean<- sqrt(2/(pi*k^2))*t(i)%*%D%*%i
asy_var<- t(i)%*%D%*%R%*%D%*%i/(k^2)
z_c<- (sqrt(nc)*MAD_c -asy_mean)/sqrt(asy_var)
z_c
pval_c<- 1-pnorm(z_c, mean=0, sd=1)
pval_c
z_s<- (sqrt(ns)*MAD_s -asy_mean)/sqrt(asy_var)
z_s
pval_s<- 1-pnorm(z_s, mean=0, sd=1)
pval_s
z_t<- (sqrt(nt)*MAD_t -asy_mean)/sqrt(asy_var)
z_t
pval_t<- 1-pnorm(z_t, mean=0, sd=1)
pval_t

# "Conformity tests":
# error in their eq (15) - "n" missing inside denominator of square root in numerator
# ------------------
delta_tilde_c<- k*sqrt(nc)* (MAD_c-sqrt( 2/(pi*nc*k^2))*t(i)%*%D%*%i )/sqrt(t(i)%*%D%*%R%*%D%*%i)
delta_tilde_c
delta_tilde_s<- k*sqrt(ns)* (MAD_s-sqrt( 2/(pi*ns*k^2))*t(i)%*%D%*%i )/sqrt(t(i)%*%D%*%R%*%D%*%i)
delta_tilde_s
delta_tilde_t<- k*sqrt(nt)* (MAD_t-sqrt( 2/(pi*nt*k^2))*t(i)%*%D%*%i )/sqrt(t(i)%*%D%*%R%*%D%*%i)
delta_tilde_t
# These are the above t-stats
# let delta_star be the discrepancy
#H0: delta=0 vs. H1: delta=delta_star (>0)   1-sided z-test based on above statistics
# under H1: distribution of statistic is N[mean, 1], where mean = k*delta_star*sqrt(n)/sqrt(i'DRDi)

# Now the chi-square: for first-digit, asy distribution is n.c. Chi-Square, ncp=0.0238*n
# ( See authors on p.13)
ncp_c<- nc*0.0238
ncp_s<- ns*0.0238
ncp_t<- nt*0.0238
Q_c1<- nc*sum((rN_c1-p1)^2/p1) 
Q_c1
ncp_c
pval_c<- 1-pchisq(Q_c1, df=k-1, ncp=ncp_c)
pval_c
Q_s1<- ns*sum((rN_s1-p1)^2/p1) 
Q_s1
ncp_s
pval_s<- 1-pchisq(Q_s1, df=k-1, ncp=ncp_s)
pval_s
Q_t1<- nt*sum((rN_t1-p1)^2/p1) 
Q_t1
ncp_t
pval_t<- 1-pchisq(Q_t1, df=k-1, ncp=ncp_t)
pval_t

# cannot reject Benford first digit
# Central chi-square p-values:
pval_c_nc<- 1-pchisq(Q_c1, df=k-1, ncp=0)
pval_s_nc<- 1-pchisq(Q_s1, df=k-1, ncp=0)
pval_t_nc<- 1-pchisq(Q_t1, df=k-1, ncp=0)
pval_c_nc
pval_s_nc
pval_t_nc


# Second digit analysis
# ---------------------
coef<- coef   #sample(coef,1000)
se<- se   #sample(se,1000)
tstat<- tstat    #sample(tstat,1000)

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

xc<- coef_2   #sample(coef_2,1000)
xs<- se_2   #sample(se_2, 1000)
xt<- tstat_2   #sample(tstat_2,1000)
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

k<- 10
MAD_c<- sum(abs(p2-rN_c2))/k
MAD_c
MAD_s<- sum(abs(p2-rN_s2))/k
MAD_s
MAD_t<- sum(abs(p2-rN_t2))/k
MAD_t

# These MAD values are below the conformity criterion of 0.018 for second digit case 
# suggested by Nigrini (2102, p.160)

# Obtain asy. distribution

i<- rep(1,k)
D<- matrix(0,nrow=k,ncol=k)
R<- matrix(nrow=k,ncol=k)
rho<- matrix(nrow=k,ncol=k)
for (ii in 1:k) {
D[ii,ii]<- sqrt(p2[ii]*(1-p2[ii]))
for (jj in 1:k){
rho[ii,jj]<- -sqrt(p2[ii]*p2[jj]/((1-p2[ii])*(1-p2[jj])))
rho[ii,ii]<- 1
R[ii,jj]<- (2/pi)*(rho[ii,jj]*asin(rho[ii,jj])+sqrt(1-rho[ii,jj]^2))-(2/pi)
}
}
asy_mean<- sqrt(2/(pi*k^2))*t(i)%*%D%*%i
asy_var<- t(i)%*%D%*%R%*%D%*%i/(k^2)
z_c<- (sqrt(nc)*MAD_c -asy_mean)/sqrt(asy_var)
z_c
pval_c<- 1-pnorm(z_c, mean=0, sd=1)
pval_c
z_s<- (sqrt(ns)*MAD_s -asy_mean)/sqrt(asy_var)
z_s
pval_s<- 1-pnorm(z_s, mean=0, sd=1)
pval_s
z_t<- (sqrt(nt)*MAD_t -asy_mean)/sqrt(asy_var)
z_t
pval_t<- 1-pnorm(z_t, mean=0, sd=1)
pval_t

# "Conformity tests":
# error in their eq (15) - "n" missing inside denominator of square root in numerator
# ------------------
delta_tilde_c<- k*sqrt(nc)* (MAD_c-sqrt( 2/(pi*nc*k^2))*t(i)%*%D%*%i )/sqrt(t(i)%*%D%*%R%*%D%*%i)
delta_tilde_c
delta_tilde_s<- k*sqrt(ns)* (MAD_s-sqrt( 2/(pi*ns*k^2))*t(i)%*%D%*%i )/sqrt(t(i)%*%D%*%R%*%D%*%i)
delta_tilde_s
delta_tilde_t<- k*sqrt(nt)* (MAD_t-sqrt( 2/(pi*nt*k^2))*t(i)%*%D%*%i )/sqrt(t(i)%*%D%*%R%*%D%*%i)
delta_tilde_t
# These are the above t-stats
# let delta_star be the discrepancy
#H0: delta=0 vs. H1: delta=delta_star (>0)   1-sided z-test based on above statistics
# under H1: distribution of statistic is N[mean, 1], where mean = k*delta_star*sqrt(n)/sqrt(i'DRDi)

# Now the chi-square: for first-digit, asy distribution is n.n. Chi-Square, ncp=0.0474*n
# ( See authors on p.13)

ncp_c<- nc*0.0474
ncp_s<- ns*0.0474
ncp_t<- nt*0.0474
Q_c2<- nc*sum((rN_c2-p2)^2/p2) 
Q_c2
ncp_c
pval_c<- 1-pchisq(Q_c2, df=k-1, ncp=ncp_c)
pval_c
Q_s2<- ns*sum((rN_s2-p2)^2/p2) 
Q_s2
ncp_s
pval_s<- 1-pchisq(Q_s2, df=k-1, ncp=ncp_s)
pval_s
Q_t2<- nt*sum((rN_t2-p2)^2/p2) 
Q_t2
ncp_t
pval_t<- 1-pchisq(Q_t2, df=k-1, ncp=ncp_t)
pval_t

# cannot reject Benford second digit
# Non-central p-values:
pval_c_nc<- 1-pchisq(Q_c2, df=k-1, ncp=0)
pval_s_nc<- 1-pchisq(Q_s2, df=k-1, ncp=0)
pval_t_nc<- 1-pchisq(Q_t2, df=k-1, ncp=0)
pval_c_nc
pval_s_nc
pval_t_nc

