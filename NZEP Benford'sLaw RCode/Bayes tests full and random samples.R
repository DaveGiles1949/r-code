# Bayesian Test - Fonesca
# -----------------------

library(benford.analysis)
library(stringr)
library(readxl)

coef<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/Digits_All.xlsx",sheet="coef")
coef<- coef[!is.na(coef)]
summary(coef)
se<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/Digits_All.xlsx",sheet="se")
se<- se[!is.na(se)]
summary(se)
tstat<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/Digits_All.xlsx",sheet="tstat")
tstat<- tstat[!is.na(tstat)]
summary(tstat)

xc<- coef     
xs<- se       
xt<- tstat    
nc<- length(xc)
ns<- length(xs)
nt<- length(xt)
#######################
# First-digit analysis
# --------------------

# Full sample
# -----------

pi0<- 0.5        # prior probability for null model (Benford's Law)

bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p1<- bfdc$bfd$benford.dist     # Benford first-digit distribution
rN_c1<- bfdc$bfd$data.dist     # these values are "r/N", or "phat"
rN_s1<- bfds$bfd$data.dist
rN_t1<- bfdt$bfd$data.dist

n<- nc          # Choose between coef, se, and t-stats
#n<- ns
#n<- nt
x<- rN_c1*n     # Number of DIGIT COUNTS. 
#x<- rN_s1*n
#x<- rN_t1*n

theta<- p1
alpha_1<- rep(1,9)
alpha_2<- rep(1/9,9)
alpha_3<-22*theta       # See Fonesca, p.23 
alpha_4<- p1      
alpha<- alpha_3         # decide which alpha vector to use in the prior pdf

# Compute the log of the Bayes factor based on Multinomial-Dirichlet prior density

term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))		# the log-gamma function uses natural logs
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01<- term1+term2+term3-term4-term5
logB01<- logB01/2.303       # convert to base 10 log
logB01
post_prob<- 1/(1+(1-pi0)/(pi0*10^logB01))    # Posterior probability in favour of B-L (H0)
post_prob

# Random samples:
#---------------

set.seed(321)

pi0<- 0.5        # prior probability for null model (Benford's Law)

rn1<- 500
rn2<- 1000		# sample sizes for random samples
rn3<- 2000        # Note: coefs, se's and tstats are sampled independently of each other
rand_xc_1<- sample(xc,rn1)
rand_xs_1<- sample(xs,rn1)
rand_xt_1<- sample(xt,rn1)
rand_xc_2<- sample(xc,rn2)
rand_xs_2<- sample(xs,rn2)
rand_xt_2<- sample(xt,rn2)
rand_xc_3<- sample(xc,rn3)
rand_xs_3<- sample(xs,rn3)
rand_xt_3<- sample(xt,rn3)

bfdc1<- benford(rand_xc_1,1,sign="both")
bfds1<- benford(rand_xs_1,1, sign="both")
bfdt1<- benford(rand_xt_1,1, sign="both")
bfdc2<- benford(rand_xc_2,1,sign="both")
bfds2<- benford(rand_xs_2,1, sign="both")
bfdt2<- benford(rand_xt_2,1, sign="both")
bfdc3<- benford(rand_xc_3,1,sign="both")
bfds3<- benford(rand_xs_3,1, sign="both")
bfdt3<- benford(rand_xt_3,1, sign="both")

rN_c1<- bfdc1$bfd$data.dist   # these values are "r/N", or "phat"
rN_s1<- bfds1$bfd$data.dist
rN_t1<- bfdt1$bfd$data.dist
rN_c2<- bfdc2$bfd$data.dist   
rN_s2<- bfds2$bfd$data.dist
rN_t2<- bfdt2$bfd$data.dist
rN_c3<- bfdc3$bfd$data.dist  
rN_s3<- bfds3$bfd$data.dist
rN_t3<- bfdt3$bfd$data.dist

n<- rn1           # choose size of random sample
#n<- rn2
#n<- rn3
x<- rN_c1*n		# choose random data to be analysed - coef, se or tstst
#x<- rN_s1*n      # digit in name refers to sample size above
#x<- rN_t1*n
#x<- rN_c2*n
#x<- rN_s2*n
#x<- rN_t2*n
#x<- rN_c3*n
#x<- rN_s3*n
#x<- rN_t3*n

alpha<- alpha_3        # decide which alpha vector to use in the prior pdf
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01<- term1+term2+term3-term4-term5
logB01<- logB01/2.303   # convert to base 10 log
logB01
post_prob<- 1/(1+(1-pi0)/(pi0*10^(logB01)))    # Posterior probability in favour of B-L (H0)
post_prob

######################################
# Second-digit analysis
# --------------------
coef_n<- c()    # number of second digits
se_n<- c()
tstat_n<- c()
x<- c()

p2<- c(0.1197,0.1139,0.1088,0.1043,0.1,0.0967, 0.0934,0.0904,0.0876,0.0850)

pi0<- 0.5     # prior probability for null model (B-L)
theta<- p2
alpha_1<- rep(1,10)
alpha_2<- rep(0.1,10)
alpha_3<-12*theta       # See Fonesca, p.23
alpha_4<- p2

# Full sample of data:
# -------------------
coef_ch <- str_sub(as.character(abs(xc)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
se_ch<- str_sub(as.character(xs), 2, 2)
tstat_ch<- str_sub(as.character(abs(xt)), 2, 2)
# Convert back to numeric
coef_f <- as.numeric(coef_ch)
se_f<- as.numeric(se_ch)
tstat_f<- as.numeric(tstat_ch)
# Remove the "NAs"
coef_f <- coef_f[!is.na(coef_f)]
se_f<- se_f[!is.na(se_f)]
tstat_f<- tstat_f[!is.na(tstat_f)]

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_n[ii] <- sum(coef_f == i1)
se_n[ii] <- sum(se_f == i1)
tstat_n[ii] <- sum(tstat_f == i1)
}

pi0<- 0.5        # prior probability for null model (Benford's Law)

#n<- length(coef_f)
#n<- length(se_f)
n<- length(tstat_f)
#x<- coef_n     # choose between coefs, se's and tstats
#x<- se_n
x<- tstat_n

# Compute the log of the Bayes factor based on Multinomial-Dirichlet prior density

alpha<- alpha_4		# decide which alpha vector to use in the prior pdf
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01<- term1+term2+term3-term4-term5
logB01<- logB01/2.303   # convert to base 10 log
logB01
post_prob<- 1/(1+(1-pi0)/(pi0*10^(logB01)))    # Posterior probability in favour of B-L (H0)
post_prob

##############

# 1st random sample                # sampling from the second digits
# -----------------
coef_1<- sample(coef_f,rn1)
se_1<- sample(se_f,rn1)
tstat_1<- sample(tstat_f,rn1)

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_n[ii] <- sum(coef_1 == i1)
se_n[ii] <- sum(se_1 == i1)
tstat_n[ii] <- sum(tstat_1 == i1)
}

pi0<- 0.5        # prior probability for null model (Benford's Law)
n<- rn1
x<- coef_n		# choose random data to be analysed
x<- se_n
x<- tstat_n

alpha<- alpha_3        # decide which alpha vector to use in the prior pdf
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01<- term1+term2+term3-term4-term5
logB01<- logB01/2.303
logB01
post_prob<- 1/(1+(1-pi0)/(pi0*10^(logB01)))    # Posterior probability in favour of B-L (H0)
post_prob

# 2nd random sample
# -----------------

coef_2<- sample(coef_f,rn2)
se_2<- sample(se_f,rn2)
tstat_2<- sample(tstat_f,rn2)

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_n[ii] <- sum(coef_2 == i1)
se_n[ii] <- sum(se_2 == i1)
tstat_n[ii] <- sum(tstat_2 == i1)
}

pi0<- 0.5        # prior probability for null model (Benford's Law)
n<- rn2
#x<- coef_n		# choose random data to be analysed
#x<- se_n
x<- tstat_n

alpha<- alpha_3        # decide which alpha vector to use in the prior pdf
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01<- term1+term2+term3-term4-term5
logB01<- logB01/2.303
logB01
post_prob<- 1/(1+(1-pi0)/(pi0*10^(logB01)))    # Posterior probability in favour of B-L (H0)
post_prob

####################################################