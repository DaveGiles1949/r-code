# Freedman's Test for Benford Distribution - from Chilean Journal of Statistics
#
# Reference: Giles, D.E.A., 2013. Exact asymptotic goodness-of-fit testing 
# for discrete circular data, with applications, Chilean Journal of Statistics,
# 4, 19-34.
############################################################################

library(benford.analysis)
library(BenfordTests)
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
bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p1<- bfdc$bfd$benford.dist     
rN_c1<- bfdc$bfd$data.dist   # these values are "r/N"
rN_s1<- bfds$bfd$data.dist
rN_t1<- bfdt$bfd$data.dist

#  First digit tests
####################

# Coefficients
# ------------

sig1<-signifd(x = xc, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_c1[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_c<- (nc/9)*(sum1-sum2^2/9)

# Std. Errors
# -----------

sig1<-signifd(x = xs, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_s1[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_s<- (ns/9)*(sum1-sum2^2/9)

# t-Statistics
# ------------

sig1<-signifd(x = xt, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_t1[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_t<- (nt/9)*(sum1-sum2^2/9)

Unsq_c
Unsq_s
Unsq_t

# 1st-digit Test Critical Values:  90%   95% 97.5%  99%    (from Giles, Table 2)

#                                 0.143 0.178 0.214  0.263
   
# Random sample tests
# -------------------

set.seed(123456)
n<- 500       # choose a sample size
#n<- 1000 
#n<- 2000

xc<- sample(coef,n)
xs<- sample(se,n)
xt<- sample(tstat,n)
bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p1<- bfdc$bfd$benford.dist     
rN_c1<- bfdc$bfd$data.dist   # these values are "r/N"
rN_s1<- bfds$bfd$data.dist
rN_t1<- bfdt$bfd$data.dist

# Coefficients
# ------------

sig1<-signifd(x = xc, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_c1[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_c<- (nc/9)*(sum1-sum2^2/9)

# Std. Errors
# -----------

sig1<-signifd(x = xs, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_s1[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_s<- (ns/9)*(sum1-sum2^2/9)

# t-Statistics
# ------------

sig1<-signifd(x = xt, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_t1[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_t<- (nt/9)*(sum1-sum2^2/9)

Unsq_c
Unsq_s
Unsq_t

# 1st-digit Test Critical Values:  90%   95% 97.5%  99%    (from Giles, Table 2)

#                                 0.143 0.178 0.214  0.263


####################################################
# Second digit tests
####################

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

# Coefficients
# ------------

S<- rep(0,10)
term<- c()
for (i in 1:10) {
term[i]<- rN_c2[i]-p2[i]
}

S<- cumsum(term[1:10])
sum1<- sum((S[1])^2:(S[9])^2) 
sum2<- sum(S[1]:S[9])
Unsq_c<- (nc/10)*(sum1-sum2^2/10)

# Std. Errors
# -----------

S<- rep(0,10)
term<- c()
for (i in 1:10) {
term[i]<- rN_s2[i]-p2[i]
}

S<- cumsum(term[1:10])
sum1<- sum((S[1])^2:(S[9])^2) 
sum2<- sum(S[1]:S[9])
Unsq_s<- (ns/10)*(sum1-sum2^2/10)

# t-Statistics
# ------------

S<- rep(0,10)
term<- c()
for (i in 1:10) {
term[i]<- rN_t2[i]-p2[i]
}

S<- cumsum(term[1:10])
sum1<- sum((S[1])^2:(S[9])^2) 
sum2<- sum(S[1]:S[9])
Unsq_t<- (nt/10)*(sum1-sum2^2/10)

Unsq_c
Unsq_s
Unsq_t

# 2nd-digit Test Critical Values:  90%   95% 97.5%  99%   (from Giles, Table 2)
#                                 0.154 0.190 0.226 0.274

# Random samples
# --------------

set.seed(123456)
n<- 500
#n<- 1000
#n<- 2000
xc<- sample(coef_2,n)
xs<- sample(se_2,n)
xt<- sample(tstat_2,n)
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
rN_c2[ii] <- sum(xc == i1)/nc
rN_s2[ii] <- sum(xs == i1)/ns
rN_t2[ii] <- sum(xt == i1)/nt
}

# Coefficients
# ------------

S<- rep(0,10)
term<- c()
for (i in 1:10) {
term[i]<- rN_c2[i]-p2[i]
}

S<- cumsum(term[1:10])
sum1<- sum((S[1])^2:(S[9])^2) 
sum2<- sum(S[1]:S[9])
Unsq_c<- (nc/10)*(sum1-sum2^2/10)

# Std. Errors
# -----------

S<- rep(0,10)
term<- c()
for (i in 1:10) {
term[i]<- rN_s2[i]-p2[i]
}

S<- cumsum(term[1:10])
sum1<- sum((S[1])^2:(S[9])^2) 
sum2<- sum(S[1]:S[9])
Unsq_s<- (ns/10)*(sum1-sum2^2/10)

# t-Statistics
# ------------

S<- rep(0,10)
term<- c()
for (i in 1:10) {
term[i]<- rN_t2[i]-p2[i]
}

S<- cumsum(term[1:10])
sum1<- sum((S[1])^2:(S[9])^2) 
sum2<- sum(S[1]:S[9])
Unsq_t<- (nt/10)*(sum1-sum2^2/10)

Unsq_c
Unsq_s
Unsq_t

# 2nd-digit Test Critical Values:  90%   95% 97.5%  99%   (from Giles, Table 2)
#                                 0.154 0.190 0.226 0.274
	