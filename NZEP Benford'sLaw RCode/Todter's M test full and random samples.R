# Construct Todter's M-test (p.342): asymptotically chi-square with 1 dof.
###########################
#############################################
#
# Code written by David Giles for the paper "Benford’s Law and Regression Results Published in
# Articles in New Zealand Economic Papers", last updated June 2026.
#
# Contact: David Giles; davegiles1949@gmail.com; davegiles.ca
# -------
############################################################

library(benford.analysis)
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


#FIRST - Tests based on the full sample
########################################

xc<- coef
xs<- se
xt<- tstat
nc<- length(xc)
ns<- length(xs)
nt<-length(xt)

bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p<- bfdc$bfd$benford.dist
phat_c1<- bfdc$bfd$data.dist
phat_s1<- bfds$bfd$data.dist
phat_t1<- bfdt$bfd$data.dist

# First digit tests
###################

d1<- 1:9

dbar<- sum(phat_c1*d1)   # phat's are RELATIVE frequencies
M_c<- nc*(mean(dbar)-3.440)^2/6.057
pval<- pchisq(M_c,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_c
pval
dbar<- sum(phat_s1*d1)
M_s<- ns*(mean(dbar)-3.440)^2/6.057
pval<- pchisq(M_s,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_s
pval
dbar<- sum(phat_t1*d1)
M_t<- nt*(mean(dbar)-3.440)^2/6.057
pval<- pchisq(M_t,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_t
pval

# Second digit tests
###############33333

d2<- 0:9
p2<- benford_2<- c(0.1197,0.1139,0.1088,0.1043,0.1,0.0967, 0.0934,0.0904,0.0876,0.0850)

coef_2_ch <- str_sub(as.character(abs(xc)), 2, 2)   # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(xs), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(xt)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
# Convert back to numeric
coef_2 <- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]
nc<- length(coef_2)
ns<- length(se_2)
nt<- length(tstat_2)
phat_c2<- c()
phat_s2<- c()
phat_t2<- c()

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
phat_c2[ii] <- sum(coef_2 == i1)/nc
phat_s2[ii] <- sum(se_2 == i1)/ns
phat_t2[ii] <- sum(tstat_2 == i1)/nt
}

dbar<- sum(phat_c2*d2)
M_c<- nc*(mean(dbar)-4.1874)^2/8.2538
pval<- pchisq(M_c,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_c
pval
M_s<- ns*(mean(dbar)-4.1874)^2/8.2538
pval<- pchisq(M_s,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_s
pval
M_t<- nt*(mean(dbar)-4.1874)^2/8.2538
pval<- pchisq(M_t,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_t
pval

# SECOND - Tests based on random samples
#######################################

set.seed(123)
randn<- 500                # choose samle sixe

#randn<- 1000
#randn<- 2000

xc<- sample(coef,randn)
xs<- sample(se,randn)
xt<- sample(tstat,randn)
nc<- length(xc)
ns<- length(xs)
nt<-length(xt)

bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p<- bfdc$bfd$benford.dist
phat_c1<- bfdc$bfd$data.dist
phat_s1<- bfds$bfd$data.dist
phat_t1<- bfdt$bfd$data.dist

# First digit tests
###################

d1<- 1:9

dbar<- sum(phat_c1*d1)
M_c<- nc*(mean(dbar)-3.440)^2/6.057
pval<- pchisq(M_c,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_c
pval
dbar<- sum(phat_s1*d1)
M_s<- ns*(mean(dbar)-3.440)^2/6.057
pval<- pchisq(M_s,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_s
pval
dbar<- sum(phat_t1*d1)
M_t<- nt*(mean(dbar)-3.440)^2/6.057
pval<- pchisq(M_t,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_t
pval

# Second digit tests
###############33333

d2<- 0:9
p2<- benford_2<- c(0.1197,0.1139,0.1088,0.1043,0.1,0.0967, 0.0934,0.0904,0.0876,0.0850)

coef_2_ch <- str_sub(as.character(abs(xc)), 2, 2)   # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(xs), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(xt)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
# Convert back to numeric
coef_2 <- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]
nc<- length(coef_2)
ns<- length(se_2)
nt<- length(tstat_2)
phat_c2<- c()
phat_s2<- c()
phat_t2<- c()

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
phat_c2[ii] <- sum(coef_2 == i1)/nc
phat_s2[ii] <- sum(se_2 == i1)/ns
phat_t2[ii] <- sum(tstat_2 == i1)/nt
}

dbar<- sum(phat_c2*d2)
M_c<- nc*(mean(dbar)-4.1874)^2/8.2538
pval<- pchisq(M_c,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_c
pval
M_s<- ns*(mean(dbar)-4.1874)^2/8.2538
pval<- pchisq(M_s,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_s
pval
M_t<- nt*(mean(dbar)-4.1874)^2/8.2538
pval<- pchisq(M_t,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_t
pval


