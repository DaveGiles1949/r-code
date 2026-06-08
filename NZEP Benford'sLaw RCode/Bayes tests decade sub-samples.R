# Bayes tests for sub-samples based on decades
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

coef_6070<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="coef", range=cell_cols("A:M"))
coef_6070<- coef_6070[!is.na(coef_6070)]
se_6070<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="se",range=cell_cols("A:J") )
se_6070<- se_6070[!is.na(se_6070)]
tstat_6070<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="tstat",range=cell_cols("A:I") )
tstat_6070<- tstat_6070[!is.na(tstat_6070)]
coef_80<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="coef", range=cell_cols("N:U"))
coef_80<- coef_80[!is.na(coef_80)]
se_80<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="se",range=cell_cols("K:N") )
se_80<- se_80[!is.na(se_80)]
tstat_80<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="tstat",range=cell_cols("J:O") )
tstat_80<- tstat_80[!is.na(tstat_80)]
coef_90<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="coef", range=cell_cols("V:AE"))
coef_90<- coef_90[!is.na(coef_90)]
se_90<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="se",range=cell_cols("O:U") )
se_90<- se_90[!is.na(se_90)]
tstat_90<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="tstat",range=cell_cols("P:Y") )
tstat_90<- tstat_90[!is.na(tstat_90)]
coef_00<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="coef", range=cell_cols("AF:AO"))
coef_00<- coef_00[!is.na(coef_00)]
se_00<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="se",range=cell_cols("V:AD") )
se_00<- se_00[!is.na(se_00)]
tstat_00<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="tstat",range=cell_cols("Z:AF") )
tstat_00<- tstat_00[!is.na(tstat_00)]
coef_10<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="coef", range=cell_cols("AP:AY"))
coef_10<- coef_10[!is.na(coef_10)]
se_10<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="se",range=cell_cols("AE:AN") )
se_10<- se_10[!is.na(se_10)]
tstat_10<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="tstat",range=cell_cols("AG:AO") )
tstat_10<- tstat_10[!is.na(tstat_10)]
coef_20<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="coef", range=cell_cols("AZ:BE"))
coef_20<- coef_20[!is.na(coef_20)]
se_20<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="se",range=cell_cols("AO:AT") )
se_20<- se_20[!is.na(se_20)]
tstat_20<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/NZEP.xlsx",sheet="tstat",range=cell_cols("AP:AR") )
tstat_20<- tstat_20[!is.na(tstat_20)]


#  First digit tests
####################
pi0<- 0.5        # prior probability for null model (B-L)

p1<- c(0.30103000,0.17609126,0.12493874,0.09691001,0.07918125,0.06694679,0.05799195,0.05115252,0.04575749)
theta<- p1
alpha_1<- rep(1,9)
alpha_2<- rep(1/9,9)
alpha_3<-22*theta       # See Fonesca, p.23 
alpha_4<- p1      

alpha<- alpha_3         # decide which alpha vector to use in the prior pdf

bfdc_6070<- benford(coef_6070,1,sign="both")  
bfds_6070<- benford(se_6070,1, sign="both")
bfdt_6070<- benford(tstat_6070,1, sign="both")
                                    
rN_c1_6070<- bfdc_6070$bfd$data.dist   # these values are "r/N"
rN_s1_6070<- bfds_6070$bfd$data.dist
rN_t1_6070<- bfdt_6070$bfd$data.dist

bfdc_80<- benford(coef_80,1,sign="both")
bfds_80<- benford(se_80,1, sign="both")
bfdt_80<- benford(tstat_80,1, sign="both")
rN_c1_80<- bfdc_80$bfd$data.dist   # these values are "r/N"
rN_s1_80<- bfds_80$bfd$data.dist
rN_t1_80<- bfdt_80$bfd$data.dist

bfdc_90<- benford(coef_90,1,sign="both")
bfds_90<- benford(se_90,1, sign="both")
bfdt_90<- benford(tstat_90,1, sign="both")
rN_c1_90<- bfdc_90$bfd$data.dist   # these values are "r/N"
rN_s1_90<- bfds_90$bfd$data.dist
rN_t1_90<- bfdt_90$bfd$data.dist

bfdc_00<- benford(coef_00,1,sign="both")
bfds_00<- benford(se_00,1, sign="both")
bfdt_00<- benford(tstat_00,1, sign="both")
rN_c1_00<- bfdc_00$bfd$data.dist   # these values are "r/N"
rN_s1_00<- bfds_00$bfd$data.dist
rN_t1_00<- bfdt_00$bfd$data.dist

bfdc_10<- benford(coef_10,1,sign="both")
bfds_10<- benford(se_10,1, sign="both")
bfdt_10<- benford(tstat_10,1, sign="both")
rN_c1_10<- bfdc_10$bfd$data.dist   # these values are "r/N"
rN_s1_10<- bfds_10$bfd$data.dist
rN_t1_10<- bfdt_10$bfd$data.dist

bfdc_20<- benford(coef_20,1,sign="both")
bfds_20<- benford(se_20,1, sign="both")
bfdt_20<- benford(tstat_20,1, sign="both")
rN_c1_20<- bfdc_20$bfd$data.dist   # these values are "r/N"
rN_s1_20<- bfds_20$bfd$data.dist
rN_t1_20<- bfdt_20$bfd$data.dist

# Coefficients
# ------------

# 1966-1979
nc<- length(coef_6070)
x<- rN_c1_6070*nc
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))     # log-gamma uses natural logs
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_6070<- (term1+term2+term3-term4-term5)/2.303     # convrt back to base 10 logs
post_prob_6070<- 1/(1+(1-pi0)/(pi0*10^(logB01_6070)))    # Posterior probability in favour of Benford's Law (H0)

# 1980-1989
nc<- length(coef_80)
x<- rN_c1_80*nc
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_80<- (term1+term2+term3-term4-term5)/2.303
post_prob_80<- 1/(1+(1-pi0)/(pi0*10^(logB01_80)))    # Posterior probability in favour of B-L (H0)

# 1990-1999
nc<- length(coef_90)
x<- rN_c1_90*nc
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_90<- (term1+term2+term3-term4-term5)/2.303
post_prob_90<- 1/(1+(1-pi0)/(pi0*10^(logB01_90)))    # Posterior probability in favour of B-L (H0)

# 2000-2009
nc<- length(coef_00)
x<- rN_c1_00*nc
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_00<- (term1+term2+term3-term4-term5)/2.303
post_prob_00<- 1/(1+(1-pi0)/(pi0*10^(logB01_00)))    # Posterior probability in favour of B-L (H0)

# 2010-2019
nc<- length(coef_10)
x<- rN_c1_10*nc
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_10<- (term1+term2+term3-term4-term5)/2.303
post_prob_10<- 1/(1+(1-pi0)/(pi0*10^(logB01_10)))    # Posterior probability in favour of B-L (H0)

# 2020-2025
nc<- length(coef_20)
x<- rN_c1_20*nc
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_20<- (term1+term2+term3-term4-term5)/2.303
post_prob_20<- 1/(1+(1-pi0)/(pi0*10^(logB01_20)))    # Posterior probability in favour of B-L (H0)

c(logB01_6070, logB01_80,logB01_90,logB01_00,logB01_10,logB01_20)
c(post_prob_6070, post_prob_80,post_prob_90,post_prob_00,post_prob_10,post_prob_20)


# Std. Errors
# -----------

# 1966-1979
ns<- length(se_6070)
x<- rN_s1_6070*ns
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_6070<- (term1+term2+term3-term4-term5)/2.303
post_prob_6070<- 1/(1+(1-pi0)/(pi0*10^(logB01_6070)))    # Posterior probability in favour of B-L (H0)

# 1980-1989
ns<- length(se_80)
x<- rN_s1_80*ns
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_80<- (term1+term2+term3-term4-term5)/2.303
post_prob_80<- 1/(1+(1-pi0)/(pi0*10^(logB01_80)))    # Posterior probability in favour of B-L (H0)

# 1990-1999
ns<- length(se_90)
x<- rN_s1_90*ns
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_90<- (term1+term2+term3-term4-term5)/2.303
post_prob_90<- 1/(1+(1-pi0)/(pi0*10^(logB01_90)))    # Posterior probability in favour of B-L (H0)

# 2000-2009
ns<- length(se_00)
x<- rN_s1_00*ns
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_00<- (term1+term2+term3-term4-term5)/2.303
post_prob_00<- 1/(1+(1-pi0)/(pi0*10^(logB01_00)))    # Posterior probability in favour of B-L (H0)

# 2010-2019
ns<- length(se_10)
x<- rN_s1_10*ns
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_10<- (term1+term2+term3-term4-term5)/2.303
post_prob_10<- 1/(1+(1-pi0)/(pi0*10^(logB01_10)))    # Posterior probability in favour of B-L (H0)

# 2020-2025
ns<- length(se_20)
x<- rN_s1_20*ns
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_20<- (term1+term2+term3-term4-term5)/2.303
post_prob_20<- 1/(1+(1-pi0)/(pi0*10^(logB01_20)))    # Posterior probability in favour of B-L (H0)

c(logB01_6070, logB01_80,logB01_90,logB01_00,logB01_10,logB01_20)
c(post_prob_6070, post_prob_80,post_prob_90,post_prob_00,post_prob_10, post_prob_20)


# t-Statistics
# ------------

# 1966-1979
nt<- length(tstat_6070)
x<- rN_t1_6070*nt
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_6070<- (term1+term2+term3-term4-term5)/2.303
post_prob_6070<- 1/(1+(1-pi0)/(pi0*10^(logB01_6070)))    # Posterior probability in favour of B-L (H0)

# 1980-1989
nt<- length(tstat_80)
x<- rN_t1_80*nt
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_80<- (term1+term2+term3-term4-term5)/2.303
post_prob_80<- 1/(1+(1-pi0)/(pi0*10^(logB01_80)))    # Posterior probability in favour of B-L (H0)

# 1990-1999
nt<- length(tstat_90)
x<- rN_t1_90*nt
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_90<- (term1+term2+term3-term4-term5)/2.303
post_prob_90<- 1/(1+(1-pi0)/(pi0*10^(logB01_90)))    # Posterior probability in favour of B-L (H0)

# 2000-2009
nt<- length(tstat_00)
x<- rN_t1_00*nt
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_00<- (term1+term2+term3-term4-term5)/2.303
post_prob_00<- 1/(1+(1-pi0)/(pi0*10^(logB01_00)))    # Posterior probability in favour of B-L (H0)

# 2010-2019
nt<- length(tstat_10)
x<- rN_t1_10*nt
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_10<- (term1+term2+term3-term4-term5)/2.303
post_prob_10<- 1/(1+(1-pi0)/(pi0*10^(logB01_10)))    # Posterior probability in favour of B-L (H0)

# 2020-2025
nt<- length(tstat_20)
x<- rN_t1_20*nt
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_20<- (term1+term2+term3-term4-term5)/2.303
post_prob_20<- 1/(1+(1-pi0)/(pi0*10^(logB01_20)))    # Posterior probability in favour of B-L (H0)

c(logB01_6070, logB01_80,logB01_90,logB01_00,logB01_10,logB01_20)
c(post_prob_6070, post_prob_80,post_prob_90,post_prob_00,post_prob_10, post_prob_20)


# Second digit tests
####################

pi0<- 0.5        # prior probability for null model (B-L)
p2<- c(0.1197,0.1139,0.1088,0.1043,0.1,0.0967, 0.0934,0.0904,0.0876,0.0850)
theta<- p2
alpha_1<- rep(1,10)
alpha_2<- rep(0.1,10)
alpha_3<-12*theta       # See Fonesca, p.23 
alpha_4<- p2      
alpha<- alpha_4         # decide which alpha vector to use in the prior pdf


# Benford 2nd-digit probs. (for i = 0 to 9)

N_c2<- c()
N_s2<- c()
N_t2<- c()

# 1966-1979

coef_2_ch <- str_sub(as.character(abs(coef_6070)), 2, 2)   # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_6070), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_6070)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
# Convert back to numeric
coef_2 <- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9. These are the digit counts
N_c2[ii] <- sum(coef_2 == i1)
N_s2[ii] <- sum(se_2 == i1)
N_t2[ii] <- sum(tstat_2 == i1)
}

# Coefficients
# ------------

x<- N_c2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_c_6070<- (term1+term2+term3-term4-term5)/2.303
post_prob_c_6070<- 1/(1+(1-pi0)/(pi0*10^(logB01_c_6070)))    # Posterior probability in favour of B-L (H0)

# Std. Errors
# -----------

x<- N_s2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_s_6070<- (term1+term2+term3-term4-term5)/2.303
post_prob_s_6070<- 1/(1+(1-pi0)/(pi0*10^(logB01_s_6070)))    # Posterior probability in favour of B-L (H0)

# t-Statistics
# ------------

x<- N_t2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_t_6070<- (term1+term2+term3-term4-term5)/2.303
post_prob_t_6070<- 1/(1+(1-pi0)/(pi0*10^(logB01_t_6070)))    # Posterior probability in favour of B-L (H0)

# 1980-1989

coef_2_ch <- str_sub(as.character(abs(coef_80)), 2, 2)   # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_80), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_80)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
# Convert back to numeric
coef_2 <- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
N_c2[ii] <- sum(coef_2 == i1)
N_s2[ii] <- sum(se_2 == i1)
N_t2[ii] <- sum(tstat_2 == i1)
}

# Coefficients
# ------------

x<- N_c2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_c_80<- (term1+term2+term3-term4-term5)/2.303
post_prob_c_80<- 1/(1+(1-pi0)/(pi0*10^(logB01_c_80)))    # Posterior probability in favour of B-L (H0)

# Std. Errors
# -----------

x<- N_s2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_s_80<- (term1+term2+term3-term4-term5)/2.303
post_prob_s_80<- 1/(1+(1-pi0)/(pi0*10^(logB01_s_80)))    # Posterior probability in favour of B-L (H0)

# t-Statistics
# ------------

x<- N_t2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_t_80<- (term1+term2+term3-term4-term5)/2.303
post_prob_t_80<- 1/(1+(1-pi0)/(pi0*10^(logB01_t_80)))    # Posterior probability in favour of B-L (H0)

# 1990-1999

coef_2_ch <- str_sub(as.character(abs(coef_90)), 2, 2)   # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_90), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_90)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
# Convert back to numeric
coef_2 <- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]
for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
N_c2[ii] <- sum(coef_2 == i1)
N_s2[ii] <- sum(se_2 == i1)
N_t2[ii] <- sum(tstat_2 == i1)
}
# Coefficients
# ------------

x<- N_c2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_c_90<- (term1+term2+term3-term4-term5)/2.303
post_prob_c_90<- 1/(1+(1-pi0)/(pi0*10^(logB01_c_90)))    # Posterior probability in favour of B-L (H0)

# Std. Errors
# -----------

x<- N_s2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_s_90<- (term1+term2+term3-term4-term5)/2.303
post_prob_s_90<- 1/(1+(1-pi0)/(pi0*10^(logB01_s_90)))    # Posterior probability in favour of B-L (H0)

# t-Statistics
# ------------

x<- N_t2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_t_90<- (term1+term2+term3-term4-term5)/2.303
post_prob_t_90<- 1/(1+(1-pi0)/(pi0*10^(logB01_t_90)))    # Posterior probability in favour of B-L (H0)

# 2000-2009

coef_2_ch <- str_sub(as.character(abs(coef_00)), 2, 2)   # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_00), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_00)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
# Convert back to numeric
coef_2 <- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]
for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
N_c2[ii] <- sum(coef_2 == i1)
N_s2[ii] <- sum(se_2 == i1)
N_t2[ii] <- sum(tstat_2 == i1)
}
# Coefficients
# ------------

x<- N_c2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_c_00<- (term1+term2+term3-term4-term5)/2.303
post_prob_c_00<- 1/(1+(1-pi0)/(pi0*10^(logB01_c_00)))    # Posterior probability in favour of B-L (H0)

# Std. Errors
# -----------

x<- N_s2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_s_00<- (term1+term2+term3-term4-term5)/2.303
post_prob_s_00<- 1/(1+(1-pi0)/(pi0*10^(logB01_s_00)))    # Posterior probability in favour of B-L (H0)

# t-Statistics
# ------------

x<- N_t2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_t_00<- (term1+term2+term3-term4-term5)/2.303
post_prob_t_00<- 1/(1+(1-pi0)/(pi0*10^(logB01_t_00)))    # Posterior probability in favour of B-L (H0)

# 2010-2019

coef_2_ch <- str_sub(as.character(abs(coef_10)), 2, 2)   # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_10), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_10)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
# Convert back to numeric
coef_2 <- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]
for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
N_c2[ii] <- sum(coef_2 == i1)
N_s2[ii] <- sum(se_2 == i1)
N_t2[ii] <- sum(tstat_2 == i1)
}
# Coefficients
# ------------

x<- N_c2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_c_10<- (term1+term2+term3-term4-term5)/2.303
post_prob_c_10<- 1/(1+(1-pi0)/(pi0*10^(logB01_c_10)))    # Posterior probability in favour of B-L (H0)

# Std. Errors
# -----------

x<- N_s2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_s_10<- (term1+term2+term3-term4-term5)/2.303
post_prob_s_10<- 1/(1+(1-pi0)/(pi0*10^(logB01_s_10)))    # Posterior probability in favour of B-L (H0)

# t-Statistics
# ------------

x<- N_t2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_t_10<- (term1+term2+term3-term4-term5)/2.303
post_prob_t_10<- 1/(1+(1-pi0)/(pi0*10^(logB01_t_10)))    # Posterior probability in favour of B-L (H0)

# 2020-2025

coef_2_ch <- str_sub(as.character(abs(coef_20)), 2, 2)   # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_20), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_20)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
# Convert back to numeric
coef_2 <- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]
for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
N_c2[ii] <- sum(coef_2 == i1)
N_s2[ii] <- sum(se_2 == i1)
N_t2[ii] <- sum(tstat_2 == i1)
}
# Coefficients
# ------------

x<- N_c2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_c_20<- (term1+term2+term3-term4-term5)/2.303
post_prob_c_20<- 1/(1+(1-pi0)/(pi0*10^(logB01_c_20)))    # Posterior probability in favour of B-L (H0)

# Std. Errors
# -----------

x<- N_s2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_s_20<- (term1+term2+term3-term4-term5)/2.303
post_prob_s_20<- 1/(1+(1-pi0)/(pi0*10^(logB01_s_20)))    # Posterior probability in favour of B-L (H0)

# t-Statistics
# ------------

x<- N_t2
term1<- sum(x*log(theta))
term2<- sum(lgamma(alpha))
term3<- lgamma(sum(alpha+x))
term4<- lgamma(sum(alpha))
term5<- sum(lgamma(alpha+x))
logB01_t_20<- (term1+term2+term3-term4-term5)/2.303
post_prob_t_20<- 1/(1+(1-pi0)/(pi0*10^(logB01_t_20)))    # Posterior probability in favour of B-L (H0)


c(logB01_c_6070,logB01_c_80,logB01_c_90,logB01_c_00,logB01_c_10,logB01_c_20)
c(post_prob_c_6070,post_prob_c_80,post_prob_c_90,post_prob_c_00,post_prob_c_10,post_prob_c_20)
c(logB01_s_6070,logB01_s_80,logB01_s_90,logB01_s_00,logB01_s_10,logB01_s_20)
c(post_prob_s_6070,post_prob_s_80,post_prob_s_90,post_prob_s_00,post_prob_s_10,post_prob_s_20)
c(logB01_t_6070,logB01_t_80,logB01_t_90,logB01_t_00,logB01_t_10,logB01_t_20)
c(post_prob_t_6070,post_prob_t_80,post_prob_t_90,post_prob_t_00,post_prob_t_10,post_prob_t_20)

