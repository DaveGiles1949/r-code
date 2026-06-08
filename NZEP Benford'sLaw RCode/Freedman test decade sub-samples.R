# Freedman's test for sub-samples based on decades
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
d1<- 1:9
bfdc_6070<- benford(coef_6070,1,sign="both")  # The chi-square test is included in the results of this command
bfds_6070<- benford(se_6070,1, sign="both")
bfdt_6070<- benford(tstat_6070,1, sign="both")
p1<- bfdc_6070$bfd$benford.dist                                     # Benford probabilities for 1st digit
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
sig1<-signifd(x = coef_6070, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_c1_6070[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_c_6070<- (nc/9)*(sum1-sum2^2/9)

# 1980-1989
nc<- length(coef_80)
sig1<-signifd(x = coef_80, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_c1_80[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_c_80<- (nc/9)*(sum1-sum2^2/9)

# 1990-1999
nc<- length(coef_90)
sig1<-signifd(x = coef_90, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_c1_90[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_c_90<- (nc/9)*(sum1-sum2^2/9)

# 2000-2009
nc<- length(coef_00)
sig1<-signifd(x = coef_00, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_c1_00[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_c_00<- (nc/9)*(sum1-sum2^2/9)

# 2010-2019
nc<- length(coef_10)
sig1<-signifd(x = coef_10, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_c1_10[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_c_10<- (nc/9)*(sum1-sum2^2/9)

# 2020-2025
nc<- length(coef_20)
sig1<-signifd(x = coef_20, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_c1_20[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_c_20<- (nc/9)*(sum1-sum2^2/9)

c(Unsq_c_6070, Unsq_c_80,Unsq_c_90,Unsq_c_00,Unsq_c_10,Unsq_c_20)

bfdc_6070$stats$chisq[[1]]       # chi-square stat
pchisq(bfdc_6070$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value
bfdc_80$stats$chisq[[1]]
pchisq(bfdc_80$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value
bfdc_90$stats$chisq[[1]]         # chi-square stat
pchisq(bfdc_90$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value                        
bfdc_00$stats$chisq[[1]]         # chi-square stat   
pchisq(bfdc_00$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value
bfdc_10$stats$chisq[[1]]         # chi-square stat
pchisq(bfdc_10$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value
bfdc_20$stats$chisq[[1]]         # chi-square stat
pchisq(bfdc_20$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value

# Std. Errors
# -----------

# 1966-1979
ns<- length(se_6070)
sig1<-signifd(x = se_6070, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_s1_6070[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_s_6070<- (ns/9)*(sum1-sum2^2/9)

# 1980-1989
ns<- length(se_80)
sig1<-signifd(x = se_80, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_s1_80[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_s_80<- (ns/9)*(sum1-sum2^2/9)

# 1990-1999
ns<- length(se_90)
sig1<-signifd(x = se_90, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_s1_90[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_s_90<- (ns/9)*(sum1-sum2^2/9)

# 2000-2009
ns<- length(se_00)
sig1<-signifd(x = se_00, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_s1_00[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_s_00<- (ns/9)*(sum1-sum2^2/9)

# 2010-2019
ns<- length(se_10)
sig1<-signifd(x = se_10, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_s1_10[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_s_10<- (ns/9)*(sum1-sum2^2/9)

# 2020-2025
ns<- length(se_20)
sig1<-signifd(x = se_20, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_s1_20[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_s_20<- (ns/9)*(sum1-sum2^2/9)

c(Unsq_s_6070, Unsq_s_80,Unsq_s_90,Unsq_s_00,Unsq_s_10,Unsq_s_20)

bfds_6070$stats$chisq[[1]]       # chi-square stat
pchisq(bfds_6070$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value
bfds_80$stats$chisq[[1]]
pchisq(bfds_80$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value
bfds_90$stats$chisq[[1]]         # chi-square stat
pchisq(bfds_90$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value                        
bfds_00$stats$chisq[[1]]         # chi-square stat   
pchisq(bfds_00$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value
bfds_10$stats$chisq[[1]]         # chi-square stat
pchisq(bfds_10$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value
bfds_20$stats$chisq[[1]]         # chi-square stat
pchisq(bfds_20$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value

# t-Statistics
# ------------
# 1966-1979
nt<- length(tstat_6070)
sig1<-signifd(x = tstat_6070, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_t1_6070[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_t_6070<- (nt/9)*(sum1-sum2^2/9)

# 1980-1989
nt<- length(tstat_80)
sig1<-signifd(x = tstat_80, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_t1_80[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_t_80<- (nt/9)*(sum1-sum2^2/9)

# 1990-1999
nt<- length(tstat_90)
sig1<-signifd(x = tstat_90, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_t1_90[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_t_90<- (nt/9)*(sum1-sum2^2/9)

# 2000-2009
nt<- length(tstat_00)
sig1<-signifd(x = tstat_00, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_t1_00[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_t_00<- (nt/9)*(sum1-sum2^2/9)

# 2010-2019
nt<- length(tstat_10)
sig1<-signifd(x = tstat_10, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_t1_10[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_t_10<- (nt/9)*(sum1-sum2^2/9)

# 2020-2025
nt<- length(tstat_20)
sig1<-signifd(x = tstat_20, digits = 1)
S<- rep(0,9)
term<- c()
for (i in 1:9) {
term[i]<- rN_t1_20[i]-p1[i]
}

S<- cumsum(term[1:9])
sum1<- sum((S[1])^2:(S[8])^2) 
sum2<- sum(S[1]:S[8])
Unsq_t_20<- (nt/9)*(sum1-sum2^2/9)

c(Unsq_t_6070, Unsq_t_80,Unsq_t_90,Unsq_t_00,Unsq_t_10,Unsq_t_20)

bfdt_6070$stats$chisq[[1]]       # chi-square stat
pchisq(bfdt_6070$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value
bfdt_80$stats$chisq[[1]]
pchisq(bfdt_80$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value
bfdt_90$stats$chisq[[1]]         # chi-square stat
pchisq(bfdt_90$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value                        
bfdt_00$stats$chisq[[1]]         # chi-square stat   
pchisq(bfdt_00$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value
bfdt_10$stats$chisq[[1]]         # chi-square stat
pchisq(bfdt_10$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value
bfdt_20$stats$chisq[[1]]         # chi-square stat
pchisq(bfdt_20$stats$chisq[[1]],8, ncp = 0, lower.tail = FALSE, log.p = FALSE)   # p-value


# Second digit tests
####################

# Benford 2nd-digit probs. (for i = 0 to 9)
p2<- c(0.1197,0.1139,0.1088,0.1043,0.1,0.0967, 0.0934,0.0904,0.0876,0.0850)
rN_c2<- c()
rN_s2<- c()
rN_t2<- c()

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
nc<- length(coef_6070)
ns<- length(se_6070)
nt<- length(tstat_6070)

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
Unsq_c_6070<- (nc/10)*(sum1-sum2^2/10)

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
Unsq_s_6070<- (ns/10)*(sum1-sum2^2/10)

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
Unsq_t_6070<- (nt/10)*(sum1-sum2^2/10)

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
nc<- length(coef_80)
ns<- length(se_80)
nt<- length(tstat_80)

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
Unsq_c_80<- (nc/10)*(sum1-sum2^2/10)

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
Unsq_s_80<- (ns/10)*(sum1-sum2^2/10)

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
Unsq_t_80<- (nt/10)*(sum1-sum2^2/10)


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
nc<- length(coef_90)
ns<- length(se_90)
nt<- length(tstat_90)

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
Unsq_c_90<- (nc/10)*(sum1-sum2^2/10)

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
Unsq_s_90<- (ns/10)*(sum1-sum2^2/10)

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
Unsq_t_90<- (nt/10)*(sum1-sum2^2/10)


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
nc<- length(coef_00)
ns<- length(se_00)
nt<- length(tstat_00)

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
Unsq_c_00<- (nc/10)*(sum1-sum2^2/10)

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
Unsq_s_00<- (ns/10)*(sum1-sum2^2/10)

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
Unsq_t_00<- (nt/10)*(sum1-sum2^2/10)


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
nc<- length(coef_10)
ns<- length(se_10)
nt<- length(tstat_10)

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
Unsq_c_10<- (nc/10)*(sum1-sum2^2/10)

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
Unsq_s_10<- (ns/10)*(sum1-sum2^2/10)

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
Unsq_t_10<- (nt/10)*(sum1-sum2^2/10)


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
nc<- length(coef_20)
ns<- length(se_20)
nt<- length(tstat_20)

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
Unsq_c_20<- (nc/10)*(sum1-sum2^2/10)

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
Unsq_s_20<- (ns/10)*(sum1-sum2^2/10)

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
Unsq_t_20<- (nt/10)*(sum1-sum2^2/10)

c(Unsq_c_6070,Unsq_c_80,Unsq_c_90,Unsq_c_00,Unsq_c_10,Unsq_c_20)
c(Unsq_s_6070,Unsq_s_80,Unsq_s_90,Unsq_s_00,Unsq_s_10,Unsq_s_20)
c(Unsq_t_6070,Unsq_t_80,Unsq_t_90,Unsq_t_00,Unsq_t_10,Unsq_t_20)

