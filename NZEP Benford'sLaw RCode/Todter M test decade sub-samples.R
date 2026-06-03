# Todter M test for decade sub-samples
#
#############################################
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


# Start computing for various sub-samples:
########################################

# Main Benford analysis for 1st. and 2nd. digits
# ----------------------------------------------

# 1966-1979
#----------
xc<- coef_6070
xs<- se_6070
xt<- tstat_6070
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

# 1980-1989
#----------
xc<- coef_80
xs<- se_80
xt<- tstat_80
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

# 1990-1999
#----------
xc<- coef_90
xs<- se_90
xt<- tstat_90
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

# 2000-2010
#----------
xc<- coef_00
xs<- se_00
xt<- tstat_00
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

# 2011-2019
#----------
xc<- coef_10
xs<- se_10
xt<- tstat_10
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

# 2020-2025
#----------
xc<- coef_20
xs<- se_20
xt<- tstat_20
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

