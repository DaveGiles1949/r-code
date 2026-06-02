# Todter M test for decade sub-samples
#
#############################################
library(benford.analysis)
library(stringr)

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

coef_20<- c(coef_2020,coef_2021,coef_2022,coef_2023,coef_2024,coef_2025)
coef_10<- c(coef_2010,coef_2011,coef_2012,coef_2013,coef_2014,coef_2015,coef_2016,coef_2017,coef_2018,coef_2019)
coef_00<- c(coef_2000,coef_2001,coef_2002,coef_2003,coef_2004,coef_2005,coef_2006,coef_2007,coef_2008,coef_2009)
coef_90<- c(coef_1990,coef_1991,coef_1992,coef_1993,coef_1994,coef_1995,coef_1996,coef_1997,coef_1998,coef_1999)
coef_80<- c(coef_1980,coef_1981,coef_1983,coef_1985,coef_1986,coef_1987,coef_1988,coef_1989)
coef_6070<- c(coef_1966,coef_1968,coef_1969,coef_1970,coef_1971,coef_1972,coef_1973,coef_1974,coef_1975,coef_1976, coef_1977,coef_1978,coef_1979)

se_20<- c(se_2020,se_2021,se_2022,se_2023,se_2024,se_2025)
se_10<- c(se_2010,se_2011,se_2012,se_2013,se_2014,se_2015,se_2016,se_2017,se_2018,se_2019)
se_00<- c(se_2000,se_2001,se_2003,se_2004,se_2005,se_2006,se_2007,se_2008,se_2009)
se_90<- c(se_1990,se_1992,se_1995,se_1996,se_1997,se_1998,se_1999)
se_80<- c(se_1980,se_1981,se_1983,se_1988)
se_6070<- c(se_1966,se_1968,se_1969,se_1971,se_1972,se_1973,se_1974,se_1975,se_1977,se_1979)

tstat_20<- c(tstat_2021,tstat_2022,tstat_2024)
tstat_10<- c(tstat_2010,tstat_2011,tstat_2012,tstat_2014,tstat_2015,tstat_2016,tstat_2017,tstat_2018,tstat_2019)
tstat_00<- c(tstat_2000,tstat_2001,tstat_2002,tstat_2003,tstat_2004,tstat_2006,tstat_2008)
tstat_90<- c(tstat_1990,tstat_1991,tstat_1992,tstat_1993,tstat_1994,tstat_1995,tstat_1996,tstat_1997,tstat_1998,tstat_1999)
tstat_80<- c(tstat_1983,tstat_1985,tstat_1986,tstat_1987,tstat_1988,tstat_1989)
tstat_6070<- c(tstat_1970,tstat_1971,tstat_1973,tstat_1974,tstat_1975,tstat_1976,tstat_1977,tstat_1978,tstat_1979)


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

