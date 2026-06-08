# Chi-square tests for sub-samples based on decades
##################################################
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

# Start computing for various sub-samples:
########################################

# Main Benford analysis for 1st. and 2nd. digits
# ----------------------------------------------

#  First digit tests
####################

bfdc_6070<- benford(coef_6070,1,sign="both")  # The chi-square test is included in the results of this command
chisqc_6070<- chisq(bfdc_6070)[[1]]
pvalc_6070<- chisq(bfdc_6070)[3]
lamc_6070<- 0.0238*length(coef_6070)
plamc_6070<- 1-pchisq(as.numeric(chisqc_6070), df=8,ncp=lamc_6070)
bfds_6070<- benford(se_6070,1, sign="both")
chisqs_6070<- chisq(bfds_6070)[[1]]
pvals_6070<- chisq(bfds_6070)[3]
lams_6070<- 0.0238*length(se_6070)
plams_6070<- 1-pchisq(as.numeric(chisqs_6070), df=8,ncp=lams_6070)
bfdt_6070<- benford(tstat_6070,1, sign="both")
chisqt_6070<- chisq(bfdt_6070)[1]
pvalt_6070<- chisq(bfdt_6070)[3]
lamt_6070<- 0.0238*length(tstat_6070)
plamt_6070<- 1-pchisq(as.numeric(chisqt_6070), df=8,ncp=lamt_6070)

bfdc_80<- benford(coef_80,1,sign="both")
chisqc_80<- chisq(bfdc_80)[1]
pvalc_80<- chisq(bfdc_80)[3]
lamc_80<- 0.0238*length(coef_80)
plamc_80<- 1-pchisq(as.numeric(chisqc_80), df=8,ncp=lamc_80)
bfds_80<- benford(se_80,1, sign="both")
chisqs_80<- chisq(bfds_80)[1]
pvals_80<- chisq(bfds_80)[3]
lams_80<- 0.0238*length(se_80)
plams_80<- 1-pchisq(as.numeric(chisqs_80), df=8,ncp=lams_80)
bfdt_80<- benford(tstat_80,1, sign="both")
chisqt_80<- chisq(bfdt_80)[1]
pvalt_80<- chisq(bfdt_80)[3]
lamt_80<- 0.0238*length(tstat_80)
plamt_80<- 1-pchisq(as.numeric(chisqt_80), df=8,ncp=lamt_80)

bfdc_90<- benford(coef_90,1,sign="both")
chisqc_90<- chisq(bfdc_90)[1]
pvalc_90<- chisq(bfdc_90)[3]
lamc_90<- 0.0238*length(coef_90)
plamc_90<- 1-pchisq(as.numeric(chisqc_90), df=8,ncp=lamc_90)
bfds_90<- benford(se_90,1, sign="both")
chisqs_90<- chisq(bfds_90)[1]
pvals_90<- chisq(bfds_90)[3]
lams_90<- 0.0238*length(se_90)
plams_90<- 1-pchisq(as.numeric(chisqs_90), df=8,ncp=lams_90)
bfdt_90<- benford(tstat_90,1, sign="both")
chisqt_90<- chisq(bfdt_90)[1]
pvalt_90<- chisq(bfdt_90)[3]
lamt_90<- 0.0238*length(tstat_90)
plamt_90<- 1-pchisq(as.numeric(chisqt_90), df=8,ncp=lamt_90)

bfdc_00<- benford(coef_00,1,sign="both")
chisqc_00<- chisq(bfdc_00)[1]
pvalc_00<- chisq(bfdc_00)[3]
lamc_00<- 0.0238*length(coef_00)
plamc_00<- 1-pchisq(as.numeric(chisqc_00), df=8,ncp=lamc_00)
bfds_00<- benford(se_00,1, sign="both")
chisqs_00<- chisq(bfds_00)[1]
pvals_00<- chisq(bfds_00)[3]
lams_00<- 0.0238*length(se_00)
plams_00<- 1-pchisq(as.numeric(chisqs_00), df=8,ncp=lams_00)
bfdt_00<- benford(tstat_00,1, sign="both")
chisqt_00<- chisq(bfdt_00)[1]
pvalt_00<- chisq(bfdt_00)[3]
lamt_00<- 0.0238*length(tstat_00)
plamt_00<- 1-pchisq(as.numeric(chisqt_00), df=8,ncp=lamt_00)

bfdc_10<- benford(coef_10,1,sign="both")
chisqc_10<- chisq(bfdc_10)[1]
pvalc_10<- chisq(bfdc_10)[3]
lamc_10<- 0.0238*length(coef_10)
plamc_10<- 1-pchisq(as.numeric(chisqc_10), df=8,ncp=lamc_10)
bfds_10<- benford(se_10,1, sign="both")
chisqs_10<- chisq(bfds_10)[1]
pvals_10<- chisq(bfds_10)[3]
lams_10<- 0.0238*length(se_10)
plams_10<- 1-pchisq(as.numeric(chisqs_10), df=8,ncp=lams_10)
bfdt_10<- benford(tstat_10,1, sign="both")
chisqt_10<- chisq(bfdt_10)[1]
pvalt_10<- chisq(bfdt_10)[3]
lamt_10<- 0.0238*length(tstat_10)
plamt_10<- 1-pchisq(as.numeric(chisqt_10), df=8,ncp=lamt_10)

bfdc_20<- benford(coef_20,1,sign="both")
chisqc_20<- chisq(bfdc_20)[1]
pvalc_20<- chisq(bfdc_20)[3]
lamc_20<- 0.0238*length(coef_20)
plamc_20<- 1-pchisq(as.numeric(chisqc_20), df=8,ncp=lamc_20)
bfds_20<- benford(se_20,1, sign="both")
chisqs_20<- chisq(bfds_20)[1]
pvals_20<- chisq(bfds_20)[3]
lams_20<- 0.0238*length(se_20)
plams_20<- 1-pchisq(as.numeric(chisqs_20), df=8,ncp=lams_20)
bfdt_20<- benford(tstat_20,1, sign="both")
chisqt_20<- chisq(bfdt_20)[1]
pvalt_20<- chisq(bfdt_20)[3]
lamt_20<- 0.0238*length(tstat_20)
plamt_20<- 1-pchisq(as.numeric(chisqt_20), df=8,ncp=lamt_20)

as.numeric(c(chisqc_6070, chisqc_80, chisqc_90,chisqc_00,chisqc_10,chisqc_20))
as.numeric(c(pvalc_6070, pvalc_80, pvalc_90,pvalc_00,pvalc_10,pvalc_20))
c(plamc_6070,plamc_80,plamc_90,plamc_00,plamc_10,plamc_20)

as.numeric(c(chisqs_6070, chisqs_80, chisqs_90,chisqs_00,chisqs_10,chisqs_20))
as.numeric(c(pvals_6070, pvals_80, pvals_90,pvals_00,pvals_10,pvals_20))
c(plams_6070,plams_80,plams_90,plams_00,plams_10,plams_20)

as.numeric(c(chisqt_6070, chisqt_80, chisqt_90,chisqt_00,chisqt_10,chisqt_20))
as.numeric(c(pvalt_6070, pvalt_80, pvalt_90,pvalt_00,pvalt_10,pvalt_20))
c(plamt_6070,plamt_80,plamt_90,plamt_00,plamt_10,plamt_20)

#########################
# Second Digit Tests
# ------------------

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
nc<- length(coef_2)
ns<- length(se_2)
nt<- length(tstat_2)

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
rN_c2[ii] <- sum(coef_2 == i1)/nc
rN_s2[ii] <- sum(se_2 == i1)/ns
rN_t2[ii] <- sum(tstat_2 == i1)/nt
}

chisqc_6070<- nc*sum((rN_c2-p2)^2/p2)
pc_6070<- 1-pchisq(chisqc_6070, df=9)
lamc_6070<- 0.0474*nc
plamc_6070<- 1-pchisq(chisqc_6070,df=9,ncp=lamc_6070)
chisqs_6070<- ns*sum((rN_s2-p2)^2/p2)
ps_6070<- 1-pchisq(chisqs_6070, df=9)
lams_6070<- 0.0474*ns
plams_6070<- 1-pchisq(chisqs_6070,df=9,ncp=lams_6070)
chisqt_6070<- nt*sum((rN_t2-p2)^2/p2)
pt_6070<- 1-pchisq(chisqt_6070, df=9)
lamt_6070<- 0.0474*nt
plamt_6070<- 1-pchisq(chisqt_6070,df=9,ncp=lams_6070)

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
nc<- length(coef_2)
ns<- length(se_2)
nt<- length(tstat_2)

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
rN_c2[ii] <- sum(coef_2 == i1)/nc
rN_s2[ii] <- sum(se_2 == i1)/ns
rN_t2[ii] <- sum(tstat_2 == i1)/nt
}

chisqc_80<- nc*sum((rN_c2-p2)^2/p2)
pc_80<- 1-pchisq(chisqc_80, df=9)
lamc_80<- 0.0474*nc
plamc_80<- 1-pchisq(chisqc_80,df=9,ncp=lamc_80)
chisqs_80<- ns*sum((rN_s2-p2)^2/p2)
ps_80<- 1-pchisq(chisqs_80, df=9)
lams_80<- 0.0474*ns
plams_80<- 1-pchisq(chisqs_80,df=9,ncp=lams_80)
chisqt_80<- nt*sum((rN_t2-p2)^2/p2)
pt_80<- 1-pchisq(chisqt_80, df=9)
lamt_80<- 0.0474*nt
plamt_80<- 1-pchisq(chisqt_80,df=9,ncp=lams_80)

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
nc<- length(coef_2)
ns<- length(se_2)
nt<- length(tstat_2)

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
rN_c2[ii] <- sum(coef_2 == i1)/nc
rN_s2[ii] <- sum(se_2 == i1)/ns
rN_t2[ii] <- sum(tstat_2 == i1)/nt
}

chisqc_90<- nc*sum((rN_c2-p2)^2/p2)
pc_90<- 1-pchisq(chisqc_90, df=9)
lamc_90<- 0.0474*nc
plamc_90<- 1-pchisq(chisqc_90,df=9,ncp=lamc_90)
chisqs_90<- ns*sum((rN_s2-p2)^2/p2)
ps_90<- 1-pchisq(chisqs_90, df=9)
lams_90<- 0.0474*ns
plams_90<- 1-pchisq(chisqs_90,df=9,ncp=lams_90)
chisqt_90<- nt*sum((rN_t2-p2)^2/p2)
pt_90<- 1-pchisq(chisqt_90, df=9)
lamt_90<- 0.0474*nt
plamt_90<- 1-pchisq(chisqt_90,df=9,ncp=lams_90)

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
nc<- length(coef_2)
ns<- length(se_2)
nt<- length(tstat_2)

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
rN_c2[ii] <- sum(coef_2 == i1)/nc
rN_s2[ii] <- sum(se_2 == i1)/ns
rN_t2[ii] <- sum(tstat_2 == i1)/nt
}

chisqc_00<- nc*sum((rN_c2-p2)^2/p2)
pc_00<- 1-pchisq(chisqc_00, df=9)
lamc_00<- 0.0474*nc
plamc_00<- 1-pchisq(chisqc_00,df=9,ncp=lamc_00)
chisqs_00<- ns*sum((rN_s2-p2)^2/p2)
ps_00<- 1-pchisq(chisqs_00, df=9)
lams_00<- 0.0474*ns
plams_00<- 1-pchisq(chisqs_00,df=9,ncp=lams_00)
chisqt_00<- nt*sum((rN_t2-p2)^2/p2)
pt_00<- 1-pchisq(chisqt_00, df=9)
lamt_00<- 0.0474*nt
plamt_00<- 1-pchisq(chisqt_00,df=9,ncp=lams_00)

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
nc<- length(coef_2)
ns<- length(se_2)
nt<- length(tstat_2)

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
rN_c2[ii] <- sum(coef_2 == i1)/nc
rN_s2[ii] <- sum(se_2 == i1)/ns
rN_t2[ii] <- sum(tstat_2 == i1)/nt
}

chisqc_10<- nc*sum((rN_c2-p2)^2/p2)
pc_10<- 1-pchisq(chisqc_10, df=9)
lamc_10<- 0.0474*nc
plamc_10<- 1-pchisq(chisqc_10,df=9,ncp=lamc_10)
chisqs_10<- ns*sum((rN_s2-p2)^2/p2)
ps_10<- 1-pchisq(chisqs_10, df=9)
lams_10<- 0.0474*ns
plams_10<- 1-pchisq(chisqs_10,df=9,ncp=lams_10)
chisqt_10<- nt*sum((rN_t2-p2)^2/p2)
pt_10<- 1-pchisq(chisqt_10, df=9)
lamt_10<- 0.0474*nt
plamt_10<- 1-pchisq(chisqt_10,df=9,ncp=lams_10)

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
nc<- length(coef_2)
ns<- length(se_2)
nt<- length(tstat_2)

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
rN_c2[ii] <- sum(coef_2 == i1)/nc
rN_s2[ii] <- sum(se_2 == i1)/ns
rN_t2[ii] <- sum(tstat_2 == i1)/nt
}

chisqc_20<- nc*sum((rN_c2-p2)^2/p2)
pc_20<- 1-pchisq(chisqc_20, df=9)
lamc_20<- 0.0474*nc
plamc_20<- 1-pchisq(chisqc_20,df=9,ncp=lamc_20)
chisqs_20<- ns*sum((rN_s2-p2)^2/p2)
ps_20<- 1-pchisq(chisqs_20, df=9)
lams_20<- 0.0474*ns
plams_20<- 1-pchisq(chisqs_20,df=9,ncp=lams_20)
chisqt_20<- nt*sum((rN_t2-p2)^2/p2)
pt_20<- 1-pchisq(chisqt_20, df=9)
lamt_20<- 0.0474*nt
plamt_20<- 1-pchisq(chisqt_20,df=9,ncp=lams_20)


as.numeric(c(chisqc_6070, chisqc_80, chisqc_90,chisqc_00,chisqc_10,chisqc_20))
as.numeric(c(pc_6070, pc_80, pc_90,pc_00,pc_10,pc_20))
c(plamc_6070,plamc_80,plamc_90,plamc_00,plamc_10,plamc_20)

as.numeric(c(chisqs_6070, chisqs_80, chisqs_90,chisqs_00,chisqs_10,chisqs_20))
as.numeric(c(ps_6070, ps_80, ps_90,ps_00,ps_10,ps_20))
c(plams_6070,plams_80,plams_90,plams_00,plams_10,plams_20)

as.numeric(c(chisqt_6070, chisqt_80, chisqt_90,chisqt_00,chisqt_10,chisqt_20))
as.numeric(c(pt_6070, pt_80, pt_90,pt_00,pt_10,pt_20))
c(plamt_6070,plamt_80,plamt_90,plamt_00,plamt_10,plamt_20)

