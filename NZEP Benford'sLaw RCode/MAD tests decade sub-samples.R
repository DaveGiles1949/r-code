# MAD tests for sub-samples based on decades
#############################################
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
k<- 9
bfdc_6070<- benford(coef_6070,1,sign="both")  
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

# ------------
# See Cerqueti and Lupi, p.11 - delta_star = 0.01 for first-digit case
# The MAD values
# Obtain asy. distribution for via z-statistic

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
# Use thes values to standardize Normal-statistics based on MAD_N:

# 1966-1979
nc<- length(coef_6070)
MAD_K_c<- sqrt(nc)*sum(abs(p1-rN_c1_6070))
MAD_N_c<- MAD_K_c/(k*sqrt(nc))
ns<- length(se_6070)
MAD_K_s<- sqrt(ns)*sum(abs(p1-rN_s1_6070))
MAD_N_s<- MAD_K_s/(k*sqrt(ns))
nt<- length(tstat_6070)
MAD_K_t<- sqrt(nt)*sum(abs(p1-rN_t1_6070))
MAD_N_t<- MAD_K_t/(k*sqrt(nt))
z_c<- (sqrt(nc)*MAD_N_c -asy_mean)/sqrt(asy_var)
pval_c<- 1-pnorm(z_c, mean=0, sd=1)
z_s<- (sqrt(ns)*MAD_N_s -asy_mean)/sqrt(asy_var)
pval_s<- 1-pnorm(z_s, mean=0, sd=1)
z_t<- (sqrt(nt)*MAD_N_t -asy_mean)/sqrt(asy_var)
pval_t<- 1-pnorm(z_t, mean=0, sd=1)
c(MAD_N_c,MAD_N_s,MAD_N_t)
c(pval_c,pval_s,pval_t)
c(MAD_K_c,MAD_K_s,MAD_K_t)
# Asy. critical values for MAD_K tests: 2.869 3.084 3.485   (90, 95, 99) # Table 4 Kossler et al. (2024)

# 1980-1989
nc<- length(coef_80)
MAD_K_c<- sqrt(n)*sum(abs(p1-rN_c1_80))
MAD_N_c<- MAD_K_c/(k*sqrt(n))
ns<- length(se_80)
MAD_K_s<- sqrt(ns)*sum(abs(p1-rN_s1_80))
MAD_N_s<- MAD_K_s/(k*sqrt(ns))
nt<- length(tstat_80)
MAD_K_t<- sqrt(nt)*sum(abs(p1-rN_t1_80))
MAD_N_t<- MAD_K_t/(k*sqrt(nt))
z_c<- (sqrt(nc)*MAD_N_c -asy_mean)/sqrt(asy_var)
pval_c<- 1-pnorm(z_c, mean=0, sd=1)
z_s<- (sqrt(ns)*MAD_N_s -asy_mean)/sqrt(asy_var)
pval_s<- 1-pnorm(z_s, mean=0, sd=1)
z_t<- (sqrt(nt)*MAD_N_t -asy_mean)/sqrt(asy_var)
pval_t<- 1-pnorm(z_t, mean=0, sd=1)
c(MAD_N_c,MAD_N_s,MAD_N_t)
c(pval_c,pval_s,pval_t)
c(MAD_K_c,MAD_K_s,MAD_K_t)
# Asy. critical values for MAD_K tests: 2.869 3.084 3.485   (90, 95, 99) # Table 4 Kossler et al. (2024)

# 1990-1999
nc<- length(coef_90)
MAD_K_c<- sqrt(n)*sum(abs(p1-rN_c1_90))
MAD_N_c<- MAD_K_c/(k*sqrt(n))
ns<- length(se_90)
MAD_K_s<- sqrt(ns)*sum(abs(p1-rN_s1_90))
MAD_N_s<- MAD_K_s/(k*sqrt(ns))
nt<- length(tstat_90)
MAD_K_t<- sqrt(nt)*sum(abs(p1-rN_t1_90))
MAD_N_t<- MAD_K_t/(k*sqrt(nt))
z_c<- (sqrt(nc)*MAD_N_c -asy_mean)/sqrt(asy_var)
pval_c<- 1-pnorm(z_c, mean=0, sd=1)
z_s<- (sqrt(ns)*MAD_N_s -asy_mean)/sqrt(asy_var)
pval_s<- 1-pnorm(z_s, mean=0, sd=1)
z_t<- (sqrt(nt)*MAD_N_t -asy_mean)/sqrt(asy_var)
pval_t<- 1-pnorm(z_t, mean=0, sd=1)
c(MAD_N_c,MAD_N_s,MAD_N_t)
c(pval_c,pval_s,pval_t)
c(MAD_K_c,MAD_K_s,MAD_K_t)
# Asy. critical values for MAD_K tests: 2.869 3.084 3.485   (90, 95, 99) # Table 4 Kossler et al. (2024)

# 2000-2009
nc<- length(coef_00)
MAD_K_c<- sqrt(n)*sum(abs(p1-rN_c1_00))
MAD_N_c<- MAD_K_c/(k*sqrt(n))
ns<- length(se_00)
MAD_K_s<- sqrt(ns)*sum(abs(p1-rN_s1_00))
MAD_N_s<- MAD_K_s/(k*sqrt(ns))
nt<- length(tstat_00)
MAD_K_t<- sqrt(nt)*sum(abs(p1-rN_t1_00))
MAD_N_t<- MAD_K_t/(k*sqrt(nt))
z_c<- (sqrt(nc)*MAD_N_c -asy_mean)/sqrt(asy_var)
pval_c<- 1-pnorm(z_c, mean=0, sd=1)
z_s<- (sqrt(ns)*MAD_N_s -asy_mean)/sqrt(asy_var)
pval_s<- 1-pnorm(z_s, mean=0, sd=1)
z_t<- (sqrt(nt)*MAD_N_t -asy_mean)/sqrt(asy_var)
pval_t<- 1-pnorm(z_t, mean=0, sd=1)
c(MAD_N_c,MAD_N_s,MAD_N_t)
c(pval_c,pval_s,pval_t)
c(MAD_K_c,MAD_K_s,MAD_K_t)
# Asy. critical values for MAD_K tests: 2.869 3.084 3.485   (90, 95, 99) # Table 4 Kossler et al. (2024

# 2010-2019
nc<- length(coef_10)
MAD_K_c<- sqrt(n)*sum(abs(p1-rN_c1_10))
MAD_N_c<- MAD_K_c/(k*sqrt(n))
ns<- length(se_10)
MAD_K_s<- sqrt(ns)*sum(abs(p1-rN_s1_10))
MAD_N_s<- MAD_K_s/(k*sqrt(ns))
nt<- length(tstat_10)
MAD_K_t<- sqrt(nt)*sum(abs(p1-rN_t1_10))
MAD_N_t<- MAD_K_t/(k*sqrt(nt))
z_c<- (sqrt(nc)*MAD_N_c -asy_mean)/sqrt(asy_var)
pval_c<- 1-pnorm(z_c, mean=0, sd=1)
z_s<- (sqrt(ns)*MAD_N_s -asy_mean)/sqrt(asy_var)
pval_s<- 1-pnorm(z_s, mean=0, sd=1)
z_t<- (sqrt(nt)*MAD_N_t -asy_mean)/sqrt(asy_var)
pval_t<- 1-pnorm(z_t, mean=0, sd=1)
c(MAD_N_c,MAD_N_s,MAD_N_t)
c(pval_c,pval_s,pval_t)
c(MAD_K_c,MAD_K_s,MAD_K_t)
# Asy. critical values for MAD_K tests: 2.869 3.084 3.485   (90, 95, 99) # Table 4 Kossler et al. (2024)

# 2020-2025
nc<- length(coef_20)
MAD_K_c<- sqrt(n)*sum(abs(p1-rN_c1_20))
MAD_N_c<- MAD_K_c/(k*sqrt(n))
ns<- length(se_20)
MAD_K_s<- sqrt(ns)*sum(abs(p1-rN_s1_20))
MAD_N_s<- MAD_K_s/(k*sqrt(ns))
nt<- length(tstat_20)
MAD_K_t<- sqrt(nt)*sum(abs(p1-rN_t1_20))
MAD_N_t<- MAD_K_t/(k*sqrt(nt))
z_c<- (sqrt(nc)*MAD_N_c -asy_mean)/sqrt(asy_var)
pval_c<- 1-pnorm(z_c, mean=0, sd=1)
z_s<- (sqrt(ns)*MAD_N_s -asy_mean)/sqrt(asy_var)
pval_s<- 1-pnorm(z_s, mean=0, sd=1)
z_t<- (sqrt(nt)*MAD_N_t -asy_mean)/sqrt(asy_var)
pval_t<- 1-pnorm(z_t, mean=0, sd=1)
c(MAD_N_c,MAD_N_s,MAD_N_t)
c(pval_c,pval_s,pval_t)
c(MAD_K_c,MAD_K_s,MAD_K_t)
# Asy. critical values for MAD_K tests: 2.869 3.084 3.485   (90, 95, 99) # Table 4 Kossler et al. (2024)

