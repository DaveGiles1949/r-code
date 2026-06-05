# Chi-square tests, random samples
#--------------------------------

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

set.seed(111)
n<- 500          # choose sample size
#n<- 1000
#n<- 2000
xc<- sample(coef,n)     
xs<- sample(se,n)       
xt<- sample(tstat,n)    
nc<- length(xc)
ns<- length(xs)
nt<- length(xt)

bfdc<- benford(xc,1,sign="both")  # The chi-square test is included in the results of this command
chisqc<- chisq(bfdc)[[1]]
pvalc<- chisq(bfdc)[3]
lamc<- 0.0238*nc
plamc<- 1-pchisq(as.numeric(chisqc), df=8,ncp=lamc)
bfds<- benford(xs,1, sign="both")
chisqs<- chisq(bfds)[[1]]
pvals<- chisq(bfds)[3]
lams<- 0.0238*ns
plams<- 1-pchisq(as.numeric(chisqs), df=8,ncp=lams)
bfdt<- benford(xt,1, sign="both")
chisqt<- chisq(bfdt)[1]
pvalt<- chisq(bfdt)[3]
lamt<- 0.0238*nt
plamt<- 1-pchisq(as.numeric(chisqt), df=8,ncp=lamt)

as.numeric(c(chisqc,chisqs,chisqt))
as.numeric(c(pvalc,pvals,pvalt))

c(plamc,plams,plamt)

##########################
# Second digit analysis
#######################

set.seed(111)
n<- 2000          # choose sample size
#n<- 1000
#n<- 2000
xc<- sample(coef,n)     
xs<- sample(se,n)       
xt<- sample(tstat,n)    

# Benford 2nd-digit probs. (for i = 0 to 9)
p2<- c(0.1197,0.1139,0.1088,0.1043,0.1,0.0967, 0.0934,0.0904,0.0876,0.0850)
rN_c2<- c()
rN_s2<- c()
rN_t2<- c()

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

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
rN_c2[ii] <- sum(coef_2 == i1)/nc
rN_s2[ii] <- sum(se_2 == i1)/ns
rN_t2[ii] <- sum(tstat_2 == i1)/nt
}

chisqc<- nc*sum((rN_c2-p2)^2/p2)
pc<- 1-pchisq(chisqc, df=9)
lamc<- 0.0474*nc
plamc<- 1-pchisq(chisqc,df=9,ncp=lamc)
chisqs<- ns*sum((rN_s2-p2)^2/p2)
ps<- 1-pchisq(chisqs, df=9)
lams<- 0.0474*ns
plams<- 1-pchisq(chisqs,df=9,ncp=lams)
chisqt<- nt*sum((rN_t2-p2)^2/p2)
pt<- 1-pchisq(chisqt, df=9)
lamt<- 0.0474*nt
plamt<- 1-pchisq(chisqt,df=9,ncp=lams)

as.numeric(c(chisqc,chisqs,chisqt))
as.numeric(c(pc,ps,pt))

c(plamc,plams,plamt)

