# LOADING OF DATA BY YEAR,AND BASIC PLOTS FOR SUB-SAMPLES BASED ON DECADES
# SPECIFICALLY, PLOTS OF THE BENFORD ANALYSIS OF THE 1ST. AND 2ND. DIGITS 
# OF THE COEFFICIENTS, STD. ERRORS, AND T-STATISTICS
# 
# IN, ADDITION, RANDOM SAMPLES ARE ALSO DRAWN AND ANALYZED
# =======================================================================================

library(benford.analysis)
library(stringr)
set.seed(123)
library(readxl)

c<- c()
s<- c()
t<- c()
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

coef_total<- c(coef_6070, coef_80, coef_90, coef_00, coef_10, coef_20)
se_total<- c(se_6070, se_80, se_90, se_00, se_10, se_20)
tstat_total<- c(tstat_6070, tstat_80, tstat_90, tstat_00, tstat_10, tstat_20)

# First Digits
##############
d1<- 1:9
bfdc_6070<- benford(coef_6070,1,sign="both")
bfds_6070<- benford(se_6070,1, sign="both")
bfdt_6070<- benford(tstat_6070,1, sign="both")

bfdc_80<- benford(coef_80,1,sign="both")
bfds_80<- benford(se_80,1, sign="both")
bfdt_80<- benford(tstat_80,1, sign="both")

bfdc_90<- benford(coef_90,1,sign="both")
bfds_90<- benford(se_90,1, sign="both")
bfdt_90<- benford(tstat_90,1, sign="both")

bfdc_00<- benford(coef_00,1,sign="both")
bfds_00<- benford(se_00,1, sign="both")
bfdt_00<- benford(tstat_00,1, sign="both")

bfdc_10<- benford(coef_10,1,sign="both")
bfds_10<- benford(se_10,1, sign="both")
bfdt_10<- benford(tstat_10,1, sign="both")

bfdc_20<- benford(coef_20,1,sign="both")
bfds_20<- benford(se_20,1, sign="both")
bfdt_20<- benford(tstat_20,1, sign="both")

par(mfrow=c(1,1))

plot(d1,bfdc_6070$bfd$data.dist, col="red", xlab="1st. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(1,9,8), ylim=c(0,0.4), main="Figure 9 (a): 1966-1979
Relative Frequency of 1st. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
points(d1,bfds_6070$bfd$data.dist, col="blue", pch=8)
points(d1,bfdt_6070$bfd$data.dist, col="black", pch=17)
lines(d1, bfdc_6070$bfd$benford.dist)
text(5, 0.025, "n = 3,013", cex=0.9)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d1,bfdc_80$bfd$data.dist, col="red", xlab="1st. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(1,9,8), ylim=c(0,0.4), main="Figure 9 (b): 1980-1989
Relative Frequency of 1st. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
points(d1,bfds_80$bfd$data.dist, col="blue", pch=8)
points(d1,bfdt_80$bfd$data.dist, col="black", pch=17)
lines(d1, bfdc_80$bfd$benford.dist)
text(5, 0.025, "n = 950", cex=0.9)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d1,bfdc_90$bfd$data.dist, col="red", xlab="1st. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(1,9,8), ylim=c(0,0.4), main="Figure 9 (c): 1990-1999
Relative Frequency of 1st. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
points(d1,bfds_90$bfd$data.dist, col="blue", pch=8)
points(d1,bfdt_90$bfd$data.dist, col="black", pch=17)
lines(d1, bfdc_90$bfd$benford.dist)
text(5, 0.025, "n = 2,996", cex=0.9)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d1,bfdc_00$bfd$data.dist, col="red", xlab="1st. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(1,9,8), ylim=c(0,0.4), main="Figure 9 (d): 2000-2009
Relative Frequency of 1st. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
points(d1,bfds_00$bfd$data.dist, col="blue", pch=8)
points(d1,bfdt_00$bfd$data.dist, col="black", pch=17)
lines(d1, bfdc_00$bfd$benford.dist)
text(5, 0.025, "n = 4,264", cex=0.9)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d1,bfdc_10$bfd$data.dist, col="red", xlab="1st. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(1,9,8), ylim=c(0,0.4), main="Figure 9 (e): 2010-2019
Relative Frequency of 1st. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
points(d1,bfds_10$bfd$data.dist, col="blue", pch=8)
points(d1,bfdt_10$bfd$data.dist, col="black", pch=17)
lines(d1, bfdc_10$bfd$benford.dist)
text(5, 0.025, "n = 11,744", cex=0.9)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d1,bfdc_20$bfd$data.dist, col="red", xlab="1st. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(1,9,8), ylim=c(0,0.4), main="Figure 9 (f): 2020-2025
Relative Frequency of 1st. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
points(d1,bfds_20$bfd$data.dist, col="blue", pch=8)
points(d1,bfdt_20$bfd$data.dist, col="black", pch=17)
lines(d1, bfdc_20$bfd$benford.dist)
text(5, 0.025, "n = 7,192", cex=0.9)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

#####################
#
#SECOND DIGITS
#
d2<- 0:9
benford_2<- c(0.1197,0.1139,0.1088,0.1043,0.1,0.0967, 0.0934,0.0904,0.0876,0.0850)
coef_share_6070<- c()
se_share_6070<- c()
tstat_share_6070<- c()
coef_share_80<- c()
se_share_80<- c()
tstat_share_80<- c()
coef_share_90<- c()
se_share_90<- c()
tstat_share_90<- c()
coef_share_00<- c()
se_share_00<- c()
tstat_share_00<- c()
coef_share_10<- c()
se_share_10<- c()
tstat_share_10<- c()
coef_share_20<- c()
se_share_20<- c()
tstat_share_20<- c()

coef_2_ch <- str_sub(as.character(abs(coef_6070)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_6070), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_6070)), 2, 2)
# Convert back to numeric
coef_2<- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_share_6070[ii] <- sum(coef_2 == i1)/length(coef_2)
se_share_6070[ii] <- sum(se_2 == i1)/length(se_2)
tstat_share_6070[ii] <- sum(tstat_2 == i1)/length(tstat_2)
}

coef_2_ch <- str_sub(as.character(abs(coef_80)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_80), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_80)), 2, 2)
# Convert back to numeric
coef_2<- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_share_80[ii] <- sum(coef_2 == i1)/length(coef_2)
se_share_80[ii] <- sum(se_2 == i1)/length(se_2)
tstat_share_80[ii] <- sum(tstat_2 == i1)/length(tstat_2)
}

coef_2_ch <- str_sub(as.character(abs(coef_90)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_90), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_90)), 2, 2)
# Convert back to numeric
coef_2<- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_share_90[ii] <- sum(coef_2 == i1)/length(coef_2)
se_share_90[ii] <- sum(se_2 == i1)/length(se_2)
tstat_share_90[ii] <- sum(tstat_2 == i1)/length(tstat_2)
}

coef_2_ch <- str_sub(as.character(abs(coef_00)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_00), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_00)), 2, 2)
# Convert back to numeric
coef_2<- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_share_00[ii] <- sum(coef_2 == i1)/length(coef_2)
se_share_00[ii] <- sum(se_2 == i1)/length(se_2)
tstat_share_00[ii] <- sum(tstat_2 == i1)/length(tstat_2)
}

coef_2_ch <- str_sub(as.character(abs(coef_10)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_10), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_10)), 2, 2)
# Convert back to numeric
coef_2<- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_share_10[ii] <- sum(coef_2 == i1)/length(coef_2)
se_share_10[ii] <- sum(se_2 == i1)/length(se_2)
tstat_share_10[ii] <- sum(tstat_2 == i1)/length(tstat_2)
}

coef_2_ch <- str_sub(as.character(abs(coef_20)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se_20), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat_20)), 2, 2)
# Convert back to numeric
coef_2<- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_share_20[ii] <- sum(coef_2 == i1)/length(coef_2)
se_share_20[ii] <- sum(se_2 == i1)/length(se_2)
tstat_share_20[ii] <- sum(tstat_2 == i1)/length(tstat_2)
}

par(mfrow=c(1,1))

plot(d2,coef_share_6070, col="red", xlab="2nd. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(0,9,9), ylim=c(0.0,0.4), main="Figure 10 (a): 1966-1979
Relative Frequency of 2nd. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
lines(d2, benford_2)
text(2, 0.05, "n = 2,890", cex=0.9)
points(d2,se_share_6070, col="blue", pch=8)
points(d2,tstat_share_6070, col="black", pch=17)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d2,coef_share_80, col="red", xlab="2nd. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(0,9,9), ylim=c(0.0,0.4), main="Figure 10 (b): 1980-1989
Relative Frequency of 2nd. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
lines(d2, benford_2)
text(2, 0.05, "n = 894", cex=0.9)
points(d2,se_share_80, col="blue", pch=8)
points(d2,tstat_share_80, col="black", pch=17)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d2,coef_share_90, col="red", xlab="2nd. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(0,9,9), ylim=c(0.0,0.4), main="Figure 10 (c): 1990-1999
Relative Frequency of 2nd. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
lines(d2, benford_2)
text(2,0.05, "n = 2,740", cex=0.9)
points(d2,se_share_90, col="blue", pch=8)
points(d2,tstat_share_90, col="black", pch=17)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d2,coef_share_00, col="red", xlab="2nd. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(0,9,9), ylim=c(0.0,0.4), main="Figure 10 (d): 2000-2009
Relative Frequency of 2nd. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
lines(d2, benford_2)
text(2,0.05, "n = 3,952", cex=0.9)
points(d2,se_share_00, col="blue", pch=8)
points(d2,tstat_share_00, col="black", pch=17)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d2,coef_share_10, col="red", xlab="2nd. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(0,9,9), ylim=c(0.0,0.4), main="Figure 10 (e): 2010-2019
Relative Frequency of 2nd. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
lines(d2, benford_2)
text(2,0.05, "n = 10,661", cex=0.9)
points(d2,se_share_10, col="blue", pch=8)
points(d2,tstat_share_10, col="black", pch=17)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d2,coef_share_20, col="red", xlab="2nd. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(0,9,9), ylim=c(0.0,0.4), main="Figure 10 (f): 2020-2025
Relative Frequency of 2nd. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
lines(d2, benford_2)
text(2, 0.05, "n = 6,242", cex=0.9)
points(d2,se_share_20, col="blue", pch=8)
points(d2,tstat_share_20, col="black", pch=17)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

#######################################################

# NOW SELECT RANDOM SAMPLES OF SIZES 500, 1,000, 2,000 
#######
######################################################%
# First Digits
##############
par(mfrow=c(1,1))
set.seed(123)

rand_c1000_1<- sample(coef_total,1000)
rand_s1000_1<- sample(se_total,1000)
rand_t1000_1<- sample(tstat_total,1000)
rand_c2000_1<- sample(coef_total,2000)
rand_s2000_1<- sample(se_total,2000)
rand_t2000_1<- sample(tstat_total,2000)
rand_c500_1<- sample(coef_total,500)
rand_s500_1<- sample(se_total,500)
rand_t500_1<- sample(tstat_total,500)

bfd_c1000_1<- benford(rand_c1000_1,1,sign="both")
bfd_s1000_1<- benford(rand_s1000_1,1,sign="both")
bfd_t1000_1<- benford(rand_t1000_1,1,sign="both")
bfd_c2000_1<- benford(rand_c2000_1,1,sign="both")
bfd_s2000_1<- benford(rand_s2000_1,1,sign="both")
bfd_t2000_1<- benford(rand_t2000_1,1,sign="both")
bfd_c500_1<- benford(rand_c500_1,1,sign="both")
bfd_s500_1<- benford(rand_s500_1,1,sign="both")
bfd_t500_1<- benford(rand_t500_1,1,sign="both")

plot(d1,bfd_c500_1$bfd$data.dist, col="red", xlab="1st. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(1,9,8), ylim=c(0,0.35), main="Figure B.1 (a): Random sample, n = 500
Relative Frequency of 1st. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
points(d1,bfd_s500_1$bfd$data.dist, col="blue", pch=8)
points(d1,bfd_t500_1$bfd$data.dist, col="black", pch=17)
lines(d1, bfd_c500_1$bfd$benford.dist)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d1,bfd_c1000_1$bfd$data.dist, col="red", xlab="1st. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(1,9,8), ylim=c(0,0.35), main="Figure B.2 (a): Random sample, n = 1,000
Relative Frequency of 1st. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
points(d1,bfd_s1000_1$bfd$data.dist, col="blue", pch=8)
points(d1,bfd_t1000_1$bfd$data.dist, col="black", pch=17)
lines(d1, bfd_c1000_1$bfd$benford.dist)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d1,bfd_c2000_1$bfd$data.dist, col="red", xlab="1st. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(1,9,8), ylim=c(0,0.35), main="Figure B.3 (a): Random sample, n = 2000
Relative Frequency of 1st. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
points(d1,bfd_s2000_1$bfd$data.dist, col="blue", pch=8)
points(d1,bfd_t2000_1$bfd$data.dist, col="black", pch=17)
lines(d1, bfd_c2000_1$bfd$benford.dist)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)


#####################
#
#SECOND DIGITS
#
d2<- 0:9
benford_2<- c(0.1197,0.1139,0.1088,0.1043,0.1,0.0967, 0.0934,0.0904,0.0876,0.0850)
coef_share1000_1<- c()
se_share1000_1<- c()
tstat_share1000_1<- c()
coef_share2000_1<- c()
se_share2000_1<- c()
tstat_share2000_1<- c()
coef_share500_1<- c()
se_share500_1<- c()
tstat_share500_1<- c()

coef_2_ch <- str_sub(as.character(abs(rand_c1000_1)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(rand_s1000_1), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(rand_t1000_1)), 2, 2)
# Convert back to numeric
coef_2<- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_share1000_1[ii] <- sum(coef_2 == i1)/length(coef_2)
se_share1000_1[ii] <- sum(se_2 == i1)/length(se_2)
tstat_share1000_1[ii] <- sum(tstat_2 == i1)/length(tstat_2)
}

coef_2_ch <- str_sub(as.character(abs(rand_c2000_1)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(rand_s2000_1), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(rand_t2000_1)), 2, 2)
# Convert back to numeric
coef_2<- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_share2000_1[ii] <- sum(coef_2 == i1)/length(coef_2)
se_share2000_1[ii] <- sum(se_2 == i1)/length(se_2)
tstat_share2000_1[ii] <- sum(tstat_2 == i1)/length(tstat_2)
}


coef_2_ch <- str_sub(as.character(abs(rand_c500_1)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(rand_s500_1), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(rand_t500_1)), 2, 2)
# Convert back to numeric
coef_2<- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
coef_share500_1[ii] <- sum(coef_2 == i1)/length(coef_2)
se_share500_1[ii] <- sum(se_2 == i1)/length(se_2)
tstat_share500_1[ii] <- sum(tstat_2 == i1)/length(tstat_2)
}

par(mfrow=c(1,1))

plot(d2,coef_share500_1, col="red", xlab="2nd. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(0,9,9), ylim=c(0.0,0.4), main="Figure B.1 (b): Random sample, n = 500
Relative Frequency of 2nd. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
lines(d2, benford_2)
points(d2,se_share500_1, col="blue", pch=8)
points(d2,tstat_share500_1, col="black", pch=17)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d2,coef_share1000_1, col="red", xlab="2nd. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(0,9,9), ylim=c(0.0,0.4), main="Figure 14 (b): Random sample, n = 1,000
Relative Frequency of 2nd. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
lines(d2, benford_2)
points(d2,se_share1000_1, col="blue", pch=8)
points(d2,tstat_share1000_1, col="black", pch=17)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)

plot(d2,coef_share2000_1, col="red", xlab="2nd. Digit", ylab = "Relative frequency",
pch=19, xaxp=c(0,9,9), ylim=c(0.0,0.4), main="Figure B.3 (b): Random sample, n = 2,000
Relative Frequency of 2nd. Digits of Coefficients,
Standard Errors, & t-Statistics",cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
lines(d2, benford_2)
points(d2,se_share2000_1, col="blue", pch=8)
points(d2,tstat_share2000_1, col="black", pch=17)
legend(6.5,0.35, legend=c("Benford's law", "coefficients","standard errors", "t-statistics"),
       col=c("black", "red", "blue","black"), lty=c(1,NA,NA,NA),pch=c(NA,19,8,17), cex=0.8)
########################################################################