# Confidence Intervals for ALL decade sub-samples
# for the first and second digit cases
#############################################
#
# Code written by David Giles for the paper "Benford’s Law and Regression Results Published in
# Articles in New Zealand Economic Papers", last updated June 2026.
#
# Contact: David Giles; davegiles1949@gmail.com; davegiles.ca
# -------
############################################################

library(benford.analysis)
library(ggplot2)
library(ggpubr)
library(gridExtra)
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


# set up data for the ALL of the decade sub-samples
#
# DECADE 1966-1979
#================= 

# First digits
# ------------

xc<- abs(coef_6070)
xs<- se_6070
xt<- abs(tstat_6070)
nc<- length(xc)
ns<- length(xs)
nt<- length(xt)
bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p<- bfdc$bfd$benford.dist      # Benford's first digit distribution
phat_c1<- bfdc$bfd$data.dist
phat_s1<- bfds$bfd$data.dist
phat_t1<- bfdt$bfd$data.dist

# 95% Confidence Intervals (Goodman)

B95<- qchisq((0.05/9),1,lower.tail=FALSE) # 95%
f_c<- nc*phat_c1
LB95_c1<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c1<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
f_s<- ns*phat_s1
LB95_s1<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s1<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
f_t<- nt*phat_t1
LB95_t1<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t1<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

# Second digits: 
###############

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
# 95% Confidence Intervals (Goodman)

B95<- qchisq((0.05/10),1,lower.tail=FALSE) # 95%
f_c<- nc*phat_c2
LB95_c2<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c2<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
f_s<- ns*phat_s2
LB95_s2<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s2<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))          
f_t<- nt*phat_t2
LB95_t2<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t2<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

#########################
# PLOTS:
#######

d1<- 1:9
d2<-0:9
df1<- data.frame(d1,p)
df2<- data.frame(d2,p2)

# Digits for coefficients
# -----------------------

graph_c1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.1(a): 95% Goodman Confidence Intervals
for Coefficient 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_c1, ymax=UB95_c1),width=.2, col="red")+
  scale_y_continuous(limits = c(0, 0.37))+ 
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1966-1979 data
(n = 3,013)")

###                                                
graph_c2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.1(b): 95% Goodman Confidence Intervals
for Coefficient 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_c2, ymax=UB95_c2), width=.2, col="red")+
scale_y_continuous(limits = c(0, 0.35))+ 
scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1966-1979 data
(n = 2,890)")

ggarrange(graph_c1, graph_c2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
 
# Digits for std. errors
# ----------------------

graph_s1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.1(c): 95% Goodman Confidence Intervals
for Std. Error 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
  plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_s1, ymax=UB95_s1), width=.2,col="blue")+
scale_y_continuous(limits = c(0, 0.43))+
  scale_x_continuous(breaks=seq(1,9,1))+
 scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1966-1979 data
(n = 1,258)")

###     
graph_s2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.1(d): 95% Goodman Confidence Intervals 
for Std. Error 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_s2, ymax=UB95_s2), width=.2, col="blue")+
scale_y_continuous(limits = c(0, 0.35))+
 scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1966-1979 data
(n = 1,177)")
ggarrange(graph_s1, graph_s2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# Digits for t-stats.
# ------------------

graph_t1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.1(e): 95% Goodman Confidence Intervals 
for t-Stat. 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_t1, ymax=UB95_t1), width=.2)+
scale_y_continuous(limits = c(0, 0.35))+
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1966-1979 data
(n = 1,232)")

###             
graph_t2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.1(f): 95% Goodman Confidence Intervals 
for t-Stat. 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_t2, ymax=UB95_t2), width=.2)+
scale_y_continuous(limits = c(0, 0.35))+
scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1966-1979 data
(n = 1,211)")
ggarrange(graph_t1, graph_t2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# DECADE 1980-1989
#================= 

# First digits
# ------------

xc<- abs(coef_80)
xs<- se_80
xt<- abs(tstat_80)
nc<- length(xc)
ns<- length(xs)
nt<- length(xt)
bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p<- bfdc$bfd$benford.dist      # Benford's first digit distribution
phat_c1<- bfdc$bfd$data.dist
phat_s1<- bfds$bfd$data.dist
phat_t1<- bfdt$bfd$data.dist

# 95% Confidence Intervals (Goodman)

B95<- qchisq((0.05/9),1,lower.tail=FALSE) # 95%
f_c<- nc*phat_c1
LB95_c1<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c1<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
f_s<- ns*phat_s1
LB95_s1<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s1<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
f_t<- nt*phat_t1
LB95_t1<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t1<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

# Second digits: 
###############

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
# 95% Confidence Intervals (Goodman)

B95<- qchisq((0.05/10),1,lower.tail=FALSE) # 95%
f_c<- nc*phat_c2
LB95_c2<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c2<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
f_s<- ns*phat_s2
LB95_s2<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s2<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))          
f_t<- nt*phat_t2
LB95_t2<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t2<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

#########################
# PLOTS:
#######

d1<- 1:9
d2<-0:9
df1<- data.frame(d1,p)
df2<- data.frame(d2,p2)

# Digits for coefficients
# -----------------------

graph_c1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.2(a): 95% Goodman Confidence Intervals
for Coefficient 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_c1, ymax=UB95_c1),width=.2, col="red")+
  scale_y_continuous(limits = c(0, 0.37))+ 
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1980-1989 data
(n = 950)")

###                                                
graph_c2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.2(b): 95% Goodman Confidence Intervals
for Coefficient 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_c2, ymax=UB95_c2), width=.2, col="red")+
scale_y_continuous(limits = c(0, 0.35))+ 
scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1980-1989 data
(n = 894)")

ggarrange(graph_c1, graph_c2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
 
# Digits for std. errors
# ----------------------

graph_s1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.2(c): 95% Goodman Confidence Intervals
for Std. Error 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
  plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_s1, ymax=UB95_s1), width=.2,col="blue")+
scale_y_continuous(limits = c(0, 0.43))+
  scale_x_continuous(breaks=seq(1,9,1))+
 scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1980-1989 data
(n = 449)")

###     
graph_s2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.2(d): 95% Goodman Confidence Intervals 
for Std. Error 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_s2, ymax=UB95_s2), width=.2, col="blue")+
scale_y_continuous(limits = c(0, 0.35))+
 scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1980-1989 data
(n = 389)")
ggarrange(graph_s1, graph_s2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# Digits for t-stats.
# ------------------

graph_t1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.2(e): 95% Goodman Confidence Intervals 
for t-Stat. 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_t1, ymax=UB95_t1), width=.2)+
scale_y_continuous(limits = c(0, 0.4))+
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1980-1989 data
(n = 261)")

###             
graph_t2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.2(f): 95% Goodman Confidence Intervals 
for t-Stat. 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_t2, ymax=UB95_t2), width=.2)+
scale_y_continuous(limits = c(0, 0.35))+
scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1980-1989 data
(n = 251)")

ggarrange(graph_t1, graph_t2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
#

# DECADE 1990-1999
#================= 

# First digits
# ------------

xc<- abs(coef_90)
xs<- se_90
xt<- abs(tstat_90)
nc<- length(xc)
ns<- length(xs)
nt<- length(xt)
bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p<- bfdc$bfd$benford.dist      # Benford's first digit distribution
phat_c1<- bfdc$bfd$data.dist
phat_s1<- bfds$bfd$data.dist
phat_t1<- bfdt$bfd$data.dist

# 95% Confidence Intervals (Goodman)

B95<- qchisq((0.05/9),1,lower.tail=FALSE) # 95%
f_c<- nc*phat_c1
LB95_c1<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c1<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
f_s<- ns*phat_s1
LB95_s1<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s1<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
f_t<- nt*phat_t1
LB95_t1<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t1<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

# Second digits: 
###############

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
# 95% Confidence Intervals (Goodman)

B95<- qchisq((0.05/10),1,lower.tail=FALSE) # 95%
f_c<- nc*phat_c2
LB95_c2<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c2<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
f_s<- ns*phat_s2
LB95_s2<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s2<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))          
f_t<- nt*phat_t2
LB95_t2<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t2<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

#########################
# PLOTS:
#######

d1<- 1:9
d2<-0:9
df1<- data.frame(d1,p)
df2<- data.frame(d2,p2)

# Digits for coefficients
# -----------------------

graph_c1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.3(a): 95% Goodman Confidence Intervals
for Coefficient 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_c1, ymax=UB95_c1),width=.2, col="red")+
  scale_y_continuous(limits = c(0, 0.37))+ 
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1990-1999 data
(n = 2,996)")

###                                                
graph_c2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.3(b): 95% Goodman Confidence Intervals
for Coefficient 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_c2, ymax=UB95_c2), width=.2, col="red")+
scale_y_continuous(limits = c(0, 0.35))+ 
scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1990-1999 data
(n = 2,740)")

ggarrange(graph_c1, graph_c2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
 
# Digits for std. errors
# ----------------------

graph_s1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.3(c): 95% Goodman Confidence Intervals
for Std. Error 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
  plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_s1, ymax=UB95_s1), width=.2,col="blue")+
scale_y_continuous(limits = c(0, 0.43))+
  scale_x_continuous(breaks=seq(1,9,1))+
 scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1990-1999 data
(n = 910)")

###     
graph_s2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.3(d): 95% Goodman Confidence Intervals 
for Std. Error 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_s2, ymax=UB95_s2), width=.2, col="blue")+
scale_y_continuous(limits = c(0, 0.35))+
 scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1990-1999 data
(n = 864)")
ggarrange(graph_s1, graph_s2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# Digits for t-stats.
# ------------------

graph_t1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.3(e): 95% Goodman Confidence Intervals 
for t-Stat. 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_t1, ymax=UB95_t1), width=.2)+
scale_y_continuous(limits = c(0, 0.35))+
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1990-1999 data
(n = 1,685)")

###             
graph_t2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.3(f): 95% Goodman Confidence Intervals 
for t-Stat. 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_t2, ymax=UB95_t2), width=.2)+
scale_y_continuous(limits = c(0, 0.35))+
scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "1990-1999 data
(n = 1,656)")

ggarrange(graph_t1, graph_t2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
#

# DECADE 2000-2009
#================= 

# First digits
# ------------

xc<- abs(coef_00)
xs<- se_00
xt<- abs(tstat_00)
nc<- length(xc)
ns<- length(xs)
nt<- length(xt)
bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p<- bfdc$bfd$benford.dist      # Benford's first digit distribution
phat_c1<- bfdc$bfd$data.dist
phat_s1<- bfds$bfd$data.dist
phat_t1<- bfdt$bfd$data.dist

# 95% Confidence Intervals (Goodman)

B95<- qchisq((0.05/9),1,lower.tail=FALSE) # 95%
f_c<- nc*phat_c1
LB95_c1<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c1<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
f_s<- ns*phat_s1
LB95_s1<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s1<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
f_t<- nt*phat_t1
LB95_t1<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t1<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

# Second digits: 
###############

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
# 95% Confidence Intervals (Goodman)

B95<- qchisq((0.05/10),1,lower.tail=FALSE) # 95%
f_c<- nc*phat_c2
LB95_c2<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c2<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
f_s<- ns*phat_s2
LB95_s2<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s2<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))          
f_t<- nt*phat_t2
LB95_t2<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t2<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

#########################
# PLOTS:
#######

d1<- 1:9
d2<-0:9
df1<- data.frame(d1,p)
df2<- data.frame(d2,p2)

# Digits for coefficients
# -----------------------

graph_c1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.4(a): 95% Goodman Confidence Intervals
for Coefficient 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_c1, ymax=UB95_c1),width=.2, col="red")+
  scale_y_continuous(limits = c(0, 0.37))+ 
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2000-2009 data
(n = 4,264)")

###                                                
graph_c2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.4(b): 95% Goodman Confidence Intervals
for Coefficient 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_c2, ymax=UB95_c2), width=.2, col="red")+
scale_y_continuous(limits = c(0, 0.35))+ 
scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2000-2009 data
(n = 3,952)")

ggarrange(graph_c1, graph_c2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
 
# Digits for std. errors
# ----------------------

graph_s1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.4(c): 95% Goodman Confidence Intervals
for Std. Error 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
  plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_s1, ymax=UB95_s1), width=.2,col="blue")+
scale_y_continuous(limits = c(0, 0.43))+
  scale_x_continuous(breaks=seq(1,9,1))+
 scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2000-2009 data
(n = 2,487)")

###     
graph_s2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.4(d): 95% Goodman Confidence Intervals 
for Std. Error 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_s2, ymax=UB95_s2), width=.2, col="blue")+
scale_y_continuous(limits = c(0, 0.35))+
 scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2000-2009 data
(n = 2,212)")
ggarrange(graph_s1, graph_s2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# Digits for t-stats.
# ------------------

graph_t1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.4(e): 95% Goodman Confidence Intervals 
for t-Stat. 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_t1, ymax=UB95_t1), width=.2)+
scale_y_continuous(limits = c(0, 0.35))+
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2000-2009 data
(n = 1,402)")

###             
graph_t2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.4(f): 95% Goodman Confidence Intervals 
for t-Stat. 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_t2, ymax=UB95_t2), width=.2)+
scale_y_continuous(limits = c(0, 0.35))+
scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2000-2009 data
(n = 1,377)")

ggarrange(graph_t1, graph_t2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
#

# DECADE 2010-2019
#================= 

# First digits
# ------------

xc<- abs(coef_10)
xs<- se_10
xt<- abs(tstat_10)
nc<- length(xc)
ns<- length(xs)
nt<- length(xt)
bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p<- bfdc$bfd$benford.dist      # Benford's first digit distribution
phat_c1<- bfdc$bfd$data.dist
phat_s1<- bfds$bfd$data.dist
phat_t1<- bfdt$bfd$data.dist

# 95% Confidence Intervals (Goodman)

B95<- qchisq((0.05/9),1,lower.tail=FALSE) # 95%
f_c<- nc*phat_c1
LB95_c1<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c1<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
f_s<- ns*phat_s1
LB95_s1<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s1<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
f_t<- nt*phat_t1
LB95_t1<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t1<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

# Second digits: 
###############

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
# 95% Confidence Intervals (Goodman)

B95<- qchisq((0.05/10),1,lower.tail=FALSE) # 95%
f_c<- nc*phat_c2
LB95_c2<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c2<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
f_s<- ns*phat_s2
LB95_s2<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s2<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))          
f_t<- nt*phat_t2
LB95_t2<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t2<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

#########################
# PLOTS:
#######

d1<- 1:9
d2<-0:9
df1<- data.frame(d1,p)
df2<- data.frame(d2,p2)

# Digits for coefficients
# -----------------------

graph_c1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.5(a): 95% Goodman Confidence Intervals
for Coefficient 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_c1, ymax=UB95_c1),width=.2, col="red")+
  scale_y_continuous(limits = c(0, 0.37))+ 
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2010-2019 data
(n = 11,744)")

###                                                
graph_c2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.5(b): 95% Goodman Confidence Intervals
for Coefficient 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_c2, ymax=UB95_c2), width=.2, col="red")+
scale_y_continuous(limits = c(0, 0.35))+ 
scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2010-2019 data
(n = 10,661)")

ggarrange(graph_c1, graph_c2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
 
# Digits for std. errors
# ----------------------

graph_s1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.5(c): 95% Goodman Confidence Intervals
for Std. Error 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
  plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_s1, ymax=UB95_s1), width=.2,col="blue")+
scale_y_continuous(limits = c(0, 0.43))+
  scale_x_continuous(breaks=seq(1,9,1))+
 scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2010-2019 data
(n = 5,224)")

###     
graph_s2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.5(d): 95% Goodman Confidence Intervals 
for Std. Error 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_s2, ymax=UB95_s2), width=.2, col="blue")+
scale_y_continuous(limits = c(0, 0.35))+
 scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2010-2019 data
(n = 4,710)")
ggarrange(graph_s1, graph_s2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# Digits for t-stats.
# ------------------

graph_t1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.5(e): 95% Goodman Confidence Intervals 
for t-Stat. 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_t1, ymax=UB95_t1), width=.2)+
scale_y_continuous(limits = c(0, 0.35))+
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2010-2019 data
(n = 5,225)")

###             
graph_t2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.5(f): 95% Goodman Confidence Intervals 
for t-Stat. 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_t2, ymax=UB95_t2), width=.2)+
scale_y_continuous(limits = c(0, 0.35))+
scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2010-2019 data
(n = 5.091)")

ggarrange(graph_t1, graph_t2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
#
# DECADE 2020-2025
#================= 

# First digits
# ------------

xc<- abs(coef_20)
xs<- se_20
xt<- abs(tstat_20)
nc<- length(xc)
ns<- length(xs)
nt<- length(xt)
bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p<- bfdc$bfd$benford.dist      # Benford's first digit distribution
phat_c1<- bfdc$bfd$data.dist
phat_s1<- bfds$bfd$data.dist
phat_t1<- bfdt$bfd$data.dist

# 95% Confidence Intervals (Goodman)

B95<- qchisq((0.05/9),1,lower.tail=FALSE) # 95%
f_c<- nc*phat_c1
LB95_c1<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c1<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
f_s<- ns*phat_s1
LB95_s1<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s1<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
f_t<- nt*phat_t1
LB95_t1<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t1<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

# Second digits: 
###############

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
# 95% Confidence Intervals (Goodman)

B95<- qchisq((0.05/10),1,lower.tail=FALSE) # 95%
f_c<- nc*phat_c2
LB95_c2<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c2<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
f_s<- ns*phat_s2
LB95_s2<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s2<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))          
f_t<- nt*phat_t2
LB95_t2<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t2<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

#########################
# PLOTS:
#######

d1<- 1:9
d2<-0:9
df1<- data.frame(d1,p)
df2<- data.frame(d2,p2)

# Digits for coefficients
# -----------------------

graph_c1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.6(a): 95% Goodman Confidence Intervals
for Coefficient 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_c1, ymax=UB95_c1),width=.2, col="red")+
  scale_y_continuous(limits = c(0, 0.37))+ 
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2020-2025 data
(n = 7,192)")

###                                                
graph_c2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.6(b): 95% Goodman Confidence Intervals
for Coefficient 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_c2, ymax=UB95_c2), width=.2, col="red")+
scale_y_continuous(limits = c(0, 0.35))+ 
scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2020-2025 data
(n = 6,242)")

ggarrange(graph_c1, graph_c2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
 
# Digits for std. errors
# ----------------------

graph_s1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.6(c): 95% Goodman Confidence Intervals
for Std. Error 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
  plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_s1, ymax=UB95_s1), width=.2,col="blue")+
scale_y_continuous(limits = c(0, 0.43))+
  scale_x_continuous(breaks=seq(1,9,1))+
 scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2020-2025 data
(n = 4,191)")

###     
graph_s2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.6(d): 95% Goodman Confidence Intervals 
for Std. Error 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_s2, ymax=UB95_s2), width=.2, col="blue")+
scale_y_continuous(limits = c(0, 0.35))+
 scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2020-2025 data
(n = 3,533)")
ggarrange(graph_s1, graph_s2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# Digits for t-stats.
# ------------------

graph_t1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure A.6(e): 95% Goodman Confidence Intervals 
for t-Stat. 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_t1, ymax=UB95_t1), width=.2)+
scale_y_continuous(limits = c(0, 0.4))+
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2020-2025 data
(n = 1,036)")

###             
graph_t2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure A.6(f): 95% Goodman Confidence Intervals 
for t-Stat. 2nd. Digits") +
  xlab("2nd. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=LB95_t2, ymax=UB95_t2), width=.2)+
scale_y_continuous(limits = c(0, 0.35))+
scale_x_continuous(breaks=seq(0,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(linetype = c(0), color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "2020-2025 data
(n = 1,020)")

ggarrange(graph_t1, graph_t2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
