# Confidence Intervals for Random sample data
#############################################
library(benford.analysis)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(stringr)
library(readxl)

set.seed(123)

coef<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/Digits_All.xlsx",sheet="coef")
coef<- coef[!is.na(coef)]
summary(coef)
se<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/Digits_All.xlsx",sheet="se")
se<- se[!is.na(se)]
summary(se)
tstat<- read_excel("C:/Users/OEM/Sync/Benford Law/NZEP/Digits_All.xlsx",sheet="tstat")
tstat<- tstat[!is.na(tstat)]
summary(tstat)

# Sample size = 500
# =================

xc<- sample(coef,500)
xs<- sample(se,500)
xt<- sample(tstat,500)
nc<- ns<- nt<- 500
bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p<- bfdc$bfd$benford.dist      # Benford's first digit distribution
phat_c1<- bfdc$bfd$data.dist
phat_s1<- bfds$bfd$data.dist
phat_t1<- bfdt$bfd$data.dist

# Confidence Intervals (Goodman)

B95<- qchisq((0.05/9),1,lower.tail=FALSE) # 95%
B90<- qchisq((0.1/9), 1, lower.tail=FALSE) # 90%
B99<- qchisq((0.01/9), 1, lower.tail=FALSE) # 99%
f_c<- nc*phat_c1
LB95_c1<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c1<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))

LB95_c1
UB95_c1
p


LB90_c1<- (  B90+2*f_c - (B90 * (B90+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B90))
UB90_c1<-  (  B90+2*f_c + (B90 * (B90+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B90))

LB90_c1
UB90_c1
p

LB99_c1<- (  B99+2*f_c - (B99 * (B99+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B99))
UB99_c1<-  (  B99+2*f_c + (B99 * (B99+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B99))

LB99_c1
UB99_c1
p

f_s<- ns*phat_s1
LB95_s1<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s1<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))

LB95_s1
UB95_s1
p            


LB90_s1<- (  B90+2*f_s - (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))
UB90_s1<- (  B90+2*f_s + (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))

LB90_s1
UB90_s1
p              

LB99_s1<- (  B99+2*f_s - (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))
UB99_s1<- (  B99+2*f_s + (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))

LB99_s1
UB99_s1
p            

f_t<- nt*phat_t1
LB95_t1<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t1<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

LB95_t1
UB95_t1
p            


LB90_t1<- (  B90+2*f_t - (B90 * (B90+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B90))
UB90_t1<- (  B90+2*f_t + (B90 * (B90+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B90))

LB90_t1
UB90_t1
p              

LB99_t1<- (  B99+2*f_t - (B99 * (B99+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B99))
UB99_t1<- (  B99+2*f_t + (B99 * (B99+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B99))

LB99_t1
UB99_t1
p      
# Second digits
###############

digits<- c(0,1,2,3,4,5,6,7,8,9)
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
# Confidence Intervals (Goodman)

B95<- qchisq((0.05/10),1,lower.tail=FALSE) # 95%
B90<- qchisq((0.1/10), 1, lower.tail=FALSE) # 90%
B99<- qchisq((0.01/10), 1, lower.tail=FALSE) # 99%

f_c<- nc*phat_c2
LB95_c2<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c2<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))

LB95_c2
UB95_c2
p2


LB90_c2<- (  B90+2*f_c - (B90 * (B90+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B90))
UB90_c2<-  (  B90+2*f_c + (B90 * (B90+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B90))

LB90_c2
UB90_c2
p2

LB99_c2<- (  B99+2*f_c - (B99 * (B99+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B99))
UB99_c2<-  (  B99+2*f_c + (B99 * (B99+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B99))

LB99_c2
UB99_c2
p2

f_s<- ns*phat_s2
LB95_s2<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s2<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))

LB95_s2
UB95_s2
p2           


LB90_s2<- (  B90+2*f_s - (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))
UB90_s2<- (  B90+2*f_s + (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))

LB90_s2
UB90_s2
p2              

LB99_s2<- (  B99+2*f_s - (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))
UB99_s2<- (  B99+2*f_s + (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))

LB99_s2
UB99_s2
p2
             
f_t<- nt*phat_t2
LB95_t2<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t2<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

LB95_t2
UB95_t2
p2            


LB90_t2<- (  B90+2*f_t - (B90 * (B90+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B90))
UB90_t2<- (  B90+2*f_t + (B90 * (B90+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B90))

LB90_t2
UB90_t2
p2              

LB99_t2<- (  B99+2*f_t - (B99 * (B99+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B99))
UB99_t2<- (  B99+2*f_t + (B99 * (B99+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B99))

LB99_t2
UB99_t2
p2        

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
ggtitle("Figure C.1 (a): 95% Goodman Confidence Intervals
for Coefficient 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_c1, ymax=UB95_c1),width=.2, col="red")+
  scale_y_continuous(limits = c(0, 0.4))+ 
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 500)")

###                                                
graph_c2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure C.1 (b): 95% Goodman Confidence Intervals
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 500)")

ggarrange(graph_c1, graph_c2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
 
# Digits for std. errors
# ----------------------

graph_s1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure C.2 (a): 95% Goodman Confidence Intervals
for Std. Error 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
  plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_s1, ymax=UB95_s1), width=.2,col="blue")+
scale_y_continuous(limits = c(0, 0.4))+
  scale_x_continuous(breaks=seq(1,9,1))+
 scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 500)")

###     
graph_s2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure C.2 (b): 95% Goodman Confidence Intervals 
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 500)")
ggarrange(graph_s1, graph_s2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# Digits for t-stats.
# ------------------

graph_t1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure C.3 (a): 95% Goodman Confidence Intervals 
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 500)")

###             
graph_t2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure C.3 (b): 95% Goodman Confidence Intervals 
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 500)")

ggarrange(graph_t1, graph_t2, ncol=1, nrow=2, common.legend = TRUE, legend="right")


# Sample size = 1,00
# =================

xc<- sample(coef,1000)
xs<- sample(se,1000)
xt<- sample(tstat,1000)
nc<- ns<- nt<- 1000
bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p<- bfdc$bfd$benford.dist      # Benford's first digit distribution
phat_c1<- bfdc$bfd$data.dist
phat_s1<- bfds$bfd$data.dist
phat_t1<- bfdt$bfd$data.dist

# Confidence Intervals (Goodman)

B95<- qchisq((0.05/9),1,lower.tail=FALSE) # 95%
B90<- qchisq((0.1/9), 1, lower.tail=FALSE) # 90%
B99<- qchisq((0.01/9), 1, lower.tail=FALSE) # 99%
f_c<- nc*phat_c1
LB95_c1<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c1<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))

LB95_c1
UB95_c1
p


LB90_c1<- (  B90+2*f_c - (B90 * (B90+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B90))
UB90_c1<-  (  B90+2*f_c + (B90 * (B90+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B90))

LB90_c1
UB90_c1
p

LB99_c1<- (  B99+2*f_c - (B99 * (B99+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B99))
UB99_c1<-  (  B99+2*f_c + (B99 * (B99+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B99))

LB99_c1
UB99_c1
p

f_s<- ns*phat_s1
LB95_s1<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s1<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))

LB95_s1
UB95_s1
p            


LB90_s1<- (  B90+2*f_s - (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))
UB90_s1<- (  B90+2*f_s + (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))

LB90_s1
UB90_s1
p              

LB99_s1<- (  B99+2*f_s - (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))
UB99_s1<- (  B99+2*f_s + (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))

LB99_s1
UB99_s1
p            

f_t<- nt*phat_t1
LB95_t1<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t1<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

LB95_t1
UB95_t1
p            


LB90_t1<- (  B90+2*f_t - (B90 * (B90+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B90))
UB90_t1<- (  B90+2*f_t + (B90 * (B90+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B90))

LB90_t1
UB90_t1
p              

LB99_t1<- (  B99+2*f_t - (B99 * (B99+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B99))
UB99_t1<- (  B99+2*f_t + (B99 * (B99+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B99))

LB99_t1
UB99_t1
p      
# Second digits
###############

digits<- c(0,1,2,3,4,5,6,7,8,9)
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
# Confidence Intervals (Goodman)

B95<- qchisq((0.05/10),1,lower.tail=FALSE) # 95%
B90<- qchisq((0.1/10), 1, lower.tail=FALSE) # 90%
B99<- qchisq((0.01/10), 1, lower.tail=FALSE) # 99%

f_c<- nc*phat_c2
LB95_c2<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c2<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))

LB95_c2
UB95_c2
p2


LB90_c2<- (  B90+2*f_c - (B90 * (B90+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B90))
UB90_c2<-  (  B90+2*f_c + (B90 * (B90+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B90))

LB90_c2
UB90_c2
p2

LB99_c2<- (  B99+2*f_c - (B99 * (B99+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B99))
UB99_c2<-  (  B99+2*f_c + (B99 * (B99+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B99))

LB99_c2
UB99_c2
p2

f_s<- ns*phat_s2
LB95_s2<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s2<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))

LB95_s2
UB95_s2
p2           


LB90_s2<- (  B90+2*f_s - (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))
UB90_s2<- (  B90+2*f_s + (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))

LB90_s2
UB90_s2
p2              

LB99_s2<- (  B99+2*f_s - (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))
UB99_s2<- (  B99+2*f_s + (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))

LB99_s2
UB99_s2
p2
             
f_t<- nt*phat_t2
LB95_t2<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t2<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

LB95_t2
UB95_t2
p2            


LB90_t2<- (  B90+2*f_t - (B90 * (B90+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B90))
UB90_t2<- (  B90+2*f_t + (B90 * (B90+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B90))

LB90_t2
UB90_t2
p2              

LB99_t2<- (  B99+2*f_t - (B99 * (B99+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B99))
UB99_t2<- (  B99+2*f_t + (B99 * (B99+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B99))

LB99_t2
UB99_t2
p2        

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
ggtitle("Figure C.4 (a): 95% Goodman Confidence Intervals
for Coefficient 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_c1, ymax=UB95_c1),width=.2, col="red")+
  scale_y_continuous(limits = c(0, 0.4))+ 
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 1,000)")

###                                                
graph_c2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure C.4 (b): 95% Goodman Confidence Intervals
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 1,000)")

ggarrange(graph_c1, graph_c2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
 
# Digits for std. errors
# ----------------------

graph_s1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure C.5 (a): 95% Goodman Confidence Intervals
for Std. Error 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
  plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_s1, ymax=UB95_s1), width=.2,col="blue")+
scale_y_continuous(limits = c(0, 0.4))+
  scale_x_continuous(breaks=seq(1,9,1))+
 scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 1,000)")

###     
graph_s2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure C.5 (b): 95% Goodman Confidence Intervals 
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 1,000)")
ggarrange(graph_s1, graph_s2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# Digits for t-stats.
# ------------------

graph_t1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure C.6 (a): 95% Goodman Confidence Intervals 
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 1,000)")

###             
graph_t2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure C.6 (b): 95% Goodman Confidence Intervals 
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 1,000)")

ggarrange(graph_t1, graph_t2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# Sample size = 2,000
# =================

xc<- sample(coef,2000)
xs<- sample(se,2000)
xt<- sample(tstat,2000)
nc<- ns<- nt<- 2000
bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p<- bfdc$bfd$benford.dist      # Benford's first digit distribution
phat_c1<- bfdc$bfd$data.dist
phat_s1<- bfds$bfd$data.dist
phat_t1<- bfdt$bfd$data.dist

# Confidence Intervals (Goodman)

B95<- qchisq((0.05/9),1,lower.tail=FALSE) # 95%
B90<- qchisq((0.1/9), 1, lower.tail=FALSE) # 90%
B99<- qchisq((0.01/9), 1, lower.tail=FALSE) # 99%
f_c<- nc*phat_c1
LB95_c1<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c1<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))

LB95_c1
UB95_c1
p


LB90_c1<- (  B90+2*f_c - (B90 * (B90+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B90))
UB90_c1<-  (  B90+2*f_c + (B90 * (B90+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B90))

LB90_c1
UB90_c1
p

LB99_c1<- (  B99+2*f_c - (B99 * (B99+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B99))
UB99_c1<-  (  B99+2*f_c + (B99 * (B99+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B99))

LB99_c1
UB99_c1
p

f_s<- ns*phat_s1
LB95_s1<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s1<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))

LB95_s1
UB95_s1
p            


LB90_s1<- (  B90+2*f_s - (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))
UB90_s1<- (  B90+2*f_s + (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))

LB90_s1
UB90_s1
p              

LB99_s1<- (  B99+2*f_s - (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))
UB99_s1<- (  B99+2*f_s + (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))

LB99_s1
UB99_s1
p            

f_t<- nt*phat_t1
LB95_t1<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t1<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

LB95_t1
UB95_t1
p            


LB90_t1<- (  B90+2*f_t - (B90 * (B90+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B90))
UB90_t1<- (  B90+2*f_t + (B90 * (B90+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B90))

LB90_t1
UB90_t1
p              

LB99_t1<- (  B99+2*f_t - (B99 * (B99+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B99))
UB99_t1<- (  B99+2*f_t + (B99 * (B99+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B99))

LB99_t1
UB99_t1
p      
# Second digits
###############

digits<- c(0,1,2,3,4,5,6,7,8,9)
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
# Confidence Intervals (Goodman)

B95<- qchisq((0.05/10),1,lower.tail=FALSE) # 95%
B90<- qchisq((0.1/10), 1, lower.tail=FALSE) # 90%
B99<- qchisq((0.01/10), 1, lower.tail=FALSE) # 99%

f_c<- nc*phat_c2
LB95_c2<- (  B95+2*f_c - (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))
UB95_c2<-  (  B95+2*f_c + (B95 * (B95+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B95))

LB95_c2
UB95_c2
p2


LB90_c2<- (  B90+2*f_c - (B90 * (B90+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B90))
UB90_c2<-  (  B90+2*f_c + (B90 * (B90+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B90))

LB90_c2
UB90_c2
p2

LB99_c2<- (  B99+2*f_c - (B99 * (B99+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B99))
UB99_c2<-  (  B99+2*f_c + (B99 * (B99+4*f_c*(nc-f_c)/nc  ))^0.5)/(2*(nc+B99))

LB99_c2
UB99_c2
p2

f_s<- ns*phat_s2
LB95_s2<- (  B95+2*f_s - (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))
UB95_s2<- (  B95+2*f_s + (B95 * (B95+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B95))

LB95_s2
UB95_s2
p2           


LB90_s2<- (  B90+2*f_s - (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))
UB90_s2<- (  B90+2*f_s + (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))

LB90_s2
UB90_s2
p2              

LB99_s2<- (  B99+2*f_s - (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))
UB99_s2<- (  B99+2*f_s + (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))

LB99_s2
UB99_s2
p2
             
f_t<- nt*phat_t2
LB95_t2<- (  B95+2*f_t - (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t2<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B95))

LB95_t2
UB95_t2
p2            


LB90_t2<- (  B90+2*f_t - (B90 * (B90+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B90))
UB90_t2<- (  B90+2*f_t + (B90 * (B90+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B90))

LB90_t2
UB90_t2
p2              

LB99_t2<- (  B99+2*f_t - (B99 * (B99+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B99))
UB99_t2<- (  B99+2*f_t + (B99 * (B99+4*f_t*(nt-f_t)/nt  ))^0.5)/(2*(nt+B99))

LB99_t2
UB99_t2
p2        

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
ggtitle("Figure C.7 (a): 95% Goodman Confidence Intervals
for Coefficient 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_c1, ymax=UB95_c1),width=.2, col="red")+
  scale_y_continuous(limits = c(0, 0.4))+ 
scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 2,000)")

###                                                
graph_c2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure C.7 (b): 95% Goodman Confidence Intervals
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n =2,000)")

ggarrange(graph_c1, graph_c2, ncol=1, nrow=2, common.legend = TRUE, legend="right")
 
# Digits for std. errors
# ----------------------

graph_s1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure C.8 (a): 95% Goodman Confidence Intervals
for Std. Error 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
  plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_s1, ymax=UB95_s1), width=.2,col="blue")+
scale_y_continuous(limits = c(0, 0.4))+
  scale_x_continuous(breaks=seq(1,9,1))+
 scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 2,000)")

###     
graph_s2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure C.8 (b): 95% Goodman Confidence Intervals 
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 2,000)")
ggarrange(graph_s1, graph_s2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# Digits for t-stats.
# ------------------

graph_t1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure C.9 (a): 95% Goodman Confidence Intervals 
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 2,000)")

###             
graph_t2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure C.9 (b): 95% Goodman Confidence Intervals 
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Random sample
(n = 2,000)")

ggarrange(graph_t1, graph_t2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

##############################################



