# Tests and Goodman confidence intervals proposed by Lesperance et al.
# --------------------------------------------------------------------

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

# Construct Todter's M-test (p.342): asymptotically chi-square with 1 dof.
###########################

xc<- coef
xs<- se
xt<- tstat
nc<- length(xc)
ns<- length(xs)
nt<-length(xt)

sig1<-signifd(x = xc, digits = 1)
M_c<- nc*(mean(sig1)-3.440)^2/6.057
pval<- pchisq(M_c,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_c
pval
sig1<-signifd(x = xs, digits = 1)
M_s<- ns*(mean(sig1)-3.440)^2/6.057
pval<- pchisq(M_s,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_s
pval
sig1<-signifd(x = xt, digits = 1)
M_t<- nt*(mean(sig1)-3.440)^2/6.057
pval<- pchisq(M_t,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_t
pval

# Construct L'Esperance et al. tests
####################################

bfdc<- benford(xc,1,sign="both")
bfds<- benford(xs,1, sign="both")
bfdt<- benford(xt,1, sign="both")
p<- bfdc$bfd$benford.dist
phat_c1<- bfdc$bfd$data.dist
phat_s1<- bfds$bfd$data.dist
phat_t1<- bfdt$bfd$data.dist

#Chi-square Tests:
##################

chisq_c<- nc*sum((phat_c1-p)^2/p) 
chisq_s<- ns*sum((phat_s1-p)^2/p) 
chisq_t<- nt*sum((phat_t1-p)^2/p) 

chisq_c
pchisq(chisq_c,8, ncp = 0, lower.tail = FALSE, log.p = FALSE)
chisq_s
pchisq(chisq_s,8, ncp = 0, lower.tail = FALSE, log.p = FALSE)
chisq_t
pchisq(chisq_t,8, ncp = 0, lower.tail = FALSE, log.p = FALSE)

# U2d Tests:
###########

T<- cumsum(p)
S_c<- cumsum(phat_c1)
S_s<- cumsum(phat_s1)
S_t<- cumsum(phat_t1)
Z_c<- S_c-T
Z_s<- S_s-T
Z_t<- S_t-T
t<- c()
for (ii in 1:8) {
t[ii]<- (p[ii]+p[ii+1])/2
}
t[9]<- (p[9]+p[1])/2
Zbar_c<- sum(t*Z_c)
Zbar_s<- sum(t*Z_s)
Zbar_t<- sum(t*Z_t)

U2d_c<- nc*sum((Z_c-Zbar_c)^2*t)
U2d_s<- ns*sum((Z_s-Zbar_s)^2*t)
U2d_t<- nt*sum((Z_t-Zbar_t)^2*t)

U2d_c
U2d_s
U2d_t

# Asy. crit. vals.:        10%     5%       2.5%    1%
#                         0.163   0.205    0.247   0.304

W2d_c<- nc*sum(Z_c^2*t)
W2d_s<- ns*sum(Z_s^2*t)
W2d_t<- nt*sum(Z_t^2*t)

W2d_c
W2d_s
W2d_t

# Asy. crit. vals.:        10%     5%       2.5%    1%
#                         0.351   0.471    0.597   0.768


term_c<- (Z_c^2*t/(T-T^2))    # Last element in sum =0/0, so remove NA (i.e., set to 0)"
term_c<- term_c[!is.na(term_c)]
A2d_c<- nc*sum(term_c)
term_s<- (Z_s^2*t/(T-T^2))
term_s<- term_s[!is.na(term_s)]
A2d_s<- ns*sum(term_s)
term_t<- (Z_t^2*t/(T-T^2))
term_t<- term_t[!is.na(term_t)]
A2d_t<- nt*sum(term_t)

A2d_c
A2d_s
A2d_t

# Asy. crit. vals.:        10%     5%       2.5%    1%
#                         1.743   2.304    2.890   3.688



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
UB95_s1<- (  B95+2*f_s + (B95 * (B95+4*f_s*(nc-f_s)/ns  ))^0.5)/(2*(ns+B95))

LB95_s1
UB95_s1
p            


LB90_s1<- (  B90+2*f_s - (B90 * (B90+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B90))
UB90_s1<- (  B90+2*f_s + (B90 * (B90+4*f_s*(nc-f_s)/ns  ))^0.5)/(2*(ns+B90))

LB90_s1
UB90_s1
p              

LB99_s1<- (  B99+2*f_s - (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))
UB99_s1<- (  B99+2*f_s + (B99 * (B99+4*f_s*(nc-f_s)/ns  ))^0.5)/(2*(ns+B99))

LB99_s1
UB99_s1
p            

f_t<- nt*phat_t1
LB95_t1<- (  B95+2*f_t - (B95 * (B95+4*f_t*(ns-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t1<- (  B95+2*f_t + (B95 * (B95+4*f_t*(nc-f_t)/nt  ))^0.5)/(2*(nt+B95))

LB95_t1
UB95_t1
p            


LB90_t1<- (  B90+2*f_t - (B90 * (B90+4*f_t*(ns-f_t)/nt  ))^0.5)/(2*(nt+B90))
UB90_t1<- (  B90+2*f_t + (B90 * (B90+4*f_t*(nc-f_t)/nt  ))^0.5)/(2*(nt+B90))

LB90_t1
UB90_t1
p              

LB99_t1<- (  B99+2*f_t - (B99 * (B99+4*f_t*(ns-f_t)/nt  ))^0.5)/(2*(nt+B99))
UB99_t1<- (  B99+2*f_t + (B99 * (B99+4*f_t*(nc-f_t)/nt  ))^0.5)/(2*(nt+B99))

LB99_t1
UB99_t1
p      

##########################################
# SECOND DIGITS
#################

digits<- c(0,1,2,3,4,5,6,7,8,9)
benford_2<- c(0.1197,0.1139,0.1088,0.1043,0.1,0.0967, 0.0934,0.0904,0.0876,0.0850)
coef_2_ch <- str_sub(as.character(abs(coef)), 2, 2)   # Need to take abs. value to eliminate the "-" CHARACTERS
se_2_ch<- str_sub(as.character(se), 2, 2)
tstat_2_ch<- str_sub(as.character(abs(tstat)), 2, 2)    # Need to take abs. value to eliminate the "-" CHARACTERS
# Convert back to numeric
coef_2 <- as.numeric(coef_2_ch)
se_2<- as.numeric(se_2_ch)
tstat_2<- as.numeric(tstat_2_ch)
# Remove the "NAs"
coef_2 <- coef_2[!is.na(coef_2)]
se_2<- se_2[!is.na(se_2)]
tstat_2<- tstat_2[!is.na(tstat_2)]

# Construct Todter's M-test (p.342)
###########################

M_c<- length(coef_2)*(mean(coef_2)-4.187)^2/8.254
pval<- pchisq(M_c,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_c
pval
M_s<- length(se_2)*(mean(se_2)-4.187)^2/8.254
pval<- pchisq(M_s,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_s
pval
M_t<- length(tstat_2)*(mean(tstat_2)-4.187)^2/8.254
pval<- pchisq(M_t,1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
M_t
pval

# Construct L'Esperance et al. tests
####################################

xc<- coef_2
xs<- se_2
xt<- tstat_2
nc<- length(xc)
ns<- length(xs)
nt<-length(xt)
# Benford 2nd digit probs (for i = 0 to 9)
p2<- c(0.1197,0.1139,0.1088,0.1043,0.1,0.0967, 0.0934,0.0904,0.0876,0.0850)

phat_c2<- c()
phat_s2<- c()
phat_t2<- c()

for (ii in 1:10) {                                  
i1<- ii-1                      # run from zero to 9
phat_c2[ii] <- sum(coef_2 == i1)/nc
phat_s2[ii] <- sum(se_2 == i1)/ns
phat_t2[ii] <- sum(tstat_2 == i1)/nt
}


#Chi-square Tests:
##################

chisq_c<- nc*sum((phat_c2-p2)^2/p2) 
chisq_s<- ns*sum((phat_s2-p2)^2/p2) 
chisq_t<- nt*sum((phat_t2-p2)^2/p2) 

chisq_c
pchisq(chisq_c,9, ncp = 0, lower.tail = FALSE, log.p = FALSE)
chisq_s
pchisq(chisq_s,9, ncp = 0, lower.tail = FALSE, log.p = FALSE)
chisq_t
pchisq(chisq_t,9, ncp = 0, lower.tail = FALSE, log.p = FALSE)

# U2d Tests:
###########

T<- cumsum(p2)
S_c<- cumsum(phat_c2)
S_s<- cumsum(phat_s2)
S_t<- cumsum(phat_t2)
Z_c<- S_c-T
Z_s<- S_s-T
Z_t<- S_t-T
t<- c()
for (ii in 1:9) {
t[ii]<- (p2[ii]+p2[ii+1])/2
}
t[10]<- (p2[9]+p2[1])/2
Zbar_c<- sum(t*Z_c)
Zbar_s<- sum(t*Z_s)
Zbar_t<- sum(t*Z_t)

U2d_c<- nc*sum((Z_c-Zbar_c)^2*t)
U2d_s<- ns*sum((Z_s-Zbar_s)^2*t)
U2d_t<- nt*sum((Z_t-Zbar_t)^2*t)

U2d_c
U2d_s
U2d_t

# Asy. crit. vals.:        10%     5%       2.5%    1%
#                         0.163   0.205    0.247   0.304

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
UB90_s2<- (  B90+2*f_s + (B90 * (B90+4*f_s*(nc-f_s)/ns  ))^0.5)/(2*(ns+B90))

LB90_s2
UB90_s2
p2              

LB99_s2<- (  B99+2*f_s - (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))
UB99_s2<- (  B99+2*f_s + (B99 * (B99+4*f_s*(ns-f_s)/ns  ))^0.5)/(2*(ns+B99))

LB99_s2
UB99_s2
p2
             
f_t<- nt*phat_t2
LB95_t2<- (  B95+2*f_t - (B95 * (B95+4*f_t*(ns-f_t)/nt  ))^0.5)/(2*(nt+B95))
UB95_t2<- (  B95+2*f_t + (B95 * (B95+4*f_t*(ns-f_t)/nt  ))^0.5)/(2*(nt+B95))

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
ggtitle("Figure 6(a): 95% Goodman Confidence Intervals
for Coefficient 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_c1, ymax=UB95_c1),width=.2, col="red")+
scale_y_continuous(limits = c(0, 0.35))+
  scale_x_continuous(breaks=seq(1,9,1))+
scale_alpha_manual(name = NULL,
   values = c(1,1),
   breaks = c("Benford"),
  guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Full sample
(n = 30,159)")

###                                                
graph_c2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure 6(b): 95% Goodman Confidence Intervals
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Full sample
(n = 27,379)")

ggarrange(graph_c1, graph_c2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# Digits for std. errors
# ----------------------

graph_s1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure 7(a): 95% Goodman Confidence Intervals
for Std. Error 1st. Digits") +
  xlab("1st. Digit") + ylab("Relative Frequency")+
theme(plot.title.position = 'plot', 
  plot.title = element_text(hjust = 0.5))+
  geom_point(aes(alpha="Benford"),size=2)+
  geom_errorbar(aes(ymin=LB95_s1, ymax=UB95_s1), width=.2,col="blue")+
scale_y_continuous(limits = c(0, 0.35))+
scale_x_continuous(breaks=seq(1,9,1))+
 scale_alpha_manual(name = NULL,
   values = c(1, 1),
   breaks = c("Benford"),
   guide = guide_legend(override.aes = list(color = "black", size=2) ) )+
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Full sample
(n = 14,519)")

###     
graph_s2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure 7(b): 95% Goodman Confidence Intervals 
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Full sample
(n = 12,885)")
ggarrange(graph_s1, graph_s2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

# Digits for t-stats.
# ------------------

graph_t1<- ggplot(df1, aes(x=d1, y=p))+
ggtitle("Figure 8(a): 95% Goodman Confidence Intervals 
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Full sample
(n = 10,841)")

###             
graph_t2<- ggplot(df2, aes(x=d2, y=p2))+
ggtitle("Figure 8(b): 95% Goodman Confidence Intervals 
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
annotate("text", x = I(0.85), y = I(0.75), size=3, label = "Full sample
(n = 10,605)")

ggarrange(graph_t1, graph_t2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

##############################################
