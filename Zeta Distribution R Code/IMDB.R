require(igraph)
require(data.table)
require(dplyr)
library(VGAM)            # For computing Riemann zeta function & its first 2 derivatives
library (tolerance)      # For computing the MLE of the zeta distribution parameter using the 'zm.ll' function
set.seed(1234)

IMDB<- read.table(file = 'C:/Users/David/Dropbox/Power Laws/Networks/out_actor2.txt',header=TRUE)

samp_size<-50
IMDB_samp<- head(IMDB, samp_size)
el<- as.matrix(IMDB_samp)
#el<- as.matrix(IMDB)   # Use this for the full IMDB database

#IMDB_samp<- IMDB[sample(1:nrow(IMDB),size=samp_size,replace=TRUE),]
#from<- as.vector(IMDB_samp$from)
#to<- as.vector(IMDB_samp$to)
#EL<- cbind(from,to)

g<- graph_from_edgelist(el, directed = FALSE)
degree<- degree(g)
summary(degree)
N<- length(degree)
plot(g, main="Network of 50 Actors")        

# Gamma constants -  from Choudhury (1995), Table 5 (Note the shift in decimal places, as indicated in that table.)
gamma1<-  -0.072815845483676724861
gamma2<-  -0.009690036031902870231084845
gamma3<-   0.0020538344203033458662
gamma4<-   0.0023253700654673000575
gamma5<-   0.00079332381730106270175
gamma6<-  -0.00023876934543019960987
gamma7<-  -0.00052728956705775104607
gamma8<-  -0.00035212335380303950960
gamma9<-  -0.000034394774418088048178
gamma10<-  0.00020533281490906479468
gamma11<-  0.00027018443954390352667
gamma12<-  0.00016727291210514019335
gamma13<- -0.000027463806603760158860
gamma14<- -0.00020920926205929994584
gamma15<- -0.00028346865532024144664
gamma16<- -0.00019969685830896977471
gamma17<-  0.000026277037109918336699
gamma18<-  0.00030736840814925282659
gamma19<-  0.00050360545304735562606
gamma20<-  0.00046634356151155944940

predtilde<- c()
predhat<- c()
predcup<- c()
zeta.data<- degree

sumlogdata<- sum(log(zeta.data))     # We need this for Firth's estimator

# Obtain the MLE of "s":

out.zeta<- zm.ll(zeta.data, N = Inf, dist = "Zeta")
stilde<- stats4::coef(out.zeta)[[1]]           # Keep value only, not name ("s") as well. This speeds up the subsequent summing.
se_tilde<- sqrt(stats4::vcov(out.zeta)[[1]])   # Asymptotic std. error of the MLE

zetas<- zeta(stilde,0,1)  # using the VGAM package in this & the next 2 lines
zeta1<- zeta(stilde,1,1)
zeta2<- zeta(stilde,2,1)

# Use Choudhury's equation 20 to get the third derivative:
zeta3<- -(6/(stilde-1)^4+gamma3-gamma4*(stilde-1)+(gamma5*(stilde-1)^2)/2-(gamma6*(stilde-1)^3)/6+(gamma7*(stilde-1)^4)/24-(gamma8*(stilde-1)^5)/120
+(gamma9*(stilde-1)^6)/720-(gamma10*(stilde-1)^7)/5040+(gamma11*(stilde-1)^8)/40320-(gamma12*(stilde-1)^9)/362880+(gamma13*(stilde-1)^10)/3628800
- (gamma14*(stilde-1)^11)/39916800 + (gamma15*(stilde-1)^12)/479001600 - (gamma16*(stilde-1)^13)/6227020800 + (gamma17*(stilde-1)^14)/87178291200 
- (gamma18*(stilde-1)^15)/1307674368000 )
bias<- (zetas/(2*N))*(3*zetas*zeta1*zeta2 - 2*zeta1^3 - zetas^2*zeta3)/((zeta1^2 - zetas*zeta2)^2)   # This is just the second-order bias!
shat<- stilde-bias                     # Compute the Cox-Snell bias-corrected estimate

#Taking account of the fact that we have a canonical linear exponential model -

f<- function (x) -sumlogdata -(N+1)*zeta(x,1,1)/zeta(x,0,1) +
(zeta(x,0,1)*
( -(6/(x-1)^4+gamma3-gamma4*(x-1)+(gamma5*(x-1)^2)/2-(gamma6*(x-1)^3)/6+(gamma7*(x-1)^4)/24-(gamma8*(x-1)^5)/120
+(gamma9*(x-1)^6)/720-(gamma10*(x-1)^7)/5040+(gamma11*(x-1)^8)/40320-(gamma12*(x-1)^9)/362880+(gamma13*(x-1)^10)/3628800
- (gamma14*(x-1)^11)/39916800 + (gamma15*(x-1)^12)/479001600 - (gamma16*(x-1)^13)/6227020800 + (gamma17*(x-1)^14)/87178291200 
- (gamma18*(x-1)^15)/1307674368000 )) 

-zeta(x,1,1)*zeta(x,2,1)) /(2*(zeta(x,0,1)*zeta(x,2,1) - zeta(x,1,1)^2))

scup<- uniroot(f , c(1.1,5), tol = 1e-6)$root
c(stilde, shat, scup, se_tilde)

summary(degree)
LIM<- 100
x<- c(1:LIM)
pred<- dzeta(x, scup)
hist(degree, prob=TRUE,breaks=240, xlim=c(0,100), xlab="Degree",main= "",sub="Figure 1: Distribution of degrees for full IMDB network")
lines(x, pred,col="red")
legend(40,0.1,legend = c(expression("Fitted zeta: Firth estimate")), 
                         box.col="white", col="red",lty=1,ncol=1)
 
#distances(g)
#d<- distances(g)
#d[1,]<- 0
#d[,1]<- 0
#max(d)
#mean(d)
#median(d)
###############################################