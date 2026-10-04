set.seed(12345)
library(VGAM)           # Need this for the rtopple random number generator
library(bootstrap)      # Need this for the Jackknife( ) command
n<- 10
nrep<- 50000 #    THE RESULTS STABILIZE WITH THIS MANY MONTE CARLO REPLICATIONS
sumvhata<- sumvhata2<- 0
sumvtilda<- sumvtilda2<- 0

vvv<- 0.1  # shape parameter "nu" 

# START THE SIMULATION LOOP
# ========================
for(jj in 1:nrep) {

mlsumjack<- momsumjack<- 0

# GENERATE THE TOPP-LEONE VARIATES 

x<- rtopple(n, shape=vvv)      
xbar<- mean(x)

# OBTAIN THE ML AND MOM ESTIMATORS OF "NU" USING THE FULL SAMPLE
# THE MAXIMUM LIKELIHOOD ESTIMATOR OF NU HAS A CLOSED-FORM SOLUTION. THE MOM ESTIMATOR DOES NOT.

vhat<- -n/(sum(log(x))+sum(log(2-x)))   # MLE
vtild<- uniroot(function(z) gamma(2+2*z)*(xbar-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root          # MOM

# Start of DELETE-1 Jackknife
#========================
#Note - we must NOT use "n" in the next line, as sample size changes during jackknifing. Use "length(x)" instead.

vhat_jack<- function(x) { -length(x)/(sum(log(x))+sum(log(2-x)))}
jack_vhat<- n*vhat-(n-1)*mean(jackknife(x,vhat_jack)$jack.values)

vtild_obj<- function(x) {uniroot(function(z) gamma(2+2*z)*(mean(x)-1)+4^z*(gamma(1+z))^2, lower=0,upper=70, tol=0.0001)$root}
temp<- jackknife(x,vtild_obj)
jack_vtild<-  n*vtild-(n-1)*mean(temp$jack.values)

# End of Jackknife
#========================

sumvhata<- sumvhata+jack_vhat
sumvhata2<- sumvhata2+jack_vhat^2
sumvtilda<- sumvtilda+jack_vtild
sumvtilda2<- sumvtilda2+jack_vtild^2
}

# END OF THE MONTE CARLO SIMULATION LOOP
# ==========================

# CALCULATE THE BIAS AND THE % BIAS OF THE JACKKNIFE-BIAS-ADJUSTED ESTIMATORS
biasvhata<- (sumvhata/nrep)-vvv
biasvtilda<- (sumvtilda/nrep)-vvv
pbiasvhata<- 100*biasvhata/vvv
pbiasvtilda<- 100*biasvtilda/vvv

# CALCULATE THE MSE AND THE % MSE OF THE JACKKNIFE-BIAS-ADJUSTED ESTIMATORS
msevhata<- biasvhata^2+((sumvhata2/nrep)-(sumvhata/nrep)^2)
pmsevhata<- 100*msevhata/vvv^2
msevtilda<- biasvtilda^2+ ((sumvtilda2/nrep)-(sumvtilda/nrep)^2)
pmsevtilda<- 100*msevtilda/vvv^2
#
#
cat("Nu =", vvv, ", Sample Size =", n, ", Repetitions =", nrep, "\n")
#
cat("Simulated % Bias of Jackknife-Adjusted MLE =", pbiasvhata, "\n")
#
cat("Simulated % MSE of Jackknife-Adjusted MLE =", pmsevhata, "\n")
#
cat("Simulated % Bias of Jackknife-Adjusted MOM =", pbiasvtilda, "\n")
#
cat("Simulated % MSE of Jackknife-Adjusted MOM =", pmsevtilda, "\n")
