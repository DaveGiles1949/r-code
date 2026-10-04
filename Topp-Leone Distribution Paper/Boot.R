set.seed(12345)
library(VGAM)      # Need this for the rtopple random number generator
n<- 500
nrep<- 50000 #    THE RESULTS SEEM TO STABILIZE WITH THIS MANY M.C. REPLICATIONS
nboot<- 1000
sumvhata<- sumvhata2<- 0
sumvtilda<- sumvtilda2<- 0

vvv<- 0.1  # shape parameter "nu" 

# NOTE THAT WE ARE SETTING B=1 THROUGHOUT, SO THERE IS JUST ONE PARAMETER TO ESTIMATE
# AND THIS ESTIMATOR OF NU HAS A CLOSED-FORM SOLUTION

# START THE SIMULATION LOOP
# ========================
for(jj in 1:nrep) {

mlsumboot<- momsumboot<- 0

# GENERATE THE TOPP-LEONE VARIATES 

x<- rtopple(n, shape=vvv)
xbar<- mean(x)
# OBTAIN THE ML AND MOM ESTIMATORS OF "NU"

vhat<- -n/(sum(log(x))+sum(log(2-x)))      # MLE
vtild<- uniroot(function(z) gamma(2+2*z)*(xbar-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root             # MOM

# Start of Bootstrap Loop
#========================

for(ii in 1:nboot) {

xx<- sample(x,n,replace=TRUE)
xxbar<- mean(xx)
vhatboot<- -n/(sum(log(xx))+sum(log(2-xx)))
vtildboot<- uniroot(function(z) gamma(2+2*z)*(xxbar-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root
mlsumboot<- mlsumboot+vhatboot
momsumboot<- momsumboot+vtildboot

}

# End of Bootstrap Loop
#========================

# CALCULATE THE BOOTSTRAP BIASES
mlbootbias<- (mlsumboot/nboot)-vhat
mombootbias<- (momsumboot/nboot)-vtild

# COMPUTE THE BIAS-ADJUSTED ESTIMATORS
vmla<- vhat-mlbootbias
sumvhata<- sumvhata+vmla
sumvhata2<- sumvhata2+vmla^2
vmoma<- vtild-mombootbias
sumvtilda<- sumvtilda+vmoma
sumvtilda2<- sumvtilda2+vmoma^2

}

# END OF THE MONTE CARLO SIMULATION LOOP
# ==========================

# CALCULATE THE BIAS AND % BIAS OF THE BOOTSTRAP-BIAS-ADJUSTED ESTIMATORS
biasvhata<- (sumvhata/nrep)-vvv
biasvtilda<- (sumvtilda/nrep)-vvv
pbiasvhata<- 100*biasvhata/vvv
pbiasvtilda<- 100*biasvtilda/vvv

# CALCULATE THE MSE AND THE % MSE OF THE BIAS-ADJUSTED ESTIMATOR
msevhata<- biasvhata^2+((sumvhata2/nrep)-(sumvhata/nrep)^2)
pmsevhata<- 100*msevhata/vvv^2
msevtilda<- biasvtilda^2+ ((sumvtilda2/nrep)-(sumvtilda/nrep)^2)
pmsevtilda<- 100*msevtilda/vvv^2
#
#
cat("Nu =", vvv, ", Sample Size =", n, ", Repetitions =", nrep, ", Bootstraps =", nboot, "\n")
#
cat("Simulated % Bias of Bootstrap-Adjusted MLE =", pbiasvhata, "\n")
#
cat("Simulated % MSE of Bootstrap-Adjusted MLE =", pmsevhata, "\n")
#
cat("Simulated % Bias of Bootstrap-Adjusted MOM =", pbiasvtilda, "\n")
#
cat("Simulated % MSE of Bootstrap-Adjusted MOM =", pmsevtilda, "\n")
