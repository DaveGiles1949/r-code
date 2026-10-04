set.seed(123)
library(VGAM)
n<-10
nrep<- 50000 #    THE RESULTS SEEM TO STABILIZE WITH 25K M.C. REPLICATIONS
vvv<- 0.1 # shape parameter "nu" - should be positive

# NOTE THAT WE ARE SETTING B=1 THROUGHOUT, SO THERE IS JUST ONE PARAMETER TO ESTIMATE
# AND THE MLE FOR NU HAS A CLOSED-FORM SOLUTION

vhat<- c()
vmom<- c()
vcs<- c()
lvhat<- length(vhat)
lvcs<- length(vcs)
lvmom<- length(vmom)

# START THE SIMULATION LOOP
# ========================
while(lvhat<nrep || lvcs<nrep || lvmom<nrep) {

# GENERATE THE TOPP-LEONE VARIATES (SEE NADARAJAH & KOTZ, J. APPLIED STATS, 2003, 30, P.317)
#u<- runif(min=0.001, max=1,n)   #  If u = 0 then x = 0, and this isn't allowed
#x<- 1-sqrt(1-u^(1/vvv))
x<- rtopple(n, shape=vvv)
#mean(x)
#1-4^vvv*gamma(1+vvv)^2/gamma(2+2*vvv)
#mean(x^2)
#(2+vvv)/(1+vvv)-2^(1+2*vvv)*(gamma(1+vvv))^2/gamma(2+2*vvv)
#mean(x^3)
#(4+vvv)/(1+vvv)-3*4^(1+vvv)*gamma(1+vvv)*gamma(vvv+3)/gamma(4+2*vvv)
# If we set n = 25,000, then the first 3 moments seem to check out O.K. for all values of "nu"

xbar<- mean(x)
# OBTAIN THE MLE AND MOM ESTIMATOR OF "NU"
ml<- -n/(sum(log(x))+sum(log(2-x)))
#if(ml<1 && ml>0){
vhat<- c(vhat, ml)
#}
lvhat<- length(vhat)
mom<- uniroot(function(z) gamma(2+2*z)*(xbar-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root
#if(mom>0 && mom<1) {
vmom<- c(vmom,mom)
#}
lvmom<- length(vmom)

# NOW SET THINGS UP FOR ANALYTICAL BIAS CALCULATIONS, BASED ON THE MLE OF THE PARAMETER:

# COMPUTE THE BIAS-ADJUSTED ESTIMATOR:
#if(ml<1 && ml>0){
ba<- (n-1)*ml/n

#if(ba>0 && ba<1) {
vcs<- c(vcs,ba)
#}
#}
lvcs<- length(vcs)

}
# END OF THE SIMULATION LOOP
# ==========================
summary(vhat)
length(vhat)
summary(vcs)
length(vcs)
summary(vmom)
length(vmom)

# WE NEED TO USE ONLY THE FIRST "NREP" ELEMENTS FOR CONSTRUCTING BIASES & MSE'S


# CALCULATE THE BIASES AND MSE'S OF THE 3 ESTIMATORS
# BASED ON THE MONTE CARLO SIMULATION

biasvhat<- mean(vhat[1:nrep])- vvv
biasvcs<- mean(vcs[1:nrep])- vvv
biasvmom<- mean(vmom[1:nrep]) - vvv
msevhat<- var(vhat[1:nrep])+biasvhat^2
msevcs<- var(vcs[1:nrep])+biasvcs^2
msevmom<- var(vmom[1:nrep])+biasvmom^2

# CALCULATE THE % BIASES & MSE's:
pbiasvhat<- 100*biasvhat/vvv
pbiasvcs<- 100*biasvcs/vvv
pbiasvmom<- 100*biasvmom/vvv
pmsevhat<- 100*msevhat/vvv^2
pmsevcs<- 100*msevcs/vvv^2
pmsevmom<- 100*msevmom/vvv^2

cat("Nu =", vvv, ", Sample Size =", n, ", Repetitions =", nrep, "\n")
cat("Simulated % Bias of Nu-hat =", pbiasvhat, "Simulated % Bias of Adjusted Nu-hat =", pbiasvcs, "Simulated % Bias of Nu-MOM =", pbiasvmom, "\n")
cat("Simulated % MSE of Nu-hat =", pmsevhat, "Simulated % MSE of Adjusted Nu-hat =", pmsevcs, "Simulated % MSE of Nu-MOM = " , pmsevmom, "\n")