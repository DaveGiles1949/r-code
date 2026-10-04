library(bootstrap)
library(dgof)
library(stats)
nboot<- 1000
set.seed(12345)

# Need to generate random T-L values even if estimator of v is greater than 1 (in which case "rtopple" in VGAM cannot be used)
# n is sample size; v is shape parameter
rtopp <- function(n, v) {
  if (v <= 0) stop("v must be greater than 0")
  u <- runif(n)
  x <- 1 - sqrt(1 - u^(1 / v))
  return(x)
}

# Functions for Kurtosis calculations (Ref., Nadarajah & Kotz, 2003):
mu2<- function(x) { 1/(1+x)-16^x*(gamma(1+x))^4/(gamma(2+2*x))^2 }
mu4<- function(x) { 2/((1+x)*(2+x))+3*2^(3+4*x)*(1+x)*(3+2*x)*(gamma(1+x))^4/(gamma(4+2*x))^2-3*16^(1+2*x)*(1+x)^4*(gamma(1+x))^8/(gamma(3+2*x))^4  }

#Application 1:
###############

df1<- read.table("C:/Users/OEM/Sync/Bias correction/Topp-Leone/2026/Data/SC16.txt", header=TRUE)
x<- df1$X
x
n<- length(x)
xbar<- mean(x)
# OBTAIN THE MLE AND MOM ESTIMATOR OF "NU"
ml<- -n/(sum(log(x))+sum(log(2-x)))
mom<- uniroot(function(z) gamma(2+2*z)*(xbar-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root
cs<- (n-1)*ml/n
ml
cs
mom
n
# Start of DELETE-1 Jackknife
#========================
#Note - we must NOT use "n" in the next line, as sample size changes during jackknifing. Use "length(x)" instead.

vhat_jack<- function(x) { -length(x)/(sum(log(x))+sum(log(2-x)))}
jack_vhat<- n*ml-(n-1)*mean(jackknife(x,vhat_jack)$jack.values)

vtild_obj<- function(x) {uniroot(function(z) gamma(2+2*z)*(mean(x)-1)+4^z*(gamma(1+z))^2, lower=0,upper=70, tol=0.0001)$root}
temp<- jackknife(x,vtild_obj)
jack_vtild<-  n*mom-(n-1)*mean(temp$jack.values)
jack_vhat
jack_vtild

# End of Jackknife
#========================

# BOOTSTRAP THE STD. ERRORS FOR MLE, CS-MLE, MOM & Jackknife-corrected estimators
# -------------------------------------------------------------------------------
mlboot<- c()
momboot<- c()
csboot<- c()
jack_vhatboot<- c()
jack_vtildboot<- c()

for(ii in 1:nboot) {

xx<- sample(x,n,replace=TRUE)
xxbar<- mean(xx)
vhatboot<- -n/(sum(log(xx))+sum(log(2-xx)))
vhatcsboot<- (n-1)*vhatboot/n
vtildboot<- uniroot(function(z) gamma(2+2*z)*(mean(xx)-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root
vhat_jack<- function(xx) { -length(xx)/(sum(log(xx))+sum(log(2-xx)))}
jack_vhats<- n*vhatboot-(n-1)*mean(jackknife(xx,vhat_jack)$jack.values)

vtild_obj<- function(xx) {uniroot(function(z) gamma(2+2*z)*(xxbar-1)+4^z*(gamma(1+z))^2, lower=0,upper=70, tol=0.0001)$root}
temp<- jackknife(xx,vtild_obj)
jack_vtilds<-  n*vtildboot-(n-1)*mean(temp$jack.values)

mlboot<- c(mlboot,vhatboot)
momboot<- c(momboot,vtildboot)
csboot<- c(csboot, vhatcsboot)
jack_vhatboot<- c(jack_vhatboot,jack_vhats)
jack_vtildboot<- c(jack_vtildboot,jack_vtilds)

}
# End of Bootstrap Loop

seml<- sd(mlboot)
semom<- sd(momboot)
secs<- sd(csboot)
sejhat<- sd(jack_vhatboot)
sejtild<- sd(jack_vtildboot)
c(seml,secs, semom, sejhat, sejtild) # Bootstrapped std. errors of MLE's, MOM & jack estimators

sum(mlboot < 1)/nboot    # Comparing shape parameter estimate with "1"
sum(csboot < 1)/nboot
sum(jack_vhatboot < 1)/nboot
sum(momboot < 1)/nboot
sum(jack_vtildboot > 1)/nboot

# CALCULATE THE BOOTSTRAP BIASES

mlbootbias<- mean(mlboot)-ml
mombootbias<- mean(momboot)-mom

# COMPUTE THE BOOTSTRAP-BIAS-ADJUSTED ESTIMATORS
mlba<- ml-mlbootbias
momba<- mom-mombootbias

c(mlba,momba)     # The bootstrap bias-adjusted MLE and MOM point estimates

# AND THEIR BOOTSTRAPPED STD. ERRORS:

# Set up the "sample" of bootstrap-bias-adjusted estimates
mlbootad<- mlboot-mlbootbias
mombootad<- momboot-mombootbias

semlboot<- sd(mlbootad)
semomboot<- sd(mombootad)
c(semlboot, semomboot)
sum(mlbootad < 1)/nboot    # Comparing shape parameter estimate with "1"
sum(mombootad < 1)/nboot

# Plot the data and fitted models:
denml <-  function(z){ 2*ml*(1-z)*(z*(2-z))^(ml-1)
}
dencs <-  function(z){ 2*cs*(1-z)*(z*(2-z))^(cs-1)
}
denmom <-  function(z){ 2*mom*(1-z)*(z*(2-z))^(mom-1)
}
denmlba <-  function(z){ 2*mlba*(1-z)*(z*(2-z))^(mlba-1)
}
denmomba <-  function(z){ 2*momba*(1-z)*(z*(2-z))^(momba-1)
}

hist(x, prob=TRUE, main="Fig. 1: SC16 Data", xlab="Capacity Factor", breaks=12)  # Altun & Hamedani
curve(denml(x), add=TRUE, col="red",lwd=2)
curve(dencs(x), add=TRUE, col="blue",lwd=2)
curve(denmom(x), add=TRUE, col="green",lwd=2)
curve(denmlba(x), add=TRUE, col="orange",lwd=2)
curve(denmomba(x), add=TRUE, col="purple",lwd=2)

legend(0.6, 4, legend=c("Actual (Hist)","ML", "C-S", "MOM","ML-Boot", "MOM-Boot"),
       col=c("black","red", "blue", "green", "orange", "purple"), lty=c(1,1,1,1,1,1), cex=0.8)


#========================

# Goodness-of-fit tests:
# ---------------------

#Use Dn test - it has the best power

xml<- rtopp(n, ml)
xcs<- rtopp(n,cs)
xmom<- rtopp(n,mom)
xmlba<- rtopp(n, mlba)
xmomba<- rtopp(n,momba)
xjackvhat<- rtopp(n, jack_vhat)
xjackvtild<- rtopp(n,jack_vtild)
# change xml to xcs, etc. in next line for various tests
Dn<- ks.test(x,xml)
Dn
# See Table 1 of Al-Zahrani (2012)for appropriate small-sample p-value ranges for the T-L distribution case

# Kurtosis calculations:
# ---------------------
est<- c(ml,cs,mlba,jack_vhat,mom,momba,jack_vtild)
kurtosis<- mu4(est)/(mu2(est)^2)
kurtosis

# QQ plot
# -------
p <- (1:n - 0.5) / n
sorted_data <- sort(x)
v<- cs
# Calculate theoretical quantiles using the formula
theoretical_quantiles <- 1 - sqrt(1 - p^(1/v))

# Create the Q-Q plot
plot(theoretical_quantiles, sorted_data,
     main = "Q-Q Plot for Topp-Leone Distribution",
     xlab = "Theoretical Quantiles",
     ylab = "Sample Quantiles",
     pch = 19, col = "blue")

# Add a reference line
abline(0, 1, col = "red", lwd = 2)

####################################################################################

# Application 2:
# ============

df2<- read.table("C:/Users/OEM/Sync/Bias correction/Topp-Leone/2026/Data/electronics.txt", header=TRUE)
x<- df2$X/1000
x
n<- length(x)
xbar<- mean(x)
# OBTAIN THE MLE AND MOM ESTIMATOR OF "NU"
ml<- -n/(sum(log(x))+sum(log(2-x)))
mom<- uniroot(function(z) gamma(2+2*z)*(xbar-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root
cs<- (n-1)*ml/n
ml
cs
mom
n
# Start of DELETE-1 Jackknife
#========================
#Note - we must NOT use "n" in the next line, as sample size changes during jackknifing. Use "length(x)" instead.

vhat_jack<- function(x) { -length(x)/(sum(log(x))+sum(log(2-x)))}
jack_vhat<- n*ml-(n-1)*mean(jackknife(x,vhat_jack)$jack.values)

vtild_obj<- function(x) {uniroot(function(z) gamma(2+2*z)*(mean(x)-1)+4^z*(gamma(1+z))^2, lower=0,upper=70, tol=0.0001)$root}
temp<- jackknife(x,vtild_obj)
jack_vtild<-  n*mom-(n-1)*mean(temp$jack.values)
jack_vhat
jack_vtild

# End of Jackknife
#========================

# BOOTSTRAP THE STD. ERRORS FOR MLE, CS-MLE, MOM & Jackknife-corrected estimators
# -------------------------------------------------------------------------------
mlboot<- c()
momboot<- c()
csboot<- c()
jack_vhatboot<- c()
jack_vtildboot<- c()

for(ii in 1:nboot) {

xx<- sample(x,n,replace=TRUE)
xxbar<- mean(xx)
vhatboot<- -n/(sum(log(xx))+sum(log(2-xx)))
vhatcsboot<- (n-1)*vhatboot/n
vtildboot<- uniroot(function(z) gamma(2+2*z)*(xxbar-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root
vhat_jack<- function(xx) { -length(xx)/(sum(log(xx))+sum(log(2-xx)))}
jack_vhats<- n*vhatboot-(n-1)*mean(jackknife(xx,vhat_jack)$jack.values)

vtild_obj<- function(xx) {uniroot(function(z) gamma(2+2*z)*(mean(xx)-1)+4^z*(gamma(1+z))^2, lower=0,upper=70, tol=0.0001)$root}
temp<- jackknife(xx,vtild_obj)
jack_vtilds<-  n*vtildboot-(n-1)*mean(temp$jack.values)

mlboot<- c(mlboot,vhatboot)
momboot<- c(momboot,vtildboot)
csboot<- c(csboot, vhatcsboot)
jack_vhatboot<- c(jack_vhatboot,jack_vhats)
jack_vtildboot<- c(jack_vtildboot,jack_vtilds)

}
# End of Bootstrap Loop

seml<- sd(mlboot)
semom<- sd(momboot)
secs<- sd(csboot)
sejhat<- sd(jack_vhatboot)
sejtild<- sd(jack_vtildboot)
c(seml,secs, semom, sejhat, sejtild) # Bootstrapped std. errors of MLE's, MOM & jack estimators
sum(mlboot < 1)/nboot    # Comparing shape parameter estimate with "1"
sum(csboot < 1)/nboot
sum(jack_vhatboot < 1)/nboot
sum(momboot < 1)/nboot
sum(jack_vtildboot < 1)/nboot

# CALCULATE THE BOOTSTRAP BIASES

mlbootbias<- mean(mlboot)-ml
mombootbias<- mean(momboot)-mom

# COMPUTE THE BOOTSTRAP-BIAS-ADJUSTED ESTIMATORS
mlba<- ml-mlbootbias
momba<- mom-mombootbias

c(mlba,momba)     # The bootstrap bias-adjusted MLE and MOM point estimates

# AND THEIR BOOTSTRAPPED STD. ERRORS:

# Set up the "sample" of bootstrap-bias-adjusted estimates
mlbootad<- mlboot-mlbootbias
mombootad<- momboot-mombootbias

semlboot<- sd(mlbootad)
semomboot<- sd(mombootad)
c(semlboot, semomboot)
sum(mlbootad < 1)/nboot    # Comparing shape parameter estimate with "1"
sum(mombootad < 1)/nboot

# Plot the data and fitted models:
denml <-  function(z){ 2*ml*(1-z)*(z*(2-z))^(ml-1)
}
dencs <-  function(z){ 2*cs*(1-z)*(z*(2-z))^(cs-1)
}
denmom <-  function(z){ 2*mom*(1-z)*(z*(2-z))^(mom-1)
}
denmlba <-  function(z){ 2*mlba*(1-z)*(z*(2-z))^(mlba-1)
}
denmomba <-  function(z){ 2*momba*(1-z)*(z*(2-z))^(momba-1)
}

hist(x, prob=TRUE, main="Fig. 1: SC16 Data", xlab="Capacity Factor", breaks=12)  # Altun & Hamedani
curve(denml(x), add=TRUE, col="red",lwd=2)
curve(dencs(x), add=TRUE, col="blue",lwd=2)
curve(denmom(x), add=TRUE, col="green",lwd=2)
curve(denmlba(x), add=TRUE, col="orange",lwd=2)
curve(denmomba(x), add=TRUE, col="purple",lwd=2)

legend(0.6, 4, legend=c("Actual (Hist)","ML", "C-S", "MOM","ML-Boot", "MOM-Boot"),
       col=c("black","red", "blue", "green", "orange", "purple"), lty=c(1,1,1,1,1,1), cex=0.8)


#========================

# Goodness-of-fit tests:
# ---------------------
#Use Dn test - it has the best power

xml<- rtopp(n,ml)
xcs<- rtopp(n,cs)
xmom<- rtopp(n,mom)
xmlba<- rtopp(n, mlba)
xmomba<- rtopp(n,momba)
xjackvhat<- rtopp(n, jack_vhat)
xjackvtild<- rtopp(n,jack_vtild)
# change xml to xcs, etc in next line for various tests
Dn<- ks.test(x,xml)
Dn
# See Table of Al-Zahrani (2012)for appropriate small-sample p-value ranges for the T-L distribution case
# Kurtosis calculations:
# ---------------------
est<- c(ml,cs,mlba,jack_vhat,mom,momba,jack_vtild)
kurtosis<- mu4(est)/(mu2(est)^2)
kurtosis

# QQ plot
# -------
p <- (1:n - 0.5) / n
sorted_data <- sort(x)
v<- cs
# Calculate theoretical quantiles using the formula
theoretical_quantiles <- 1 - sqrt(1 - p^(1/v))

# Create the Q-Q plot
plot(theoretical_quantiles, sorted_data,
     main = "Q-Q Plot for Topp-Leone Distribution",
     xlab = "Theoretical Quantiles",
     ylab = "Sample Quantiles",
     pch = 19, col = "blue")

# Add a reference line
abline(0, 1, col = "red", lwd = 2)
##################################

# Application 3 - Leukemia data (Feigl & Zelen, 1965)
# =========================

df3<- read.table("C:/Users/OEM/Sync/Bias correction/Topp-Leone/2026/Data/Leukemia.txt", header=TRUE)
x<- df3$X
x<- x/(max(x))    # 0 < x < 1
n<- length(x)
xbar<- mean(x)
# OBTAIN THE MLE AND MOM ESTIMATOR OF "NU"
ml<- -n/(sum(log(x))+sum(log(2-x)))
mom<- uniroot(function(z) gamma(2+2*z)*(xbar-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root
cs<- (n-1)*ml/n
ml
cs
mom
n
# Start of DELETE-1 Jackknife
#========================
#Note - we must NOT use "n" in the next line, as sample size changes during jackknifing. Use "length(x)" instead.

vhat_jack<- function(x) { -length(x)/(sum(log(x))+sum(log(2-x)))}
jack_vhat<- n*ml-(n-1)*mean(jackknife(x,vhat_jack)$jack.values)

vtild_obj<- function(x) {uniroot(function(z) gamma(2+2*z)*(mean(x)-1)+4^z*(gamma(1+z))^2, lower=0,upper=70, tol=0.0001)$root}
temp<- jackknife(x,vtild_obj)
jack_vtild<-  n*mom-(n-1)*mean(temp$jack.values)
jack_vhat
jack_vtild

# End of Jackknife
#========================

# BOOTSTRAP THE STD. ERRORS FOR MLE, CS-MLE, MOM & Jackknife-corrected estimators
# -------------------------------------------------------------------------------
mlboot<- c()
momboot<- c()
csboot<- c()
jack_vhatboot<- c()
jack_vtildboot<- c()

for(ii in 1:nboot) {

xx<- sample(x,n,replace=TRUE)
xxbar<- mean(xx)
vhatboot<- -n/(sum(log(xx))+sum(log(2-xx)))
vhatcsboot<- (n-1)*vhatboot/n
vtildboot<- uniroot(function(z) gamma(2+2*z)*(mean(xx)-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root
vhat_jack<- function(xx) { -length(xx)/(sum(log(xx))+sum(log(2-xx)))}
jack_vhats<- n*vhatboot-(n-1)*mean(jackknife(xx,vhat_jack)$jack.values)

vtild_obj<- function(xx) {uniroot(function(z) gamma(2+2*z)*(xxbar-1)+4^z*(gamma(1+z))^2, lower=0,upper=70, tol=0.0001)$root}
temp<- jackknife(xx,vtild_obj)
jack_vtilds<-  n*vtildboot-(n-1)*mean(temp$jack.values)

mlboot<- c(mlboot,vhatboot)
momboot<- c(momboot,vtildboot)
csboot<- c(csboot, vhatcsboot)
jack_vhatboot<- c(jack_vhatboot,jack_vhats)
jack_vtildboot<- c(jack_vtildboot,jack_vtilds)

}
# End of Bootstrap Loop

seml<- sd(mlboot)
semom<- sd(momboot)
secs<- sd(csboot)
sejhat<- sd(jack_vhatboot)
sejtild<- sd(jack_vtildboot)
c(seml,secs, semom, sejhat, sejtild) # Bootstrapped std. errors of MLE's, MOM & jack estimators

sum(mlboot < 1)/nboot    # Comparing shape parameter estimate with "1"
sum(csboot < 1)/nboot
sum(jack_vhatboot < 1)/nboot
sum(momboot < 1)/nboot
sum(jack_vtildboot < 1)/nboot

# CALCULATE THE BOOTSTRAP BIASES

mlbootbias<- mean(mlboot)-ml
mombootbias<- mean(momboot)-mom

# COMPUTE THE BOOTSTRAP-BIAS-ADJUSTED ESTIMATORS
mlba<- ml-mlbootbias
momba<- mom-mombootbias

c(mlba,momba)     # The bootstrap bias-adjusted MLE and MOM point estimates

# AND THEIR BOOTSTRAPPED STD. ERRORS:

# Set up the "sample" of bootstrap-bias-adjusted estimates
mlbootad<- mlboot-mlbootbias
mombootad<- momboot-mombootbias

semlboot<- sd(mlbootad)
semomboot<- sd(mombootad)
c(semlboot, semomboot)
sum(mlbootad < 1)/nboot    # Comparing shape parameter estimate with "1"
sum(mombootad < 1)/nboot

# Plot the data and fitted models:
denml <-  function(z){ 2*ml*(1-z)*(z*(2-z))^(ml-1)
}
dencs <-  function(z){ 2*cs*(1-z)*(z*(2-z))^(cs-1)
}
denmom <-  function(z){ 2*mom*(1-z)*(z*(2-z))^(mom-1)
}
denmlba <-  function(z){ 2*mlba*(1-z)*(z*(2-z))^(mlba-1)
}
denmomba <-  function(z){ 2*momba*(1-z)*(z*(2-z))^(momba-1)
}

hist(x, prob=TRUE, main="Fig. 1: SC16 Data", xlab="Capacity Factor", breaks=12)  # Altun & Hamedani
curve(denml(x), add=TRUE, col="red",lwd=2)
curve(dencs(x), add=TRUE, col="blue",lwd=2)
curve(denmom(x), add=TRUE, col="green",lwd=2)
curve(denmlba(x), add=TRUE, col="orange",lwd=2)
curve(denmomba(x), add=TRUE, col="purple",lwd=2)

legend(0.6, 4, legend=c("Actual (Hist)","ML", "C-S", "MOM","ML-Boot", "MOM-Boot"),
       col=c("black","red", "blue", "green", "orange", "purple"), lty=c(1,1,1,1,1,1), cex=0.8)


#========================

# Goodness-of-fit tests:
# ---------------------

#Use Dn test - it has the best power

xml<- rtopp(n, ml)
xcs<- rtopp(n,cs)
xmom<- rtopp(n,mom)
xmlba<- rtopp(n, mlba)
xmomba<- rtopp(n,momba)
xjackvhat<- rtopp(n, jack_vhat)
xjackvtild<- rtopp(n,jack_vtild)
# change xml to xcs, etc. in next line for various tests
Dn<- ks.test(x,xml)
Dn
# See Table of Al-Zahrani (2012)for appropriate small-sample p-value ranges for the T-L distribution case
# Kurtosis calculations:
# ---------------------
est<- c(ml,cs,mlba,jack_vhat,mom,momba,jack_vtild)
kurtosis<- mu4(est)/(mu2(est)^2)
kurtosis

# QQ plot
# -------
p <- (1:n - 0.5) / n
sorted_data <- sort(x)
v<- cs
# Calculate theoretical quantiles using the formula
theoretical_quantiles <- 1 - sqrt(1 - p^(1/v))

# Create the Q-Q plot
plot(theoretical_quantiles, sorted_data,
     main = "Q-Q Plot for Topp-Leone Distribution",
     xlab = "Theoretical Quantiles",
     ylab = "Sample Quantiles",
     pch = 19, col = "blue")

# Add a reference line
abline(0, 1, col = "red", lwd = 2)

####################################################################################

