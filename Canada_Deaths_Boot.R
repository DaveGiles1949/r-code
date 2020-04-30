# Last update: 30 April, 2020
# ---------------------------

library(growthcurver)

set.seed(1234)
pred<- c()
pre<- c()
f_star<- c()
forc<- c()
lower<- c()
upper<- c()
today <- Sys.Date()
n_boot<-1000

#cases<- read.csv("https://raw.githubusercontent.com/DaveGiles1949/r-code/master/Canadian_Covid-19_Deaths.txt", header=TRUE)
file_name <- "C:/Users/David Giles/Desktop/Virus/Canadian_Covid-19_Deaths.txt"
cases <- read.table(file_name, header = TRUE)
dead<- cases$DEATHS
deaths<- dead
n_max<- length(dead)
# NOTE: "Day 1" for this analysis is 14 March, 2020

  t<- seq(1:length(dead))
  t7<- seq(1: (length(dead)+7))
  
 
  # Now use "Growthcurver" to fit Logistic Growth Model using latest sample
  
  gc_fit <- SummarizeGrowth(t, deaths)               # Use these estimated coefficients as basis for point predections 
                                                     # & for setting up the bootstrapping
  
  fit <- gc_fit$vals$n0*gc_fit$vals$k / (gc_fit$vals$n0 +(gc_fit$vals$k - gc_fit$vals$n0)*exp(-t*gc_fit$vals$r))
  pred<- round(c(pred, fit),0) 
  
   #  Set up for sequential predictions for up to 7 days beyond latest full sample
  n_pred<- n_max + 7
  start<- n_max + 1
  
  for (jj in start:n_pred) {
  
    # Point predictions, beyond sample:
  f <- gc_fit$vals$n0*gc_fit$vals$k / (gc_fit$vals$n0 +(gc_fit$vals$k - gc_fit$vals$n0)*exp(-jj*gc_fit$vals$r))
  pre<- round(c(pre, f),0) 
  
  # Now start the Bootstrap loop:
 
   
  for (ii in 1:n_boot) {
    resid <- deaths - (gc_fit$vals$n0*gc_fit$vals$k / (gc_fit$vals$n0 +(gc_fit$vals$k - gc_fit$vals$n0)*exp(-t*gc_fit$vals$r)))
    boot_res<- sample(resid, size=n_max, replace=TRUE)
    
    deaths_star<- round(gc_fit$vals$n0*gc_fit$vals$k / (gc_fit$vals$n0 +(gc_fit$vals$k - gc_fit$vals$n0)*exp(-t*gc_fit$vals$r)) 
  + boot_res,0)
    gc_fit_star <- SummarizeGrowth(t, deaths_star)
    
    f_star[ii] <- gc_fit_star$vals$n0*gc_fit_star$vals$k / (gc_fit_star$vals$n0 +(gc_fit_star$vals$k - gc_fit_star$vals$n0)*exp(-jj*gc_fit_star$vals$r))
  
  # Now we need to create the Bootstrap the prediction interval, for the given forecast period

   }

ci<- quantile(f_star, c(0.025, 0.975))
lower[jj]<- round(ci[[1]],0)
upper[jj]<- round(ci[[2]],0)

forc[jj]<- round(mean(f_star),0)

}


# Get rid of the "NAs" and shorten the vector lengths
fore<- subset(forc, forc>0)
low<- subset(lower, lower>0)
up<- subset(upper, upper>0)


fitted<- c(fit,fore)

# In R, the calendar starts at 1970-01-01

poi1<- round(18334+n_max+7,0)    # Day 18334 = 2020-03-13. This is Day zero
class(poi1)<- "Date"
poi2<- round(18334+n_max,0)    
class(poi2)<- "Date" 

par(mfrow=c(1,1))
plot(t7,fitted, type="l",
     main="Predicted Covid-19 Deaths in Canada", 
     ylab="Deaths", xlab="Days", ylim=c(0,max(up)), col="blue")
lines(t7,upper, col="red", lty=2)
lines(t7,lower, col="red", lty=2)
lines(t,deaths, type="p",lty=0, cex=0.5, pch=16)
abline(v=n_max, col="purple")
legend("topleft",inset=0.025,
       c("Actual","Logistic Model Prediction", "lower 95% C.I.", "upper 95% C.I."),
       col=c("black","blue", "red","red"), lty=c(0,1,2,2), pch=c(16,46,46,46), box.lty=0)
text(10, 400, cex=0.9, col="blue", paste0("Sample up to end of day ", poi2))
text(33,25, cex=0.8,font=3, paste0("Produced on ", today))
text(n_pred,deaths[n_max]-400, cex=0.7, col="red", max(low))
text(n_pred,deaths[n_max], cex=0.7, col="red" , max(up))
text(n_pred,deaths[n_max]-200, cex=0.7, col="blue", fitted[n_pred])
text(n_pred-1,deaths[n_max]-600, cex=0.7, poi1)
low
fore
up

# END OF FILE  ###################################################################
