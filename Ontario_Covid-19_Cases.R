# CODE FOR MODELLING THE GROWTH CURVE OF TOTAL CONFIRMED COVID-19 CASES IN ONTARIO

# (AVAILABLE AT:  https://raw.githubusercontent.com/DaveGiles1949/r-code/master/Ontario_Covid-19_Cases.R )
# ---------------------------------------------------------------------------------------------------------

# DATA SOURCE:    https://www.covid-19canada.com/

# AUTHOR:         DAVID GILES (davegiles1949@gmail.com)

# LAST UPDATED:   25 April, 2020
# ---------------------------------------------------------------------------------------------------------

library(growthcurver)  

today <- Sys.Date()

inflection<- c()
doub_time<- c()
est_doub_time<- c()
gof<- c()
pred<- c()

#cases<- read.csv("https://raw.githubusercontent.com/DaveGiles1949/r-code/master/Ontario_Covid-19_Cases.txt", header=TRUE)

file_name <- "C:/Users/David Giles/Desktop/Virus/Ontario_Covid-19_Cases.txt"
cases <- read.table(file_name, header = TRUE)
tot<- cases$TOTAL_CASES
n_min<- 26             # Smallest number of days to include in the sequential anlaysis
n_max<- length(tot)    # Largest number of days to include in the sequential analysis
                       # NOTE: "Day 1" for this analysis is 2 March, 2020
# Run the model over various samples of data, beginning on 2 March

for (ii in n_min:n_max) {
t<- seq(1:ii)
total<- tot[1:ii]
# Now use "Growthcurver" to fit Logistic Growth Model (Ref: https://cran.r-project.org/web/packages/growthcurver/growthcurver.pdf)
gc_fit <- SummarizeGrowth(t, total)
inflection<- c(inflection,gc_fit$vals$t_mid)
doub_time<- c(doub_time,log(2)/log(tot[ii]/tot[ii-1]))
est_doub_time<- c(est_doub_time,gc_fit$vals$t_gen )
gof<- c(gof,gc_fit$vals$auc_l/gc_fit$vals$auc_e)

}

# Get the 1 week-ahead prediction path from the end of the longest available sample, for plotting later:
n_pred<- n_max+7

for (kk in 1:n_pred) {
        f <- gc_fit$vals$n0*gc_fit$vals$k / (gc_fit$vals$n0 +(gc_fit$vals$k - gc_fit$vals$n0)*exp(-kk*gc_fit$vals$r))
        pred<- round(c(pred, f),0)  
}

# "DT" is the (max.)doubling time
# auc_l The area under the curve of the fitted logistic equation from time 0 to time t
# auc_e The area under the curve of the measurements.
# t_mid is the time of the point of inflection (half way to max. capacity)
gc_fit$vals
gc_fit$vals$note            # Are there any warning notes?

# Simple Exponential Growth Model using latest sample
egm<- lm(log(total) ~ t)
summary(egm)
egm_pred<- exp(predict(egm))               # Convert prediction of "log(total)" to one for "total" itself
Duan<- mean(exp(log(total)-predict(egm)))  
egm_fit<- egm_pred*Duan
                                           # Reduce bias of the "raw" prediction by using Duan's non-parametric smoother
                                           # (See https://davegiles.blogspot.com/2014/12/s.html)

# Plot the results based on the latest data:
# In R, the calendar starts at 1970-01-01

poi1<- round(18322+gc_fit$vals$t_mid,0)    # Day 18322 is 2020-03-01. This is Day zero
class(poi1)<- "Date"   
poi2<- round(18322+n_max,0)    
class(poi2)<- "Date" 
par(mfrow=c(1,1))
plot(gc_fit,
     main="Covid-19 Confirmed Cases in Ontario", 
          ylab="Cases", xlab="Day", 
     type="l")
lines(t,egm_fit, col="blue")
legend("topleft",inset=0.025,
       c("Confirmed Cases","Naive Exponential Prediction", "Logistic Growth Prediction"),
       col=c("black","blue","red"), lty=c(5,1,1), pch=c(16,46,46), box.lty=0)
text(15,4000, cex=0.9, col="blue", paste0("Sample up to end of Day ", poi2))
text(15,3500, cex=0.9, col="red",paste0("Logistic doubling time = ", round(gc_fit$vals$t_gen,1), " days") )
text(15,3000, cex=0.9,col="red", paste0("Median date (point of inflection) = ", poi1))
text(15,2500, cex=0.9,col="red", paste0("Area under logistic / Area under observed = ", round(gc_fit$vals$auc_l/gc_fit$vals$auc_e,4)))
text(5,500, cex=0.8,"(t = 0 is 1 March, 2020)")
text(40,400, cex=0.8,font=3, paste0("Produced on ", today))

# Plot the predicted time-path for the confirmed cases, up to 1 week ahead:
# In R, the calendar starts at 1970-01-01

poi3<- round(18322+n_pred,0)    # Day 18322 is 2020-03-01
class(poi3)<- "Date"   
par(mfrow=c(1,1))
plot(pred, col="red",
     main="Projected Ontario Cases, Up to 1 Week Ahead", 
     ylab="Cases", xlab="Day", 
     type="l")
lines(t,total,type="b")
legend("topleft",inset=0.025,
       c("Confirmed Cases","Logistic Growth Prediction"),
       col=c("black","red"), lty=c(NA,1), pch=c(1,46), box.lty=0)
text(10, 4000, cex=0.9, col="blue", paste0("Sample up to end of  ", poi2))
text(n_max-2, pred[n_pred], cex=0.8, col="red", paste0(poi3, " = ", pred[n_pred], " cases"))
text(10,700, cex=0.8,"(t = 0 is 1 March, 2020)")
text(45,500, cex=0.8,font=3, paste0("Produced on ", today ))

# Now plot a summary of the results using the successive sample periods:

Obs<- seq(n_min:n_max)+n_min-1
par(mfrow=c(2,2))
plot(Obs,inflection, main="Estimated Inflection Date",
        ylab= "Inflection (Day)", xlab="Sample Size (Days)", type="b", col="red")
plot(Obs,gof, main="Area Under Logistic / Area Under Actual",
     ylab= "Ratio", xlab="Sample Size (Days)", type="b", col="red")
abline(h=1, col="purple")
plot(Obs,est_doub_time, main="Logistic Doubling Time",
     ylab= "Doubling Time (Days)", xlab="Sample Size (Days)", col="red",type="b", ylim = rev(range(est_doub_time)))
text(35,4.5, cex=0.8,col="blue", "Note: Reversed y-axis")
plot(Obs,doub_time, main="Actual Doubling Time",
     ylab= "Doubling Time (Days)", xlab="Sample Size (Days)", type="b", col="red", ylim = rev(range(doub_time)))
text(35,12, cex=0.8,col="blue", "Note: Reversed y-axis")

# END OF FILE  ###################################################################
