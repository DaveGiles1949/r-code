# CODE FOR MODELLING THE GROWTH CURVE OF TOTAL CONFIRMED COVID-19 CASES IN CANADA

# (AVAILABLE AT:  https://raw.githubusercontent.com/DaveGiles1949/r-code/master/Canadian_Covid-19_Deaths.R )
# ---------------------------------------------------------------------------------------------------------

# DATA SOURCE:    https://www.covid-19canada.com/

# AUTHOR:         DAVID GILES (davegiles1949@gmail.com)

# LAST UPDATED:   11 April, 2020
# ---------------------------------------------------------------------------------------------------------

library(growthcurver)

today <- Sys.Date()

inflection<- c()
doub_time<- c()
pred<- c()
est_doub_time<- c()
gof<- c()

#cases<- read.csv("https://raw.githubusercontent.com/DaveGiles1949/r-code/master/Canadian_Covid-19_Deaths.txt", header=TRUE)
file_name <- "C:/Users/David Giles/Desktop/Virus/Canadian_Covid-19_Deaths.txt"
cases <- read.table(file_name, header = TRUE)
dead<- cases$DEATHS
n_min<- 10             # Smallest number of days to include in the sequential anlaysis
n_max<- length(dead)   # Largest number of days to include in the sequential analysis
                       # NOTE: "Day 1" for this analysis is 14 March, 2020

# Run the model over various samples of data, beginning on 14 March
for (ii in n_min:n_max) {
t<- seq(1:ii)
deaths<- dead[1:ii]

# Now use "Growthcurver" to fit Logistic Growth Model
gc_fit <- SummarizeGrowth(t, deaths)
gc_fit    # "DT" is the (max.)doubling time
          # auc_l The area under the curve of the fitted logistic equation from time 0 to time t
          # auc_e The area under the curve of the measurements.
          # t_mid is the time of the point of inflection (half way to max. capacity)
gc_fit$vals
gc_fit$vals$note            # Are there any warning notes?
inflection<- c(inflection,gc_fit$vals$t_mid)
doub_time<- c(doub_time,log(2)/log(deaths[ii]/deaths[ii-1]))
est_doub_time<- c(est_doub_time,gc_fit$vals$t_gen )
gof<- c(gof,gc_fit$vals$auc_l/gc_fit$vals$auc_e)

}

# Simple Exponential Growth Model using latest sample
egm<- lm(log(deaths) ~ t)
summary(egm)
egm_pred<- exp(predict(egm))               # Convert prediction of "log(deaths)" to one for "deaths" itself
Duan<- mean(exp(log(deaths)-predict(egm)))  
egm_fit<- egm_pred*Duan
                                           # Reduce bias of the "raw" prediction by using Duan's non-parametric smoother
                                           # (See https://davegiles.blogspot.com/2014/12/s.html)

# Use the Logistic model to predict beyond the end of the full sample by 1 week:

n_pred<- n_max + 7

for (jj in 1:n_pred) {
        f <- gc_fit$vals$n0*gc_fit$vals$k / (gc_fit$vals$n0 +(gc_fit$vals$k - gc_fit$vals$n0)*exp(-jj*gc_fit$vals$r))
        pred<- round(c(pred, f),0)                
}

# In R, the calendar starts at 1970-01-01

poi1<- round(gc_fit$vals$t_mid+18334,0)    # Day 18334 = 2020-03-13. This is Day zero
class(poi1)<- "Date"                   
poi2<- round(18334+n_max,0)    
class(poi2)<- "Date" 

# Plot the results based on the latest data:
par(mfrow=c(1,1))
plot(gc_fit,
     main="Covid-19 Deaths in Canada", 
          ylab="Deaths", xlab="Day",
     type="l")
lines(t,egm_fit, col="blue")
legend("topleft",inset=0.025,
       c("Deaths","Naive Exponential Prediction", "Logistic Growth Prediction"),
       col=c("black","blue","red"), lty=c(1,1,1), pch=c(16,46,46), box.lty=0)
text(13,1, cex=0.8,"(t = 0 is 13 March, 2020)")
text(7, 230, cex=0.9, col="blue", paste0("Sample up to end of Day ", poi2))
text(7,190, cex=0.9, col="red",paste0("Logistic doubling time = ", round(gc_fit$vals$t_gen,1), " days") )
text(7,140, cex=0.9,col="red", paste0("Median date (point of inflection) = ", poi1))
text(7,100, cex=0.9,col="red", paste0("Area under logistic / Area under actual = ", round(gc_fit$vals$auc_l/gc_fit$vals$auc_e,4)))
text(24,1, cex=0.8,font=3, paste0("Produced on ", today))

# Plot the predicted time-path for deaths, up to 1 week ahead:
poi3<- round(18334+n_pred,0)    # Day 18334 is 2020-03-13
class(poi3)<- "Date"
par(mfrow=c(1,1))
plot(pred, col="red",
     main="Projected Deaths, Up to 1 Week Ahead", 
     ylab="Deaths", xlab="Day", 
     type="l")
lines(t,deaths,type="o")
legend("topleft",inset=0.025,
       c("Deaths","Logistic Growth Prediction"),
       col=c("black","red"), lty=c(1,1), pch=c(1,46), box.lty=0)
text(10, 400, cex=0.9, col="blue", paste0("Sample up to end of  ", poi2))
text(n_max+1, pred[n_pred], cex=0.8, col="red", paste0(poi3, " = ", pred[n_pred], " deaths"))
text(5,80, cex=0.8,"(t = 0 is 13 March, 2020)")
text(27,1, cex=0.8,font=3, paste0("Produced on ", today))

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
text(20,1.5, cex=0.8,col="blue", "Note: Reversed y-axis")
plot(Obs,doub_time, main="Actual Doubling Time",
     ylab= "Doubling Time (Days)", xlab="Sample Size (Days)", type="b", col="red", ylim = rev(range(doub_time)))
text(21,8.5, cex=0.8,col="blue", "Note: Reversed y-axis")
##########################################################
# Alternative plot

#plot(Obs,doub_time, main="Actual Doubling Time & Inflection Day", col.axis="red",
#     ylab= "Doubling Time (Days)", xlab="Sample Size (Days)", type="l", lty=2,col="red", ylim = rev(range(doub_time)))
#par(new = TRUE)
#plot(Obs, inflection, type = "l",lty=1, axes = FALSE, bty = "n", xlab = "", ylab = "", col="blue")
#axis(side=4, at = pretty(range(inflection)), col.axis="blue")
#mtext("Inflection (Day)", side=4, line=3)

# END OF FILE  ###################################################################
