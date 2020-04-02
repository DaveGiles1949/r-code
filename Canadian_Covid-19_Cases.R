# CODE FOR MODELLING THE GROWTH CURVE OF TOTAL COVID-19 CASES IN CANADA
# ---------------------------------------------------------------------

# DATA SOURCE: https://www.covid-19canada.com/

# AUTHOR: DAVID GILES (davegiles1949@gmail.com)

# LAST UPDATED: 1 April, 2020
# ------------------------------------------------------------------------------------------------

library(growthcurver)

today <- Sys.Date()
n_min<- 26             # Smallest number of days to include in the sequential anlaysis
n_max<- 31             # Largest number of days to include in the sequential analysis
                       # NOTE: "Day 1" is 2 March, 2020

upper<- c()
lower<- c()
saturation<- c()
inflection<- c()
max_doub_time<- c()
exp_growth<- c()

cases<- read.csv("https://raw.githubusercontent.com/DaveGiles1949/r-code/master/Canadian_Covid-19_Cases.txt", header=TRUE)
tot<- cases$TOTAL_CASES

# Run the model over various samples of data, beginning on 2 March
for (ii in n_min:n_max) {
t<- seq(1:ii)
total<- tot[1:ii]
egm<- lm(log(total) ~ t)                           # Simple Exponential Growth Model
exp_growth<- c(exp_growth,coefficients(egm)[[2]])

# Now use "Growthcurver" to fit Logistic Growth Model
gc_fit <- SummarizeGrowth(t, total)
gc_fit    # "DT" is the (max.)doubling time
          # auc_l The area under the curve of the fitted logistic equation from time 0 to time t
          # auc_e The area under the curve of the measurements.
          # t_mid is the time of the point of inflection (half way to max. capacity)
gc_fit$vals
gc_fit$vals$note            # Are there any warning notes?
inflection<- c(inflection,gc_fit$vals$t_mid)
saturation<-c(saturation,gc_fit$vals$k)
max_doub_time<- c(max_doub_time,log(2)/gc_fit$vals$r)
upper<- c(upper,gc_fit$vals$k+1.96*gc_fit$vals$k_se)   # 95% Confidence Interval
lower<- c(lower,gc_fit$vals$k-1.96*gc_fit$vals$k_se)
}

# Simple Exponential Growth Model using latest sample
egm<- lm(log(total) ~ t)
summary(egm)
egm_pred<- exp(predict(egm))               # Convert prediction of "log(total)" to one for "total" itself
Duan<- mean(exp(log(total)-predict(egm)))  
egm_fit<- egm_pred*Duan
                                           # Reduce bias of the "raw" prediction by using Duan's non-parametric smoother
                                           # (See https://davegiles.blogspot.com/2014/12/s.html)

# In R, the calendar starts at 1970-01-01
#day<- t+18322
#class(day) <- "Date"
poi<- round(gc_fit$vals$t_mid+18322,0)
class(poi)<- "Date"

# Plot the results based on the latest data:
par(mfrow=c(1,1))
plot(gc_fit,
     main="Covid-19 Confirmed Cases in Canada", 
          ylab="Cases", xlab="Day",
     type="l")
lines(t,egm_fit, col="blue")
legend("topleft",inset=0.025,
       c("Total Cases","Naive Exponential Model", "Logistic Growth Model"),
       col=c("black","blue","red"), lty=c(5,1,1), pch=c(16,46,46), box.lty=0)
text(5,800, cex=0.8,"(t = 0 is 1 March, 2020)")
text(15,6000, cex=0.9, col="red", paste0("Max. doubling time = ", round(log(2)/gc_fit$vals$r,2), " days"))
text(15,5600, cex=0.9,col="red", paste0("Saturation level = ", round(gc_fit$vals$k,0) , " cases"))
text(15,5200, cex=0.9,col="red", paste0("(95% C.I. = [", round(gc_fit$vals$k-1.96*gc_fit$vals$k_se,0), " , ",round(gc_fit$vals$k+1.96*gc_fit$vals$k_se,0), "] cases)"))
text(15,4000, cex=0.9,col="red", paste0("Median date (point of inflection) = ", poi))
text(26,800, cex=0.8,font=3, paste0("As at: ", today))

# Now plot a summary of the results using the successive sample periods:
Obs<- seq(n_min:n_max)+n_min-1

par(mfrow=c(2,2))
plot(Obs,inflection, main="Estimated Inflection Date",
        ylab= "Inflection (Days)", xlab="Sample Size (Days)", type="b", col="red")
plot(Obs,saturation, main="Estimated Saturation Level",
     ylab= "Saturation (Cases)", xlab="Sample Size (Days)", type="b", ylim=c(min(lower), max(upper)))
lines(Obs, lower, col="red",lty=2)
lines(Obs, upper, col="red",lty=2)
text(28,17000, cex=0.8,col="red", "95% Confidence Band")
plot(Obs,max_doub_time, main="Estimated Max. Doubling Time",
     ylab= "Max. Doubling Time (Days)", xlab="Sample Size (Days)", type="b", col="red")
# Just for interest - what did the naive Exponential have to say?
plot(Obs,exp_growth, main="Estimated Exponential Growth Rate",
     ylab= "Exponential Growth Rate", xlab="Sample Size (Days)", type="b", col="blue")
text(28,0.215, cex=0.8, col="blue", "From Naive Exponential Model")
# END OF FILE  ###################################################################
