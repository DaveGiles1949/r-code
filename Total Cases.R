# CODE FOR MODELLING THE GROWTH CURVE OF TOTAL COVID-19 CASES IN CANADA
# ---------------------------------------------------------------------

# DATA SOURCE: https://www.covid-19canada.com/

# AUTHOR: DAVID GILES (davegiles1949@gmail.com)

# ------------------------------------------------------------------------------------------------

library(growthcurver)
library(curl)
library(data.table)    # For the "fread" command to read data file from webpage
n_min<- 26             # Smallest number of days to include in the sequential anlaysis
n_max<- 30             # Largest number of days to include in the sequential analysis
                       # NOTE: "Day 1" is 2 March, 2020

upper<- c()
lower<- c()
saturation<- c()
inflection<- c()
max_doub_time<- c()
exp_growth<- c()

file_name <- "C:/Users/David Giles/Desktop/Virus/total_cases.txt"
cases <- read.table(file_name, header = TRUE)
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
egm_pred<- exp(predict(egm))                 # Convert prediction of "log(total)" to one for "total" itself
Duan<- mean(exp(log(total)-predict(egm)))  
egm_fit<- egm_pred*Duan
# Adjust the "raw" prediction using Duan's non-parametric smoother
# (See https://davegiles.blogspot.com/2014/12/s.html)

# Plot the results based on the latest data:
par(mfrow=c(1,1))
plot(gc_fit,
     main="Covid-19 Cases in Canada",
     ylab="Cases",
     type="l")
lines(t,egm_fit, col="blue")
legend(1,8500,
       c("Total Cases","Naive Exponential Model", "Logistic Growth Model"),
       col=c("black","blue","red"), lty=c(5,1,1), pch=c(16,46,46), box.lty=0)
text(5,800, "(t = 0 is 1 March, 2020)")
text(15,6000, col="red", "Max. doubling time = 2.7 days")
text(15,5600, col="red", "Saturation level = 15,182 cases")
text(15,5200, col="red", "(95% C.I. = [13,305 , 16,960] cases)")
text(15,4000, col="red", "Median date (point of inflection) = Day 29 (i.e., March 30)")

# Now plot a summary of the results using the successive sample periods
Obs<- seq(n_min:n_max)+n_min-1

par(mfrow=c(2,2))
plot(Obs,inflection, main="Estimated Inflection Date",
        ylab= "Inflection (Days)", xlab="Sample Size (Days)", type="b", col="red")

plot(Obs,saturation, main="Estimated Saturation Level",
     ylab= "Saturation (Cases)", xlab="Sample Size (Days)", type="b", ylim=c(min(lower), max(upper)))
lines(Obs, lower, col="red",lty=2)
lines(Obs, upper, col="red",lty=2)
text(28,30000, col="red", "95% Confidence Band")
plot(Obs,max_doub_time, main="Estimated Max. Doubling Time",
     ylab= "Max. Doubling Time (Days)", xlab="Sample Size (Days)", type="b", col="red")
# Just for interest - what did the naive Exponential have to say?
plot(Obs,exp_growth, main="Estimated Exponential Growth Rate",
     ylab= "Exponential Growth Rate", xlab="Sample Size (Days)", type="b", col="blue")

# END OF FILE  ###################################################################
