# "JUNK" Regression
#
# David Giles; davegiles.ca; davegiles1949@gmail.com; March 2026
################################################################

x<- c(1.5,2,2.5,2,2,2,1.7,1.5,1.3,3.5,3.5,3.5,3.7,3.9,4.1,4.3,4.5,4.5,
4.5,5.5,5.5,5.5,5.5,6.5,6.5,6.5,6.5,6.0,5.75,6.25,7.5,7.5,7.5,7.5,7.75,8.0,8.25,8.0,8.25,8.5)
y<- c(5,5,5,4,3,2,1.7,1.8,2.0,8.0,7.0,6.0,5.5,5.3,5.3,5.5,6.0,7.0,8.0,
12.0,11.0,10.0,9.0,12.0,11.0,10.0,9.0,10.5,11.5,9.5,15.5,14.5,13.5,12.5,14.0,14.5,15.0,13.5,13.0,12.5)

model<- lm(y ~ x)
summary(model)
plot(x, y, pch=25, col="blue",xlim=c(0,20), ylim=c(0,10))
abline(model, col = "red", lwd = 3, lty=2)

plot(x, y, pch = 16,xlim=c(1,10), ylim=c(0,20), main = "Spurious Regression !")
abline(model, col = "red", lwd = 3)
coef <- round(coef(model), 3)
se<- round(sqrt(diag(vcov(model))),3)
text(3, 15,  paste("y = ", coef [1], "+", coef [2], "x"),col="red")
text(3.2,14, paste("(",se [1],")","(",se [2],")"),cex=0.8)
rsq<- round(summary(model)$r.squared,3)
text(3,13, paste("R-squared = ",rsq))

