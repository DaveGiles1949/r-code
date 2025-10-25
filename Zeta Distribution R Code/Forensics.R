library(fitPS)
data(Psurveys) 
roux<- as.data.frame(Psurveys$roux)
fit.z<- fitDist(Psurveys$roux)
fit.z
write.table(roux,file="C:/Users/David/Sync/Power Laws/Bias Correction/Comm Stats Submission/roux.txt", row.names=FALSE)
