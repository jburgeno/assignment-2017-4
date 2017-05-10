
## Montgomery, Hollendbach, and Ward
## Replication files necessary for replicating
## Table 3,  Footnote 21, Table 4, and  Figure 3
## Last edited: Jacob M. Montgomery
## Date: 12-19-2011

# You will want to replace this line with the appropriate directory
setwd("/Users/jmontgomery/Dropbox/EBMA/PRESIDENTIAL DATA for Jacob/Final submission")
library(foreign)
library(nls2)
library(ensembleBMA)

hibbspreds <- read.csv("Prediction_hibbs.csv") # Hibbs predictive intervals generated in Stata.  See *.do file.
master.data<-read.csv("presdata.csv") # read in data



# Function to fit each model on the appropriate training years, fit the EBMA model on the validation period, and make out-of-sample forecasts
my.eBMA <- function(tyn = 15, a = 1, train.years =  10, hibbpreds=hibbspreds){
   master.years <- seq(1952, 2008, by=4)
   in.data <- matrix(NA, (tyn-1), 6)
   out.data <- matrix(NA, 1, 6)
   rownames(in.data) <- master.years[1:(tyn-1)]
   rownames(out.data) <- master.years[tyn]
   colnames(out.data) <- colnames(in.data) <- c("Campbell", "Lewis-Beck", "EWT2C2", "Fair", "Hibbs", "Abramowitz")
   in.master <- master.data[master.data$year >=1952 & master.data$year<master.years[tyn],] #Validation period
   out.master <- master.data[master.data$year == master.years[tyn],] # Test period

   # Fit each of the models on the training data, and estimate y.hat for the validation period
   c.model <- lm(dv~septpoll+gdpqtr2half, data=master.data[master.data$year<master.years[tyn],]) #Training period
   in.data[,1] <- predict(c.model, newdata=in.master) # Validation Period
   out.data[,1] <- predict(c.model, newdata=out.master) #Test Period
    
   lb.model <- lm(dv ~ julypop + incxgnp + jobhousu + closeinc, data=master.data[master.data$year<master.years[tyn],])
   in.data[,2] <- predict(lb.model, newdata=in.master)
   out.data[,2] <- predict(lb.model, newdata=out.master)
   
   ew.model <- lm(dv ~ l1cumleigrowth + incumbentpoll, data=master.data[master.data$year<master.years[tyn],])
   in.data[,3] <- predict(ew.model, newdata=in.master)
   out.data[,3] <- predict(ew.model, newdata=out.master)

   fair.model <- lm(dv ~ G + P + Z + adper + adur + war + I, data=master.data[master.data$year<master.years[tyn],])
   in.data[,4] <- predict(fair.model, newdata=in.master)
   out.data[,4] <- predict(fair.model, newdata=out.master)

   hibbs.model <- nls(dv ~ beta0 + bdlnr* (((1.0*wtq16*dnlr) + (g*dnlr.L1) + ((g^2)*dnlr.L2) + ((g^3)*dnlr.L3) +
                                             ((g^4)*dnlr.L4) + ((g^5)*dnlr.L5) + ((g^6)*dnlr.L6) + ((g^7)*dnlr.L7) +
                                             ((g^8)*dnlr.L8) + ((g^9)*dnlr.L9) + ((g^10)*dnlr.L10) + ((g^11)*dnlr.L11) +
                                             ((g^12)*dnlr.L12)+ ((g^13)*dnlr.L13)+ ((g^14)*dnlr.L14))/
                                            (1.0*wtq16 + g + g^2 + g^3 + g^3 + g^4 + g^5 + g^6 + g^7 + g^8 + g^9 +
                                             g^10 + g^11 + g^12 + g^13 + g^14)) + bkia*fatalities,
                      start=list(beta0=45, g=0.95, bdlnr=4, bkia=-0.1),
                      data=master.data[master.data$year<master.years[tyn],])
   in.data[,5] <- predict(hibbs.model, newdata=in.master)
   out.data[,5] <- predict(hibbs.model, newdata=out.master)
   
   ab.model <- lm(dv~q2gdp+term+juneapp,data=master.data[master.data$year<master.years[tyn],])
   in.data[,6] <- predict(ab.model, newdata=in.master)
   out.data[,6] <- predict(ab.model, newdata=out.master)

   # record the predictive invervals for each individual function
   sixseven <- matrix(nrow=6, ncol=2)
   sixseven[1,]=predict(c.model, newdata=out.master, interval="prediction", level=.67)[2:3]
   sixseven[2,]=predict(lb.model, newdata=out.master, interval="prediction", level=.67)[2:3]
   sixseven[3,]=predict(ew.model, newdata=out.master, interval="prediction", level=.67)[2:3]
   sixseven[4,]=predict(fair.model, newdata=out.master, interval="prediction", level=.67)[2:3]
   sixseven[5,]=c(hibbspreds$X66_low[hibbspreds$Year==master.years[tyn]], hibbspreds$X66_high[hibbspreds$Year==master.years[tyn]]) # Read in from external file
   sixseven[6,]=predict(ab.model, newdata=out.master, interval="prediction", level=.67)[2:3]

   ninezero <- matrix(nrow=6, ncol=2)
   ninezero[1,]=predict(c.model, newdata=out.master, interval="prediction", level=.90)[2:3]
   ninezero[2,]=predict(lb.model, newdata=out.master, interval="prediction", level=.90)[2:3]
   ninezero[3,]=predict(ew.model, newdata=out.master, interval="prediction", level=.90)[2:3]
   ninezero[4,]=predict(fair.model, newdata=out.master, interval="prediction", level=.90)[2:3]
   ninezero[5,]=c(hibbspreds$X90_low[hibbspreds$Year==master.years[tyn]], hibbspreds$X90_high[hibbspreds$Year==master.years[tyn]]) 
   ninezero[6,]=predict(ab.model, newdata=out.master, interval="prediction", level=.90)[2:3]
   
   # Now fit the ebma model
    full.forecasts <- rbind(in.data, out.data) 
    full.observed <- c(master.data$dv[10:(9+tyn)]) # Reducing the validation period to begin in 1952

   # Stupid thing to make data work with ensembleBMA function
   dates <- rep(NA, tyn)
   for (i in 1:tyn){
     dates[i] <- paste("2011", "01", 10+i, "01", sep="")
    }

   pred.date <- dates[tyn]
   my.E.data <- ensembleData(forecasts=(full.forecasts)^(1/a), dates=dates, observations=full.observed,
                             initializationTime=1, forecastHour=1) #Make a dataset of the appropriate format for the ensembleBMA package
   fit.eBMA <- ensembleBMAnormal(my.E.data, trainingDays=train.years, dates=pred.date, minCRPS=TRUE,
                              control=controlBMAnormal(biasCorrection="none")) # Fit the EBMA models
   conf.int <- quantileForecast(fit.eBMA, my.E.data, c(.025, .05, .166666, .5, .833333, .95, .975)) #Make needed confidence intervals for the test-period year
   
   err <- cbind(out.data-full.observed[tyn], conf.int[4]-full.observed[tyn]) #Errors for the EBMA model
   colnames(err) <- c(colnames(in.data), "EBMA")

   observed <- full.observed[tyn]

   # Code whether or not, the test-period observations fell within the predictive intervals
   cov.90 <- c(observed > ninezero[,1] & observed < ninezero[,2],
               observed<=conf.int[6] & observed>=conf.int[2])
   names(cov.90) <- c(colnames(in.data), "EBMA")
   
   cov.67 <- c(observed > sixseven[,1] & observed < sixseven[,2],
               observed<=conf.int[5] & observed>=conf.int[3])
   names(cov.67) <- c(colnames(in.data), "EBMA")
 
   out <- list(in.data=in.data, out.data=out.data, full.observed=full.observed, pred.date=pred.date,
               E.data=my.E.data, fit.eBMA=fit.eBMA, conf.int=conf.int, err=err, observed=full.observed[tyn],
               cov.67=cov.67, cov.90=cov.90, sixseven=sixseven, ninezero=ninezero)
}


######################################################
# Conduct the analysis################################
######################################################

# Fit the EBMA models for each year
fit.2008 <- my.eBMA(tyn=15, train.years=14)
fit.2004 <- my.eBMA(tyn=14, train.years=13)
fit.2000 <- my.eBMA(tyn=13, train.years=12)
fit.1996 <- my.eBMA(tyn=12, train.years=11)
fit.1992 <- my.eBMA(tyn=11, train.years=10)
fit.1988 <- my.eBMA(tyn=10, train.years=9)
fit.1984 <- my.eBMA(tyn=9, train.years=8)
fit.1980 <- my.eBMA(tyn=8, train.years=7)
fit.1976 <- my.eBMA(tyn=7, train.years=6)

# Table 3
# Results for 2004
t3.res.04<-matrix(NA, nrow=6, ncol=4)
colnames(t3.res.04)<-c("Weights", "RMSE", "MAE", "Pred.Error")
rownames(t3.res.04)<-c(colnames(fit.2004$in.data))
model.errors.2004 <- fit.2004$in.data-fit.2004$full.observed[1:13]
t3.res.04[,"RMSE"]<-sqrt(colMeans((model.errors.2004^2)))
t3.res.04[,"MAE"]<-colMeans(abs(model.errors.2004))
t3.res.04[,"Pred.Error"]<-fit.2004$err[1:6]
t3.res.04[,"Weights"]<-fit.2004$fit.eBMA$weights
EBMA.errors <- ((fit.2004$in.data)%*%fit.2004$fit.eBMA$weights)-(fit.2004$full.observed[1:13])
EBMA<-c(NA, sqrt(mean(EBMA.errors^2)), mean(abs(EBMA.errors)), fit.2004$err[7])
t3.res.04<-rbind(t3.res.04, EBMA)
t3.res.04<-round(t3.res.04, digits=3)
print(t3.res.04[c(1,6,5,4,2,3,7),], digits=3, na.print="", zero.print=2)

#Table 3 Results for 2008
t3.res.08<-matrix(NA, nrow=6, ncol=4)
colnames(t3.res.08)<-c("Weights", "RMSE", "MAE", "Pred.Error")
rownames(t3.res.08)<-c(colnames(fit.2008$in.data))
model.errors.2008 <- fit.2008$in.data-fit.2008$full.observed[1:14]
t3.res.08[,"RMSE"]<-sqrt(colMeans((model.errors.2008^2)))
t3.res.08[,"MAE"]<-colMeans(abs(model.errors.2008))
t3.res.08[,"Pred.Error"]<-fit.2008$err[1:6]
t3.res.08[,"Weights"]<-fit.2008$fit.eBMA$weights
EBMA.errors <- ((fit.2008$in.data)%*%fit.2008$fit.eBMA$weights)-(fit.2008$full.observed[1:14])
EBMA<-c(NA, sqrt(mean(EBMA.errors^2)), mean(abs(EBMA.errors)), fit.2008$err[7])
t3.res.08<-rbind(t3.res.08, EBMA)
t3.res.08<-round(t3.res.08, digits=3)
print(t3.res.08[c(1,6,5,4,2,3,7),], digits=3, na.print="", zero.print=2)

# Footnote 21
round(cor(fit.2008$in.data[,c(1,6,5,4,2,3)]), digits=2)


# Table 4
t4<-matrix(NA, nrow=7, ncol=4)
rownames(t4)<-colnames(fit.1976$err)
colnames(t4)<- c("RMSE", "MAE", "67", "90")
errors <- rbind(fit.1976$err, fit.1980$err, fit.1984$err, fit.1988$err, fit.1992$err, fit.1996$err ,fit.2000$err, fit.2004$err, fit.2008$err)
t4[,"RMSE"]<-sqrt(colMeans(errors^2))
t4[,"MAE"] <- colMeans(abs(errors))
coverage.67 <- rbind(fit.1976$cov.67, fit.1980$cov.67, fit.1984$cov.67, fit.1988$cov.67, fit.1992$cov.67, fit.1996$cov.67 ,fit.2000$cov.67, fit.2004$cov.67, fit.2008$cov.67)
t4[,"67"] <- colMeans(coverage.67)
coverage.90 <- rbind(fit.1976$cov.90, fit.1980$cov.90, fit.1984$cov.90, fit.1988$cov.90, fit.1992$cov.90, fit.1996$cov.90 ,fit.2000$cov.90, fit.2004$cov.90, fit.2008$cov.90)
t4[,"90"]<-colMeans(coverage.90)
round(t4[c(1,6,5,4,2,3,7),], digits=2)



#Figure 3

#A function to make each sub-plots for the main figure.  Some graphic paramters to control the size of each
plot.fun <- function(obj=fit.2008, year="2008", left=35, right=65, right.tex=5){
  order.vec <- c(6+2,2,1,3+1,4+2,5+2)
  ebma.vec <- 7+3
  plot(obj$out.data, order.vec, xlim=c(left, right), ylim=c(0, 12.5), xlab="", ylab="", cex=1.25, pch=19)
  segments(obj$sixseven[,1], order.vec, obj$sixseven[,2], order.vec, lwd=3)
  segments(obj$ninezero[,1], order.vec, obj$ninezero[,2], order.vec, lwd=1)
  abline(v=obj$observed, lty=2)
  points(obj$conf.int[4], ebma.vec,  pch=19, cex=1.25)
  segments(obj$conf.int[3], ebma.vec, obj$conf.int[5], ebma.vec, lwd=3)
  segments(obj$conf.int[2], ebma.vec, obj$conf.int[6], ebma.vec, lwd=1)
  text(left+right.tex, 12 ,year, cex=2, col="gray")
  text(rep(left+right.tex, 6), order.vec+.35, colnames(obj$in.data))
  text(left+right.tex, ebma.vec+.35, "EBMA")
}

# Actually make the plots
par(mfrow=c(3,3), mar=c(1.25,.5,.5,.5), las=1, tcl=.5, mgp=c(1.5,.25,0), yaxt="n")
plot.fun(obj=fit.2008, year="2008", left=35, right=57, right.tex=2.9)
plot.fun(obj=fit.2004, year="2004", left=44, right=63, right.tex=2.5)
plot.fun(obj=fit.2000, year="2000", left=42, right=60, right.tex=2.2)
plot.fun(obj=fit.1996, year="1996", left=47, right=62, right.tex=1.85)
plot.fun(obj=fit.1992, year="1992", right=60, right.tex=3.1)
plot.fun(obj=fit.1988, year="1988", right=61, left=45, right.tex=2)
plot.fun(obj=fit.1984, year="1984", left=47, right=66, right.tex=2.35)
plot.fun(obj=fit.1980, year="1980", left=10, right=85, right.tex=8.7)
plot.fun(obj=fit.1976, year="1976", left=24, right=77, right.tex=6.65)



