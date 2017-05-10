## Montgomery, Hollendbach, and Ward
## Replication files necessary for replicating
## Table 5
## Last edited: Jacob M. Montgomery
## Date: 12-19-2011

# You will want to replace this line with the appropriate directory
setwd("/Users/jmontgomery/Dropbox/EBMA/PRESIDENTIAL DATA for Jacob/Final submission")
library(foreign)
library(ensembleBMA)
library(EBMAforecast) # This package will be distribute on CRAN

# Read in data.  
courts <- read.csv("CourtsData.csv")

# Dividing the data into validation/test sets based on docket number
in.sample.docket <- unique(courts$docket)[1:44]
out.sample.docket <- unique(courts$docket)[45:68]

courts$in.samp <- 1

for (i in out.sample.docket){
  courts$in.samp[courts$docket==i] <- 0
}

final <- courts[!is.na(courts$act),]#Removing missing
in.sample <- final[final$in.samp==1,]
out.sample <- final[final$in.samp==0,]


### Specifying the two "model's" predictions
in.pred.mod <- in.sample$mod # The Model's prediciton
in.pred.exp <- in.sample$exp.tot # the average prediction of the experts

# Test-period forecasts
out.pred.mod <- out.sample$mod
out.pred.exp <- out.sample$exp.tot

in.data <- cbind(in.pred.mod, in.pred.exp)
out.data <- cbind(out.pred.mod, out.pred.exp)

# Necessary to keep the bivariate regression from exploding
in.data[in.data==1] <- .999
in.data[in.data==0] <- .001
out.data[out.data==1] <- .999
out.data[out.data==0] <- .001

colnames(in.data) <- colnames(out.data) <- c("Model", "Experts")

#How each justice actually voted on each case
y.in <- in.sample$act
y.out <- out.sample$act

#Validation-sample fit 
BMA.fit <- Ensemble.logit(y=y.in, pp.raw=in.data, tol=.0001, exp=4)

#Table 5
BMA.pred <- predict.Ensemble.logit(obj=BMA.fit, newdata=out.data, y.out=y.out)
print.Ensemble.logit(BMA.pred)



