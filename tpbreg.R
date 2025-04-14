rm(list=ls())


library(betareg)
library(ggplot2)
library(numDeriv)
library(twopiece)


source("routines_tpbetareg.R")

####################################################################################
# Example
####################################################################################

# Using betareg
library(betareg)
data("GasolineYield", package = "betareg") 
head(GasolineYield)
n <- nrow(GasolineYield)

## Table 1 
gy <- betareg(yield ~ batch + temp, data = GasolineYield, link = "logit") 
summary(gy)

# Direct implementation
Xb <- matrix(0, ncol = 10, nrow = n)
Xb[cbind(1:n,GasolineYield$batc)] <- 1

X <- cbind(1,Xb[,1:9],GasolineYield$temp)
y <- GasolineYield$yield

# It is recommended to use different methods and initial points
OPT <- MLEBR(y = y, des = X, init = rep(0,12), method = "nlminb", control = list( iter.max = 10000))

OPT2 <- MLETPBR(y = y, des = X, init = rep(0,13), method = "nlminb", control = list( iter.max = 10000))

# Comparison
cbind(coef(gy),OPT$MLE)

c(gy$loglik, - OPT$OPT$objective)

OPT$loglik(c(coef(gy)[1:11],log(coef(gy)[12])))

OPT$loglik(OPT$OPT$par)

# Model comparison

c(OPT$AIC, OPT2$AIC)

c(OPT$BIC, OPT2$BIC)





