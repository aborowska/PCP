# setwd('../../../New Projects/Censored Posterior/PCP/include/data/')

data <- read.csv('Perc_Rets_GSPC_IBM_AAPL_MSFT_JPM_GE.csv')
data2 <- read.csv('Perc_Rets_KO_T_WMT_XOM.csv')
data3 <- read.csv('Perc_Rets_FORD_BAC_SBUX_PVH.csv')

tickers <- c("GSPC","IBM","AAPL","MSFT","JPM","GE")
tickers2 <- c("KO","T","WMT","XOM")
tickers3 <- c('FORD','BAC','SBUX','PVH')

tmin <- as.Date("1998-01-03")
tmax <- as.Date("2017-08-31")
T = nrow(data3)
time <- seq(tmin,tmax, length.out = T)

plot_on = FALSE

require("rugarch")
newspec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                   mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                   distribution.model = 'std')

param_names <- c("mu", "omega", "alpha", "beta", "nu")

# install.packages("moments")
require("moments")
mom_names <- c( "Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","Skew","Kurt")

if (plot_on==TRUE){
  par(mfrow=c(3,2))
  for (i in 1:6){
    plot(as.double(data[,i])~time,xlab="",ylab="",type='l',main=tickers[i])
  }
  par(mfrow=c(2,2))
  for (i in 1:4){
    plot(as.double(data2[,i])~time,xlab="",ylab="",type='l',main=tickers2[i])
  }
  par(mfrow=c(2,2))
  for (i in 1:4){
    plot(as.double(data3[,i])~time,xlab="",ylab="",type='l',main=tickers3[i])
  }  
}

##### DATA 1: FULL SAMPLE
resid_mom <- matrix(rep(NaN,30),ncol=8,nrow=6)
colnames(resid_mom) <- mom_names
rownames(resid_mom) <- tickers

coefs <- matrix(rep(NaN,30),ncol=5,nrow=6)
colnames(coefs) <- param_names
rownames(coefs) <- tickers

par(mfrow=c(3,2))
for (i in 1:6){
  fitted <- ugarchfit(newspec, data[,i])
  coefs[i,] <- fitted@fit$coef
  # resid <- fitted@fit$residuals
  resid <- as.vector(residuals(fitted,standardize=TRUE))
  resid_mom[i,] <- c(summary(resid),skewness(resid),kurtosis(resid))
  if (plot_on==TRUE){
    hist(resid,breaks = 40,main=tickers[i])
  }
}

##### DATA 1: IN-SAMPLE
resid_mom_short <- matrix(rep(NaN,30),ncol=8,nrow=6)
colnames(resid_mom_short) <- mom_names
rownames(resid_mom_short) <- tickers

coefs_short <- matrix(rep(NaN,30),ncol=5,nrow=6)
colnames(coefs_short) <- param_names
rownames(coefs_short) <- tickers

par(mfrow=c(3,2))
for (i in 1:6){
  fitted <- ugarchfit(newspec, data[1:(T-2000),i])
  # resid <- fitted@fit$residuals
  coefs_short[i,] <- fitted@fit$coef
  resid <- as.vector(residuals(fitted,standardize=TRUE))
  resid_mom_short[i,] <- c(summary(resid),skewness(resid),kurtosis(resid))
  if (plot_on==TRUE){
    hist(resid,breaks = 40,main=tickers[i])
  }
}

##### DATA 2: FULL SAMPLE
resid_mom2 <- matrix(rep(NaN,30),ncol=8,nrow=4)
colnames(resid_mom2) <- mom_names
rownames(resid_mom2) <- tickers2

coefs2 <- matrix(rep(NaN,30),ncol=5,nrow=4)
colnames(coefs2) <- param_names
rownames(coefs2) <- tickers2

par(mfrow=c(2,2))
for (i in 1:4){
  fitted <- ugarchfit(newspec, data2[,i])
  # resid <- fitted@fit$residuals
  coefs2[i,] <- fitted@fit$coef
  resid <- as.vector(residuals(fitted,standardize=TRUE))
  resid_mom2[i,] <- c(summary(resid),skewness(resid),kurtosis(resid))
  if (plot_on==TRUE){
    hist(resid,breaks = 40,main=tickers2[i])
  }
}

##### DATA 2: IN-SAMPLE
resid_mom_short2 <- matrix(rep(NaN,30),ncol=8,nrow=4)
colnames(resid_mom_short2) <- mom_names
rownames(resid_mom_short2) <- tickers2

coefs_short2 <- matrix(rep(NaN,30),ncol=5,nrow=4)
colnames(coefs_short2) <- param_names
rownames(coefs_short2) <- tickers2

par(mfrow=c(2,2))
for (i in 1:4){
  fitted <- ugarchfit(newspec, data2[1:(T-2000),i])
  # resid <- fitted@fit$residuals
  coefs_short2[i,] <- fitted@fit$coef
  resid <- as.vector(residuals(fitted,standardize=TRUE))
  resid_mom_short2[i,] <- c(summary(resid),skewness(resid),kurtosis(resid))
  if (plot_on==TRUE){
    hist(resid,breaks = 40,main=tickers2[i])
  }
}



##### Bine and save
coefs_all <- rbind(coefs,coefs2)
coefs_short_all <- rbind(coefs_short,coefs_short2)

resid_mom_all <- rbind(resid_mom,resid_mom2)
resid_mom_short_all <- rbind(resid_mom_short,resid_mom_short2)

save(file="results_rugarch_all_data.RData", time, newspec,
     resid_mom_all, resid_mom_short_all,
     coefs_all, coefs_short_all)