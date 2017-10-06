require(quantmod)
require(PerformanceAnalytics)


Symbols<-c  ("XOM","MSFT","JNJ","GE","CVX","WFC","PG","JPM","VZ","PFE","T","IBM","MRK","BAC","DIS","ORCL","PM","INTC","SLB")
length(Symbols)

#Set start date
start_date=as.Date("1998-01-01")
end_date=as.Date("2017-08-31")

#Create New environment to contain stock price data
dataEnv<-new.env()

getSymbols("GOOG", src = "yahoo",from=start_date, end=end_date)

#download data          
getSymbols("MSFT", src = "google",from=start_date, end=end_date)
getSymbols("GOOG", src = "google",from=start_date, end=end_date)
getSymbols("GE", src = "google",from=start_date, end=end_date)
getSymbols("IBM", src = "google",from=start_date, end=end_date)
getSymbols("AAPL", src = "google",from=start_date, end=end_date)
getSymbols("JPM", src = "google",from=start_date, end=end_date)
getSymbols("SPX", src = "google",from=start_date, end=end_date)


GOOG_cl <- GOOG[,4]
GE_cl <- GE[,4]
MSFT_cl <- MSFT[,4]
IBM_cl <- IBM[,4]
AAPL_cl <- AAPL[,4]
JPM_cl <- JPM[,4]
GSPC_cl <- GSPC[,4]



GE_r <- diff(log(GE_cl))[-1]
GOOG_r <- diff(log(GOOG_cl))[-1]
MSFT_r <- diff(log(MSFT_cl))[-1]
IBM_r <- diff(log(IBM_cl))[-1]
JPM_r <- diff(log(JPM_cl))[-1]
GSPC_r <- diff(log(GSPC_cl))[-1]


par(mfrow=c(3,2))
plot(as.double(GE_r),xlab="",ylab="",type='l',main='GE')
plot(as.double(GOOG_r),xlab="",ylab="",type='l',main='GOOG')
plot(as.double(IBM_r),xlab="",ylab="",type='l',main='IBM')
plot(as.double(MSFT_r),xlab="",ylab="",type='l',main='MSFT')
plot(as.double(JPM_r),xlab="",ylab="",type='l',main='JPM')
plot(as.double(GSPC_r),xlab="",ylab="",type='l',main='GSPC')


