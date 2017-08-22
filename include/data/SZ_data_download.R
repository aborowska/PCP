install.packages("quantmod")
require("quantmod")
sp500 <- new.env()
AAPL <- new.env()
getSymbols("AAPL", env = AAPL, src = "google", from = as.Date("1996-04-15"), to = as.Date("2015-10-05"))

getSymbols("INDEXSP:.INX", env = sp500, src = "google", from = as.Date("1996-04-15"), to = as.Date("2015-10-05"))
GSPC <- sp500$GSPC
Price <- GSPC$GSPC.Adjusted
ret <- diff(log(Price))
setwd("C:/Users/aba228/Dropbox/New Projects/Censored Posterior/PCP/include/data") 
write.csv(ret,file="SZ_snp_data.csv")
