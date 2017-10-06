% s(1).x = 1;
% s(1).y = 2;
% s(2).x = 3;
% s(2).y = 4;
% s(3).x = 5;
% s(3).y = 6;
% [s([1 2]).x] % ans =  1 3


start_date = '01011998';
end_date = '31082017';

tickers = {'^GSPC' 'IBM' 'AAPL' 'MSFT' 'JPM' 'GE'};
tickers2 = {'KO', 'T', 'WMT', 'XOM'}; % NAIS: GE JP Morgan Coca-Cola AT&T Wal-Mart Exxon
tickers3 ={'FORD','BAC','SBUX','PVH'};


AAPL = hist_stock_data(start_date, end_date, 'AAPL');
date_aapl = AAPL.Date;
AAPL2 = AAPL.AdjClose;
AAPL = AAPL.Close;

IBM = hist_stock_data(start_date, end_date, 'IBM');
date_ibm = IBM.Date;
IBM2 = IBM.AdjClose;
IBM = IBM.Close;

data = hist_stock_data(start_date, end_date, tickers);
data_AdjClose = [data([1:end]).AdjClose];
data_Returns = 100*diff(log(data_AdjClose));
csvwrite('Perc_Rets_GSPC_IBM_AAPL_MSFT_JPM_GE.csv',data_Returns);




data2 = hist_stock_data(start_date, end_date, tickers2);
data_AdjClose2 = [data2([1:end]).AdjClose];
data_Returns2 = 100*diff(log(data_AdjClose2));
csvwrite('Perc_Rets_KO_T_WMT_XOM.csv',data_Returns2);


data3 = hist_stock_data(start_date, end_date, tickers3);
data_AdjClose3 = [data3([1:end]).AdjClose];
data_Returns3 = 100*diff(log(data_AdjClose3));
csvwrite('Perc_Rets_FORD_BAC_SBUX_PVH.csv',data_Returns3);



