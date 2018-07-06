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



% shorter insample?
start_date = '01012002';
end_date = '31122015';
tickers = {'^GSPC' 'IBM' 'AAPL' 'MSFT' 'JPM' 'GE'};

data = hist_stock_data(start_date, end_date, tickers);
data_AdjClose = [data([1:end]).AdjClose];
data_Returns = 100*diff(log(data_AdjClose));
csvwrite('Perc_Rets_GSPC_IBM_AAPL_MSFT_JPM_GE_short.csv',data_Returns);




% shorter insample? 1 year ~ 250
% % % IN: 2004 2005 2006 2007 
% % % OUT: 2008 2009 2010 2011
start_date = '01012004'; 
end_date = '31122011';
tickers = {'^GSPC' 'IBM' 'MSFT'};

data = hist_stock_data(start_date, end_date, tickers);
data_AdjClose = [data([1:end]).AdjClose];
data_Returns = 100*diff(log(data_AdjClose));
% % % % % csvwrite('Perc_Rets_GSPC_IBM_MSFT_T2000_crisis.csv',data_Returns(:,1:3));


%% IN: 2010 2011 2012 2013 
% OUT: 2014 2015 2016 2017

start_date = '01012010'; 
end_date = '31122017';
tickers = {'^GSPC' 'IBM' 'MSFT'};

data = hist_stock_data(start_date, end_date, tickers);
data_AdjClose = [data([1:end]).AdjClose];
data_Returns = 100*diff(log(data_AdjClose));
csvwrite('Perc_Rets_IBM_MSFT_T2000.csv',data_Returns(:,2:3));

T = size(data_Returns,1);
AX = [0 T min(min(data_Returns)) max(max(data_Returns))];
for ii = 1:3
    subplot(1,3,ii)
    plot(data_Returns(:,ii))
    axis(AX)
    line([T-1000 T-1000], AX(3:4),'Color',[1 0 0])
end 



figure(2)
subplot(2,2,1)
histogram(data_Returns(1:(T-1000),1))
subplot(2,2,2)
histogram(data_Returns(((T-1000)+1):end,1))

subplot(2,2,3)
histogram(data_Returns(1:(T-1000),2),'FaceColor',[0.7 0.2 0.1])
subplot(2,2,4)
histogram(data_Returns(((T-1000)+1):end,2),'FaceColor',[0.7 0.2 0.1])





% shorter insample? 1 year ~ 250
% % % IN: 2007 2008 2009 2010
% % % OUT: 2011 2012 2013 2014
start_date = '01012007'; 
end_date = '31122014';
tickers = {'^GSPC' 'IBM' 'MSFT'};

data = hist_stock_data(start_date, end_date, tickers);
data_AdjClose = [data([1:end]).AdjClose];
data_Returns = 100*diff(log(data_AdjClose));
csvwrite('Perc_Rets_GSPC_IBM_MSFT_T2000_crisis.csv',data_Returns(:,1:3));


%%% 
% start_date = '01012015'; 
% end_date = '31122017';
% tickers = {'^GSPC' 'IBM' 'MSFT'};
% 
% data = hist_stock_data(start_date, end_date, tickers);
% data_AdjClose = [data([1:end]).AdjClose];
% data_Returns = 100*diff(log(data_AdjClose));

start_date = '01012007'; 
end_date = '23122016';  % T = 1012 + H = 1500
tickers = {'^GSPC' 'IBM' 'MSFT'};

data = hist_stock_data(start_date, end_date, tickers);
data_AdjClose = [data([1:end]).AdjClose];
data_Returns = 100*diff(log(data_AdjClose));
csvwrite('Perc_Rets_GSPC_IBM_MSFT_T2000_crisis2.csv',data_Returns(:,1:3));
