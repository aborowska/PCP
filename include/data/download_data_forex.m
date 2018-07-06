% 04-01-1999 : 31-08-2017
% GBP JPY USD / EUR
M = csvread('forex.csv',1,2);
M = M(:,1:3);

M = 100*diff(log(M));

GBP = M(:,1);
JPY = M(:,2);
USD = M(:,3);

plot(M)

histogram(GBP)
hold all
histogram(JPY)
histogram(USD)


load Data_MarkPound
r  = price2ret(Data);
pR = 100*r;
T  = length(r);