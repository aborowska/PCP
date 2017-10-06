function [y, T, y_plot, data_name, time] = Load_data_empirical(data, H)   
    
    switch data
        case 1 % SZ SnP500 data
            y = load('GSPC.txt'); % 15-04-1996 : 05-10-2015
            y = 100*diff(log10(y));
            TT = length(y);
            T = TT - H; 
            y_plot = y;
            y = y(1:(T+H));
    %         time = linspace((1996 + 3.5/12),(2015 + (9 + 1/6)/12),TT);
            time = [(1996 + 3.5/12),(2015 + (9 + 1/6)/12)];
            data_name = 'SnP500_SZ';
        case 2 % updated SnP500
            y = load('^GSPC_new.csv'); % 02-01-1998 : 31-08-2017 T=0
            y = 100*diff(log(y));
            time = [1998,(2017 + 8/12)];        
            TT = length(y);
%             T = TT - H;
            T = 2500;
            t2 = time(1) +(T+H)*(time(2)-time(1))/TT;
            time(2) = t2;% roughly: 20-11-2015
%             y_plot = y;
            y = y(1:(T+H));
            y_plot = y;
            data_name = 'SnP500';
        case 3 % updated IBM
            y = load('IBM_new.csv'); % 02-01-1998 : 31-08-2017 T=0
            y = 100*diff(log(y));
            time = [1998,(2017 + 8/12)];        
            TT = length(y);
            T = TT - H;
            y_plot = y;
            y = y(1:(T+H));          
            data_name = 'IBM';
        case 4 % updated MSFT
            y = load('Perc_Rets_GSPC_IBM_AAPL_MSFT_JPM_GE.csv');
            y = y(:,4);
            time = [1998,(2017 + 8/12)];       
            TT = length(y);
            T = TT - H;
            y_plot = y;
            data_name = 'MSFT';            
    end
end    