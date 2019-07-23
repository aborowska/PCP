addpath(genpath('include/'));

clear all

save_on = false; %true;
plot_on = true; %false;

sdd = 1;        
s = RandStream('mt19937ar','Seed',sdd);
RandStream.setGlobalStream(s); 

data = 103; %104; %44; % MSFT 4; % GSPC 2;
H = 1500; %1000; %2275 +248;

arg0 = 1;
arg1 = 1;%1;
arg2 = 1;
arg3 = 1;
arg4 = 0;
arg5 = 0; % grid
arg6 = 1;
arg7 = 1;
arg8 = 0;

PCP_t_garch11_empirical;

% load(name,'-regexp','^mit','^CV','^mu','^Sigma')