addpath(genpath('include/'));

clear all

save_on = false ;
plot_on = true;

sdd = 1;        
s = RandStream('mt19937ar','Seed',sdd);
RandStream.setGlobalStream(s); 

data = 102; % MSFT 4; % GSPC 2;
H = 1000;

arg0 = 1;
arg1 = 1;%1;
arg2 = 1;
arg3 = 1;
arg4 = 1;
arg5 = 0; % grid
arg6 = 1;
arg7 = 1;
arg8 = 1;

PCP_t_gas_empirical;

% load(name,'-regexp','^mit','^CV','^mu','^Sigma')