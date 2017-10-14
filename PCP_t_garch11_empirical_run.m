addpath(genpath('include/'));

clear all

save_on = false; %true;
plot_on = false;

sdd = 1;        
s = RandStream('mt19937ar','Seed',sdd);
RandStream.setGlobalStream(s); 

data = 3; % MSFT 4; % GSPC 2;
H = 2000;

arg0 = 1;
arg1 = 0;%1;
arg2 = 0;
arg3 = 0;
arg4 = 0;
arg5 = 0;
arg6 = 0;
arg7 = 0;
arg8 = 0;

PCP_t_garch11_empirical;

load(name,'-regexp','^mit','^CV','^mu','^Sigma')