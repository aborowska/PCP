clear all

% cd ../
addpath(genpath('include/'));

data = 22;
H = 2275 +248;
model = 'tgarch11';
[y, T, ~, data_name, time] = Load_data_empirical(data, H);   
fprintf('Data loaded: %s\n',data_name);

y_S = var(y(1:T));
%  theta[i] <= nu
%  theta[i+N] <= mu
%  theta[i+2*N] <= omega
%  theta[i+3*N] <= alpha
%  theta[i+4*N] <= beta 
mu_init = [8,0,1,0.1,0.8];
hyper = 0.01;

options = optimset('Display','off');
% w = warning('query','last');
% id = w.identifier;
id = 'optim:fminunc:SwitchingMethod';
warning('off',id);

kernel_init = @(xx) - posterior_t_garch11(xx, y(1:T), y_S, hyper)/T;
[mu_MLE1,~,~,~,~,Sigma1] = fminunc(kernel_init,mu_init,options);
[mu_MLE,~,exitflag,output,~,Sigma] = fminunc(kernel_init,mu_MLE1,options);

% name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'.mat'];
% % load(name,'-regexp','^mu','^Sigma','^draw','^accept','^mit','^CV','\w*_MLE','^THR');
% load(name,'-regexp','^mit');

M = 10000;
BurnIn = 1000;

BurnIn_PCP = BurnIn/10; 
II = 10;
thinning = 1;


p_bar0 = 0.005;
p_bar1 = 0.01;
p_bar = 0.05; 
P_bars = [p_bar0, p_bar1, p_bar];
QUANT_MLE = tinv(P_bars, mu_MLE(1))';
y_sort = sort(y(1:T));
THR_emp = y_sort(round(P_bars*T));

threshold = y_sort(round(2*p_bar*T));
        
x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);

cont = MitISEM_Control;
cont.mit.CV_max = 1; %2?
cont.mit.iter_max = 10;
cont.mit.Hmax = 3; %6;
cont.mit.dfnc = 3;%5;
%     df = 5; % default df for a mit
df = 3; % default df for a mit
cont.disp = true;