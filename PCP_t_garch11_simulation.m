clear all
addpath(genpath('include/'));

model = 't_garch';
plot_on = false;
sdd=1;
s = RandStream('mt19937ar','Seed',sdd);
RandStream.setGlobalStream(s); 
p_bar = 0.1;

%  theta[i] <= nu
%  theta[i+N] <= mu
%  theta[i+2*N] <= omega
%  theta[i+3*N] <= alpha
%  theta[i+4*N] <= beta 
 
T = 2000;
H = 2000;

sigma1 = 1;
sigma2 = 4; 
c = (sigma2 - sigma1)/sqrt(2*pi); % mean of eps
kappa = 0.5*(sigma1^2 + sigma2^2 - ((sigma2-sigma1)^2)/pi); % var of eps
sigma1_k = sigma1/sqrt(kappa);
sigma2_k = sigma2/sqrt(kappa);

omega = 1;
alpha = 0.1;
beta = 0.8;
%mu_true = [0, omega, alpha, beta];
param_true = [c,sigma2,omega,alpha,beta];
mu_init = [8,0,1,0.1,0.8];

x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);
    
M = 10000; % number of draws 
BurnIn = 1000;

cont = MitISEM_Control;
cont.mit.CV_max = 1; %2?
cont.mit.iter_max = 10;
cont.mit.Hmax = 3; %6;
cont.mit.dfnc = 3;%5;
df = 3; % default df for a mit
cont.disp = true;


thinning = 10;
BurnIn_PCP = BurnIn/10;

options = optimset('Display','off');
% w = warning('query','last');
% id = w.identifier;
id = 'optim:fminunc:SwitchingMethod';
warning('off',id);

%% GARCH(1,1)   
eps = randn(T+H,1);
ind = (eps>0);
eps(ind) = c + sigma1.*eps(ind);
eps(~ind) = c + sigma2.*eps(~ind);
eps = eps/sqrt(kappa);  
y = zeros(T+H,1);
h_true = zeros(T+H,1);

for ii = 1:T+H
    if (ii == 1)
        h_true(ii,1) = omega;
    else
        h_true(ii,1) = omega*(1-alpha-beta) + alpha*(y(ii-1,1))^2 + beta*h_true(ii-1,1);
    end
    y(ii,1) = sqrt(h_true(ii,1))*eps(ii,1);
end

%% Uncensored Posterior
fprintf('*** Uncensored Posterior ***\n');
y_S = var(y(1:T));

hyper = 0.01; 
kernel_init = @(xx) - posterior_t_garch11(xx, y(1:T), y_S, hyper)/T;
kernel = @(xx) posterior_t_garch11_mex(xx, y(1:T), y_S, GamMat, hyper);
[mu_MLE,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init,options);
Sigma = inv(T*Sigma);    
%  sim:      11.2279    0.0842    1.1092    0.1387    0.7268


%     kernel_init = @(xx) -posterior_t_garch11_mex(xx, y(1:T), y_S)/T;
%     kernel = @(xx) posterior_garch11_mex(xx, y(1:T), y_S);
try
    [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
    [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
    lnd = dmvgt(draw, mit, true, GamMat); 
catch
    draw = rmvt(mu_MLE,Sigma,df,M+BurnIn);
    mit = struct('mu',mu_MLE,'Sigma',reshape(Sigma,1,length(mu_MLE)^2),'df', df, 'p', 1);
    [mit, CV] = MitISEM_new(mit, kernel, mu_init, cont, GamMat);            
    [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
    lnd = dmvgt(draw, mit, true, GamMat); 
end

lnw = lnk - lnd;
lnw = lnw - max(lnw);
[ind, a] = fn_MH(lnw);
draw = draw(ind,:);
accept = a/(M+BurnIn);
draw = draw(BurnIn+1:BurnIn+M,:);    

% compute the implied volatility for the last in-sample period
hT = volatility_t_garch11(draw,y(1:T),y_S,0);  
 
mean_draw = mean(draw);
median_draw = median(draw);
std_draw = std(draw);

%% CENSORED: Threshold = 10% perscentile of the data sample
fprintf('*** Censored Posterior, threshold 10%% ***\n');
threshold = sort(y(1:T));
threshold = threshold(round(2*p_bar*T));

kernel_init = @(xx) - C_posterior_t_garch11(xx, y(1:T,1), threshold, y_S, hyper)/T;    
[mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
Sigma_C = inv(T*Sigma_C);  
%    47.3460    0.2598    2.3792    0.2650    0.4609

%%  PARTIALLY CENSORED based on a grid
grid = 2.1:0.1:20;
% N = 180;
% grid = linspace(-1.5,1.5,180)
kernel = @(xx) C_posterior_t_garch11_2_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);
draw_PC = sim_cond_inv_trans(draw, grid, kernel);

mean_draw_PC = mean(draw_PC);
median_draw_PC = median(draw_PC);
std_draw_PC = std(draw_PC);

% compute the implied volatility for the last in-sample period
hT_PC = volatility_t_garch11(draw_PC,y(1:T),y_S,0);  

%% CENSORED MLE PARAMETERS    
threshold_m = 0.1; %<---------- HiGhER?
quantile = tinv(threshold_m, mu_MLE(1));
fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold_m);
kernel_init = @(xx) - C_posterior_t_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S, GamMat, hyper)/T;    
kernel = @(xx) C_posterior_t_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S, GamMat, hyper);

kernel(mu_C)
draw_PCm = sim_cond_inv_trans(draw, grid, kernel);


mean_draw_PCm = mean(draw_PCm);
median_draw_PCm = median(draw_PCm);
std_draw_PCm = std(draw_PCm);


% compute the implied volatility for the last in-sample period
hT_PCm = volatility_t_garch11(draw_PCm,y(1:T),y_S,0);  



if plot_on
    figure(1)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');

    hold on
    histogram(draw(:,1))
    histogram(draw_PC(:,1))
    histogram(draw_PCm(:,1))
    hold off 
    
    
    hist([draw(:,1),draw_PC(:,1),draw_PCm(:,1)],20)
    
    leg = legend('Posterior','PC 10%','PC var MLE');  
    set(leg,'Interpreter','latex','FontSize',10)
    title('Draws of \nu in GARCH-t (sim. data)')
end


name = ['PCP_sim_',model,'_sigma2_',num2str(sigma2),'_T',num2str(T),'.mat'];
save(name,...
    'y','draw','draw_PC','draw_PCm',...
    'mean_draw','mean_draw_PC','mean_draw_PCm',...
    'median_draw','median_draw_PC','median_draw_PCm',...
    'std_draw','std_draw_PC','std_draw_PCm',...
    'accept','mit','CV');