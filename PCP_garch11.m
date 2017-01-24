clear all
close all

addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

model = 'garch11';
parameters = {'$\\mu$','$\\sigma$','$\\phi$'};

sigma1 = 1;
sigma2 = 1;
c = (sigma2 - sigma1)/sqrt(2*pi); % mean of eps
kappa = 0.5*(sigma1^2 + sigma2^2 - ((sigma2-sigma1)^2)/pi); % var of eps
sigma1_k = sigma1/sqrt(kappa);
sigma2_k = sigma2/sqrt(kappa);


T = 1000; % time series length
p_bar1 = 0.01;
p_bar = 0.05;
% Metropolis-Hastings for the parameters
M = 10000; % number of draws 
BurnIn = 1000;

x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);

cont = MitISEM_Control;
cont.mit.iter_max = 10;

df = 5; % default df for a mit

% various display options
cont.disp = true;%false;

v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014a)')
    fn_hist = @(xx) hist(xx,20);
else
    fn_hist = @(xx) histogram(xx,20);
end

plot_on = false;
save_on = false; %true;

options = optimset('Display','off');
% w = warning('query','last');
% id = w.identifier;
id = 'optim:fminunc:SwitchingMethod';
warning('off',id);


%% GARCH(1,1)
eps = randn(T,1);
ind = (eps>0);
eps(ind) = c + sigma1.*eps(ind);
eps(~ind) = c + sigma2.*eps(~ind);
eps = eps/sqrt(kappa);  

omega = 0.01;
% omega = 0.1;
alpha = 0.1;
beta = 0.8;
param_true = [c,sigma2,omega,alpha,beta];
y = zeros(T,1);
h_true = zeros(T,1);
% eps = randn(T,1);
for ii = 1:T
    if (ii == 1)
        h_true(ii,1) = omega/(1-alpha-beta);
    else
        h_true(ii,1) = omega + alpha*(y(ii-1,1))^2 + beta*h_true(ii-1,1);
    end
    y(ii,1) = sqrt(h_true(ii,1))*eps(ii,1);
end

if plot_on
    figure(1)
    hold on
    plot(eps,'Color',[0.64 0.64 0.74])
    plot(y,'k')
    plot(h_true,'r')
    hold off
    leg = legend('$\varepsilon$','$y$','$h$');
    set(leg,'Interpreter','latex','FontSize',11)
    plotTickLatex2D('FontSize',12);
end

% MC VaRs under the true model
eps_sort = randn(M,1);
ind = (eps_sort>0);
eps_sort(ind) = c + sigma1.*eps_sort(ind);
eps_sort(~ind) = c + sigma2.*eps_sort(~ind);
eps_sort = eps_sort/sqrt(kappa); 

y_sort = sqrt(h_true(T,1))*eps_sort;
y_sort = sort(y_sort);
VaR_1 = y_sort(p_bar1*M); 
VaR_5 = y_sort(p_bar*M); 

%% Misspecified model: GARCH(1,1) normal 
%% Uncensored Posterior
% theta  = [mu, omega, alpha, beta]
% mu_init = [0, 0.01, 0.15, 0.8];
% fn_trans_param = @(xx,mm) transform_param_garch(xx, mm);
% fn_jacobian = @(xx) jacobian_garch(xx);
% kernel_init = @(xx) -posterior_garch11(fn_trans_param(xx,'back'),y)/T;
% options = optimset('display','off','TolFun',1e-5,'LargeScale','off','TolX',1e-5,'HessUpdate','bfgs','FinDiffType','central',...
%      'maxiter',5000,'MaxFunEvals',5000);
% [mu, ~, Hessian] = estimate(kernel_init,mu_init,fn_trans_param,fn_jacobian,options);
% Sigma = inv(T*Hessian);


% theta  = [mu, omega, alpha, beta]
mu_init = [0, 0.1, 0.05, 0.85];
y_S = var(y);
kernel_init = @(xx) -posterior_garch11(xx,y,y_S)/T;
kernel = @(xx) posterior_garch11(xx,y, y_S);
[mu,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init);
Sigma = inv(T*Sigma);
% df = 5;
% draw = rmvt(mu,Sigma,df,M+BurnIn);
% lnk = kernel(draw);
% lnd = dmvgt_mex(draw, mu, Sigma, df, 1, GamMat, double(1));
[mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
[draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
lnd = dmvgt(draw, mit, true, GamMat); 
    
lnw = lnk - lnd;
lnw = lnw - max(lnw);
[ind, a] = fn_MH(lnw);
draw = draw(ind,:);
lnw = lnw(ind);
accept = a/(M+BurnIn);
draw = draw(BurnIn+1:BurnIn+M,:);    
lnw = lnw(BurnIn+1:BurnIn+M,:);    
 
h_post = volatility_garch11(draw,y,y_S);
y_post = draw(:,1) + sqrt(h_post).*randn(M,1);
y_post = sort(y_post);
VaR_1_post = y_post(p_bar1*M); 
VaR_5_post = y_post(p_bar*M); 

%% Threshold = 10% perscentile of the data sample
threshold = sort(y);
threshold = threshold(2*p_bar*T);
%% CENSORED
kernel_init = @(xx) - C_posterior_garch11(xx, y, threshold, y_S)/T;    
kernel = @(xx) C_posterior_garch11(xx, y, threshold, y_S);
kernel_mex = @(xx) C_posterior_garch11_mex(xx, y, threshold, y_S);

lnk_C = kernel(draw_C(1:100,:));
h_C = volatility_garch11(draw_C(1:100,:),y,y_S);
[lnk_mex,~,h_mex] = kernel_mex(draw_C(1:100,:));


try
    [mit_C, CV_C] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);   
catch
    [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options)
%     [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_C,options)
    Sigma_C = inv(T*Sigma_C);
    draw_C = rmvt(mu_C,Sigma_C,df,M+BurnIn);
    lnk_C = kernel(draw_C);
    lnd_C = dmvgt_mex(draw_C, mu_C, Sigma_C, df, 1, GamMat, double(1));

    mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
    CV = cont.mit.CV_old;
end
[draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
lnd_C = dmvgt(draw_C, mit_C, true, GamMat);     
lnw_C = lnk_C - lnd_C;
lnw_C = lnw_C - max(lnw_C);
[ind, a] = fn_MH(lnw_C);
draw_C = draw_C(ind,:);
accept_C(s,1) = a/(M+BurnIn);
draw_C = draw_C(BurnIn+1:BurnIn+M,:);

h_post_C = volatility_garch11(draw_C,y,y_S);
y_post_C = draw_C(:,1) + sqrt(h_post_C).*randn(M,1);
y_post_C = sort(y_post_C);
VaR_1_post_C = y_post_C(p_bar1*M); 
VaR_5_post_C = y_post+C(p_bar*M); 




%% Threshold = 0
threshold0 = 0;
%% CENSORED
kernel_init = @(xx) - C_posterior_garch11(xx, y, threshold0, y_S)/T;    
kernel = @(xx) C_posterior_garch11(xx, y, threshold0, y_S);
% mu_init = [-0.5, 0.01, 0.05, 0.7];
try
    [mit_C0, CV_C0] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
catch
    [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu_init,options);
    Sigma_C0 = inv(T*Sigma_C0);
    mit_C0 = struct('mu',mu_C0,'Sigma',reshape(Sigma_C0,1,length(mu_C0)^2),'df', df, 'p', 1);
    CV_0 = cont.mit.CV_old;
end
[draw_C0, lnk_C0] = fn_rmvgt_robust(M+BurnIn, mit_C0, kernel, false);
lnd_C0 = dmvgt(draw_C0, mit_C0, true, GamMat);    
lnw_C0 = lnk_C0 - lnd_C0;
lnw_C0 = lnw_C0 - max(lnw_C0);
[ind, a] = fn_MH(lnw_C0);
draw_C0 = draw_C0(ind,:);
accept_C0(s,1) = a/(M+BurnIn);
draw_C0 = draw_C0(BurnIn+1:BurnIn+M,:);

h_post_C0 = volatility_garch11(draw_C0,y,y_S);
y_post_C0 = draw_C0(:,1) + sqrt(h_post_C0).*randn(M,1);
VaR_1_post_C0(s,1) = y_post_C0(p_bar1*M); 
VaR_5_post_C0(s,1) = y_post_C0(p_bar*M); 