clear all
close all

addpath(genpath('include/'));

% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s); 

model = 'ar1';
parameters = {'$\\mu$','$\\sigma$','$\\phi$'};

sigma1 = 1;
sigma2 = 2;
c = (sigma2 - sigma1)/sqrt(2*pi);

T = 1000; % time series length
p_bar1 = 0.01;
p_bar = 0.05;
M = 10000; % number of draws 

x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);

v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014a)')
    fn_hist = @(xx) hist(xx,20);
else
    fn_hist = @(xx) histogram(xx,20);
end

plot_on = false;
save_on = true;

%% simple AR(1)
eps = randn(T,1);
ind = (eps>0);
eps(ind) = c + sigma1.*eps(ind);
eps(~ind) = c + sigma2.*eps(~ind);

rho = 0.8;
y = zeros(T,1);

y(1,1) = eps(1,1);
for ii = 2:T
    y(ii,1) = rho*y(ii-1,1) + eps(ii,1);
end
threshold = sort(y);
threshold = threshold(2*p_bar*T);

% MC VaRs under the true model
eps_sort = randn(M,1);
ind = (eps_sort>0);
eps_sort(ind) = c + sigma1.*eps_sort(ind);
eps_sort(~ind) = c + sigma2.*eps_sort(~ind);

y_sort = rho*y(T,1) + eps_sort;
y_sort = sort(y_sort);
VaR_1 = y_sort(p_bar1*M); 
VaR_5 = y_sort(p_bar*M); 


%% Uncensored Posterior
% Misspecified model: AR1 normal with unknown mu and sigma
% Metropolis-Hastings for the parameters
BurnIn = 1000;

% Uncensored likelihood
kernel_init = @(xx) -loglik_ar1(xx,y);
[mu,~,~,~,~,Sigma] = fminunc(kernel_init,[0,1,0.9]);
Sigma = inv(T*Sigma);
df = 5;
draw = rmvt(mu,Sigma,df,M+BurnIn);
kernel = @(xx) posterior_ar1(xx,y);
lnk = kernel(draw);

lnd = dmvgt_mex(draw, mu, Sigma, df, 1, GamMat, double(1));
lnw = lnk - lnd;
lnw = lnw - max(lnw);
[ind, a] = fn_MH(lnw);
draw = draw(ind,:);
accept = a/(M+BurnIn);
draw = draw(BurnIn+1:BurnIn+M,:);    
 
y_post = draw(:,1) + draw(:,3).*y(T,1) + draw(:,2).*randn(M,1);
y_post = sort(y_post);
VaR_1_post = y_post(p_bar1*M); 
VaR_5_post = y_post(p_bar*M); 

%% PARTIAL CENSORING: keep rho uncensored, then censor mu i sigma

draw_mix = draw;
partition = 3;
% conditional candidate 
% joint cnadidate for the joint censored posterior
%     kernel = @(a) posterior_arch(a, data, S, true);
%     [mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
kernel_init = @(xx) - C_posterior_ar1(xx, y, threshold);
kernel = @(xx) C_posterior_ar1(xx, y, threshold);
mu_init = [0, 1, 0.5];
cont = MitISEM_Control;
[mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);

cont.mit.iter_max = 0;
[mit1, CV1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);



% 1. Threshold = 10% perscentile of the data sample
% profile on
kernel_init = @(xx) - C_posterior_ar1(xx, y, threshold)/T;
[mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu);
% profile off
% profile viewer
Sigma_C = inv(T*Sigma_C);
draw_C = rmvt(mu_C,Sigma_C,df,M+BurnIn);
kernel = @(ss) C_posterior_ar1(ss, y, threshold);
lnk_C = kernel(draw_C);
lnd_C = dmvgt_mex(draw_C, mu_C, Sigma_C, df, 1, GamMat, double(1));
lnw_C = lnk_C - lnd_C;
lnw_C = lnw_C - max(lnw_C);
[ind, a] = fn_MH(lnw_C);
draw_C = draw_C(ind,:);
accept_C = a/(M+BurnIn);
draw_C = draw_C(BurnIn+1:BurnIn+M,:);
 
y_post_C = draw_C(:,1) + draw_C(:,3).*y(T,1) + draw_C(:,2).*randn(M,1);
y_post_C = sort(y_post_C);
VaR_1_post_C = y_post_C(p_bar1*M); 
VaR_5_post_C = y_post_C(p_bar*M); 

% 2. Threshold = 0             
threshold0 = 0;
kernel_init = @(xx) - C_posterior_ar1(xx, y, threshold0)/T;
[mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu_C);
Sigma_C0 = inv(T*Sigma_C0);
draw_C0 = rmvt(mu_C0,Sigma_C0,df,M+BurnIn);
kernel = @(ss) C_posterior_ar1(ss, y, threshold0);
lnk_C0 = kernel(draw_C0);
lnd_C0 = dmvgt_mex(draw_C0, mu_C0, Sigma_C0, df, 1, GamMat, double(1));
lnw_C0 = lnk_C0 - lnd_C0;
lnw_C0 = lnw_C0 - max(lnw_C0);
[ind, a] = fn_MH(lnw_C0);
draw_C0 = draw_C0(ind,:);
accept_C0 = a/(M+BurnIn);
draw_C0 = draw_C0(BurnIn+1:BurnIn+M,:);
 
y_post_C0 = draw_C0(:,1) + draw_C0(:,3).*y(T,1) + draw_C0(:,2).*randn(M,1);
y_post_C0 = sort(y_post_C0);
VaR_1_post_C0 = y_post_C0(p_bar1*M); 
VaR_5_post_C0 = y_post_C0(p_bar*M); 

param_true = [c,sigma1,sigma2,rho];
if save_on
    save(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'.mat'],...
    'y','draw','draw_C','draw_C0','param_true',...
    'accept','accept_C','accept_C0',...
    'VaR_1','VaR_1_post','VaR_1_post_C','VaR_1_post_C0',...
    'VaR_5','VaR_5_post','VaR_5_post_C','VaR_5_post_C0');
end
 