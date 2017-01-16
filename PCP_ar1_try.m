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
% Metropolis-Hastings for the parameters
M = 10000; % number of draws 
BurnIn = 1000;

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
lnw = lnw(ind);
accept = a/(M+BurnIn);
draw = draw(BurnIn+1:BurnIn+M,:);    
lnw = lnw(BurnIn+1:BurnIn+M,:);    
 
y_post = draw(:,1) + draw(:,3).*y(T,1) + draw(:,2).*randn(M,1);
y_post = sort(y_post);
VaR_1_post = y_post(p_bar1*M); 
VaR_5_post = y_post(p_bar*M); 

%% PARTIAL CENSORING: keep rho uncensored, then censor mu i sigma

 
partition = 3;
% conditional candidate 
% joint cnadidate for the joint censored posterior
%     kernel = @(a) posterior_arch(a, data, S, true);
%     [mit1, summary1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
% kernel_init = @(xx) - C_posterior_ar1(xx, y, threshold);
% kernel = @(xx) C_posterior_ar1(xx, y, threshold);
kernel_init = @(xx) - C_posterior_ar1_mex(xx, y, threshold);
kernel = @(xx) C_posterior_ar1_mex(xx, y, threshold);
mu_init = [0, 1, 0.5];
cont = MitISEM_Control;
[mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);

param_true = [c,sigma1,sigma2,rho];
if save_on
    save(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP.mat'],...
    'y','draw','accept','param_true',...
    'mit','CV','cont',...
    'VaR_1','VaR_1_post',...
    'VaR_5','VaR_5_post');
end

load(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP.mat'])

% cont.mit.iter_max = 0;
% [mit1, CV1] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);


% profile on
draw_partial = sim_cond_mit(mit, draw, partition, GamMat);
% profile off
% profile viewer
% % % TO DO
% % % for each rho run MH!!!

%% Highest weights rhps
draw_org = draw;
[~,ind_max] = sort(lnw);
draw_max = draw(ind_max,:);
II = 100;
draw_max = draw_max(M-II+1:M,:);

[draw_partial, accept_partial] = sim_cond_mit_MH(mit, draw_max, partition, M, kernel, GamMat);
%% Short version
draw_short = draw(1:10,:);
[draw_partial2, accept_partial2] = sim_cond_mit_MH(mit, draw_short, partition, 10*M, kernel, GamMat);
%% True rho
[draw_partial3, accept_partial3] = sim_cond_mit_MH(mit, [0,0,0.8], partition, M, kernel, GamMat);



y_post_partial = draw_partial(:,1) + draw_partial(:,3).*y(T,1) + draw_partial(:,2).*randn(100*M,1);
y_post_partial = sort(y_post_partial);
VaR_1_post_partial = y_post_partial(p_bar1*100*M); 
VaR_5_post_partial = y_post_partial(p_bar*100*M); 

if plot_on
    subplot(1,3,1)
    hold on
    fn_hist(draw(:,1))
    fn_hist(draw_partial(:,1)) 
    YL = get(gca,'YLim');
    line([c c], YL,'Color','r','LineWidth',3); 
    hold off
    xlabel('\mu')
    
    subplot(1,3,2)
    hold on
    fn_hist(draw(:,2))
    fn_hist(draw_partial(:,2)) 
    YL = get(gca,'YLim');
    line([sigma2 sigma2], YL,'Color','r','LineWidth',3); 
    hold off   
    xlabel('\sigma')
    
    subplot(1,3,3)    
    hold on
    fn_hist(draw(:,3))
    fn_hist(draw_partial(:,3)) 
    YL = get(gca,'YLim');
    line([rho rho], YL,'Color','r','LineWidth',3);     
    hold off 
    xlabel('\rho') 
end

kernel = @(xx) C_posterior_ar1(xx, y, threshold);
lnk_partial= kernel(draw_partial);
lnd_partial = dmvgt_mex(draw_partial, mu_C, Sigma_C, df, 1, GamMat, double(1));
lnw_partial = lnk_partial - lnd_partial;
lnw_partial = lnw_partial - max(lnw_partial);
[ind, a] = fn_MH(lnw_partial);
draw_partial = draw_partial(ind,:);
accept_partial = a/lenght(draw_partial);
% draw_partial = draw_partial(BurnIn+1:BurnIn+M,:);


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
 


 