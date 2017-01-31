clear all
close all

addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

model = 'ar1';
fprintf('Model: %s.\n',model)
parameters = {'$\\mu$','$\\sigma$','$\\phi$'};

sigma1 = 1;
sigma2 = 1;
c = (sigma2 - sigma1)/sqrt(2*pi);

S = 100; % number of MC replications

% theoretical quantiles
q1 = zeros(S,1);
q5 = zeros(S,1);

% simulated quantiles: 
% true model, posterior, censored posterior 10%, censored posterior at 0
VaR_1 = zeros(S,1);
VaR_1_post = zeros(S,1);
VaR_1_post_C = zeros(S,1);
VaR_1_post_C0 = zeros(S,1);

VaR_5 = zeros(S,1);
VaR_5_post = zeros(S,1);
VaR_5_post_C = zeros(S,1);
VaR_5_post_C0 = zeros(S,1);

T = 10000; %time series length
p_bar1 = 0.01;
p_bar = 0.05;
M = 10000; % number of draws 

x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);

% various display options
v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014a)')
    fn_hist = @(xx) hist(xx,20);
else
    fn_hist = @(xx) histogram(xx,20);
end

plot_on = false;
save_on = true;

options = optimset('Display','off');
% w = warning('query','last');
% id = w.identifier;
id = 'optim:fminunc:SwitchingMethod';
warning('off',id);

for s = 1:S
%     if (mod(s,10)==0)
        fprintf(['\n',model, ' simualtion no. %i\n'],s)
%     end
    
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

    % true VaRs
    q1(s,1) = norminv(p_bar1,c+rho*y(T,1),sigma2);
    q5(s,1) = norminv(p_bar,c+rho*y(T,1),sigma2); 
    
    threshold = sort(y);
    threshold = threshold(2*p_bar*T);
    
    % MC VaRs under the true model
    eps_sort = randn(M,1);
    ind = (eps_sort>0);
    eps_sort(ind) = c + sigma1.*eps_sort(ind);
    eps_sort(~ind) = c + sigma2.*eps_sort(~ind);

    y_sort = rho*y(T,1) + eps_sort;
    y_sort = sort(y_sort);
    VaR_1(s,1) = y_sort(p_bar1*M);
    VaR_5(s,1) = y_sort(p_bar*M); 

    %% Uncensored Posterior
    % Misspecified model: AR1 normal with unknown mu and sigma
    % Metropolis-Hastings for the parameters
    BurnIn = 1000;

    % Uncensored likelihood
    kernel_init = @(xx) -loglik_ar1(xx,y);
    [mu,~,~,~,~,Sigma] = fminunc(kernel_init,[0,1,0.9],options);
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
    VaR_1_post(s,1) = y_post(p_bar1*M); 
    VaR_5_post(s,1) = y_post(p_bar*M); 

    % 1. Threshold = 10% perscentile of the data sample
    kernel_init = @(xx) - C_posterior_ar1(xx, y, threshold)/T;
    [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu,options);
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
    VaR_1_post_C(s,1) = y_post_C(p_bar1*M); 
    VaR_5_post_C(s,1) = y_post_C(p_bar*M); 

    % 2. Threshold = 0             
    threshold0 = 0;
    kernel_init = @(xx) - C_posterior_ar1(xx, y, threshold0)/T;
    [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu_C,options);
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
    VaR_1_post_C0(s,1) = y_post_C0(p_bar1*M); 
    VaR_5_post_C0(s,1) = y_post_C0(p_bar*M); 
end

% MSEs
MSE_1 = mean((VaR_1 - q1).^2);
MSE_1_post = mean((VaR_1_post - q1).^2);
MSE_1_post_C = mean((VaR_1_post_C - q1).^2);
MSE_1_post_C0 = mean((VaR_1_post_C0 - q1).^2);

MSE_5 = mean((VaR_5 - q5).^2);
MSE_5_post = mean((VaR_5_post - q5).^2);
MSE_5_post_C = mean((VaR_5_post_C - q5).^2);
MSE_5_post_C0 = mean((VaR_5_post_C0 - q5).^2);

if plot_on
    figure(5)
    hold on
    plot(VaR_5,'k')
    plot(VaR_5_post,'b')
    plot(VaR_5_post_C,'c')
    plot(VaR_5_post_C0,'g')
    plot(q5,'r')
    hold off
end


param_true = [c,sigma2,rho];
if save_on
    save(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_MC.mat'],...
    'y','draw','draw_C','draw_C0','param_true','q1','q5',...
    'accept','accept_C','accept_C0',...
    'VaR_1','VaR_1_post','VaR_1_post_C','VaR_1_post_C0',...
    'VaR_5','VaR_5_post','VaR_5_post_C','VaR_5_post_C0',...
    'MSE_1','MSE_1_post','MSE_1_post_C','MSE_1_post_C0',...
    'MSE_5','MSE_5_post','MSE_5_post_C','MSE_5_post_C0')
end

% print_table_cp_mc(model,parameters,sigma1,sigma2)