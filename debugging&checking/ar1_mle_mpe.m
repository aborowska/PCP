% clear all
model = 'ar1';
fprintf('Model: %s.\n',model)
parameters = {'$\\mu$','$\\sigma$','$\\phi$'};
partition = 3;

sigma1 = 1;
sigma2 = 2;    
c = (sigma2 - sigma1)/sqrt(2*pi);
rho = 0.8;
param_true = [c, sigma2, rho];
mu_init = [0,1,0.9];

T = 100;
H = 0;

% quantiles of interest
p_bar1 = 0.01;
p_bar = 0.05;

% Metropolis-Hastings for the parameters
M = 10000; % number of draws 
BurnIn = 1000;

options = optimset('Display','off');
% w = warning('query','last');
% id = w.identifier;
id = 'optim:fminunc:SwitchingMethod';
warning('off',id);
  
%%    
S = 100;

MU_MLE = zeros(S,3);
MU_C = zeros(S,3);
MU_C0 = zeros(S,3);
MU_Cah = zeros(S,3);
MU_Cm = zeros(S,3);


for ss = 1:S
    eps = randn(T+H,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);

    y = zeros(T+H,1);
    y(1,1) = eps(1,1);
    for ii = 2:T+H
        y(ii,1) = rho*y(ii-1,1) + eps(ii,1);
    end

    kernel_init = @(xx) -posterior_ar1_mex(xx,y(1:T))/T;
    [mu_mle,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init, options);

    threshold = sort(y(1:T));
    threshold = threshold(2*p_bar*T);   
    kernel_init = @(xx) - C_posterior_ar1_mex(xx, y(1:T), threshold)/T;
    [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);

    threshold0 = 0;
    kernel_init = @(xx) - C_posterior_ar1_mex(xx, y(1:T), threshold0)/T; 
    [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu_init,options);

    threshold = 1.0;
    kernel_init = @(xx) - C_posterior_ar1_varc_noparam_mex(xx, y(1:T,1), threshold)/T;    
    [mu_Cah,~,~,~,~,Sigma_Cah] = fminunc(kernel_init,mu_init,options);

    threshold = 0.1;
    quantile = norminv(threshold);
    kernel_init = @(xx) - C_posterior_ar1_varc_mle_mex(xx, y(1:T), mu_mle, quantile)/T; 
    [mu_Cm,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_init,options);

    MU_MLE(ss,:) = mu_mle;
    MU_C(ss,:) = mu_C;
    MU_C0(ss,:) = mu_C0;
    MU_Cah(ss,:) = mu_Cah;
    MU_Cm(ss,:) = mu_Cm;
end

display(mean(MU_MLE))
display(mean(MU_C))
display(mean(MU_C0))
display(mean(MU_Cah))
display(mean(MU_Cm))