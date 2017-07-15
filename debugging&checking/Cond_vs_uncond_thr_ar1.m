% UNCONDITIONAL VS CONDITIONAL THRESHOLD
clear all
close all

addpath(genpath('include/'));
options = optimset('Display','off');
% w = warning('query','last');
% id = w.identifier;
id = 'optim:fminunc:SwitchingMethod';
warning('off',id);
% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s); 

model = 'ar1';
fprintf('Model: %s.\n',model)
parameters = {'$\\mu$','$\\sigma$','$\\phi$'};

sigma1 = 1;
rho = 0.8;
mu_init = [0,1,0.9];
S = 100;
H = 0;

MU1 = zeros(3,6,S);
MU2 = zeros(3,6,S);

for s = 1:S
    %% SIGMA2 = 1
    sigma2 = 1;    
    c = (sigma2 - sigma1)/sqrt(2*pi);

    %% T = 100
    T = 100;
    eps = randn(T+H,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);

    y = zeros(T+H,1);
    y(1,1) = eps(1,1);
    for ii = 2:T+H
        y(ii,1) = rho*y(ii-1,1) + eps(ii,1);
    end  

    fprintf('*** Censored Posterior, time varying threshold, THR = 1 ***\n');
    threshold = 1;
    kernel_init = @(xx) - C_posterior_ar1_varc_mex( transform_param_ar(xx,'back'), y(1:T), threshold)/T; 
    muuu = fminunc(kernel_init,transform_param_ar(mu_init,'opt'),options);
    MU1(1,1:3,s) = transform_param_ar(muuu,'back');

    fprintf('*** Censored Posterior, threshold = 0 ***\n');
    threshold0 = 0;
    kernel_init = @(xx) - C_posterior_ar1_mex(transform_param_ar(xx,'back'), y(1:T), threshold0)/T; 
    muuu = fminunc(kernel_init,transform_param_ar(mu_init,'opt'),options);
    MU1(1,4:6,s) = transform_param_ar(muuu,'back');

    %% T = 1000
    T = 1000;
    eps = randn(T+H,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);

    y = zeros(T+H,1);
    y(1,1) = eps(1,1);
    for ii = 2:T+H
        y(ii,1) = rho*y(ii-1,1) + eps(ii,1);
    end 

    fprintf('*** Censored Posterior, time varying threshold, THR = 1 ***\n');
    threshold = 1;
    kernel_init = @(xx) - C_posterior_ar1_varc_mex( transform_param_ar(xx,'back'), y(1:T), threshold)/T; 
    muuu = fminunc(kernel_init,transform_param_ar(mu_init,'opt'),options);
    MU1(2,1:3,s) = transform_param_ar(muuu,'back');

    fprintf('*** Censored Posterior, threshold = 0 ***\n');
    threshold0 = 0;
    kernel_init = @(xx) - C_posterior_ar1_mex(transform_param_ar(xx,'back'), y(1:T), threshold0)/T; 
    muuu = fminunc(kernel_init,transform_param_ar(mu_init,'opt'),options);
    MU1(2,4:6,s) = transform_param_ar(muuu,'back');

    %% T = 10000
    T = 10000;
    eps = randn(T+H,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);

    y = zeros(T+H,1);
    y(1,1) = eps(1,1);
    for ii = 2:T+H
        y(ii,1) = rho*y(ii-1,1) + eps(ii,1);
    end  

    fprintf('*** Censored Posterior, time varying threshold, THR = 1 ***\n');
    threshold = 1;
    kernel_init = @(xx) - C_posterior_ar1_varc_mex( transform_param_ar(xx,'back'), y(1:T), threshold)/T; 
    muuu = fminunc(kernel_init,transform_param_ar(mu_init,'opt'),options);
    MU1(3,1:3,s) = transform_param_ar(muuu,'back');

    fprintf('*** Censored Posterior, threshold = 0 ***\n');
    threshold0 = 0;
    kernel_init = @(xx) - C_posterior_ar1_mex(transform_param_ar(xx,'back'), y(1:T), threshold0)/T; 
    muuu = fminunc(kernel_init,transform_param_ar(mu_init,'opt'),options);
    MU1(3,4:6,s) = transform_param_ar(muuu,'back');

    %% SIGMA2 = 2
    sigma2 = 2;    
    c = (sigma2 - sigma1)/sqrt(2*pi);

    %% T = 100
    T = 100;
    eps = randn(T+H,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);

    y = zeros(T+H,1);
    y(1,1) = eps(1,1);
    for ii = 2:T+H
        y(ii,1) = rho*y(ii-1,1) + eps(ii,1);
    end  

    fprintf('*** Censored Posterior, time varying threshold, THR = 1 ***\n');
    threshold = 1;
    kernel_init = @(xx) - C_posterior_ar1_varc_mex( transform_param_ar(xx,'back'), y(1:T), threshold)/T; 
    muuu = fminunc(kernel_init,transform_param_ar(mu_init,'opt'),options);
    MU2(1,1:3,s) = transform_param_ar(muuu,'back');

    fprintf('*** Censored Posterior, threshold = 0 ***\n');
    threshold0 = 0;
    kernel_init = @(xx) - C_posterior_ar1_mex(transform_param_ar(xx,'back'), y(1:T), threshold0)/T; 
    muuu = fminunc(kernel_init,transform_param_ar(mu_init,'opt'),options);
    MU2(1,4:6,s) = transform_param_ar(muuu,'back');

    %% T = 1000
    T = 1000;
    eps = randn(T+H,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);

    y = zeros(T+H,1);
    y(1,1) = eps(1,1);
    for ii = 2:T+H
        y(ii,1) = rho*y(ii-1,1) + eps(ii,1);
    end  

    fprintf('*** Censored Posterior, time varying threshold, THR = 1 ***\n');
    threshold = 1;
    kernel_init = @(xx) - C_posterior_ar1_varc_mex( transform_param_ar(xx,'back'), y(1:T), threshold)/T; 
    muuu = fminunc(kernel_init,transform_param_ar(mu_init,'opt'),options);
    MU2(2,1:3,s) = transform_param_ar(muuu,'back');

    fprintf('*** Censored Posterior, threshold = 0 ***\n');
    threshold0 = 0;
    kernel_init = @(xx) - C_posterior_ar1_mex(transform_param_ar(xx,'back'), y(1:T), threshold0)/T; 
    muuu = fminunc(kernel_init,transform_param_ar(mu_init,'opt'),options);
    MU2(2,4:6,s) = transform_param_ar(muuu,'back');

    %% T = 10000
    T = 10000;
    eps = randn(T+H,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);

    y = zeros(T+H,1);
    y(1,1) = eps(1,1);
    for ii = 2:T+H
        y(ii,1) = rho*y(ii-1,1) + eps(ii,1);
    end  

    fprintf('*** Censored Posterior, time varying threshold, THR = 1 ***\n');
    threshold = 1;
    kernel_init = @(xx) - C_posterior_ar1_varc_mex( transform_param_ar(xx,'back'), y(1:T), threshold)/T; 
    muuu = fminunc(kernel_init,transform_param_ar(mu_init,'opt'),options);
    MU2(3,1:3,s) = transform_param_ar(muuu,'back');

    fprintf('*** Censored Posterior, threshold = 0 ***\n');
    threshold0 = 0;
    kernel_init = @(xx) - C_posterior_ar1_mex(transform_param_ar(xx,'back'), y(1:T), threshold0)/T; 
    muuu = fminunc(kernel_init,transform_param_ar(mu_init,'opt'),options);
    MU2(3,4:6,s) = transform_param_ar(muuu,'back');
end
 

MU1_true = repmat([0,1,0.8],3,2);
MU2_true = repmat([c,2,0.8],3,2);

MSE1 = bsxfun(@minus,MU1,MU1_true);
MSE1 = mean(MSE1.^2,3);

MSE2 = bsxfun(@minus,MU2,MU2_true);
MSE2 = mean(MSE2.^2,3);