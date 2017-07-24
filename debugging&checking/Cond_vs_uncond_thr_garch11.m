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

model = 'garch11';
fprintf('Model: %s.\n',model)
parameters = {'$\\mu$','$\\sigma$','$\\phi$'};

sigma1 = 1;
omega = 1;
alpha = 0.1;
beta = 0.8;

mu_init = [0, 1, 0.05, 0.85];
S = 1;
H = 0;

MU1 = zeros(3,8,S);
MU2 = zeros(3,8,S);

for s = 1:S
    %% SIGMA2 = 1
    sigma2 = 1;    
    c = (sigma2 - sigma1)/sqrt(2*pi);
    kappa = 0.5*(sigma1^2 + sigma2^2 - ((sigma2-sigma1)^2)/pi); % var of eps
    sigma1_k = sigma1/sqrt(kappa);
    sigma2_k = sigma2/sqrt(kappa);
    
    %% T = 1000
    T = 1000;
      
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

    fprintf('*** Censored Posterior, time varying threshold ***\n');
%     threshold = 1;
    kernel_init = @(xx) - C_posterior_garch11_varc( transform_param_garch(xx,'opt'), y(1:T))/T; 
    muuu = fminunc(kernel_init,transform_param_garch(mu_init,'back'),options);
    MU1(1,1:4,s) = transform_param_garch(muuu,'opt');

    fprintf('*** Censored Posterior, threshold = 0 ***\n');
    threshold0 = 0;
    kernel_init = @(xx) - C_posterior_garch11_mex(transform_param_garch(xx,'opt'), y(1:T), threshold0)/T; 
    muuu = fminunc(kernel_init,transform_param_garch(mu_init,'back'),options);
    MU1(1,5:8,s) = transform_param_garch(muuu,'opt');

 
    %% SIGMA2 = 2
    sigma2 = 2;    
    c = (sigma2 - sigma1)/sqrt(2*pi);
    kappa = 0.5*(sigma1^2 + sigma2^2 - ((sigma2-sigma1)^2)/pi); % var of eps
    sigma1_k = sigma1/sqrt(kappa);
    sigma2_k = sigma2/sqrt(kappa);
    %% T = 1000
    T = 1000;
      
    eps = randn(T+H,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);
    eps = eps/sqrt(kappa);  
    y = zeros(T+H,1);
    y_S = var(y(1:T));
    h_true = zeros(T+H,1);

    for ii = 1:T+H
        if (ii == 1)
            h_true(ii,1) = omega;
        else
            h_true(ii,1) = omega*(1-alpha-beta) + alpha*(y(ii-1,1))^2 + beta*h_true(ii-1,1);
        end
        y(ii,1) = sqrt(h_true(ii,1))*eps(ii,1);
    end

    fprintf('*** Censored Posterior, time varying threshold ***\n');
%     threshold = 1;
    kernel_init = @(xx) - C_posterior_garch11_varc( transform_param_garch(xx,'opt'), y(1:T))/T; 
    muuu = fminunc(kernel_init,transform_param_garch(mu_init,'back'),options);
    MU2(1,1:4) = transform_param_garch(muuu,'opt');

    kernel_init = @(xx) - C_posterior_garch11_varc(xx, y(1:T))/T; 
    muuu = fminunc(kernel_init,mu_init,options);
   
    
    fprintf('*** Censored Posterior, threshold = 0 ***\n');
    threshold0 = 0;
    kernel_init = @(xx) - C_posterior_garch11_mex(transform_param_garch(xx,'opt'), y(1:T), threshold0)/T; 
    muuu = fminunc(kernel_init,transform_param_garch(mu_init,'back'),options);
    MU2(1,5:8,s) = transform_param_garch(muuu,'opt');

end
 

MU1_true = repmat([0,1,0.8],3,2);
MU2_true = repmat([c,2,0.8],3,2);

MSE1 = bsxfun(@minus,MU1,MU1_true);
MSE1 = mean(MSE1.^2,3);

MSE2 = bsxfun(@minus,MU2,MU2_true);
MSE2 = mean(MSE2.^2,3);