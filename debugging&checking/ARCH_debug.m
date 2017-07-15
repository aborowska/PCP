close all
clear all
sdd=1
T = 1000; S = 1; II = 10; sigma2 = 2 ;
 

    addpath(genpath('include/'));

%     s = RandStream('mt19937ar','Seed',1);
%     RandStream.setGlobalStream(s); 

    model = 'arch1';    
    partition = 3;
    fprintf('Model: %s.\n',model)
    parameters = {'$\\mu$','$\\omega$','$\\mu2$','$\\alpha$'};
    fprintf('Time series length T = %d.\n',T);
    
    sigma1 = 1;
    % sigma2 = 2;
    c = (sigma2 - sigma1)/sqrt(2*pi); % mean of eps
    kappa = 0.5*(sigma1^2 + sigma2^2 - ((sigma2-sigma1)^2)/pi); % var of eps
    sigma1_k = sigma1/sqrt(kappa);
    sigma2_k = sigma2/sqrt(kappa);
    
    mu2 = 0; % gama = mu - mu2;
    omega = 1;
    alpha = 0.1;

    % theta  = [mu, omega, mu2, alpha]
    mu_true = [0, omega, 0, alpha];
    param_true = [c,sigma2,omega,mu2,alpha];
    mu_init = [0, 1, 0.1, 0.05];

    % S = 20; % number of MC replications
    H = 100;

    varc = true; % run the version with time varying threshold

    % quantiles of interest
    p_bar1 = 0.01;
    p_bar = 0.05;
    % theoretical quantiles
    q1 = zeros(S,H);
    q5 = zeros(S,H);

    % Metropolis-Hastings for the parameters
    M = 10000; % number of draws 
    BurnIn = 1000;

    x_gam = (0:0.00001:50)'+0.00001;
    GamMat = gamma(x_gam);

    cont = MitISEM_Control;
    cont.mit.CV_max = 1; %2?
    cont.mit.iter_max = 10;
    cont.mit.Hmax = 6;
    cont.mit.dfnc = 5;
    df = 5; % default df for a mit

    %% various display options
    cont.disp = false; %true; %false;

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
    
    
       
    s = RandStream('mt19937ar','Seed',sdd);
    RandStream.setGlobalStream(s); 
    
    sigma1_k = sigma1/sqrt(kappa);
    sigma2_k = sigma2/sqrt(kappa);

    %% ARCH(1,1)
    eps = randn(T+H,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);
    eps = eps/sqrt(kappa);  
    y = zeros(T+H,1);
    h_true = zeros(T+H,1);
    h_true(1,1) = omega;
    y(1,1) = sqrt(h_true(1,1)).*eps(1,1);
    for ii = 2:T+H
        h_true(ii,1) = omega*(1-alpha) + alpha*(y(ii-1,1)).^2;
        y(ii,1) = sqrt(h_true(ii,1)).*eps(ii,1);
    end
    % true VaRs
    q1 = norminv(p_bar1, c, sigma2_k*sqrt(h_true(T+1:T+H)))';
    q5 = norminv(p_bar, c, sigma2_k*sqrt(h_true(T+1:T+H)))'; 

    % MC VaRs under the true model
    eps_sort = randn(M,H);
    ind = (eps_sort>0);
    eps_sort(ind) = c + sigma1.*eps_sort(ind);
    eps_sort(~ind) = c + sigma2.*eps_sort(~ind);
    eps_sort = eps_sort/sqrt(kappa); 

    y_sort = bsxfun(@times,eps_sort,sqrt(h_true(T+1:T+H,1))');
    y_sort = sort(y_sort);
    VaR_1 = y_sort(p_bar1*M,:); 
    VaR_5 = y_sort(p_bar*M,:); 

    %% Uncensored Posterior
    fprintf('*** Uncensored Posterior ***\n');
    y_S = var(y(1:T));
  
    kernel_init = @(xx) -posterior_arch1_mex(xx, y(1:T), y_S)/T;
    kernel = @(xx) posterior_arch1_mex(xx, y(1:T), y_S);
    [mu_MLE,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init,options);
    Sigma = inv(T*Sigma);
    
    
    
    threshold = 1.0;
    fprintf('*** Censored Posterior, ad hoc time varying threshold, THR = %3.2f ***\n',threshold);
    kernel_init = @(xx) - C_posterior_arch1_varc_noparam_mex(xx, y(1:T,1), threshold, y_S)/T;    
    kernel = @(xx) C_posterior_arch1_varc_noparam_mex(xx, y(1:T,1), threshold, y_S);
    [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
    Sigma_C = inv(T*Sigma_C);
    
    
    threshold = 0.1; %<---------- HiGhER?
    quantile = norminv(threshold);
    fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold);
    kernel_init = @(xx) - C_posterior_arch1_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S)/T;    
    kernel = @(xx) C_posterior_arch1_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S);
    
    [mu_Cm,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_init,options);
    [mu_Cm2,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_C,options);
    
    [~,~,~,THR_mex] = C_posterior_arch1_varc_mle_mex(mu_Cm, y(1:T,1), mu_MLE, quantile, y_S);
    
    
cond = zeros(T,1);
h_MLE = zeros(T,1);
THR = zeros(T,1);
for jj=1:T
    if (jj==1)
        h_MLE(jj) = y_S;            
    else
        mu = y(jj-1) - mu_MLE(3);
        h_MLE(jj) = mu_MLE(2)*(1.0-mu_MLE(4)) + mu_MLE(4)*(mu*mu) ;                  
    end       
    THR(jj) = (y(jj) - mu_MLE(1))/sqrt(h_MLE(jj));
    if (THR(jj) >= quantile) % OTHER WAY ROUND AS IN MEX:COND IS WHEN TO CENSOR
        cond(jj) = 1;
    end
    THR(jj) = mu_MLE(1) + sqrt(h_MLE(jj))*quantile; % time varying threshold for a nonstandardised variable
                                                    % to put to the cdf function
end

THR0 = zeros(T,1);
cond0 = zeros(T,1);
for jj=1:T
    if (jj==1)
        if (y(jj) >= mean(y))
            THR0(jj) = mean(y);
             cond0(jj) = 1;
        else
            THR0(jj) = NaN;
        end
    else
        if (y(jj) >= y(jj-1) + mean(y))
            THR0(jj) =  y(jj-1) +mean(y);
            cond0(jj) = 1;
        else
            THR0(jj) = NaN;
        end        
    end
end
y_cond0 = y(1:T);
y_cond0(~cond0) = NaN;
y0 = y_cond0((cond0==1))

hold on
plot(y_cond0(1:T),'k')
plot(THR,'r')
plot(THR0,'b')
hold off