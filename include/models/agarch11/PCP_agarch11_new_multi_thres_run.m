function results = PCP_agarch11_new_multi_thres_run(sdd,...
    c, sigma1, sigma2, kappa, omega, alpha, beta, ...
    p_bar0, p_bar1, p_bar, T, H,...
    sampling_opt, cont, options, GamMat,...
    THRES)
    
    s = RandStream('mt19937ar','Seed',sdd);
    RandStream.setGlobalStream(s); 
        
    M = sampling_opt. M;
    BurnIn = sampling_opt.BurnIn;
    mu_init = sampling_opt.mu_init;
    df = sampling_opt.df;   

    
    sigma2_k = sigma2/sqrt(kappa);

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

    % true VaRs
    q05 = norminv(p_bar0, c, sigma2_k*sqrt(h_true(T+1:T+H)))';
    q1 = norminv(p_bar1, c, sigma2_k*sqrt(h_true(T+1:T+H)))';
    q5 = norminv(p_bar, c, sigma2_k*sqrt(h_true(T+1:T+H)))'; 
   
    % true ESs
    cdf05 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar0))/(p_bar0);
    cdf1 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar1))/(p_bar1);
    cdf5 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar))/(p_bar);    
        
%%     % MC VaRs under the true model
%     eps_sort = randn(M,H);
%     ind = (eps_sort>0);
%     eps_sort(ind) = c + sigma1.*eps_sort(ind);
%     eps_sort(~ind) = c + sigma2.*eps_sort(~ind);
%     eps_sort = eps_sort/sqrt(kappa); 
% 
%     y_sort = bsxfun(@times,eps_sort,sqrt(h_true(T+1:T+H,1))');
%     y_sort = sort(y_sort);
%     
%     VaR_1 = y_sort(p_bar1*M,:); 
%     VaR_5 = y_sort(p_bar*M,:); 
%     VaR_05 = y_sort(p_bar0*M,:); 
% 
%     ES_1 = mean(y_sort(1:p_bar1*M,:)); 
%     ES_5 = mean(y_sort(1:p_bar*M,:)); 
%     ES_05 = mean(y_sort(1:p_bar0*M,:));     
    
    %% Misspecified model: AGARCH(1,1) normal 
    %% Uncensored Posterior
    fprintf('*** Uncensored Posterior ***\n');
    y_S = var(y(1:T));
%         mu_init(1,1) = mean(y(1:T));
%         mu_init(1,2) = y_S;
    kernel_init = @(xx) -posterior_agarch11_mex(xx, y(1:T), y_S)/T;
    kernel = @(xx) posterior_agarch11_mex(xx, y(1:T), y_S);
    try
        [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
        [draw, accept] = IndMH_mit(mit, kernel,M,BurnIn,GamMat);
    catch
        [mu,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init,options);
        Sigma = inv(T*Sigma);
        draw = rmvt(mu,Sigma,df,M+BurnIn);
        mit = struct('mu',mu,'Sigma',reshape(Sigma,1,length(mu)^2),'df', df, 'p', 1);
        [mit, CV] = MitISEM_new(mit, kernel, mu_init, cont, GamMat);            
        [draw, accept] = IndMH_mit(mit, kernel,M,BurnIn,GamMat);
    end

    mu_MLE = mean(draw);
    
    h_post = volatility_agarch11(draw,y,y_S,H);   
    y_post = randn(M,H).*sqrt(h_post);
    y_post = bsxfun(@plus,y_post,draw(:,1));
    y_post = sort(y_post);
    
    VaR_1_post = y_post(p_bar1*M,:); 
    VaR_5_post = y_post(p_bar*M,:); 
    VaR_05_post = y_post(p_bar0*M,:); 

    ES_1_post = mean(y_post(1:p_bar1*M,:)); 
    ES_5_post = mean(y_post(1:p_bar*M,:)); 
    ES_05_post = mean(y_post(1:p_bar0*M,:)); 

    results_post.mit = mit;
    results_post.CV = CV;  
    results_post.draw_post = draw;        
    results_post.accept_post = accept;
    results_post.VaR_post = [VaR_05_post;VaR_1_post;VaR_5_post];
    results_post.ES_post = [ES_05_post;ES_1_post;ES_5_post];    
    
    %% Threshold = 10% perscentile of the data sample
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES(1)*T));
    fprintf('**** THRES = %4.2f ... \n',THRES(1))    
    % CP & PCP
    kernel_init = @(xx) - C_posterior_agarch11_mex(xx, y(1:T,1), threshold, y_S)/T;    
    kernel = @(xx) C_posterior_agarch11_mex(xx, y(1:T,1), threshold, y_S);
    try
        results_10 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat); 
        fprintf('OK! ****\n')
    catch
        results_10 = {};
        fprintf('crash ****\n')        
    end
    
    %% Threshold = 20% perscentile of the data sample
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES(2)*T));
    fprintf('**** THRES = %4.2f ...\n',THRES(2))   
    % CP & PCP
    kernel_init = @(xx) - C_posterior_agarch11_mex(xx, y(1:T,1), threshold, y_S)/T;    
    kernel = @(xx) C_posterior_agarch11_mex(xx, y(1:T,1), threshold, y_S);
    try
        results_20 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat);
        fprintf('OK! ****\n')
    catch
        results_20 = {};
        fprintf('crash ****\n')                
    end    
    
    
    %% Threshold = 30% perscentile of the data sample
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES(3)*T));
    fprintf('**** THRES = %4.2f ...\n',THRES(3))   
    % CP & PCP
    kernel_init = @(xx) - C_posterior_agarch11_mex(xx, y(1:T,1), threshold, y_S)/T;    
    kernel = @(xx) C_posterior_agarch11_mex(xx, y(1:T,1), threshold, y_S);
    try
        results_30 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat); 
        fprintf('OK! ****\n')
    catch
        results_30 = {};
        fprintf('crash ****\n')               
    end
    

    %% Threshold = 40% perscentile of the data sample
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES(4)*T));
    fprintf('**** THRES = %4.2f ...\n',THRES(4))
    % CP & PCP
    kernel_init = @(xx) - C_posterior_agarch11_mex(xx, y(1:T,1), threshold, y_S)/T;    
    kernel = @(xx) C_posterior_agarch11_mex(xx, y(1:T,1), threshold, y_S);
    try
        results_40 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat); 
        fprintf('OK! ****\n')    
    catch
        results_40 = {};
        fprintf('crash ****\n')               
    end
    
    
    %% Threshold = 0
    fprintf('**** ZERO ...****\n')  
    threshold0 = 0;
    kernel_init = @(xx) - C_posterior_agarch11_mex(xx, y(1:T,1), threshold0, y_S)/T;    
    kernel = @(xx) C_posterior_agarch11_mex(xx, y(1:T,1), threshold0, y_S);
    
    try
        results_0 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat); 
        fprintf('OK! ****\n')
    catch
        results_0 = {};
        fprintf('crash ****\n')               
    end
    

    
    
    
   %% Threshold = 10% perscentile of the data sample
    quantile = norminv(THRES(1));
    fprintf('**** VAR THRES = %4.2f ...\n',THRES(1))    
    % CP & PCP
    kernel_init = @(xx) - C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S)/T;    
    kernel = @(xx) C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S);
    try
        results_v10 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat); 
        fprintf('OK! ****\n')    
    catch
        results_v10 = {};
        fprintf('crash ****\n')               
    end
    
    %% Threshold = 10% perscentile of the data sample
    quantile = norminv(THRES(2));
    fprintf('**** VAR THRES = %4.2f ...\n',THRES(2))   
    % CP & PCP
    kernel_init = @(xx) - C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S)/T;    
    kernel = @(xx) C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S);
    try
        results_v20 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat); 
        fprintf('OK! ****\n')
    catch
        results_v20 = {};
        fprintf('crash ****\n')               
    end    
    
    
    %% Threshold = 30% perscentile of the data sample
    quantile = norminv(THRES(3));
    fprintf('**** VAR THRES = %4.2f ...\n',THRES(3))   
    % CP & PCP
    kernel_init = @(xx) - C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S)/T;    
    kernel = @(xx) C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S);
    try
        results_v30 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat); 
        fprintf('OK! ****\n')  
    catch
        results_v30 = {};
        fprintf('crash ****\n')               
    end
    

    %% Threshold = 40% perscentile of the data sample
    quantile = norminv(THRES(4));
    fprintf('**** VAR THRES = %4.2f ...\n',THRES(4))
    % CP & PCP
    kernel_init = @(xx) - C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S)/T;    
    kernel = @(xx) C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S);
    try
        results_v40 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat); 
        fprintf('OK! ****\n')    
    catch
        results_v40 = {};
        fprintf('crash ****\n')                
    end
        
    %% Results   
    results.q_theor = [q05;q1;q5];
    results.cdf_theor = [cdf05;cdf1;cdf5];

    results.post = results_post;
    results.thr_10 = results_10;
    results.thr_20 = results_20;
    results.thr_30 = results_30;
    results.thr_40 = results_40;
    results.thr_0 = results_0;
    results.thr_v10 = results_v10;
    results.thr_v20 = results_v20;
    results.thr_v30 = results_v30;
    results.thr_v40 = results_v40;
 end