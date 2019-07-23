function RESULTS = PCP_agarch11_new_multi_thres_run_MCLE(sdd,...
    c, sigma1, sigma2, kappa, omega, alpha, beta, ...
    p_bar0, p_bar1, p_bar, T, H,...
    sampling_opt, cont, options, GamMat,...
    THRES, ...
    RESULTS) % RESULTS = RES{1,1}
    
    s1 = RandStream('mt19937ar','Seed',sdd);
    RandStream.setGlobalStream(s1); 
        
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
    y_S = var(y(1:T));
%     % true VaRs
%     q05 = norminv(p_bar0, c, sigma2_k*sqrt(h_true(T+1:T+H)))';
%     q1 = norminv(p_bar1, c, sigma2_k*sqrt(h_true(T+1:T+H)))';
%     q5 = norminv(p_bar, c, sigma2_k*sqrt(h_true(T+1:T+H)))'; 
%    
%     % true ESs
%     cdf05 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar0))/(p_bar0);
%     cdf1 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar1))/(p_bar1);
%     cdf5 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar))/(p_bar);              
% %     all(RESULTS.q_theor(1,:) == q05)
    %% Misspecified model: AGARCH(1,1) normal 

    draw = RESULTS.post.draw_post;

%% Threshold = 10% perscentile of the data sample
    try
        mu_MCLE = mean(RESULTS.thr_10.draw_C);
   
        quantile = norminv(THRES(1));
        fprintf('**** VAR THRES = %4.2f ...\n',THRES(1))    
        % CP & PCP
        kernel_init = @(xx) - C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MCLE, quantile, y_S)/T;    
        kernel = @(xx) C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MCLE, quantile, y_S);

        results_v10 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat); 
        fprintf('OK! ****\n')    
    
    catch
        results_v10 = {};
        fprintf('crash ****\n') 
    end
    
%% Threshold = 20% perscentile of the data sample
    try
        mu_MCLE = mean(RESULTS.thr_20.draw_C);
   
        quantile = norminv(THRES(2));
        fprintf('**** VAR THRES = %4.2f ...\n',THRES(2))    
        % CP & PCP
        kernel_init = @(xx) - C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MCLE, quantile, y_S)/T;    
        kernel = @(xx) C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MCLE, quantile, y_S);

        results_v20 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat); 
        fprintf('OK! ****\n')    
    
    catch
        results_v20 = {};
        fprintf('crash ****\n') 
    end
    
 %% Threshold = 30% perscentile of the data sample
    try
        mu_MCLE = mean(RESULTS.thr_30.draw_C);
   
        quantile = norminv(THRES(3));
        fprintf('**** VAR THRES = %4.2f ...\n',THRES(3))    
        % CP & PCP
        kernel_init = @(xx) - C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MCLE, quantile, y_S)/T;    
        kernel = @(xx) C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MCLE, quantile, y_S);

        results_v30 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat); 
        fprintf('OK! ****\n')    
    
    catch
        results_v30 = {};
        fprintf('crash ****\n') 
    end
       

    
    %% Threshold = 40% perscentile of the data sample
    try
        mu_MCLE = mean(RESULTS.thr_40.draw_C);
   
        quantile = norminv(THRES(4));
        fprintf('**** VAR THRES = %4.2f ...\n',THRES(4))    
        % CP & PCP
        kernel_init = @(xx) - C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MCLE, quantile, y_S)/T;    
        kernel = @(xx) C_posterior_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MCLE, quantile, y_S);

        results_v40 = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
                     p_bar0, p_bar1, p_bar, T, H, ...
                     sampling_opt, cont, options, GamMat); 
        fprintf('OK! ****\n')    
    
    catch
        results_v40 = {};
        fprintf('crash ****\n') 
    end
    
    

        
    %% Results   
 
    RESULTS.thr_vcp10 = results_v10;
    RESULTS.thr_vcp20 = results_v20;
    RESULTS.thr_vcp30 = results_v30;
    RESULTS.thr_vcp40 = results_v40;
    
    if false
        fields = fieldnames(RESULTS);
        MSE = zeros(1,23,3);
        MSE(1,1,:) = mean((RESULTS.post.VaR_post - RESULTS.q_theor).^2,2)
        
        for oo = 4:numel(fields)
            MSE(1,2*(oo-2)-2,:) = oo;%mean((RESULTS.(fields{oo}).VaR_C - RESULTS.q_theor).^2,2);
            MSE(1,2*(oo-2)-1,:) = oo;%mean((RESULTS.(fields{oo}).VaR_PC - RESULTS.q_theor).^2,2);
        end
    end
    
 end