function RESULTS = PCP_agarch11_new_multi_thres_MCLE_forecast(sdd,...
    c, sigma1, sigma2, kappa, omega, alpha, beta, ...
    p_bar0, p_bar1, p_bar, T, H,...
    THRES, ...
    RESULTS) % RESULTS = RES{1,1}
    
    s1 = RandStream('mt19937ar','Seed',sdd);
    RandStream.setGlobalStream(s1); 
 
    M = 10000;
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
    % true VaRs
    q05 = norminv(p_bar0, c, sigma2_k*sqrt(h_true(T+1:T+H)))';
    q1 = norminv(p_bar1, c, sigma2_k*sqrt(h_true(T+1:T+H)))';
    q5 = norminv(p_bar, c, sigma2_k*sqrt(h_true(T+1:T+H)))'; 
   
    % true ESs
    cdf05 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar0))/(p_bar0);
    cdf1 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar1))/(p_bar1);
    cdf5 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar))/(p_bar);              
% %     all(RESULTS.q_theor(1,:) == q05)


%% POSTERIOR    
    draw = RESULTS.post.draw_post;
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

    results_post.VaR_post = [VaR_05_post;VaR_1_post;VaR_5_post];
    results_post.ES_post = [ES_05_post;ES_1_post;ES_5_post];    

%% FIXED 0
   try
        fprintf('**** CONST 0 \n')            
        draw_C = RESULTS.thr_0.draw_C;
        draw_PC = RESULTS.thr_0.draw_PC;
   
        results_0 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_0 = {};
        fprintf('crash ****\n') 
    end

%% CONST
%% Threshold = 10% perscentile of the data sample
    try
        fprintf('**** CONST THRES = %4.2f ...\n',THRES(1))            
        draw_C = RESULTS.thr_10.draw_C;
        draw_PC = RESULTS.thr_10.draw_PC;
   
        results_10 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_10 = {};
        fprintf('crash ****\n') 
    end
    
    
 %% Threshold = 20% perscentile of the data sample
    try
        fprintf('**** CONST THRES = %4.2f ...\n',THRES(2))            
        draw_C = RESULTS.thr_20.draw_C;
        draw_PC = RESULTS.thr_20.draw_PC;
   
        results_20 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_20 = {};
        fprintf('crash ****\n') 
    end
    
    
 %% Threshold = 30% perscentile of the data sample
    try
        fprintf('**** CONST THRES = %4.2f ...\n',THRES(3))            
        draw_C = RESULTS.thr_30.draw_C;
        draw_PC = RESULTS.thr_30.draw_PC;
   
        results_30 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_30 = {};
        fprintf('crash ****\n') 
    end
       
 %% Threshold = 40% perscentile of the data sample
    try
        fprintf('**** CONST THRES = %4.2f ...\n',THRES(4))            
        draw_C = RESULTS.thr_40.draw_C;
        draw_PC = RESULTS.thr_40.draw_PC;
   
        results_40 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_40 = {};
        fprintf('crash ****\n') 
    end
    

%% VAR MLE
 %% Threshold = 10%  MLE implied predicitve distribution
    try
        fprintf('**** VAR THRES = %4.2f ...\n',THRES(1))            
        draw_C = RESULTS.thr_v10.draw_C;
        draw_PC = RESULTS.thr_v10.draw_PC;
   
        results_v10 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_v10 = {};
        fprintf('crash ****\n') 
    end
    
 %% Threshold = 20%  MLE implied predicitve distribution
    try
        fprintf('**** VAR THRES = %4.2f ...\n',THRES(2))            
        draw_C = RESULTS.thr_v20.draw_C;
        draw_PC = RESULTS.thr_v20.draw_PC;
   
        results_v20 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_v20 = {};
        fprintf('crash ****\n') 
    end
       
 %% Threshold = 30%  MLE implied predicitve distribution
    try
        fprintf('**** VAR THRES = %4.2f ...\n',THRES(3))            
        draw_C = RESULTS.thr_v30.draw_C;
        draw_PC = RESULTS.thr_v30.draw_PC;
   
        results_v30 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_v30 = {};
        fprintf('crash ****\n') 
    end
       
 %% Threshold = 40%  MLE implied predicitve distribution
    try
        fprintf('**** VAR THRES = %4.2f ...\n',THRES(4))            
        draw_C = RESULTS.thr_v40.draw_C;
        draw_PC = RESULTS.thr_v40.draw_PC;
   
        results_v40 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_v40 = {};
        fprintf('crash ****\n') 
    end
    

  
 %% MCLE    
    
 %% Threshold = 10%  MCLE implied predicitve distribution
    try
        fprintf('**** MCLE VAR THRES = %4.2f ...\n',THRES(1))            
        draw_C = RESULTS.thr_vcp10.draw_C;
        draw_PC = RESULTS.thr_vcp10.draw_PC;
   
        results_vcp10 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_vcp10 = {};
        fprintf('crash ****\n') 
    end
    
 %% Threshold = 20%  MCLE implied predicitve distribution
    try
        fprintf('**** MCLE VAR THRES = %4.2f ...\n',THRES(2))            
        draw_C = RESULTS.thr_vcp20.draw_C;
        draw_PC = RESULTS.thr_vcp20.draw_PC;
   
        results_vcp20 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_vcp20 = {};
        fprintf('crash ****\n') 
    end
        
 %% Threshold = 30%  MCLE implied predicitve distribution
    try
        fprintf('**** MCLE VAR THRES = %4.2f ...\n',THRES(3))            
        draw_C = RESULTS.thr_vcp30.draw_C;
        draw_PC = RESULTS.thr_vcp30.draw_PC;
   
        results_vcp30 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_vcp30 = {};
        fprintf('crash ****\n') 
    end
       
 %% Threshold = 40%  MCLE implied predicitve distribution
    try
        fprintf('**** MCLE VAR THRES = %4.2f ...\n',THRES(4))            
        draw_C = RESULTS.thr_vcp40.draw_C;
        draw_PC = RESULTS.thr_vcp40.draw_PC;
   
        results_vcp40 = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
                     p_bar0, p_bar1, p_bar, T, H); 
        fprintf('OK! ****\n')    
    
    catch
        results_vcp40 = {};
        fprintf('crash ****\n') 
    end    
    

        
    %% Results
    RESULTS.q_theor = [q05;q1;q5];
    RESULTS.cdf_theor = [cdf05;cdf1;cdf5];

    RESULTS.post = results_post;
    RESULTS.thr_0 = results_0;

    RESULTS.thr_10 = results_10;
    RESULTS.thr_20 = results_20;
    RESULTS.thr_30 = results_30;
    RESULTS.thr_40 = results_40;
        
    RESULTS.thr_v10 = results_v10;
    RESULTS.thr_v20 = results_v20;
    RESULTS.thr_v30 = results_v30;
    RESULTS.thr_v40 = results_v40;
 
    RESULTS.thr_vcp10 = results_vcp10;
    RESULTS.thr_vcp20 = results_vcp20;
    RESULTS.thr_vcp30 = results_vcp30;
    RESULTS.thr_vcp40 = results_vcp40;
    
%     if false
%         fields = fieldnames(RESULTS);
%         MSE = zeros(1,23,3);
%         MSE(1,1,:) = mean((RESULTS.post.VaR_post - RESULTS.q_theor).^2,2)
%         
%         for oo = 4:numel(fields)
%             MSE(1,2*(oo-2)-2,:) = oo;%mean((RESULTS.(fields{oo}).VaR_C - RESULTS.q_theor).^2,2);
%             MSE(1,2*(oo-2)-1,:) = oo;%mean((RESULTS.(fields{oo}).VaR_PC - RESULTS.q_theor).^2,2);
%         end
%     end
    
 end