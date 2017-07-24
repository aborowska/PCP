function results = PCP_ar1_run_varc(c, sigma1, sigma2, rho, p_bar1, p_bar, T, H, M, BurnIn, mu_init, df, cont, options, partition, II, GamMat)

    %% simple AR(1)
    eps = randn(T+H,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);

    y = zeros(T+H,1);
    y(1,1) = eps(1,1);
    for ii = 2:T+H
        y(ii,1) = rho*y(ii-1,1) + eps(ii,1);
    end
    % true VaRs
    q1 = norminv(p_bar1,c+rho*y(T:(T+H-1),1),sigma2)';
    q5 = norminv(p_bar,c+rho*y(T:(T+H-1),1),sigma2)'; 
     
    % MC VaRs under the true model
    eps_sort = randn(M,H);
    ind = (eps_sort>0);
    eps_sort(ind) = c + sigma1.*eps_sort(ind);
    eps_sort(~ind) = c + sigma2.*eps_sort(~ind);

    y_sort = bsxfun(@plus,eps_sort,rho*y(T:(T+H-1),1)');
    y_sort = sort(y_sort);
    VaR_1 = y_sort(p_bar1*M,:); 
    VaR_5 = y_sort(p_bar*M,:); 

    %% Misspecified model: AR1 normal with unknown mu and sigma
    %% UNCENSORED posterior
    fprintf('*** Uncensored Posterior ***\n');
%     kernel_init = @(xx) -posterior_ar1(xx,y(1:T))/T;

    kernel_init = @(xx) -posterior_ar1_mex(xx,y(1:T))/T;
    kernel = @(xx) posterior_ar1_mex(xx,y(1:T));
    [mu_mle,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init, options);
    
    try
        [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
    catch
        Sigma = inv(T*Sigma);
        mit = struct('mu',mu_mle,'Sigma',reshape(Sigma,1,length(mu_mle)^2),'df', df, 'p', 1);
        [mit, CV] = MitISEM_new(mit, kernel, mu_init, cont, GamMat);
    end
    [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
    lnd = dmvgt(draw, mit, true, GamMat); 
    lnw = lnk - lnd;
    lnw = lnw - max(lnw);
    [ind, a] = fn_MH(lnw);
    draw = draw(ind,:);
    lnw = lnw(ind);
    accept = a/(M+BurnIn);
    draw = draw(BurnIn+1:BurnIn+M,:);    
 
    y_post = bsxfun(@times,randn(M,H),draw(:,2));
    y_post = bsxfun(@plus,y_post,draw(:,1));
    y_post = y_post + draw(:,3)*y(T:(T+H-1),1)';
    y_post = sort(y_post);
    VaR_1_post = y_post(p_bar1*M,:); 
    VaR_5_post = y_post(p_bar*M,:); 
    mean_draw = mean(draw);
    std_draw = std(draw);
    
    %% Time varying threshold, THR = 1
%     threshold = sort(y(1:T));
%     threshold = threshold(2*p_bar*T);
%     threshold = 1; %0.9;
    %% CENSORED
    fprintf('*** Censored Posterior, time varying threshold, THR = 1 ***\n');
%     kernel_init = @(xx) - C_posterior_ar1_varc_mex(xx, y(1:T), threshold)/T; 
%     kernel = @(xx) C_posterior_ar1_varc_mex(xx, y(1:T), threshold); 
%     [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
%  
    % NOPARAMETERS
    threshold = 1;
    kernel_init = @(xx) - C_posterior_ar1_varc_noparam_mex(xx, y(1:T), threshold)/T; 
    kernel = @(xx) C_posterior_ar1_varc_noparam_mex(xx, y(1:T), threshold);     
    [mu_C2,~,~,~,~,Sigma_C2] = fminunc(kernel_init,mu_init,options);
    Sigma_C2 = inv(T*Sigma_C2);
 
%     fprintf('*** Censored Posterior, threshold = 0 ***\n');
%     threshold0 = 0;
%     kernel_init = @(xx) - C_posterior_ar1(xx, y(1:T), threshold0)/T; 
%     kernel = @(xx) C_posterior_ar1(xx, y(1:T), threshold0); 
%     [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu_init,options);
    
    try
        [mit_C, CV_C] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);   
    catch
        [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
        Sigma_C = inv(T*Sigma_C);
        mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
        CV_C = cont.mit.CV_old;
    end
    [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
    lnd_C = dmvgt(draw_C, mit_C, true, GamMat);     
    lnw_C = lnk_C - lnd_C;
    lnw_C = lnw_C - max(lnw_C);
    [ind, a] = fn_MH(lnw_C);
    draw_C = draw_C(ind,:);
    accept_C = a/(M+BurnIn);
    draw_C = draw_C(BurnIn+1:BurnIn+M,:);
    mean_draw_C = mean(draw_C);
    std_draw_C = std(draw_C);
        
    y_post_C = bsxfun(@times,randn(M,H),draw_C(:,2));
    y_post_C = bsxfun(@plus,y_post_C,draw_C(:,1));
    y_post_C = y_post_C + draw_C(:,3)*y(T:(T+H-1),1)';
    y_post_C = sort(y_post_C);
    VaR_1_post_C  = y_post_C(p_bar1*M,:); 
    VaR_5_post_C  = y_post_C(p_bar*M,:); 
    
    %% PARTIAL CENSORING: keep rho uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior, time varying threshold, THR = 1 ***\n');
    % mit_C: joint candidate for the joint censored posterior
    % Short version
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos 
    [draw_PC, a_PC] = sim_cond_mit_MH_outloop(mit_C, draw_short, partition, II, BurnIn, kernel, GamMat, cont.disp);
    accept_PC = mean(a_PC);
    mean_draw_PC = mean(draw_PC);
    std_draw_PC = std(draw_PC);
  
    y_post_PC = bsxfun(@times,randn(M,H),draw_PC(:,2));
    y_post_PC = bsxfun(@plus,y_post_PC,draw_PC(:,1));
    y_post_PC = y_post_PC + draw_PC(:,3)*y(T:(T+H-1),1)';
    y_post_PC = sort(y_post_PC);
    VaR_1_post_PC  = y_post_PC(p_bar1*M,:); 
    VaR_5_post_PC  = y_post_PC(p_bar*M,:); 

    %% CENSORED MLE PARAMETERS    
    threshold = 0.1; %<---------- HiGhER?
    quantile = norminv(threshold);
    fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold);

    kernel_init = @(xx) - C_posterior_ar1_varc_mle_mex(xx, y(1:T), mu_mle, quantile)/T; 
%     kernel_init = @(xx) - C_posterior_ar1_varc_mle(xx, y(1:T), mu_mle, quantile)/T; 
    kernel = @(xx) C_posterior_ar1_varc_mle_mex(xx, y(1:T), mu_mle, quantile);     
    [mu_Cm,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_init,options);
    Sigma_Cm = inv(T*Sigma_Cm);
    
    try
        [mit_Cm, CV_Cm] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);   
    catch
        [mu_Cm,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_init,options);
        Sigma_Cm = inv(T*Sigma_Cm);
        mit_Cm = struct('mu',mu_Cm,'Sigma',reshape(Sigma_Cm,1,length(mu_Cm)^2),'df', df, 'p', 1);
        CV_Cm = cont.mit.CV_old;
    end
    [draw_Cm, lnk_Cm] = fn_rmvgt_robust(M+BurnIn, mit_Cm, kernel, false);
    lnd_Cm = dmvgt(draw_Cm, mit_Cm, true, GamMat);     
    lnw_Cm = lnk_Cm - lnd_Cm;
    lnw_Cm = lnw_Cm - max(lnw_Cm);
    [ind, a] = fn_MH(lnw_Cm);
    draw_Cm = draw_Cm(ind,:);
    accept_Cm = a/(M+BurnIn);
    draw_Cm = draw_Cm(BurnIn+1:BurnIn+M,:);
    mean_draw_Cm = mean(draw_Cm);
    std_draw_Cm = std(draw_Cm);
        
    y_post_Cm = bsxfun(@times,randn(M,H),draw_Cm(:,2));
    y_post_Cm = bsxfun(@plus,y_post_Cm,draw_Cm(:,1));
    y_post_Cm = y_post_Cm + draw_Cm(:,3)*y(T:(T+H-1),1)';
    y_post_Cm = sort(y_post_Cm);
    VaR_1_post_Cm  = y_post_Cm(p_bar1*M,:); 
    VaR_5_post_Cm  = y_post_Cm(p_bar*M,:);
    
    
    %% PARTIAL CENSORING: keep rho uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold);
    % mit_C: joint candidate for the joint censored posterior
    % Short version
%     draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos 
    [draw_PCm, a_PCm] = sim_cond_mit_MH_outloop(mit_Cm, draw_short, partition, II, BurnIn, kernel, GamMat, cont.disp);
    accept_PCm = mean(a_PCm);
    mean_draw_PCm = mean(draw_PCm);
    std_draw_PCm = std(draw_PCm);
  
    y_post_PCm = bsxfun(@times,randn(M,H),draw_PCm(:,2));
    y_post_PCm = bsxfun(@plus,y_post_PCm,draw_PCm(:,1));
    y_post_PCm = y_post_PCm + draw_PCm(:,3)*y(T:(T+H-1),1)';
    y_post_PCm = sort(y_post_PCm);
    VaR_1_post_PCm  = y_post_PCm(p_bar1*M,:); 
    VaR_5_post_PCm  = y_post_PCm(p_bar*M,:);    
    
    %% Results
    results = struct('y',y,'draw',draw,'draw_C',draw_C,'draw_PC',draw_PC,'draw_Cm',draw_Cm,'draw_PCm',draw_PCm,...
        'q1',q1,'q5',q5,...
        'mean_draw',mean_draw,'mean_draw_C',mean_draw_C,'mean_draw_PC',mean_draw_PC,'mean_draw_Cm',mean_draw_Cm,'mean_draw_PCm',mean_draw_PCm,...
        'std_draw',std_draw,'std_draw_C',std_draw_C,'std_draw_PC',std_draw_PC,'std_draw_Cm',std_draw_Cm,'std_draw_PCm',std_draw_PCm,...
        'accept',accept,'accept_C',accept_C,'accept_PC',accept_PC,'accept_Cm',accept_Cm,'accept_PCm',accept_PCm,...
        'mit',mit,'CV',CV,'mit_C',mit_C,'CV_C',CV_C,'mit_Cm',mit_C,'CV_Cm',CV_C,...
        'VaR_1',VaR_1,'VaR_1_post',VaR_1_post,'VaR_1_post_C',VaR_1_post_C,'VaR_1_post_PC',VaR_1_post_PC,'VaR_1_post_Cm',VaR_1_post_Cm,'VaR_1_post_PCm',VaR_1_post_PCm,...
        'VaR_5',VaR_5,'VaR_5_post',VaR_5_post,'VaR_5_post_C',VaR_5_post_C,'VaR_5_post_PC',VaR_5_post_PC,'VaR_5_post_Cm',VaR_5_post_Cm,'VaR_5_post_PCm',VaR_5_post_PCm);
end