function results = PCP_arch1_run_varc(sdd, c, sigma1, sigma2, kappa, omega, alpha, p_bar1, p_bar, T, H, M, BurnIn, mu_init, df, cont, options, partition, II, GamMat)
    
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

    %% Misspecified model: ARCH(1,1) normal 
    %% Uncensored Posterior
    fprintf('*** Uncensored Posterior ***\n');
    y_S = var(y(1:T));
  
    kernel_init = @(xx) -posterior_arch1_mex(xx, y(1:T), y_S)/T;
    kernel = @(xx) posterior_arch1_mex(xx, y(1:T), y_S);
    [mu_MLE,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init,options);
    Sigma = inv(T*Sigma);
    try
        [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
        [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
        lnd = dmvgt(draw, mit, true, GamMat); 
    catch
        draw = rmvt(mu_MLE,Sigma,df,M+BurnIn);
        mit = struct('mu',mu_MLE,'Sigma',reshape(Sigma,1,length(mu_MLE)^2),'df', df, 'p', 1);
        [mit, CV] = MitISEM_new(mit, kernel, mu_init, cont, GamMat);            
        [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
        lnd = dmvgt(draw, mit, true, GamMat); 
    end

    lnw = lnk - lnd;
    lnw = lnw - max(lnw);
    [ind, a] = fn_MH(lnw);
    draw = draw(ind,:);
    accept = a/(M+BurnIn);
    draw = draw(BurnIn+1:BurnIn+M,:);    

    h_post = volatility_arch1(draw,y,y_S,H);   
    y_post = randn(M,H).*sqrt(h_post);
    y_post = bsxfun(@plus,y_post,draw(:,1));
    y_post = sort(y_post);
    VaR_1_post = y_post(p_bar1*M,:); 
    VaR_5_post = y_post(p_bar*M,:); 
    mean_draw = mean(draw);
    median_draw = median(draw);
    std_draw = std(draw);


    %% CENSORED  NOPARAMETERS
    threshold = 1.0;
    fprintf('*** Censored Posterior, ad hoc time varying threshold, THR = %3.2f ***\n',threshold);
    kernel_init = @(xx) - C_posterior_arch1_varc_noparam_mex(xx, y(1:T,1), threshold, y_S)/T;    
    kernel = @(xx) C_posterior_arch1_varc_noparam_mex(xx, y(1:T,1), threshold, y_S);

%     cont.mit.CV_tol = 0.3; 
    cont.mit.CV_max = 1.9;
    CV_C = cont.mit.CV_old;
    while (CV_C(end) >= 2)
        try
            [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
            Sigma_C = inv(T*Sigma_C);
            draw_C = rmvt(mu_C,Sigma_C,df,M+BurnIn);
            mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
            [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            if CV_C(end)>2
                [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            end
            [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
            lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
        catch
            mu_C = fminunc(kernel_init,mu_init,options);
            mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma,1,length(mu_C)^2),'df', df, 'p', 1);
            [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            if CV_C(end)>2
                [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            end
            [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
            lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
        end 
    end

    lnw_C = lnk_C - lnd_C;
    lnw_C = lnw_C - max(lnw_C);
    [ind, a] = fn_MH(lnw_C);
    draw_C = draw_C(ind,:);
    accept_C = a/(M+BurnIn);
    draw_C = draw_C(BurnIn+1:BurnIn+M,:);

    h_post_C = volatility_arch1(draw_C,y,y_S,H);
    y_post_C = randn(M,H).*sqrt(h_post_C);
    y_post_C = bsxfun(@plus,y_post_C,draw_C(:,1));
    y_post_C = sort(y_post_C);
    VaR_1_post_C = y_post_C(p_bar1*M,:); 
    VaR_5_post_C = y_post_C(p_bar*M,:); 
    mean_draw_C = mean(draw_C);
    median_draw_C = median(draw_C);
    std_draw_C = std(draw_C);

    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior,  time varying threshold, THR = %3.2f ***\n',threshold);
    % mit_C: joint candidate for the joint censored posterior    
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos

    thinning = 1;
    
    [draw_PC, a_PC, lnw_PC] = sim_cond_mit_MH_outloop(mit_C, draw_short, partition, II, BurnIn, kernel, GamMat, cont.disp, thinning);
    accept_PC = mean(a_PC); 
    ind_fin = isfinite(lnw_PC);
    M_fin = sum(ind_fin);
    draw_PC = draw_PC(ind_fin,:) ;  

    mean_draw_PC = mean(draw_PC);
    median_draw_PC = median(draw_PC);
    std_draw_PC = std(draw_PC);

    h_post_PC = volatility_arch1(draw_PC,y,y_S,H);
    y_post_PC = randn(M_fin,H).*sqrt(h_post_PC);
    y_post_PC = bsxfun(@plus,y_post_PC,draw_PC(:,1));
    y_post_PC = sort(y_post_PC);
    VaR_1_post_PC = y_post_PC(round(p_bar1*M_fin),:); 
    VaR_5_post_PC = y_post_PC(round(p_bar*M_fin),:); 

    %% CENSORED MLE PARAMETERS    
    threshold = 0.1; %<---------- HiGhER?
    quantile = norminv(threshold);
    fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold);
    kernel_init = @(xx) - C_posterior_arch1_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S)/T;    
    kernel = @(xx) C_posterior_arch1_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S);
    try
        [mu_Cm,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_init,options);
    catch
        [mu_Cm,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_C,options);
    end
    Sigma_Cm = inv(T*Sigma_Cm);
      
%     cont.mit.CV_tol = 0.3; 
    cont.mit.CV_max = 1.9; %1.5;
    CV_Cm = cont.mit.CV_old;
    while (CV_Cm(end) >= 2)
        try
            draw_Cm = rmvt(mu_Cm,Sigma_Cm,df,M+BurnIn);
            mit_Cm = struct('mu',mu_Cm,'Sigma',reshape(Sigma_Cm,1,length(mu_Cm)^2),'df', df, 'p', 1);
            [mit_Cm, CV_Cm] = MitISEM_new2(mit_Cm, kernel, mu_init, cont, GamMat);   
            if CV_Cm(end)>2
                [mit_Cm, CV_Cm] = MitISEM_new2(mit_Cm, kernel, mu_init, cont, GamMat);   
            end
            [draw_Cm, lnk_Cm] = fn_rmvgt_robust(M+BurnIn, mit_Cm, kernel, false);
            lnd_Cm = dmvgt(draw_Cm, mit_Cm, true, GamMat);    
        catch
            mit_Cm = struct('mu',mu_Cm,'Sigma',reshape(Sigma,1,length(mu_Cm)^2),'df', df, 'p', 1);
            [mit_Cm, CV_Cm] = MitISEM_new2(mit_Cm, kernel, mu_init, cont, GamMat);   
            if CV_Cm(end)>2
                [mit_Cm, CV_Cm] = MitISEM_new2(mit_Cm, kernel, mu_init, cont, GamMat);   
            end
            [draw_Cm, lnk_Cm] = fn_rmvgt_robust(M+BurnIn, mit_Cm, kernel, false);
            lnd_Cm = dmvgt(draw_Cm, mit_Cm, true, GamMat);    
        end 
    end
    lnw_Cm = lnk_Cm - lnd_Cm;
    lnw_Cm = lnw_Cm - max(lnw_Cm);
    [ind, a] = fn_MH(lnw_Cm);
    draw_Cm = draw_Cm(ind,:);
    accept_Cm = a/(M+BurnIn);
    draw_Cm = draw_Cm(BurnIn+1:BurnIn+M,:);
    
    h_post_Cm = volatility_arch1(draw_Cm,y,y_S,H);
    y_post_Cm = randn(M,H).*sqrt(h_post_Cm);
    y_post_Cm = bsxfun(@plus,y_post_Cm,draw_Cm(:,1));
    y_post_Cm = sort(y_post_Cm);
    VaR_1_post_Cm = y_post_Cm(p_bar1*M,:); 
    VaR_5_post_Cm = y_post_Cm(p_bar*M,:); 
    mean_draw_Cm = mean(draw_Cm);
    median_draw_Cm = median(draw_Cm);
    std_draw_Cm = std(draw_Cm);

    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold);
    % mit_C: joint candidate for the joint censored posterior    
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
    [draw_PCm, a_PCm, lnw_PCm] = sim_cond_mit_MH_outloop(mit_Cm, draw_short, partition, II, BurnIn, kernel, GamMat, cont.disp);
    accept_PCm = mean(a_PCm); 
    ind_fin = isfinite(lnw_PCm);
    M_fin = sum(ind_fin);
    draw_PCm = draw_PCm(ind_fin,:) ;  

    mean_draw_PCm = mean(draw_PCm);
    median_draw_PCm = median(draw_PCm);
    std_draw_PCm = std(draw_PCm);

    h_post_PCm = volatility_arch1(draw_PCm,y,y_S,H);
    y_post_PCm = randn(M_fin,H).*sqrt(h_post_PCm);
    y_post_PCm = bsxfun(@plus,y_post_PCm,draw_PCm(:,1));
    y_post_PCm = sort(y_post_PCm);
    VaR_1_post_PCm = y_post_PCm(round(p_bar1*M_fin),:); 
    VaR_5_post_PCm = y_post_PCm(round(p_bar*M_fin),:); 
   
   %% Results
   results = struct('y',y,'draw',draw,'draw_C',draw_C,'draw_PC',draw_PC,'draw_Cm',draw_Cm,'draw_PCm',draw_PCm,...
        'q1',q1,'q5',q5,...
        'mean_draw',mean_draw,'mean_draw_C',mean_draw_C,'mean_draw_PC',mean_draw_PC,'mean_draw_Cm',mean_draw_Cm,'mean_draw_PCm',mean_draw_PCm,...
        'median_draw',median_draw,'median_draw_C',median_draw_C,'median_draw_PC',median_draw_PC,'median_draw_Cm',median_draw_Cm,'median_draw_PCm',median_draw_PCm,...
        'std_draw',std_draw,'std_draw_C',std_draw_C,'std_draw_PC',std_draw_PC,'std_draw_Cm',std_draw_Cm,'std_draw_PCm',std_draw_PCm,...
        'accept',accept,'accept_C',accept_C,'accept_PC',accept_PC,'accept_Cm',accept_Cm,'accept_PCm',accept_PCm,...
        'mit',mit,'CV',CV,'mit_C',mit_C,'CV_C',CV_C,'mit_Cm',mit_Cm,'CV_Cm',CV_Cm,...
        'VaR_1',VaR_1,'VaR_1_post',VaR_1_post,'VaR_1_post_C',VaR_1_post_C,'VaR_1_post_PC',VaR_1_post_PC,'VaR_1_post_Cm',VaR_1_post_Cm,'VaR_1_post_PCm',VaR_1_post_PCm,...
        'VaR_5',VaR_5,'VaR_5_post',VaR_5_post,'VaR_5_post_C',VaR_5_post_C,'VaR_5_post_PC',VaR_5_post_PC,'VaR_5_post_Cm',VaR_5_post_Cm,'VaR_5_post_PCm',VaR_5_post_PCm);

end 