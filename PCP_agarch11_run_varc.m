function results = PCP_agarch11_run_varc(sdd, c, sigma1, sigma2, kappa, omega, alpha, beta, p_bar1, p_bar, T, H, M, BurnIn, mu_init, df, cont, options, partition, II, GamMat)
    
    s = RandStream('mt19937ar','Seed',sdd);
    RandStream.setGlobalStream(s); 
    
    sigma1_k = sigma1/sqrt(kappa);
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

    %% Misspecified model: GARCH(1,1) normal 
    %% Uncensored Posterior
    fprintf('*** Uncensored Posterior ***\n');
    y_S = var(y(1:T));
    kernel_init = @(xx) -posterior_agarch11_mex(xx, y(1:T), y_S)/T;
    kernel = @(xx) posterior_agarch11_mex(xx, y(1:T), y_S);
    try
        [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
        [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
        lnd = dmvgt(draw, mit, true, GamMat); 
    catch
        [mu,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init,options);
        Sigma = inv(T*Sigma);
        draw = rmvt(mu,Sigma,df,M+BurnIn);
        mit = struct('mu',mu,'Sigma',reshape(Sigma,1,length(mu)^2),'df', df, 'p', 1);
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

    h_post = volatility_agarch11(draw,y,y_S,H);   
    y_post = randn(M,H).*sqrt(h_post);
    y_post = bsxfun(@plus,y_post,draw(:,1));
    y_post = sort(y_post);
    VaR_1_post = y_post(p_bar1*M,:); 
    VaR_5_post = y_post(p_bar*M,:); 
    mean_draw = mean(draw);
    median_draw = median(draw);
    std_draw = std(draw);

    %% Threshold = 10% perscentile of the data sample
%     threshold = sort(y(1:T));
%     threshold = threshold(2*p_bar*T);
    %% CENSORED  NOPARAMETERS
    threshold = 0.5;
    fprintf('*** Censored Posterior, time varying threshold, THR = %3.2f ***\n',threshold);
    kernel_init = @(xx) - C_posterior_agarch11_varc_noparam_mex(xx, y(1:T,1), threshold, y_S)/T;    
    kernel = @(xx) C_posterior_agarch11_varc_noparam_mex(xx, y(1:T,1), threshold, y_S);
    [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);

%     cont.mit.CV_tol = 0.3; 
    cont.mit.CV_max = 1.5;
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

    h_post_C = volatility_agarch11(draw_C,y,y_S,H);
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

    h_post_PC = volatility_agarch11(draw_PC,y,y_S,H);
    y_post_PC = randn(M_fin,H).*sqrt(h_post_PC);
    y_post_PC = bsxfun(@plus,y_post_PC,draw_PC(:,1));
    y_post_PC = sort(y_post_PC);
    VaR_1_post_PC = y_post_PC(round(p_bar1*M_fin),:); 
    VaR_5_post_PC = y_post_PC(round(p_bar*M_fin),:); 

     %% Results
    results = struct('y',y,'draw',draw,'draw_C',draw_C,'draw_PC',draw_PC,...
        'q1',q1,'q5',q5,...
        'mean_draw',mean_draw,'mean_draw_C',mean_draw_C,'mean_draw_PC',mean_draw_PC,...
        'median_draw',median_draw,'median_draw_C',median_draw_C,'median_draw_PC',median_draw_PC,...
        'std_draw',std_draw,'std_draw_C',std_draw_C,'std_draw_PC',std_draw_PC,...
        'accept',accept,'accept_C',accept_C,'accept_PC',accept_PC,...
        'mit',mit,'CV',CV,'mit_C',mit_C,'CV_C',CV_C,...
        'VaR_1',VaR_1,'VaR_1_post',VaR_1_post,'VaR_1_post_C',VaR_1_post_C,'VaR_1_post_PC',VaR_1_post_PC,...
        'VaR_5',VaR_5,'VaR_5_post',VaR_5_post,'VaR_5_post_C',VaR_5_post_C,'VaR_5_post_PC',VaR_5_post_PC);

end