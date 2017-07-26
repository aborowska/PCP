function results = PCP_garch11_run_varc(sdd, c, sigma1, sigma2, kappa, ...
    omega, alpha, beta, p_bar0, p_bar1, p_bar, T, H, M, BurnIn, mu_init,...
    df, cont, options, partition, II, GamMat)
    
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
    q05 = norminv(p_bar0, c, sigma2_k*sqrt(h_true(T+1:T+H)))';    
    q1 = norminv(p_bar1, c, sigma2_k*sqrt(h_true(T+1:T+H)))';
    q5 = norminv(p_bar, c, sigma2_k*sqrt(h_true(T+1:T+H)))'; 
    
    % true ESs
    cdf05 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar0))/(p_bar0);
    cdf1 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar1))/(p_bar1);
    cdf5 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar))/(p_bar);    
    
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
    VaR_05 = y_sort(p_bar0*M,:);

    ES_1 = mean(y_sort(1:p_bar1*M,:)); 
    ES_5 = mean(y_sort(1:p_bar*M,:)); 
    ES_05 = mean(y_sort(1:p_bar0*M,:));     

    %% Misspecified model: GARCH(1,1) normal 
    %% Uncensored Posterior
    fprintf('*** Uncensored Posterior ***\n');
    y_S = var(y(1:T)); 
    kernel_init = @(xx) -posterior_garch11_mex(xx, y(1:T), y_S)/T;
    kernel = @(xx) posterior_garch11_mex(xx, y(1:T), y_S);
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

    h_post = volatility_garch11(draw,y,y_S,H);   
    y_post = randn(M,H).*sqrt(h_post);
    y_post = bsxfun(@plus,y_post,draw(:,1));
    y_post = sort(y_post);
    
    VaR_1_post = y_post(p_bar1*M,:); 
    VaR_5_post = y_post(p_bar*M,:); 
    VaR_05_post = y_post(p_bar0*M,:); 

    ES_1_post = mean(y_post(1:p_bar1*M,:)); 
    ES_5_post = mean(y_post(1:p_bar*M,:)); 
    ES_05_post = mean(y_post(1:p_bar0*M,:)); 
    
    mean_draw = mean(draw);
    median_draw = median(draw);
    std_draw = std(draw);

    %% CENSORED  NOPARAMETERS
    threshold = 1.0;
    fprintf('*** Censored Posterior, ad hoc time varying threshold, THR = %3.2f ***\n',threshold);
    kernel_init = @(xx) - C_posterior_garch11_varc_noparam_mex(xx, y(1:T,1), threshold, y_S)/T;    
    kernel = @(xx) C_posterior_garch11_varc_noparam_mex(xx, y(1:T,1), threshold, y_S);
%     [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
      
%     cont.mit.CV_tol = 0.3; 
    cont.mit.CV_max = 1.5; %1.9;
    CV_Cah = cont.mit.CV_old;
    while (CV_Cah(end) >= 2)
        try
            [mu_Cah,~,~,~,~,Sigma_Cah] = fminunc(kernel_init,mu_init,options);
            Sigma_Cah = inv(T*Sigma_Cah);
            draw_Cah = rmvt(mu_Cah,Sigma_Cah,df,M+BurnIn);
            mit_Cah = struct('mu',mu_Cah,'Sigma',reshape(Sigma_Cah,1,length(mu_Cah)^2),'df', df, 'p', 1);
            [mit_Cah, CV_Cah] = MitISEM_new2(mit_Cah, kernel, mu_init, cont, GamMat);   
            if CV_Cah(end)>2
                [mit_Cah, CV_Cah] = MitISEM_new2(mit_Cah, kernel, mu_init, cont, GamMat);   
            end
            [draw_Cah, lnk_Cah] = fn_rmvgt_robust(M+BurnIn, mit_Cah, kernel, false);
            lnd_Cah = dmvgt(draw_Cah, mit_Cah, true, GamMat);    
        catch
            [mu_Cah,~,~,~,~,Sigma_Cah] = fminunc(kernel_init,mu_init,options);
            mit_Cah = struct('mu',mu_Cah,'Sigma',reshape(Sigma_Cah,1,length(mu_Cah)^2),'df', df, 'p', 1);
            [mit_Cah, CV_Cah] = MitISEM_new2(mit_Cah, kernel, mu_init, cont, GamMat);   
            if CV_Cah(end)>2
                [mit_Cah, CV_Cah] = MitISEM_new2(mit_Cah, kernel, mu_init, cont, GamMat);   
            end
            [draw_Cah, lnk_Cah] = fn_rmvgt_robust(M+BurnIn, mit_Cah, kernel, false);
            lnd_Cah = dmvgt(draw_Cah, mit_Cah, true, GamMat);    
        end 
    end
    lnw_Cah = lnk_Cah - lnd_Cah;
    lnw_Cah = lnw_Cah - max(lnw_Cah);
    [ind, a] = fn_MH(lnw_Cah);
    draw_Cah = draw_Cah(ind,:);
    accept_Cah = a/(M+BurnIn);
    draw_Cah = draw_Cah(BurnIn+1:BurnIn+M,:);

    h_post_Cah = volatility_garch11(draw_Cah,y,y_S,H);
    y_post_Cah = randn(M,H).*sqrt(h_post_Cah);
    y_post_Cah = bsxfun(@plus,y_post_Cah,draw_Cah(:,1));
    y_post_Cah = sort(y_post_Cah);
    
    VaR_1_post_Cah = y_post_Cah(p_bar1*M,:); 
    VaR_5_post_Cah = y_post_Cah(p_bar*M,:); 
    VaR_05_post_Cah = y_post_Cah(p_bar0*M,:); 

    ES_1_post_Cah = mean(y_post_Cah(1:p_bar1*M,:)); 
    ES_5_post_Cah = mean(y_post_Cah(1:p_bar*M,:)); 
    ES_05_post_Cah = mean(y_post_Cah(1:p_bar0*M,:)); 
    
    mean_draw_Cah = mean(draw_Cah);
    median_draw_Cah = median(draw_Cah);
    std_draw_Cah = std(draw_Cah);

    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior,  time varying threshold, THR = %3.2f ***\n',threshold);
    % mit_Cah: joint candidate for the joint censored posterior    
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
    [draw_PCah, a_PCah, lnw_PCah] = sim_cond_mit_MH_outloop(mit_Cah, draw_short, partition, II, BurnIn, kernel, GamMat, cont.disp);
    accept_PCah = mean(a_PCah); 
%     ind_fin = isfinite(lnw_PCah);
%     M_fin = sum(ind_fin);
%     draw_PCah = draw_PCah(ind_fin,:) ;  
M_fin = M;
    mean_draw_PCah = mean(draw_PCah);
    median_draw_PCah = median(draw_PCah);
    std_draw_PCah = std(draw_PCah);

    h_post_PCah = volatility_garch11(draw_PCah,y,y_S,H);
    y_post_PCah = randn(M_fin,H).*sqrt(h_post_PCah);
    y_post_PCah = bsxfun(@plus,y_post_PCah,draw_PCah(:,1));
    y_post_PCah = sort(y_post_PCah);
    
    VaR_1_post_PCah = y_post_PCah(round(p_bar1*M_fin),:); 
    VaR_5_post_PCah = y_post_PCah(round(p_bar*M_fin),:); 
    VaR_05_post_PCah = y_post_PCah(round(p_bar0*M_fin),:); 
 
    ES_1_post_PCah = mean(y_post_PCah(1:round(p_bar1*M_fin),:)); 
    ES_5_post_PCah = mean(y_post_PCah(1:round(p_bar*M_fin),:)); 
    ES_05_post_PCah = mean(y_post_PCah(1:round(p_bar0*M_fin),:));     

    %% CENSORED MLE PARAMETERS
    threshold = 0.1; %<---------- HiGhER?
    quantile = norminv(threshold);
    fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold);
    kernel_init = @(xx) - C_posterior_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S)/T;    
    kernel = @(xx) C_posterior_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S);
    [mu_Cm,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_init,options);
    Sigma_Cm = inv(T*Sigma_Cm);
      
%     cont.mit.CV_tol = 0.3; 
    cont.mit.CV_max = 1.5; %1.9;
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

    h_post_Cm = volatility_garch11(draw_Cm,y,y_S,H);
    y_post_Cm = randn(M,H).*sqrt(h_post_Cm);
    y_post_Cm = bsxfun(@plus,y_post_Cm,draw_Cm(:,1));
    y_post_Cm = sort(y_post_Cm);
    
    VaR_1_post_Cm = y_post_Cm(p_bar1*M,:); 
    VaR_5_post_Cm = y_post_Cm(p_bar*M,:); 
    VaR_05_post_Cm = y_post_Cm(p_bar0*M,:); 

    ES_1_post_Cm = mean(y_post_Cm(1:p_bar1*M,:)); 
    ES_5_post_Cm = mean(y_post_Cm(1:p_bar*M,:)); 
    ES_05_post_Cm = mean(y_post_Cm(1:p_bar0*M,:));     
     
    mean_draw_Cm = mean(draw_Cm);
    median_draw_Cm = median(draw_Cm);
    std_draw_Cm = std(draw_Cm);

    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold);
    % mit_C: joint candidate for the joint censored posterior    
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
    [draw_PCm, a_PCm, lnw_PCm] = sim_cond_mit_MH_outloop(mit_Cm, draw_short, partition, II, BurnIn, kernel, GamMat, cont.disp);
    accept_PCm = mean(a_PCm); 
    
%     ind_fin = isfinite(lnw_PCm);
%     M_fin = sum(ind_fin);
%     draw_PCm = draw_PCm(ind_fin,:) ;  
M_fin = M;
    mean_draw_PCm = mean(draw_PCm);
    median_draw_PCm = median(draw_PCm);
    std_draw_PCm = std(draw_PCm);

    h_post_PCm = volatility_garch11(draw_PCm,y,y_S,H);
    y_post_PCm = randn(M_fin,H).*sqrt(h_post_PCm);
    y_post_PCm = bsxfun(@plus,y_post_PCm,draw_PCm(:,1));
    y_post_PCm = sort(y_post_PCm);
    
    VaR_1_post_PCm = y_post_PCm(round(p_bar1*M_fin),:); 
    VaR_5_post_PCm = y_post_PCm(round(p_bar*M_fin),:); 
    VaR_05_post_PCm = y_post_PCm(round(p_bar0*M_fin),:); 

    ES_1_post_PCm = mean(y_post_PCm(1:round(p_bar1*M_fin),:)); 
    ES_5_post_PCm = mean(y_post_PCm(1:round(p_bar*M_fin),:)); 
    ES_05_post_PCm = mean(y_post_PCm(1:round(p_bar0*M_fin),:)); 
    
     %% Results
     results = struct('y',y,'draw',draw,'draw_Cah',draw_Cah,'draw_PCah',draw_PCah,'draw_Cm',draw_Cm,'draw_PCm',draw_PCm,...
        'q1',q1,'q5',q5,'q05',q05,'cdf1',cdf1,'cdf5',cdf5,'cdf05',cdf05,...
        'mean_draw',mean_draw,'mean_draw_Cah',mean_draw_Cah,'mean_draw_PCah',mean_draw_PCah,'mean_draw_Cm',mean_draw_Cm,'mean_draw_PCm',mean_draw_PCm,...
        'median_draw',median_draw,'median_draw_Cah',median_draw_Cah,'median_draw_PCah',median_draw_PCah,'median_draw_Cm',median_draw_Cm,'median_draw_PCm',median_draw_PCm,...
        'std_draw',std_draw,'std_draw_Cah',std_draw_Cah,'std_draw_PCah',std_draw_PCah,'std_draw_Cm',std_draw_Cm,'std_draw_PCm',std_draw_PCm,...
        'accept',accept,'accept_Cah',accept_Cah,'accept_PCah',accept_PCah,'accept_Cm',accept_Cm,'accept_PCm',accept_PCm,...
        'mit',mit,'CV',CV,'mit_Cah',mit_Cah,'CV_Cah',CV_Cah,'mit_Cm',mit_Cm,'CV_Cm',CV_Cm,...
        'VaR_1',VaR_1,'VaR_1_post',VaR_1_post,'VaR_1_post_Cah',VaR_1_post_Cah,'VaR_1_post_PCah',VaR_1_post_PCah,'VaR_1_post_Cm',VaR_1_post_Cm,'VaR_1_post_PCm',VaR_1_post_PCm,...
        'VaR_5',VaR_5,'VaR_5_post',VaR_5_post,'VaR_5_post_Cah',VaR_5_post_Cah,'VaR_5_post_PCah',VaR_5_post_PCah,'VaR_5_post_Cm',VaR_5_post_Cm,'VaR_5_post_PCm',VaR_5_post_PCm,...
        'VaR_05',VaR_05,'VaR_05_post',VaR_05_post,'VaR_05_post_Cah',VaR_05_post_Cah,'VaR_05_post_PCah',VaR_05_post_PCah,'VaR_05_post_Cm',VaR_05_post_Cm,'VaR_05_post_PCm',VaR_05_post_PCm,...
        'ES_1',ES_1,'ES_1_post',ES_1_post,'ES_1_post_Cah',ES_1_post_Cah,'ES_1_post_PCah',ES_1_post_PCah,'ES_1_post_Cm',ES_1_post_Cm,'ES_1_post_PCm',ES_1_post_PCm,...
        'ES_5',ES_5,'ES_5_post',ES_5_post,'ES_5_post_Cah',ES_5_post_Cah,'ES_5_post_PCah',ES_5_post_PCah,'ES_5_post_Cm',ES_5_post_Cm,'ES_5_post_PCm',ES_5_post_PCm,...
        'ES_05',ES_05,'ES_05_post',ES_05_post,'ES_05_post_Cah',ES_05_post_Cah,'ES_05_post_PCah',ES_05_post_PCah,'ES_05_post_Cm',ES_05_post_Cm,'ES_05_post_PCm',ES_05_post_PCm);


end