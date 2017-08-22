function results = PCP_ar1_run_varc(c, sigma1, sigma2, rho, ...
    p_bar0, p_bar1, p_bar, T, H, M, BurnIn, mu_init, ...
    df, cont, options, partition, II, GamMat)

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
    q05 = norminv(p_bar0,c+rho*y(T:(T+H-1),1),sigma2)';
    q1 = norminv(p_bar1,c+rho*y(T:(T+H-1),1),sigma2)';
    q5 = norminv(p_bar,c+rho*y(T:(T+H-1),1),sigma2)'; 

    % true ESs
    cdf05 = c - sigma2*normpdf(norminv(1-p_bar0))/(p_bar0);
    cdf1 = c - sigma2*normpdf(norminv(1-p_bar1))/(p_bar1);
    cdf5 = c - sigma2*normpdf(norminv(1-p_bar))/(p_bar);    
 
    % MC VaRs under the true model
    eps_sort = randn(M,H);
    ind = (eps_sort>0);
    eps_sort(ind) = c + sigma1.*eps_sort(ind);
    eps_sort(~ind) = c + sigma2.*eps_sort(~ind);

    y_sort = bsxfun(@plus,eps_sort,rho*y(T:(T+H-1),1)');
    y_sort = sort(y_sort);

    VaR_1 = y_sort(p_bar1*M,:); 
    VaR_5 = y_sort(p_bar*M,:); 
    VaR_05 = y_sort(p_bar0*M,:);

    ES_1 = mean(y_sort(1:p_bar1*M,:)); 
    ES_5 = mean(y_sort(1:p_bar*M,:)); 
    ES_05 = mean(y_sort(1:p_bar0*M,:)); 

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
    mean_draw = mean(draw);
    std_draw = std(draw);
    
    y_post = bsxfun(@times,randn(M,H),draw(:,2));
    y_post = bsxfun(@plus,y_post,draw(:,1));
    y_post = y_post + draw(:,3)*y(T:(T+H-1),1)';
    y_post = sort(y_post);

    VaR_1_post = y_post(p_bar1*M,:); 
    VaR_5_post = y_post(p_bar*M,:); 
    VaR_05_post = y_post(p_bar0*M,:); 

    ES_1_post = mean(y_post(1:p_bar1*M,:)); 
    ES_5_post = mean(y_post(1:p_bar*M,:)); 
    ES_05_post = mean(y_post(1:p_bar0*M,:));  

    
    %% Time varying threshold, THR = 1
    %% CENSORED
    fprintf('*** Censored Posterior, ad hoc time varying threshold, THR = 1 ***\n');
%     kernel_init = @(xx) - C_posterior_ar1_varc_mex(xx, y(1:T), threshold)/T; 
%     kernel = @(xx) C_posterior_ar1_varc_mex(xx, y(1:T), threshold); 
%     [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
%  
    % NOPARAMETERS
    threshold = 1;
    kernel_init = @(xx) - C_posterior_ar1_varc_noparam_mex(xx, y(1:T), threshold)/T; 
    kernel = @(xx) C_posterior_ar1_varc_noparam_mex(xx, y(1:T), threshold);     
    [mu_Cah,~,~,~,~,Sigma_Cah] = fminunc(kernel_init,mu_init,options);
    Sigma_Cah = inv(T*Sigma_Cah);
 
%     fprintf('*** Censored Posterior, threshold = 0 ***\n');
%     threshold0 = 0;
%     kernel_init = @(xx) - C_posterior_ar1(xx, y(1:T), threshold0)/T; 
%     kernel = @(xx) C_posterior_ar1(xx, y(1:T), threshold0); 
%     [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu_init,options);
    
    try
        [mit_Cah, CV_Cah] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);   
    catch
%         [mu_Cah,~,~,~,~,Sigma_Cah] = fminunc(kernel_init,mu_init,options);
%         Sigma_Cah = inv(T*Sigma_Cah);
        mit_Cah = struct('mu',mu_Cah,'Sigma',reshape(Sigma_Cah,1,length(mu_Cah)^2),'df', df, 'p', 1);
        CV_Cah = cont.mit.CV_old;
    end
    [draw_Cah, lnk_Cah] = fn_rmvgt_robust(M+BurnIn, mit_Cah, kernel, false);
    lnd_Cah = dmvgt(draw_Cah, mit_Cah, true, GamMat);     
    lnw_Cah = lnk_Cah - lnd_Cah;
    lnw_Cah = lnw_Cah - max(lnw_Cah);
    [ind, a] = fn_MH(lnw_Cah);
    draw_Cah = draw_Cah(ind,:);
    accept_Cah = a/(M+BurnIn);
    
    draw_Cah = draw_Cah(BurnIn+1:BurnIn+M,:);
    mean_draw_Cah = mean(draw_Cah);
    std_draw_Cah = std(draw_Cah);
        
    y_post_Cah = bsxfun(@times,randn(M,H),draw_Cah(:,2));
    y_post_Cah = bsxfun(@plus,y_post_Cah,draw_Cah(:,1));
    y_post_Cah = y_post_Cah + draw_Cah(:,3)*y(T:(T+H-1),1)';
    y_post_Cah = sort(y_post_Cah);

    VaR_1_post_Cah = y_post_Cah(p_bar1*M,:); 
    VaR_5_post_Cah = y_post_Cah(p_bar*M,:); 
    VaR_05_post_Cah = y_post_Cah(p_bar0*M,:); 

    ES_1_post_Cah = mean(y_post_Cah(1:p_bar1*M,:)); 
    ES_5_post_Cah = mean(y_post_Cah(1:p_bar*M,:)); 
    ES_05_post_Cah = mean(y_post_Cah(1:p_bar0*M,:));     
    
    %% PARTIAL CENSORING: keep rho uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior, ad hoc time varying threshold, THR = 1 ***\n');
    % mit_C: joint candidate for the joint censored posterior
    % Short version
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos 
    [draw_PCah, a_PCah] = sim_cond_mit_MH_outloop(mit_Cah, draw_short,...
        partition, II, BurnIn, kernel, GamMat, cont.disp);
    accept_PCah = mean(a_PCah);
    mean_draw_PCah = mean(draw_PCah);
    std_draw_PCah = std(draw_PCah);
  
    y_post_PCah = bsxfun(@times,randn(M,H),draw_PCah(:,2));
    y_post_PCah = bsxfun(@plus,y_post_PCah,draw_PCah(:,1));
    y_post_PCah = y_post_PCah + draw_PCah(:,3)*y(T:(T+H-1),1)';
    y_post_PCah = sort(y_post_PCah);

    VaR_1_post_PCah = y_post_PCah(round(p_bar1*M),:); 
    VaR_5_post_PCah = y_post_PCah(round(p_bar*M),:); 
    VaR_05_post_PCah = y_post_PCah(round(p_bar0*M),:); 
 
    ES_1_post_PCah = mean(y_post_PCah(1:round(p_bar1*M),:)); 
    ES_5_post_PCah = mean(y_post_PCah(1:round(p_bar*M),:)); 
    ES_05_post_PCah = mean(y_post_PCah(1:round(p_bar0*M),:));   

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

    VaR_1_post_Cm = y_post_Cm(p_bar1*M,:); 
    VaR_5_post_Cm = y_post_Cm(p_bar*M,:); 
    VaR_05_post_Cm = y_post_Cm(p_bar0*M,:); 

    ES_1_post_Cm = mean(y_post_Cm(1:p_bar1*M,:)); 
    ES_5_post_Cm = mean(y_post_Cm(1:p_bar*M,:)); 
    ES_05_post_Cm = mean(y_post_Cm(1:p_bar0*M,:));         
    
    
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

    VaR_1_post_PCm = y_post_PCm(round(p_bar1*M),:); 
    VaR_5_post_PCm = y_post_PCm(round(p_bar*M),:); 
    VaR_05_post_PCm = y_post_PCm(round(p_bar0*M),:); 

    ES_1_post_PCm = mean(y_post_PCm(1:round(p_bar1*M),:)); 
    ES_5_post_PCm = mean(y_post_PCm(1:round(p_bar*M),:)); 
    ES_05_post_PCm = mean(y_post_PCm(1:round(p_bar0*M),:)); 
    
    %% Results
    results = struct('y',y,'draw',draw,'draw_Cah',draw_Cah,'draw_PCah',draw_PCah,'draw_Cm',draw_Cm,'draw_PCm',draw_PCm,...
        'q1',q1,'q5',q5,'q05',q05,'cdf1',cdf1,'cdf5',cdf5,'cdf05',cdf05,...
        'mean_draw',mean_draw,'mean_draw_Cah',mean_draw_Cah,'mean_draw_PCah',mean_draw_PCah,'mean_draw_Cm',mean_draw_Cm,'mean_draw_PCm',mean_draw_PCm,...
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