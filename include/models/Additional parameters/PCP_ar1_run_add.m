function results = PCP_ar1_run_add(sdd, c, sigma1, sigma2, rho, p_bar1, p_bar, T, H, M, BurnIn, mu_init, df, cont, options, partition, II, GamMat)
   
    s = RandStream('mt19937ar','Seed',sdd);
    RandStream.setGlobalStream(s); 
    
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

    kernel_init = @(xx) -posterior_ar1_mex(xx,y(1:T))/T;
    kernel = @(xx) posterior_ar1_mex(xx,y(1:T));  
    try
        [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
    catch
        [mu,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init, options);
        Sigma = inv(T*Sigma);
        mit = struct('mu',mu,'Sigma',reshape(Sigma,1,length(mu)^2),'df', df, 'p', 1);
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
    
    %% ADDITIONAL PARAMS
    % Threshold = 10% perscentile of the data sample
    fprintf('*** Additional Parameters ***\n');    
    threshold = sort(y(1:T));
    threshold = threshold(2*p_bar*T);
    mu_add_init = [0, 1];
    mu_add = zeros(M,2);
    Sigma_add = zeros(M,2*2);    
    for ii = 1:M
        if (mod(ii,1000)==0)
            fprintf('Add param iter = %d\n',ii)
        end
        kernel_init = @(aa) - C_addparam_posterior_ar1_mex(aa, draw(ii,:), y(1:T), threshold)/T;   
        [mu_add(ii,:),~,~,~,~,S_add] = fminunc(kernel_init, mu_add_init,options);
        Sigma_add(ii,:) = reshape(inv(T*S_add),1,4);
    end
    y_post_add = bsxfun(@times, randn(M,H), mu_add(:,2).*draw(:,2));
    y_post_add = bsxfun(@plus, y_post_add, mu_add(:,1) + draw(:,1));
    y_post_add = y_post_add + draw(:,3)*y(T:(T+H-1),1)';
    y_post_add = sort(y_post_add);
    VaR_1_post_add = y_post_add(p_bar1*M,:); 
    VaR_5_post_add = y_post_add(p_bar*M,:); 

    mean_mu_add = mean(mu_add);
    median_mu_add = median(mu_add);
    std_mu_add = std(mu_add);
     
    %% Threshold = 10% perscentile of the data sample
    threshold = sort(y(1:T));
    threshold = threshold(2*p_bar*T);
    %% CENSORED
    fprintf('*** Censored Posterior, threshold 10%% ***\n');
%     kernel_init = @(xx) - C_posterior_ar1(xx, y(1:T), threshold)/T;    
%     kernel_init = @(xx) - C_loglik_ar1(xx, y, threshold)/T;
    kernel_init = @(xx) - C_posterior_ar1_mex(xx, y(1:T), threshold)/T;    
    kernel = @(xx) C_posterior_ar1_mex(xx, y(1:T), threshold);
    try
        [mit_C, CV_C] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);   
    catch
        [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
        Sigma_C = inv(T*Sigma_C);
        mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
        CV = cont.mit.CV_old;
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
    fprintf('*** Partially Censored Posterior, threshold 10%% ***\n');
    % mit_C: joint candidate for the joint censored posterior
    % Short version
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos 
%     M_short = M/II;
%     [draw_PC, a_PC] = sim_cond_mit_MH(mit_C, draw_short, partition, M_short, BurnIn, kernel, GamMat);
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
    
    %% Threshold = 0
    threshold0 = 0;
    %% CENSORED
    fprintf('*** Censored Posterior, threshold 0 ***\n');    
    kernel_init = @(xx) - C_posterior_ar1_mex(xx, y(1:T), threshold0)/T; 
    kernel = @(xx) C_posterior_ar1_mex(xx, y(1:T), threshold0);
    try
        [mit_C0, CV_C0] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
    catch
        [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu_init,options);
        Sigma_C0 = inv(T*Sigma_C0);
        mit_C0 = struct('mu',mu_C0,'Sigma',reshape(Sigma_C0,1,length(mu_C0)^2),'df', df, 'p', 1);
        [mit_C0, CV_C0] = MitISEM_new(mit_C0, kernel, mu_init, cont, GamMat);
    end
    [draw_C0, lnk_C0] = fn_rmvgt_robust(M+BurnIn, mit_C0, kernel, false);
    lnd_C0 = dmvgt(draw_C0, mit_C0, true, GamMat);    
    lnw_C0 = lnk_C0 - lnd_C0;
    lnw_C0 = lnw_C0 - max(lnw_C0);
    [ind, a] = fn_MH(lnw_C0);
    draw_C0 = draw_C0(ind,:);
    accept_C0 = a/(M+BurnIn);
    draw_C0 = draw_C0(BurnIn+1:BurnIn+M,:);
    mean_draw_C0 = mean(draw_C0);
    std_draw_C0 = std(draw_C0);

    y_post_C0 = bsxfun(@times,randn(M,H),draw_C0(:,2));
    y_post_C0 = bsxfun(@plus,y_post_C0,draw_C0(:,1));
    y_post_C0 = y_post_C0 + draw_C0(:,3)*y(T:(T+H-1),1)';
    y_post_C0 = sort(y_post_C0);
    VaR_1_post_C0 = y_post_C0(p_bar1*M,:); 
    VaR_5_post_C0 = y_post_C0(p_bar*M,:); 

    %% PARTIAL CENSORING: keep rho uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior, threshold 0 ***\n');
    % mit_C0: joint candidate for the joint censored posterior
    % Short version
%     [draw_PC0, a_PC0] = sim_cond_mit_MH(mit_C0, draw_short, partition, M_short, BurnIn, kernel, GamMat);
    [draw_PC0, a_PC0] = sim_cond_mit_MH_outloop(mit_C0, draw_short, partition, II, BurnIn, kernel, GamMat, cont.disp);
    accept_PC0 = mean(a_PC0);
    mean_draw_PC0 = mean(draw_PC0);
    std_draw_PC0 = std(draw_PC0);    

    y_post_PC0 = bsxfun(@times,randn(M,H),draw_PC0(:,2));
    y_post_PC0 = bsxfun(@plus,y_post_PC0,draw_PC0(:,1));
    y_post_PC0 = y_post_PC0 + draw_PC0(:,3)*y(T:(T+H-1),1)';
    y_post_PC0 = sort(y_post_PC0);
    VaR_1_post_PC0 = y_post_PC0(p_bar1*M,:); 
    VaR_5_post_PC0 = y_post_PC0(p_bar*M,:); 
    
  
    %% Results
    results = struct('y',y,'draw',draw,'draw_C',draw_C,'draw_PC',draw_PC,'draw_C0',draw_C0,'draw_PC0',draw_PC0,...
        'q1',q1,'q5',q5,...
        'mean_draw',mean_draw,'mean_draw_C',mean_draw_C,'mean_draw_PC',mean_draw_PC,'mean_draw_C0',mean_draw_C0,'mean_draw_PC0',mean_draw_PC0,...
        'std_draw',std_draw,'std_draw_C',std_draw_C,'std_draw_PC',std_draw_PC,'std_draw_C0',std_draw_C0,'std_draw_PC0',std_draw_PC0,...
        'accept',accept,'accept_C',accept_C,'accept_PC',accept_PC,'accept_C0',accept_C0,'accept_PC0',accept_PC0,...
        'mit',mit,'CV',CV,'mit_C',mit_C,'CV_C',CV_C,'mit_C0',mit_C0,'CV_C0',CV_C0,...
        'VaR_1',VaR_1,'VaR_1_post',VaR_1_post,'VaR_1_post_C',VaR_1_post_C,'VaR_1_post_PC',VaR_1_post_PC,'VaR_1_post_C0',VaR_1_post_C0,'VaR_1_post_PC0',VaR_1_post_PC0,...
        'VaR_5',VaR_5,'VaR_5_post',VaR_5_post,'VaR_5_post_C',VaR_5_post_C,'VaR_5_post_PC',VaR_5_post_PC,'VaR_5_post_C0',VaR_5_post_C0,'VaR_5_post_PC0',VaR_5_post_PC0,...
        'VaR_1_post_add',VaR_1_post_add,'VaR_5_post_add',VaR_5_post_add,'mean_mu_add',mean_mu_add,'median_mu_add',median_mu_add,'std_mu_add',std_mu_add);
end