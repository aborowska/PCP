function results = PCP_t_garch11_empirical(sdd, data, p_bar0, p_bar1, p_bar, T, H, M, BurnIn, ...
    mu_init, df, cont, options, partition, II, GamMat)
   
    addpath(genpath('include/'));
%  theta[i] <= nu
%  theta[i+N] <= mu
%  theta[i+2*N] <= omega
%  theta[i+3*N] <= alpha
%  theta[i+4*N] <= beta 
    mu_init = [8,0,1,0.1,0.8];
    x_gam = (0:0.00001:50)'+0.00001;
    GamMat = gamma(x_gam);
    
    M = 10000; % number of draws 
    BurnIn = 1000;
    
    cont = MitISEM_Control;
    cont.mit.CV_max = 1; %2?
    cont.mit.iter_max = 10;
    cont.mit.Hmax = 6;
    cont.mit.dfnc = 5;
    df = 5; % default df for a mit
    
    options = optimset('Display','off');
    % w = warning('query','last');
    % id = w.identifier;
    id = 'optim:fminunc:SwitchingMethod';
    warning('off',id);
    
    
    
    s = RandStream('mt19937ar','Seed',sdd);
    RandStream.setGlobalStream(s); 
        
    thinning = 10;
    BurnIn_PCP = BurnIn/10;
 
    save_on = true;
    %% Load data
    if (data == 1)
        y = load('GSPC.txt'); % 15-04-1996 : 05-10-2015
        y = 100*diff(log(y));
        TT = length(y);
        H = 2000;
        T = TT - H;
        H = 100;
        y_plot = y;
        y = y(1:(T+H));
%         time = linspace((1996 + 3.5/12),(2015 + (9 + 1/6)/12),TT);
        time = [(1996 + 3.5/12),(2015 + (9 + 1/6)/12)];
        figures_path = 'figures/tgarch11_emp/';
    else
        y = load('GSPC.txt');
        y = 100*diff(log(y));
    end
    data_stats = Plot_data(y_plot,time,save_on,figures_path);
    
    %% Uncensored Posterior
    fprintf('*** Uncensored Posterior ***\n');
    y_S = var(y(1:T));
%         mu_init(1,1) = mean(y(1:T));
%         mu_init(1,2) = y_S;
    hyper = 0.01; 
    kernel_init = @(xx) - posterior_t_garch11_mex(xx, y(1:T), y_S, GamMat, hyper)/T;
    kernel = @(xx) posterior_t_garch11_mex(xx, y(1:T), y_S, GamMat, hyper);
    [mu,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init,options);
%  7.7131    0.0625    1.0133    0.0362    0.9592


%     kernel_init = @(xx) -posterior_t_garch11_mex(xx, y(1:T), y_S)/T;
%     kernel = @(xx) posterior_garch11_mex(xx, y(1:T), y_S);
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

    h_post = volatility_t_garch11(draw,y,y_S,H);  
    rho = (draw(:,1)-2)./draw(:,1);
    eps_H = trnd(repmat(draw(:,1),1,H));
 
    h_post = bsxfun(@times,h_post,rho);
    y_post = eps_H.*sqrt(h_post);
    y_post = bsxfun(@plus,y_post,draw(:,2));
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

    %% Threshold = 10% perscentile of the data sample
    threshold = sort(y(1:T));
    threshold = threshold(round(2*p_bar*T));
    %% CENSORED
    fprintf('*** Censored Posterior, threshold 10%% ***\n');
    kernel_init = @(xx) - C_posterior_t_garch11_mex(xx, y(1:T,1), threshold, y_S, GamMat, hyper)/T;    
    kernel = @(xx) C_posterior_t_garch11_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);
    
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
            [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
            mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
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

    h_post_C = volatility_garch11(draw_C,y,y_S,H);
    y_post_C = randn(M,H).*sqrt(h_post_C);
    y_post_C = bsxfun(@plus,y_post_C,draw_C(:,1));
    y_post_C = sort(y_post_C);
    
    VaR_1_post_C = y_post_C(p_bar1*M,:); 
    VaR_5_post_C = y_post_C(p_bar*M,:); 
    VaR_05_post_C = y_post_C(p_bar0*M,:); 

    ES_1_post_C = mean(y_post_C(1:p_bar1*M,:)); 
    ES_5_post_C = mean(y_post_C(1:p_bar*M,:)); 
    ES_05_post_C = mean(y_post_C(1:p_bar0*M,:)); 
    
    mean_draw_C = mean(draw_C);
    median_draw_C = median(draw_C);
    std_draw_C = std(draw_C);

    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior, threshold 10%% ***\n');
    % mit_C: joint candidate for the joint censored posterior
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
%     M_short = M/II;
%     [draw_PC, a_PC] = sim_cond_mit_MH(mit_C, draw_short, partition, M_short, BurnIn, kernel, GamMat);
%     [draw_PC, a_PC] = sim_cond_mit_MH(mit_C, draw_short, partition, II, BurnIn, kernel, GamMat);
    [draw_PC, a_PC, lnw_PC] = sim_cond_mit_MH_outloop(mit_C, draw_short,...
        partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
    accept_PC = mean(a_PC); 
%     ind_fin = isfinite(lnw_PC);
%     M_fin = sum(ind_fin);
%     draw_PC = draw_PC(ind_fin,:) ;  
M_fin = M;
    mean_draw_PC = mean(draw_PC);
    median_draw_PC = median(draw_PC);
    std_draw_PC = std(draw_PC);

    h_post_PC = volatility_garch11(draw_PC,y,y_S,H);
%     y_post_PC = bsxfun(@times,randn(M_fin,H),sqrt(h_post_PC(T+1:T+H,1))');
    y_post_PC = randn(M_fin,H).*sqrt(h_post_PC);
    y_post_PC = bsxfun(@plus,y_post_PC,draw_PC(:,1));
    y_post_PC = sort(y_post_PC);

    VaR_1_post_PC = y_post_PC(round(p_bar1*M_fin),:); 
    VaR_5_post_PC = y_post_PC(round(p_bar*M_fin),:); 
    VaR_05_post_PC = y_post_PC(round(p_bar0*M_fin),:); 

    ES_1_post_PC = mean(y_post_PC(1:round(p_bar1*M_fin),:)); 
    ES_5_post_PC = mean(y_post_PC(1:round(p_bar*M_fin),:)); 
    ES_05_post_PC = mean(y_post_PC(1:round(p_bar0*M_fin),:));
    
    %% Threshold = 0
    threshold0 = 0;
    %% CENSORED
    fprintf('*** Censored Posterior, threshold 0 ***\n');    
    kernel_init = @(xx) - C_posterior_garch11_mex(xx, y(1:T,1), threshold0, y_S)/T;    
    kernel = @(xx) C_posterior_garch11_mex(xx, y(1:T,1), threshold0, y_S);
    
    CV_C0 = cont.mit.CV_old;
    while (CV_C0(end) >= 2)
        try
            [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu_init,options);
            Sigma_C0 = inv(T*Sigma_C0);
            mit_C0 = struct('mu',mu_C0,'Sigma',reshape(Sigma_C0,1,length(mu_C0)^2),'df', df, 'p', 1);
            draw_C0 = rmvt(mu_C0,Sigma_C0,df,M+BurnIn);
            [mit_C0, CV_C0] = MitISEM_new2(mit_C0, kernel, mu_init, cont, GamMat);   
            [draw_C0, lnk_C0] = fn_rmvgt_robust(M+BurnIn, mit_C0, kernel, false);
            lnd_C0 = dmvgt(draw_C0, mit_C0, true, GamMat);            
        catch
            mu_C0 = fminunc(kernel_init,mu_init,options);
            [~, ind_aux] = max(mit_C.p);
            Sigma_aux = mit_C.Sigma(ind_aux,:);
            mit_C0 = struct('mu',mu_C0,'Sigma',Sigma_aux,'df', df, 'p', 1);
            [mit_C0, CV_C0] = MitISEM_new2(mit_C0, kernel, mu_init, cont, GamMat);   
            if CV_C0(end)>2
                [mit_C0, CV_C0] = MitISEM_new2(mit_C0, kernel, mu_init, cont, GamMat);   
            end
            [draw_C0, lnk_C0] = fn_rmvgt_robust(M+BurnIn, mit_C0, kernel, false);
            lnd_C0 = dmvgt(draw_C0, mit_C0, true, GamMat);                
        end
    end
    
    lnw_C0 = lnk_C0 - lnd_C0;
    lnw_C0 = lnw_C0 - max(lnw_C0);
    [ind, a] = fn_MH(lnw_C0);
    draw_C0 = draw_C0(ind,:);
    accept_C0 = a/(M+BurnIn);
    draw_C0 = draw_C0(BurnIn+1:BurnIn+M,:);
    mean_draw_C0 = mean(draw_C0);
    median_draw_C0 = median(draw_C0);
    std_draw_C0 = std(draw_C0);

    h_post_C0 = volatility_garch11(draw_C0,y,y_S,H);
    y_post_C0 = randn(M,H).*sqrt(h_post_C0);
    y_post_C0 = bsxfun(@plus,y_post_C0,draw_C0(:,1));
    y_post_C0 = sort(y_post_C0);

    VaR_1_post_C0 = y_post_C0(p_bar1*M,:); 
    VaR_5_post_C0 = y_post_C0(p_bar*M,:); 
    VaR_05_post_C0 = y_post_C0(p_bar0*M,:); 

    ES_1_post_C0 = mean(y_post_C0(1:p_bar1*M,:)); 
    ES_5_post_C0 = mean(y_post_C0(1:p_bar*M,:)); 
    ES_05_post_C0 = mean(y_post_C0(1:p_bar0*M,:)); 
    
    %% PARTIAL CENSORING: keep alpha and beta uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior, threshold 0 ***\n');
    % mit_C0: joint cnadidate for the joint censored posterior
%     [draw_PC0, a_PC0] = sim_cond_mit_MH(mit_C0, draw_short, partition, M_short, BurnIn, kernel, GamMat);
    [draw_PC0, a_PC0, lnw_PC0] = sim_cond_mit_MH_outloop(mit_C0, draw_short,...
        partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
    accept_PC0 = mean(a_PC0);
%     ind_fin = isfinite(lnw_PC0);
%     M_fin = sum(ind_fin);
%     draw_PC0 = draw_PC0(ind_fin,:);
M_fin = M;    
    mean_draw_PC0 = mean(draw_PC0);
    median_draw_PC0 = median(draw_PC0);
    std_draw_PC0 = std(draw_PC0);    

    h_post_PC0 = volatility_garch11(draw_PC0,y,y_S,H);
    y_post_PC0 = randn(M_fin,H).*sqrt(h_post_PC0);
    y_post_PC0 = bsxfun(@plus,y_post_PC0,draw_PC0(:,1));
    y_post_PC0 = sort(y_post_PC0);
    
    VaR_1_post_PC0 = y_post_PC0(round(p_bar1*M_fin),:); 
    VaR_5_post_PC0 = y_post_PC0(round(p_bar*M_fin),:); 
    VaR_05_post_PC0 = y_post_PC0(round(p_bar0*M_fin),:); 
  
    ES_1_post_PC0 = mean(y_post_PC0(1:round(p_bar1*M_fin),:)); 
    ES_5_post_PC0 = mean(y_post_PC0(1:round(p_bar*M_fin),:)); 
    ES_05_post_PC0 = mean(y_post_PC0(1:round(p_bar0*M_fin),:)); 
    
    %% Results   
    results = struct('y',y,'draw',draw,'draw_C',draw_C,'draw_PC',draw_PC,'draw_C0',draw_C0,'draw_PC0',draw_PC0,...
        'q1',q1,'q5',q5,'q05',q05,'cdf1',cdf1,'cdf5',cdf5,'cdf05',cdf05,...
    'mean_draw',mean_draw,'mean_draw_C',mean_draw_C,'mean_draw_PC',mean_draw_PC,'mean_draw_C0',mean_draw_C0,'mean_draw_PC0',mean_draw_PC0,...
    'median_draw',median_draw,'median_draw_C',median_draw_C,'median_draw_PC',median_draw_PC,'median_draw_C0',median_draw_C0,'median_draw_PC0',median_draw_PC0,...
    'std_draw',std_draw,'std_draw_C',std_draw_C,'std_draw_PC',std_draw_PC,'std_draw_C0',std_draw_C0,'std_draw_PC0',std_draw_PC0,...
    'accept',accept,'accept_C',accept_C,'accept_PC',accept_PC,'accept_C0',accept_C0,'accept_PC0',accept_PC0,...
    'mit',mit,'CV',CV,'mit_C',mit_C,'CV_C',CV_C,'mit_C0',mit_C0,'CV_C0',CV_C0,...
    'VaR_1',VaR_1,'VaR_1_post',VaR_1_post,'VaR_1_post_C',VaR_1_post_C,'VaR_1_post_PC',VaR_1_post_PC,'VaR_1_post_C0',VaR_1_post_C0,'VaR_1_post_PC0',VaR_1_post_PC0,...
    'VaR_5',VaR_5,'VaR_5_post',VaR_5_post,'VaR_5_post_C',VaR_5_post_C,'VaR_5_post_PC',VaR_5_post_PC,'VaR_5_post_C0',VaR_5_post_C0,'VaR_5_post_PC0',VaR_5_post_PC0,...
    'VaR_05',VaR_05,'VaR_05_post',VaR_05_post,'VaR_05_post_C',VaR_05_post_C,'VaR_05_post_PC',VaR_05_post_PC,'VaR_05_post_C0',VaR_05_post_C0,'VaR_05_post_PC0',VaR_05_post_PC0,...
    'ES_1',ES_1,'ES_1_post',ES_1_post,'ES_1_post_C',ES_1_post_C,'ES_1_post_PC',ES_1_post_PC,'ES_1_post_C0',ES_1_post_C0,'ES_1_post_PC0',ES_1_post_PC0,...
    'ES_5',ES_5,'ES_5_post',ES_5_post,'ES_5_post_C',ES_5_post_C,'ES_5_post_PC',ES_5_post_PC,'ES_5_post_C0',ES_5_post_C0,'ES_5_post_PC0',ES_5_post_PC0,...
    'ES_05',ES_05,'ES_05_post',ES_05_post,'ES_05_post_C',ES_05_post_C,'ES_05_post_PC',ES_05_post_PC,'ES_05_post_C0',ES_05_post_C0,'ES_05_post_PC0',ES_05_post_PC0);
end