% function results = PCP_t_garch11_empirical(sdd, data, p_bar0, p_bar1, p_bar, T, H, M, BurnIn, ...
%     mu_init, df, cont, options, partition, II, GamMat)
% addpath(genpath('include/'));
% 
% clear all
% sdd = 1;
% plot_on = false;
% save_on = false; %true;

    if arg0
        model = 'tgarch11';
        parameters = {'$\\nu$','$\\mu$','$\\omega$','$\\alpha$','$\\beta$'};
        params = {'$\nu$','$\mu$','$\omega$','$\alpha$','$\beta$'};

        % quantiles of interest
        p_bar0 = 0.005;
        p_bar1 = 0.01;
        p_bar = 0.05; 
        P_bars = [p_bar0, p_bar1, p_bar];

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
        cont.mit.Hmax = 3; %6;
        cont.mit.dfnc = 3;%5;
    %     df = 5; % default df for a mit
        df = 3; % default df for a mit
        cont.disp = true;

        options = optimset('Display','off');
        % w = warning('query','last');
        % id = w.identifier;
        id = 'optim:fminunc:SwitchingMethod';
        warning('off',id);

%         thinning = 10;
        BurnIn_PCP = BurnIn/10; 
        II = 10;
        thinning = 1;
        
        %% Load data
    %     data = 2;
    %     H =2000;
        [y, T, y_plot, data_name, time] = Load_data_empirical(data, H);   
        y_S = var(y(1:T));
        
        y_sort = sort(y(1:T));
        THR_emp = y_sort(round(P_bars*T));
        
        hyper = 0.01; 
        
        fprintf('Data loaded: %s\n',data_name);

        figures_path = ['figures/',model,'/',data_name,'/']; 

        if plot_on
            data_stats = Plot_data(y_plot,H,time,save_on,figures_path);
        end

        name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'.mat'];
    end  
        
    %% Uncensored Posterior
    if arg1
        fprintf('*** Uncensored Posterior ***\n');
 
    
%     if exist(name,'file')
%         exits_on = true;
%         load(name,'mit','CV','draw','accept','mu_MLE','Sigma')
%         hT = volatility_t_garch11(draw,y(1:T),y_S,0);  
%     else
%         exits_on = false;    
    %     kernel_init = @(xx) - posterior_t_garch11_mex(xx, y(1:T), y_S, GamMat, hyper)/T;
        kernel_init = @(xx) - posterior_t_garch11(xx, y(1:T), y_S, hyper)/T;
        kernel = @(xx) posterior_t_garch11_mex(xx, y(1:T), y_S, GamMat, hyper);
        [mu_MLE1,~,exitflag,output,~,Sigma] = fminunc(kernel_init,mu_init,options);
        [mu_MLE,~,exitflag,output,~,Sigma] = fminunc(kernel_init,mu_MLE1,options);

        Sigma = inv(T*Sigma);    
        if (data == 4)
            Sigma = csvread('MSFT_Sigma_mle.csv',1,1)
        end
    % SZ SnP:  8.2280    0.0262    0.3328    0.0682    0.9270 
    %  ==> omega_SZ = 0.3328*(1 - 0.0682  -  0.9270) =    0.0016
    % new SnP: 9.1468    0.0485    3.0129    0.0732    0.9242
    % AAPL:   7.9999    0.0021    1.0009    0.1478    0.8515
    % IBM:   5.7431    0.0381    6.1932    0.0512    0.9470
    % MSFT: mu_MLE = [7.9996    0.0029    1.0028    0.1329    0.8671]
    %       mu_MLE2 = [ 5.3233    0.0164   66.2934    0.0559    0.9440]

    %     kernel_init = @(xx) -posterior_t_garch11_mex(xx, y(1:T), y_S)/T;
    %     kernel = @(xx) posterior_garch11_mex(xx, y(1:T), y_S);

    %     cont.mit.CV_tol = 0.3; 
        cont.mit.CV_max = 1.9;
        if (data == 2)
            cont.mit.CV_max = 2.1;
            cont.mit.dfnc = 4;
        end  
        
        CV = cont.mit.CV_old;
        while (CV(end) >= 2)
            try
    %             [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
                [mit, CV] = MitISEM_new2(kernel_init, kernel, mu_init, cont, GamMat);
                [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
                lnd = dmvgt(draw, mit, true, GamMat); 
            catch
                draw = rmvt(mu_MLE,Sigma,df,M+BurnIn);
                mit = struct('mu',mu_MLE,'Sigma',reshape(Sigma,1,length(mu_MLE)^2),'df', df, 'p', 1);
%                 [mit, CV] = MitISEM_new(mit, kernel, mu_init, cont, GamMat);            
                [mit, CV] = MitISEM_new2(mit, kernel, mu_init, cont, GamMat);            
                [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
                lnd = dmvgt(draw, mit, true, GamMat);
            end
        end
        
        lnw = lnk - lnd;
        lnw = lnw - max(lnw);
        [ind, a] = fn_MH(lnw);
        draw = draw(ind,:);
        accept = a/(M+BurnIn);
        draw = draw(BurnIn+1:BurnIn+M,:);    

        % compute the implied volatility for the last in-sample period
        hT = volatility_t_garch11(draw,y(1:T),y_S,0);  

        mean_draw = mean(draw);
        median_draw = median(draw);
        std_draw = std(draw);    
%     end 

        % predictive densities
        dens_post = predictive_dens_t_garch11(y(T:(T+H)), hT, draw);
         % predicitve cdfs, constant threshold for different tails
        cdf_post = predictive_cdf_t_garch11(y(T:(T+H)), hT, draw, THR_emp);
        C_score_post_05 = C_ScoringRule(dens_post, cdf_post(1,:), y((T+1):(T+H)), THR_emp(1));
        C_score_post_1 = C_ScoringRule(dens_post, cdf_post(2,:), y((T+1):(T+H)), THR_emp(2));
        C_score_post_5 = C_ScoringRule(dens_post, cdf_post(3,:), y((T+1):(T+H)), THR_emp(3));
    end
    
    %% CENSORED: Threshold = 10% perscentile of the data sample
    if arg2
        threshold = sort(y(1:T));
        threshold = threshold(round(2*p_bar*T));
        fprintf('*** Censored Posterior, threshold 10%% ***\n');
    %     kernel_init = @(xx) - C_posterior_t_garch11_mex(xx, y(1:T,1), threshold, y_S, GamMat, hyper)/T;    
        kernel_init = @(xx) - C_posterior_t_garch11_2(xx, y(1:T,1), threshold, y_S, hyper)/T;    
        kernel = @(xx) C_posterior_t_garch11_2_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);
        [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
        [mu_C2,~,~,~,~,Sigma_C2] = fminunc(kernel_init,mu_C,options);
        Sigma_C = inv(T*Sigma_C);
    %  SZ SnP:     9.5292    0.0673    1.9216    0.0583    0.9400
    % IBM:     7.9994   -0.0045    1.0022    0.1434    0.8566

        r = fn_testSigma(reshape(Sigma_C,1,5^2)); % if r==1 then there is a problem with Sigma_C
        if r
            Sigma_start = Sigma;        
        end

    %     cont.mit.CV_tol = 0.3; 
        cont.mit.CV_max = 1.9;
        CV_C = cont.mit.CV_old;
        while (CV_C(end) >= 2)
            try
                mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_start,1,length(mu_C)^2),'df', df, 'p', 1);
                [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
                if CV_C(end)>2
                    [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
                end
                [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
                lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
            catch
                mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_start,1,length(mu_C)^2),'df', df, 'p', 1);
                [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
                lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
            end 
        end   

        lnw_C = lnk_C - lnd_C;
        lnw_C = lnw_C - max(lnw_C);
        [ind, a] = fn_MH(lnw_C);
        accept_C = a/(M+BurnIn);
        draw_C = draw_C(BurnIn+1:BurnIn+M,:);

        hT_C = volatility_t_garch11(draw_C,y(1:T),y_S,0);  

        mean_draw_C = mean(draw_C);
        median_draw_C = median(draw_C);
        std_draw_C = std(draw_C);
        
        % predictive densities
        dens_C = predictive_dens_t_garch11(y(T:(T+H)), hT_C, draw_C);
        % predicitve cdfs, constant threshold for different tails
        cdf_C = predictive_cdf_t_garch11(y(T:(T+H)), hT_C, draw_C, THR_emp);       
        C_score_C_05 = C_ScoringRule(dens_C, cdf_C(1,:), y((T+1):(T+H)), THR_emp(1));
        C_score_C_1 = C_ScoringRule(dens_C, cdf_C(2,:), y((T+1):(T+H)), THR_emp(2));
        C_score_C_5 = C_ScoringRule(dens_C, cdf_C(3,:), y((T+1):(T+H)), THR_emp(3));    
    end
    
    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor nu, mu and sigma
    if arg3
        fprintf('*** Partially Censored Posterior, threshold 10%%, partition = 4 ***\n');
        % mit_C: joint candidate for the joint censored posterior
        partition = 4;
        draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
        [draw_PC, a_PC, lnw_PC] = sim_cond_mit_MH_outloop(mit_C, draw_short,...
            partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
        accept_PC = mean(a_PC); 

        hT_PC = volatility_t_garch11(draw_PC,y(1:T),y_S,0);  

        mean_draw_PC = mean(draw_PC);
        median_draw_PC = median(draw_PC);
        std_draw_PC = std(draw_PC);

        % predictive densities
        dens_PC = predictive_dens_t_garch11(y(T:(T+H)), hT_PC, draw_PC);      
        % predicitve cdfs, constant threshold for different tails
        cdf_PC = predictive_cdf_t_garch11(y(T:(T+H)), hT_PC, draw_PC, THR_emp);
        C_score_PC_05 = C_ScoringRule(dens_PC, cdf_PC(1,:), y((T+1):(T+H)), THR_emp(1));
        C_score_PC_1 = C_ScoringRule(dens_PC, cdf_PC(2,:), y((T+1):(T+H)), THR_emp(2));
        C_score_PC_5 = C_ScoringRule(dens_PC, cdf_PC(3,:), y((T+1):(T+H)), THR_emp(3));       
    end
    
    %% PARTIALLY CENSORED 2: keep omega, alpha and beta uncensored, then censor nu and mu
    if arg4
        fprintf('*** Partially Censored Posterior, threshold 10%%, partition = 3 ***\n');
        % mit_C: joint candidate for the joint censored posterior
        partition2 = 3;
        [draw_PC2, a_PC2, lnw_PC2] = sim_cond_mit_MH_outloop(mit_C, draw_short,...
            partition2, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
        accept_PC2 = mean(a_PC2); 

        hT_PC2 = volatility_t_garch11(draw_PC2,y(1:T),y_S,0);  

        mean_draw_PC2 = mean(draw_PC2);
        median_draw_PC2 = median(draw_PC2);
        std_draw_PC2 = std(draw_PC2);
        
        % predictive densities        
        dens_PC2 = predictive_dens_t_garch11(y(T:(T+H)), hT_PC2, draw_PC2);
        % predicitve cdfs, constant threshold for different tails
        cdf_PC2 = predictive_cdf_t_garch11(y(T:(T+H)), hT_PC2, draw_PC2, THR_emp);
        C_score_PC2_05 = C_ScoringRule(dens_PC2, cdf_PC2(1,:), y((T+1):(T+H)), THR_emp(1));
        C_score_PC2_1 = C_ScoringRule(dens_PC2, cdf_PC2(2,:), y((T+1):(T+H)), THR_emp(2));
        C_score_PC2_5 = C_ScoringRule(dens_PC2, cdf_PC2(3,:), y((T+1):(T+H)), THR_emp(3));                 
    end
    
    %%  PARTIALLY CENSORED based on a grid
    if arg5
        fprintf('*** Partially Censored Posterior, threshold 10%%, grid***\n');

        grid = 2.1:0.1:20;
        % N = 180;
        % grid = linspace(-1.5,1.5,180)
        kernel = @(xx) C_posterior_t_garch11_2_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);
        draw_PCg = sim_cond_inv_trans(draw, grid, kernel);

        mean_draw_PCg = mean(draw_PCg);
        median_draw_PCg = median(draw_PCg);
        std_draw_PCg = std(draw_PCg);

        % compute the implied volatility for the last in-sample period
        hT_PCg = volatility_t_garch11(draw_PCg,y(1:T),y_S,0);        
    end
    
    %% CENSORED MLE PARAMETERS
    if arg6
        if (arg1 == 0)
            load(name,'mu_MLE','Sigma');
        end
%       
        threshold_m = 0.1; %<---------- HiGhER?
        quantile = tinv(threshold_m, mu_MLE(1));
        fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold_m);
        kernel_init = @(xx) - C_posterior_t_garch11_varc_mle(xx, y(1:T,1), mu_MLE, quantile, y_S, hyper)/T;    
        kernel = @(xx) C_posterior_t_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S, GamMat, hyper);

        [mu_Cm2,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_init,options);
        [mu_Cm,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_Cm2,options);
        Sigma_Cm = inv(T*Sigma_Cm);

    %     cont.mit.CV_tol = 0.3; 
        cont.mit.CV_max = 2.1; %1.5; %1.9;
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

        
        mean_draw_Cm = mean(draw_Cm);
        median_draw_Cm = median(draw_Cm);
        std_draw_Cm = std(draw_Cm);
     
     
        % compute the implied volatility for the last in-sample period
        hT_Cm = volatility_t_garch11(draw_Cm,y(1:T),y_S,0);  

        % predictive densities        
        dens_Cm = predictive_dens_t_garch11(y(T:(T+H)), hT_Cm, draw_Cm);
        % predicitve cdfs, constant threshold for different tails
        cdf_Cm = predictive_cdf_t_garch11(y(T:(T+H)), hT_Cm, draw_Cm, THR_emp);
        C_score_Cm_05 = C_ScoringRule(dens_Cm, cdf_Cm(1,:), y((T+1):(T+H)), THR_emp(1));
        C_score_Cm_1 = C_ScoringRule(dens_Cm, cdf_Cm(2,:), y((T+1):(T+H)), THR_emp(2));
        C_score_Cm_5 = C_ScoringRule(dens_Cm, cdf_Cm(3,:), y((T+1):(T+H)), THR_emp(3));    
    end
    
    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor nu, mu and sigma
    if arg7
        if (arg1 == 0)
            load(name,'mu_MLE','draw');
        end
        if (arg6 == 0)
            load(name,'mit_Cm');
        end        
        fprintf('*** Partially Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold_m);
        % mit_C: joint candidate for the joint censored posterior    
        draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
        partition = 4; 
        [draw_PCm, a_PCm, lnw_PCm] = sim_cond_mit_MH_outloop(mit_Cm, draw_short, ...
            partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
        accept_PCm = mean(a_PCm); 
        
        mean_draw_PCm = mean(draw_PCm);
        median_draw_PCm = median(draw_PCm);
        std_draw_PCm = std(draw_PCm);
     
        % compute the implied volatility for the last in-sample period
        hT_PCm = volatility_t_garch11(draw_PCm,y(1:T),y_S,0);  

        % predictive densities        
        dens_PCm = predictive_dens_t_garch11(y(T:(T+H)), hT_PCm, draw_PCm);
        % predicitve cdfs, constant threshold for different tails
        cdf_PCm = predictive_cdf_t_garch11(y(T:(T+H)), hT_PCm, draw_PCm, THR_emp);
        C_score_PCm_05 = C_ScoringRule(dens_PCm, cdf_PCm(1,:), y((T+1):(T+H)), THR_emp(1));
        C_score_PCm_1 = C_ScoringRule(dens_PCm, cdf_PCm(2,:), y((T+1):(T+H)), THR_emp(2));
        C_score_PCm_5 = C_ScoringRule(dens_PCm, cdf_PCm(3,:), y((T+1):(T+H)), THR_emp(3));            
    end

    %% PARTIALLY CENSORED 2: keep omega, alpha and beta uncensored, then censor nu and mu
    if arg8
        if (arg1 == 0)
            load(name,'mu_MLE','draw');
        end
        if (arg6 == 0)
            load(name,'mit_Cm');
        end        
        fprintf('*** Partially Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold_m);
        % mit_C: joint candidate for the joint censored posterior    
        draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
        partition = 3; 
        [draw_PCm2, a_PCm2, lnw_PCm2] = sim_cond_mit_MH_outloop(mit_Cm, draw_short, ...
            partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
        accept_PCm2 = mean(a_PCm2); 
        
        mean_draw_PCm2 = mean(draw_PCm2);
        median_draw_PCm2 = median(draw_PCm2);
        std_draw_PCm2 = std(draw_PCm2);
     
        % compute the implied volatility for the last in-sample period
        hT_PCm2 = volatility_t_garch11(draw_PCm2,y(1:T),y_S,0);  

        % predictive densities        
        dens_PCm2 = predictive_dens_t_garch11(y(T:(T+H)), hT_PCm2, draw_PCm2);
        % predicitve cdfs, constant threshold for different tails
        cdf_PCm2 = predictive_cdf_t_garch11(y(T:(T+H)), hT_PCm2, draw_PCm2, THR_emp);
        C_score_PCm2_05 = C_ScoringRule(dens_PCm2, cdf_PCm2(1,:), y((T+1):(T+H)), THR_emp(1));
        C_score_PCm2_1 = C_ScoringRule(dens_PCm2, cdf_PCm2(2,:), y((T+1):(T+H)), THR_emp(2));
        C_score_PCm2_5 = C_ScoringRule(dens_PCm2, cdf_PCm2(3,:), y((T+1):(T+H)), THR_emp(3));   
        
        
    end
        
    %% Results
    if save_on
        name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'.mat'];
        if exist(name,'file') %exits_on
            save(name,'-regexp','^mu','^Sigma','^draw','^mean','^median','^std',...
            '^accept','^mit','^CV',...
            '^dens','^cdf','^C_score','-append');
        else
            save(name,'-regexp','^mu','^Sigma','^draw','^mean','^median','^std',...
                '^accept','^mit','^CV',...
                '^dens','^cdf','^C_score');
        end
    end
% end