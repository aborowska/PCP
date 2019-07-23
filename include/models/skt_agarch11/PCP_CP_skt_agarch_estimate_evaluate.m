function results = PCP_CP_skt_agarch_estimate_evaluate(kernel, kernel_init,...
    y, T, H, draw, df, cont, GamMat, mu_init, options,...
    M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE)  

    [mu_C,~,~,~,~,~] = fminunc(kernel_init,mu_init,options);
    [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_C,options);

    Sigma_C = inv(T*Sigma_C);

    r = fn_testSigma(reshape(Sigma_C,1,D^2)); % if r==1 then there is a problem with Sigma_C
    if ~r
        Sigma_start = Sigma_C;   
    else
        try 
            Sigma_start = nearestSPD(Sigma_C);  
        catch
%             Sigma_start = Sigma;
        end
    end

    [mean_theta_C, median_theta_C, std_theta_C, mean_accept_C, Draw_MH_C] = RWMH_skt_agarch11(kernel, ...
        @(xx)prior_skt_agarch11(xx,hyper), mu_C,[0.25, 10, 0.2, 10, 0.5, 0.007, 0.008], M, BurnIn, 'true');
% [0.3, 30, 0.3, 15, 0.5, 0.005, 0.005]

    cont.mit.N = 30000; %1e5;
    cont.mit.CV_max = 2.8; %2.1; %1.9;
    CV_C = cont.mit.CV_old;
    df = 4; %3;
    cont.mit.dfnc = 4; %5;
    while (CV_C(end) >= 2)
        try
            mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_start,1,D^2),'df', df, 'p', 1);
            [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            if CV_C(end)>2
                [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            end
            [draw_C, accept_C] = IndMH_mit(mit_C, kernel,M,BurnIn,GamMat);
        catch
            try
    %                     mit_C = struct('mu',mean_theta,'Sigma',reshape(diag(std_theta.^2),1,length(mu_C)^2),'df', df+1, 'p', 1);
                mit_C = struct('mu',median_theta_C,'Sigma',reshape(diag(std_theta_C.^2),1,length(mu_C)^2),'df', df+1, 'p', 1);
                [mit_C2, CV_C2] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
                [draw_C, accept_C] = IndMH_mit(mit_C, kernel,M,BurnIn,GamMat);
            catch
    %                     mit_C = struct('mu',mean_theta,'Sigma',reshape(diag(std_theta.^2),1,length(mu_C)^2),'df', df+1, 'p', 1);
                mit_C = struct('mu',median_theta_C,'Sigma',reshape(diag(std_theta_C.^2),1,length(mu_C)^2),'df', df+1, 'p', 1);
                [draw_C, accept_C] = IndMH_mit(mit_C, kernel,M,BurnIn,GamMat);
            end
        end 
    end   

    hT_C = volatility_skt_agarch11(draw_C,y(1:T),y_S,0);  

    % predictive densities
    dens_C = predictive_dens_skt_agarch11(y(T:(T+H)), hT_C, draw_C);
    % predicitve cdfs, constant threshold for different tails
    cdf_C = predictive_cdf_skt_agarch11(y(T:(T+H)), hT_C, draw_C, THR_emp);       
    C_score_C = C_ScoringRule(dens_C, cdf_C, y((T+1):(T+H)), THR_emp );

    % predicitve cdfs, var mle threshold for different tails  
    cdf_v_C = predictive_cdf_skt_agarch11(y(T:(T+H)), hT_C, draw_C, THR_MLE);       
    Cv_score_C = C_ScoringRule(dens_C, cdf_v_C, y((T+1):(T+H))', THR_MLE);

    results.mit_C = mit_C;
    results.CV_C = CV_C;
    results.mu_C = mu_C;
    results.Sigma_C = Sigma_C;

    results.draw_C = draw_C;        
    results.accept_C = accept_C;
    results.fT_C = hT_C;        
    results.dens_C = dens_C;
    results.cdf_C = cdf_C;
    results.C_score_C = C_score_C;
    results.cdf_v_C = cdf_v_C;
    results.Cv_score_C = Cv_score_C;


    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor nu, mu and omega

    partition = 5;
    thin = 10;
    II = 1;
    
    fprintf('*** Partially Censored Posterior, threshold 10%%, partition = %i ***\n',partition);
    % mit_C: joint candidate for the joint censored posteriors
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
    
    [draw_PC, a_PC, lnw_PC] = sim_cond_mit_MH_outloop(mit_C, draw_short,...
        partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thin);
    accept_PC = mean(a_PC); 

    hT_PC = volatility_skt_agarch11(draw_PC,y(1:T),y_S,0);  


    % predictive densities
    dens_PC = predictive_dens_skt_agarch11(y(T:(T+H)), hT_PC, draw_PC);      
    % predicitve cdfs, constant threshold for different tails
    cdf_PC = predictive_cdf_skt_agarch11(y(T:(T+H)), hT_PC, draw_PC, THR_emp);
    C_score_PC = C_ScoringRule(dens_PC, cdf_PC, y((T+1):(T+H)), THR_emp);

    cdf_v_PC = predictive_cdf_skt_agarch11(y(T:(T+H)), hT_PC, draw_PC, THR_MLE);       
    Cv_score_PC = C_ScoringRule(dens_PC, cdf_v_PC, y((T+1):(T+H))', THR_MLE);

    results.draw_PC = draw_PC;
    results.accept_PC = accept_PC;
    results.fT_PC = hT_PC;   
    results.dens_PC = dens_PC;
    results.cdf_PC = cdf_PC;
    results.C_score_PC = C_score_PC;
    results.cdf_v_PC = cdf_v_PC;
    results.Cv_score_PC = Cv_score_PC;  
    

end