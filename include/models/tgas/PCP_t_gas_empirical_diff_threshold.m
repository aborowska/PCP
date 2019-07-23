function results = t_gas_CP_PCP(THR, y, T, H, hyper, cont, GamMat, mu_init, options,...
    M, BurnIn,BurnIn_PCP, THR_emp, THR_mle)  
    
    thinning = 1;
    II = 10;
    %% CENSORED: Threshold = 20% perscentile of the data sample
        threshold2 = sort(y(1:T));
        threshold2 = threshold2(round(THR*T));
        fprintf('*** Censored Posterior, threshold 20%% ***\n');
        kernel_init = @(xx) - C_posterior_t_gas(xx, y(1:T,1), threshold2, hyper)/T;    
        kernel = @(xx) C_posterior_t_gas_mex(xx, y(1:T,1), threshold2, GamMat, hyper);
        [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options); 
        Sigma_C = inv(T*Sigma_C); 
        
        r = fn_testSigma(reshape(Sigma_C,1,5^2)); % if r==1 then there is a problem with Sigma_C
        if r
            Sigma_startb = Sigma;        
        else
            Sigma_startb = Sigma_C;
        end
        
        cont.mit.CV_max = 1.9;
        CV_C = cont.mit.CV_old;
        while (CV_C(end) >= 2)
            try
                mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_startb,1,length(mu_C)^2),'df', df, 'p', 1);
                [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
                if CV_C(end)>2
                    [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
                end
                 [draw_C, accept_C] = IndMH_mit(mit_C, kernel,M,BurnIn,GamMat);
            catch
                mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_startb,1,length(mu_C)^2),'df', df, 'p', 1);
                [draw_C, accept_C] = IndMH_mit(mit_C, kernel,M,BurnIn,GamMat);
            end 
        end  
        
        fT_C = volatility_t_gas(draw_C,y(1:T),0);  
        ind_okb = (draw_C(:,4)>0);

        % predictive densities
        dens_C = predictive_dens_t_gas(y(T:(T+H)), fT_C(ind_okb,:), draw_C(ind_okb,:));
        % predicitve cdfs, constant threshold for different tails
        cdf_C = predictive_cdf_t_gas(y(T:(T+H)), fT_C(ind_okb,:), draw_C(ind_okb,:), THR_emp);        
        C_score_C = C_ScoringRule(dens_C, cdf_C, y((T+1):(T+H)), THR_emp);    
        % predicitve cdfs, var mle threshold for different tails       
        cdf_v_C = predictive_cdf_t_gas(y(T:(T+H)), fT_C(ind_okb,:), draw_C(ind_okb,:), THR_MLE);         
        Cv_score_C = C_ScoringRule(dens_C, cdf_v_C, y((T+1):(T+H))', THR_MLE);  
      
    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor nu, mu and sigma
        fprintf('*** Partially Censored Posterior, threshold 10%%, partition = 4 ***\n');
        % mit_C: joint candidate for the joint censored posterior
        partition = 4;
        draw_short = draw((1:II:M)',:); % thinning - to get higfT quality rhos
        [draw_PC, a_PC, ~] = sim_cond_mit_MH_outloop(mit_C, draw_short,...
            partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
        accept_PC = mean(a_PC); 

        fT_PC = volatility_t_gas(draw_PC,y(1:T),0);          
        % predictive densities
        dens_PC = predictive_dens_t_gas(y(T:(T+H)), fT_PC, draw_PC);      
        % predicitve cdfs, constant threshold for different tails
        cdf_PC = predictive_cdf_t_gas(y(T:(T+H)), fT_PC, draw_PC, THR_emp);
        C_score_PC = C_ScoringRule(dens_PC, cdf_PC, y((T+1):(T+H)), THR_emp);   
        % predicitve cdfs, var mle threshold for different tails               
        cdf_v_PC = predictive_cdf_t_gas(y(T:(T+H)), fT_PC, draw_PC, THR_MLE);             
        Cv_score_PC = C_ScoringRule(dens_PC, cdf_v_PC, y((T+1):(T+H))', THR_MLE);       

    %%  
        results.draw_C = draw_C;
        results.CV_C = CV_C;
        results.accept_C = accept_C;
        results.dens_C = dens_C;
        results.cdf_C = cdf_C;
        results.cdf_v_C = cdf_v_C;
        results.C_score_C = C_score_C;
        results.Cv_score_C = Cv_score_C;
        
        results.draw_PC = draw_PC;
        results.accept_PC = accept_PC;
        results.dens_PC = dens_PC;
        results.cdf_PC = cdf_PC;
        results.cdf_v_PC = cdf_v_PC;
        results.C_score_PC = C_score_PC;
        results.Cv_score_PC = Cv_score_PC;        
end
