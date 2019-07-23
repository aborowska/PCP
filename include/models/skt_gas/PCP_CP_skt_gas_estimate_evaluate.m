function results = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
    M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE)  

% cont0 = cont; cont.mit.ISpc = 0.2; cont= cont0;
% cont.mit.N = 10000; cont.mit.ISpc = 0.1; cont.mit.dfnc = 3; cont= cont0;


    %% CENSORED: Threshold = 20% perscentile of the data sample
        %  [mean_theta, median_theta, std_theta, mean_accept, THETA] = ...
%             RWMH_skt_gas(kernel, prior, mu_init, delta, MRW, BurnInRW, plot_on)
        % mu_Cb = mean(THETA); 
        % Sigma_Cb=cov(THETA)
        % mit_C = struct('mu',mu_Cb,'Sigma',reshape(Sigma_Cb,1,length(mu_Cb)^2),'df', df, 'p', 1);

        [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options); 
        Sigma_C = inv(T*Sigma_C); 
        
        r = fn_testSigma(reshape(Sigma_C,1,length(mu_C)^2)); % if r==1 then there is a problem with Sigma_C
        if r
%             Sigma_startb = Sigma; 
            Sigma_startb = nearestSPD(Sigma_C);
        else
            Sigma_startb = Sigma_C;
        end
        
%         cont.mit.N = 30000;
        cont.mit.CV_max = 1.9;  2.2; 2.5; 1.5;1.9;
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
        
        fT_C = volatility_skt_gas(draw_C,y(1:T),0);  
        ind_okb = (draw_C(:,5)>0);

        % predictive densities
        dens_C = predictive_dens_skt_gas(y(T:(T+H)), fT_C(ind_okb,:), draw_C(ind_okb,:));
        % predicitve cdfs, constant threshold for different tails
        cdf_C = predictive_cdf_skt_gas(y(T:(T+H)), fT_C(ind_okb,:), draw_C(ind_okb,:), THR_emp);        
        C_score_C = C_ScoringRule(dens_C, cdf_C, y((T+1):(T+H)), THR_emp);    
        % predicitve cdfs, var mle threshold for different tails       
        cdf_v_C = predictive_cdf_skt_gas(y(T:(T+H)), fT_C(ind_okb,:), draw_C(ind_okb,:), THR_MLE);         
        Cv_score_C = C_ScoringRule(dens_C, cdf_v_C, y((T+1):(T+H))', THR_MLE);  
      
    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor nu, mu and sigma
        thinning = 10;
        II = 1;
        
        fprintf('*** Partially Censored Posterior, threshold 10%%, partition = 5 ***\n');
        % mit_C: joint candidate for the joint censored posterior
        partition = 5;
        draw_short = draw((1:II:M)',:); % thinning - to get higfT quality rhos
        [draw_PC, a_PC, ~] = sim_cond_mit_MH_outloop(mit_C, draw_short,...
            partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
        accept_PC = mean(a_PC); 

        fT_PC = volatility_skt_gas(draw_PC,y(1:T),0);          
        % predictive densities
        dens_PC = predictive_dens_skt_gas(y(T:(T+H)), fT_PC, draw_PC);      
        % predicitve cdfs, constant threshold for different tails
        cdf_PC = predictive_cdf_skt_gas(y(T:(T+H)), fT_PC, draw_PC, THR_emp);
        C_score_PC = C_ScoringRule(dens_PC, cdf_PC, y((T+1):(T+H)), THR_emp);   
        % predicitve cdfs, var mle threshold for different tails               
        cdf_v_PC = predictive_cdf_skt_gas(y(T:(T+H)), fT_PC, draw_PC, THR_MLE);             
        Cv_score_PC = C_ScoringRule(dens_PC, cdf_v_PC, y((T+1):(T+H))', THR_MLE);       

         
    %%  
        results.mit_C = mit_C;
        results.CV_C = CV_C;
        results.mu_C = mu_Cb;
        results.Sigma_C = Sigma_Cb;
        
        results.draw_C = draw_C;        
        results.accept_C = accept_C;
        results.dens_C = dens_C;
        results.cdf_C = cdf_C;
        results.cdf_v_C = cdf_v_C;
        results.C_score_C = C_score_C;
        results.Cv_score_C = Cv_score_C;
        results.fT_C = fT_C;
        
        results.draw_PC = draw_PC;
        results.accept_PC = accept_PC;
        results.dens_PC = dens_PC;
        results.cdf_PC = cdf_PC;
        results.cdf_v_PC = cdf_v_PC;
        results.C_score_PC = C_score_PC;
        results.Cv_score_PC = Cv_score_PC;  
        results.fT_PC = fT_PC; 
%         results.THRES = THRES;
end



if false
    for jj = 1:6
        subplot(3,2,jj)
        hold on
        histogram(draw_C(:,jj),50)
        histogram(draw_PC(:,jj),50)        
    end


end

