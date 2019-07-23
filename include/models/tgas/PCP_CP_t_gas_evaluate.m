function results = PCP_CP_t_gas_evaluate(results, y, T, H, THR_emp, THR_MLE)  


    %% CENSORED: Threshold = 20% perscentile of the data sample

%         fT_C = volatility_t_gas(draw_C,y(1:T),0);  
        ind_okb = (results.draw_C(:,4)>0);

        % predictive densities
        dens_C = predictive_dens_t_gas(y(T:(T+H)), results.fT_C(ind_okb,:), results.draw_C(ind_okb,:));
        % predicitve cdfs, constant threshold for different tails
        cdf_C = predictive_cdf_t_gas(y(T:(T+H)), results.fT_C(ind_okb,:), results.draw_C(ind_okb,:), THR_emp);        
        C_score_C = C_ScoringRule(dens_C, cdf_C, y((T+1):(T+H)), THR_emp);    
        % predicitve cdfs, var mle threshold for different tails       
        cdf_v_C = predictive_cdf_t_gas(y(T:(T+H)), results.fT_C(ind_okb,:), results.draw_C(ind_okb,:), THR_MLE);         
        Cv_score_C = C_ScoringRule(dens_C, cdf_v_C, y((T+1):(T+H))', THR_MLE);  
      
    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor nu, mu and sigma
 
%         fT_PC = volatility_t_gas(draw_PC,y(1:T),0);          

        % predictive densities
        dens_PC = predictive_dens_t_gas(y(T:(T+H)), results.fT_PC, results.draw_PC);      
        % predicitve cdfs, constant threshold for different tails
        cdf_PC = predictive_cdf_t_gas(y(T:(T+H)), results.fT_PC, results.draw_PC, THR_emp);
        C_score_PC = C_ScoringRule(dens_PC, cdf_PC, y((T+1):(T+H)), THR_emp);   
        % predicitve cdfs, var mle threshold for different tails               
        cdf_v_PC = predictive_cdf_t_gas(y(T:(T+H)), results.fT_PC, results.draw_PC, THR_MLE);             
        Cv_score_PC = C_ScoringRule(dens_PC, cdf_v_PC, y((T+1):(T+H))', THR_MLE);       

    %%  
%         results.mit_C = mit_C;
%         results.CV_C = CV_C;
%         results.mu_C = mu_C;
%         results.Sigma_C = Sigma_C;
%         
%         results.draw_C = draw_C;        
%         results.accept_C = accept_C;
        results.dens_C = dens_C;
        results.cdf_C = cdf_C;
        results.cdf_v_C = cdf_v_C;
        results.C_score_C = C_score_C;
        results.Cv_score_C = Cv_score_C;
%         results.fT_C = fT_C;
        
%         results.draw_PC = draw_PC;
%         results.accept_PC = accept_PC;
        results.dens_PC = dens_PC;
        results.cdf_PC = cdf_PC;
        results.cdf_v_PC = cdf_v_PC;
        results.C_score_PC = C_score_PC;
        results.Cv_score_PC = Cv_score_PC;  
%         results.fT_PC = fT_PC;        
end
