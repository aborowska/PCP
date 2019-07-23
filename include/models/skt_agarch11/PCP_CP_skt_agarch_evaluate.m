function results = PCP_CP_skt_agarch_evaluate(results, y, T, H, THR_emp, THR_MLE) 


    % predictive densities
    dens_C = predictive_dens_skt_agarch11(y(T:(T+H)), results.fT_C, results.draw_C);
    % predicitve cdfs, constant threshold for different tails
    cdf_C = predictive_cdf_skt_agarch11(y(T:(T+H)), results.fT_C, results.draw_C, THR_emp);       
    C_score_C = C_ScoringRule(dens_C, cdf_C, y((T+1):(T+H)), THR_emp );

    % predicitve cdfs, var mle threshold for different tails  
    cdf_v_C = predictive_cdf_skt_agarch11(y(T:(T+H)), results.fT_C, results.draw_C, THR_MLE);       
    Cv_score_C = C_ScoringRule(dens_C, cdf_v_C, y((T+1):(T+H))', THR_MLE);

%     results.mit_C = mit_C;
%     results.CV_C = CV_C;
%     results.mu_C = mu_C;
%     results.Sigma_C = Sigma_C;
% 
%     results.draw_C = draw_C;        
%     results.accept_C = accept_C;
%     results.fT_C = hT_C;        
    results.dens_C = dens_C;
    results.cdf_C = cdf_C;
    results.C_score_C = C_score_C;
    results.cdf_v_C = cdf_v_C;
    results.Cv_score_C = Cv_score_C;


    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor nu, mu and omega
 

    % predictive densities
    dens_PC = predictive_dens_skt_agarch11(y(T:(T+H)), results.fT_PC, results.draw_PC);      
    % predicitve cdfs, constant threshold for different tails
    cdf_PC = predictive_cdf_skt_agarch11(y(T:(T+H)), results.fT_PC, results.draw_PC, THR_emp);
    C_score_PC = C_ScoringRule(dens_PC, cdf_PC, y((T+1):(T+H)), THR_emp);

    cdf_v_PC = predictive_cdf_skt_agarch11(y(T:(T+H)), results.fT_PC, results.draw_PC, THR_MLE);       
    Cv_score_PC = C_ScoringRule(dens_PC, cdf_v_PC, y((T+1):(T+H))', THR_MLE);

%     results.draw_PC = draw_PC;
%     results.accept_PC = accept_PC;
%     results.fT_PC = hT_PC;   
    results.dens_PC = dens_PC;
    results.cdf_PC = cdf_PC;
    results.C_score_PC = C_score_PC;
    results.cdf_v_PC = cdf_v_PC;
    results.Cv_score_PC = Cv_score_PC;  
    

end