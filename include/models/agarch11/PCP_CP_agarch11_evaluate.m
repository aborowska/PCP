function results = PCP_CP_agarch11_evaluate(y, draw_C, draw_PC, ...
    p_bar0, p_bar1, p_bar, T, H)


    y_S = var(y(1:T));    
    M = 10000;
 
    h_post_C = volatility_agarch11(draw_C,y,y_S,H);
    y_post_C = randn(M,H).*sqrt(h_post_C);
    y_post_C = bsxfun(@plus,y_post_C,draw_C(:,1));
    y_post_C = sort(y_post_C);

    VaR_05_post_C = y_post_C(p_bar0*M,:); 
    VaR_1_post_C = y_post_C(p_bar1*M,:); 
    VaR_5_post_C = y_post_C(p_bar*M,:); 

    ES_05_post_C = mean(y_post_C(1:p_bar0*M,:)); 
    ES_1_post_C = mean(y_post_C(1:p_bar1*M,:)); 
    ES_5_post_C = mean(y_post_C(1:p_bar*M,:)); 
    

    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior ***\n');
    M_fin = size(draw_PC,1);
    h_post_PC = volatility_agarch11(draw_PC,y,y_S,H);
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
    
    results.VaR_C = [VaR_05_post_C;VaR_1_post_C;VaR_5_post_C];
    results.ES_C = [ES_05_post_C;ES_1_post_C;ES_5_post_C];   

    results.VaR_PC = [VaR_05_post_PC; VaR_1_post_PC; VaR_5_post_PC];
    results.ES_PC = [ES_05_post_PC; ES_1_post_PC; ES_5_post_PC];
end