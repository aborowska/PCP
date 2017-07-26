    %% MSE_ess
    MSE_es_1 = mean((ES_1 - cdf1).^2,2);
    MSE_es_1_post = mean((ES_1_post - cdf1).^2,2);
    if varc
        MSE_es_1_post_Cah = mean((ES_1_post_Cah - cdf1).^2,2);
        MSE_es_1_post_PCah = mean((ES_1_post_PCah - cdf1).^2,2);
        MSE_es_1_post_Cm = mean((ES_1_post_Cm - cdf1).^2,2);
        MSE_es_1_post_PCm = mean((ES_1_post_PCm - cdf1).^2,2);
    else
        MSE_es_1_post_C = mean((ES_1_post_C - cdf1).^2,2);
        MSE_es_1_post_PC = mean((ES_1_post_PC - cdf1).^2,2);
        MSE_es_1_post_C0 = mean((ES_1_post_C0 - cdf1).^2,2);
        MSE_es_1_post_PC0 = mean((ES_1_post_PC0 - cdf1).^2,2);
    end
    
    MSE_es_5 = mean((ES_5 - cdf5).^2,2);
    MSE_es_5_post = mean((ES_5_post - cdf5).^2,2);
    if varc
        MSE_es_5_post_Cah = mean((ES_5_post_Cah - cdf5).^2,2);
        MSE_es_5_post_PCah = mean((ES_5_post_PCah - cdf5).^2,2);
        MSE_es_5_post_Cm = mean((ES_5_post_Cm - cdf5).^2,2);
        MSE_es_5_post_PCm = mean((ES_5_post_PCm - cdf5).^2,2);
    else
        MSE_es_5_post_C = mean((ES_5_post_C - cdf5).^2,2);
        MSE_es_5_post_PC = mean((ES_5_post_PC - cdf5).^2,2);
        MSE_es_5_post_C0 = mean((ES_5_post_C0 - cdf5).^2,2);
        MSE_es_5_post_PC0 = mean((ES_5_post_PC0 - cdf5).^2,2);
    end
    
    MSE_es_05 = mean((ES_05 - cdf05).^2,2);
    MSE_es_05_post = mean((ES_05_post - cdf05).^2,2);
    if varc
        MSE_es_05_post_Cah = mean((ES_05_post_Cah - cdf05).^2,2);
        MSE_es_05_post_PCah = mean((ES_05_post_PCah - cdf05).^2,2);
        MSE_es_05_post_Cm = mean((ES_05_post_Cm - cdf05).^2,2);
        MSE_es_05_post_PCm = mean((ES_05_post_PCm - cdf05).^2,2);        
    else
        MSE_es_05_post_C = mean((ES_05_post_C - cdf05).^2,2);
        MSE_es_05_post_PC = mean((ES_05_post_PC - cdf05).^2,2);
        MSE_es_05_post_C0 = mean((ES_05_post_C0 - cdf05).^2,2);
        MSE_es_05_post_PC0 = mean((ES_05_post_PC0 - cdf05).^2,2);
    end
      
 