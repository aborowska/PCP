addpath(genpath('include/'));

model = 'tgarch11';
data_name = 'IBM_T2000_crisis';
name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'.mat'];
p_bar0 = 0.005;
p_bar1 = 0.01;
p_bar = 0.05; 
P_bars = [p_bar0, p_bar1, p_bar];

y = load('Perc_Rets_GSPC_IBM_MSFT_T2000_crisis2.csv');
y = y(:,2);   
time = [2007,2017];   % start_date = '01012007'; end_date = '23122016    
H = 1500;
TT = length(y);
T = TT - H;     

y_S = var(y(1:T));
y_sort = sort(y(1:T));
THR_emp = y_sort(round(P_bars*T));

 
load(name, '-regexp','^draw','^THR_emp','^mu_MLE');
            
hT = volatility_t_garch11(draw,y(1:T),y_S,0);  
hT_MLE = volatility_t_garch11(mu_MLE,y(1:T),y_S,0);  
hT_C = volatility_t_garch11(draw_C,y(1:T),y_S,0);  
hT_PC = volatility_t_garch11(draw_PC,y(1:T),y_S,0);  
hT_Cm = volatility_t_garch11(draw_Cm,y(1:T),y_S,0);  
hT_PCm = volatility_t_garch11(draw_PCm,y(1:T),y_S,0);  

QUANT_MLE = tinv(P_bars, mu_MLE(1))';
  
%% Uncensored Posterior
        % predictive densities
        dens_post = predictive_dens_t_garch11(y(T:(T+H)), hT, draw);
         % predicitve cdfs, constant threshold for different tails
        cdf_post = predictive_cdf_t_garch11(y(T:(T+H)), hT, draw, THR_emp);
        C_score_post_05 = C_ScoringRule(dens_post, cdf_post(1,:), y((T+1):(T+H)), THR_emp(1));
        C_score_post_1 = C_ScoringRule(dens_post, cdf_post(2,:), y((T+1):(T+H)), THR_emp(2));
        C_score_post_5 = C_ScoringRule(dens_post, cdf_post(3,:), y((T+1):(T+H)), THR_emp(3));

        [THR_MLE, cond_MLE] = threshold_t_garch11_varc_mle(y((T+1):(T+H)), mu_MLE, hT_MLE, QUANT_MLE);
        
        cdf_v_post = predictive_cdf_t_garch11(y(T:(T+H)), hT, draw, THR_MLE);       
        Cv_score_post_05 = C_ScoringRule(dens_post, cdf_v_post(1,:), y((T+1):(T+H)), THR_MLE(1));
        Cv_score_post_1 = C_ScoringRule(dens_post, cdf_v_post(2,:), y((T+1):(T+H)), THR_MLE(2));
        Cv_score_post_5 = C_ScoringRule(dens_post, cdf_v_post(3,:), y((T+1):(T+H)), THR_MLE(3));          
 
%% CENSORED: Threshold = 10% perscentile of the data sample
        % predictive densities
        dens_C = predictive_dens_t_garch11(y(T:(T+H)), hT_C, draw_C);
        % predicitve cdfs, constant threshold for different tails
        cdf_C = predictive_cdf_t_garch11(y(T:(T+H)), hT_C, draw_C, THR_emp);       
        C_score_C_05 = C_ScoringRule(dens_C, cdf_C(1,:), y((T+1):(T+H)), THR_emp(1));
        C_score_C_1 = C_ScoringRule(dens_C, cdf_C(2,:), y((T+1):(T+H)), THR_emp(2));
        C_score_C_5 = C_ScoringRule(dens_C, cdf_C(3,:), y((T+1):(T+H)), THR_emp(3));    
        % predicitve cdfs, var mle threshold for different tails
        
        cdf_v_C = predictive_cdf_t_garch11(y(T:(T+H)), hT_C, draw_C, THR_MLE);       
        Cv_score_C_05 = C_ScoringRule(dens_C, cdf_v_C(1,:), y((T+1):(T+H)), THR_MLE(1));
        Cv_score_C_1 = C_ScoringRule(dens_C, cdf_v_C(2,:), y((T+1):(T+H)), THR_MLE(2));
        Cv_score_C_5 = C_ScoringRule(dens_C, cdf_v_C(3,:), y((T+1):(T+H)), THR_MLE(3));  
 
%% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor nu, mu and sigma
        % predictive densities
        dens_PC = predictive_dens_t_garch11(y(T:(T+H)), hT_PC, draw_PC);      
        % predicitve cdfs, constant threshold for different tails
        cdf_PC = predictive_cdf_t_garch11(y(T:(T+H)), hT_PC, draw_PC, THR_emp);
        C_score_PC_05 = C_ScoringRule(dens_PC, cdf_PC(1,:), y((T+1):(T+H)), THR_emp(1));
        C_score_PC_1 = C_ScoringRule(dens_PC, cdf_PC(2,:), y((T+1):(T+H)), THR_emp(2));
        C_score_PC_5 = C_ScoringRule(dens_PC, cdf_PC(3,:), y((T+1):(T+H)), THR_emp(3));     
        
        cdf_v_PC = predictive_cdf_t_garch11(y(T:(T+H)), hT_PC, draw_PC, THR_MLE);       
        Cv_score_PC_05 = C_ScoringRule(dens_PC, cdf_v_PC(1,:), y((T+1):(T+H)), THR_MLE(1));
        Cv_score_PC_1 = C_ScoringRule(dens_PC, cdf_v_PC(2,:), y((T+1):(T+H)), THR_MLE(2));
        Cv_score_PC_5 = C_ScoringRule(dens_PC, cdf_v_PC(3,:), y((T+1):(T+H)), THR_MLE(3));       

    
%% CENSORED MLE PARAMETERS
        % predictive densities        
        dens_Cm = predictive_dens_t_garch11(y(T:(T+H)), hT_Cm, draw_Cm);
        % predicitve cdfs, constant threshold for different tails
        cdf_Cm = predictive_cdf_t_garch11(y(T:(T+H)), hT_Cm, draw_Cm, THR_emp);
        C_score_Cm_05 = C_ScoringRule(dens_Cm, cdf_Cm(1,:), y((T+1):(T+H)), THR_emp(1));
        C_score_Cm_1 = C_ScoringRule(dens_Cm, cdf_Cm(2,:), y((T+1):(T+H)), THR_emp(2));
        C_score_Cm_5 = C_ScoringRule(dens_Cm, cdf_Cm(3,:), y((T+1):(T+H)), THR_emp(3));    
        
        cdf_v_Cm = predictive_cdf_t_garch11(y(T:(T+H)), hT_Cm, draw_Cm, THR_MLE);       
        Cv_score_Cm_05 = C_ScoringRule(dens_Cm, cdf_v_Cm(1,:), y((T+1):(T+H)), THR_MLE(1));
        Cv_score_Cm_1 = C_ScoringRule(dens_Cm, cdf_v_Cm(2,:), y((T+1):(T+H)), THR_MLE(2));
        Cv_score_Cm_5 = C_ScoringRule(dens_Cm, cdf_v_Cm(3,:), y((T+1):(T+H)), THR_MLE(3));         
 
    
%% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor nu, mu and sigma
        % predictive densities        
        dens_PCm = predictive_dens_t_garch11(y(T:(T+H)), hT_PCm, draw_PCm);
        % predicitve cdfs, constant threshold for different tails
        cdf_PCm = predictive_cdf_t_garch11(y(T:(T+H)), hT_PCm, draw_PCm, THR_emp);
        C_score_PCm_05 = C_ScoringRule(dens_PCm, cdf_PCm(1,:), y((T+1):(T+H)), THR_emp(1));
        C_score_PCm_1 = C_ScoringRule(dens_PCm, cdf_PCm(2,:), y((T+1):(T+H)), THR_emp(2));
        C_score_PCm_5 = C_ScoringRule(dens_PCm, cdf_PCm(3,:), y((T+1):(T+H)), THR_emp(3));     
        
        cdf_v_PCm = predictive_cdf_t_garch11(y(T:(T+H)), hT_PCm, draw_PCm, THR_MLE);       
        Cv_score_PCm_05 = C_ScoringRule(dens_PCm, cdf_v_PCm(1,:), y((T+1):(T+H)), THR_MLE(1));
        Cv_score_PCm_1 = C_ScoringRule(dens_PCm, cdf_v_PCm(2,:), y((T+1):(T+H)), THR_MLE(2));
        Cv_score_PCm_5 = C_ScoringRule(dens_PCm, cdf_v_PCm(3,:), y((T+1):(T+H)), THR_MLE(3));        
        
        
%%
name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'2.mat'];
save(name,'-regexp','^mu','^draw',...
                '^dens','^cdf','^C_score',...
                '^cdf_v_','^Cv_score',...
                '\w*_MLE','^THR');