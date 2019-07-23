%% %%%% EVALUATION
%% Posterior 
hT = volatility_t_garch11(draw,y(1:T),y_S,0);
hT_MLE = volatility_t_garch11(mu_MLE,y(1:T),y_S,0);
[THR_MLE, cond_MLE] = threshold_t_garch11_varc_mle(y((T+1):(T+H)), mu_MLE, hT_MLE, QUANT_MLE);

% predictive densities
dens_post = predictive_dens_t_garch11(y(T:(T+H)), hT, draw);
% predicitve cdfs, constant threshold for different tails
cdf_post = predictive_cdf_t_garch11(y(T:(T+H)), hT, draw, THR_emp);
C_score_post_05 = C_ScoringRule(dens_post, cdf_post(1,:), y((T+1):(T+H)), THR_emp(1));
C_score_post_1 = C_ScoringRule(dens_post, cdf_post(2,:), y((T+1):(T+H)), THR_emp(2));
C_score_post_5 = C_ScoringRule(dens_post, cdf_post(3,:), y((T+1):(T+H)), THR_emp(3));

cdf_v_post = predictive_cdf_t_garch11(y(T:(T+H)), hT, draw, THR_MLE);       
Cv_score_post_05 = C_ScoringRule(dens_post, cdf_v_post(1,:), y((T+1):(T+H)), THR_MLE(1));
Cv_score_post_1 = C_ScoringRule(dens_post, cdf_v_post(2,:), y((T+1):(T+H)), THR_MLE(2));
Cv_score_post_5 = C_ScoringRule(dens_post, cdf_v_post(3,:), y((T+1):(T+H)), THR_MLE(3));          

%% Censored posterior constant threshold
hT_C = volatility_t_garch11(draw_C,y(1:T),y_S,0);  
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


%% Partially censored posterior constant threshold (censoring of omega)
hT_PC = volatility_t_garch11(draw_PC,y(1:T),y_S,0);  
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

%% Partially censored posterior constant threshold (no censoring of omega)
hT_PC2 = volatility_t_garch11(draw_PC2,y(1:T),y_S,0);  
% predictive densities        
dens_PC2 = predictive_dens_t_garch11(y(T:(T+H)), hT_PC2, draw_PC2);
% predicitve cdfs, constant threshold for different tails
cdf_PC2 = predictive_cdf_t_garch11(y(T:(T+H)), hT_PC2, draw_PC2, THR_emp);
C_score_PC2_05 = C_ScoringRule(dens_PC2, cdf_PC2(1,:), y((T+1):(T+H)), THR_emp(1));
C_score_PC2_1 = C_ScoringRule(dens_PC2, cdf_PC2(2,:), y((T+1):(T+H)), THR_emp(2));
C_score_PC2_5 = C_ScoringRule(dens_PC2, cdf_PC2(3,:), y((T+1):(T+H)), THR_emp(3));  

cdf_v_PC2 = predictive_cdf_t_garch11(y(T:(T+H)), hT_PC2, draw_PC2, THR_MLE);       
Cv_score_PC2_05 = C_ScoringRule(dens_PC2, cdf_v_PC2(1,:), y((T+1):(T+H)), THR_MLE(1));
Cv_score_PC2_1 = C_ScoringRule(dens_PC2, cdf_v_PC2(2,:), y((T+1):(T+H)), THR_MLE(2));
Cv_score_PC2_5 = C_ScoringRule(dens_PC2, cdf_v_PC2(3,:), y((T+1):(T+H)), THR_MLE(3));         

%% Censored posterior time varying MLE-based threshold
hT_Cm = volatility_t_garch11(draw_Cm,y(1:T),y_S,0);  
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


%% Partially censored posterior time varying MLE-based threshold (censoring of omega)
hT_PCm = volatility_t_garch11(draw_PCm,y(1:T),y_S,0);  
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

%% Partially censored posterior time varying MLE-based threshold (no censoring of omega)
hT_PCm2 = volatility_t_garch11(draw_PCm2,y(1:T),y_S,0);  
% predictive densities        
dens_PCm2 = predictive_dens_t_garch11(y(T:(T+H)), hT_PCm2, draw_PCm2);
% predicitve cdfs, constant threshold for different tails
cdf_PCm2 = predictive_cdf_t_garch11(y(T:(T+H)), hT_PCm2, draw_PCm2, THR_emp);
C_score_PCm2_05 = C_ScoringRule(dens_PCm2, cdf_PCm2(1,:), y((T+1):(T+H)), THR_emp(1));
C_score_PCm2_1 = C_ScoringRule(dens_PCm2, cdf_PCm2(2,:), y((T+1):(T+H)), THR_emp(2));
C_score_PCm2_5 = C_ScoringRule(dens_PCm2, cdf_PCm2(3,:), y((T+1):(T+H)), THR_emp(3));   

cdf_v_PCm2 = predictive_cdf_t_garch11(y(T:(T+H)), hT_PCm2, draw_PCm2, THR_MLE);       
Cv_score_PCm2_05 = C_ScoringRule(dens_PCm2, cdf_v_PCm2(1,:), y((T+1):(T+H)), THR_MLE(1));
Cv_score_PCm2_1 = C_ScoringRule(dens_PCm2, cdf_v_PCm2(2,:), y((T+1):(T+H)), THR_MLE(2));
Cv_score_PCm2_5 = C_ScoringRule(dens_PCm2, cdf_v_PCm2(3,:), y((T+1):(T+H)), THR_MLE(3));         


%% %%%%%% DM tests

    aaa_05 = who('-regexp','C_score\w*_05$'); 
    aaa_1 = who('-regexp','C_score\w*_1$'); 
    aaa_5 = who('-regexp','C_score\w*_5$'); 

    bbb_05 = who('-regexp','Cv_score\w*_05$'); 
    bbb_1 = who('-regexp','Cv_score\w*_1$'); 
    bbb_5 = who('-regexp','Cv_score\w*_5$'); 
    
    M = length(aaa_05);
    
    DM = NaN(M,1); % NaN will be removed while printing to tex
    DM_05 = diag(DM);
    DM_1 = diag(DM);
    DM_5 = diag(DM);
    % positive terms: the colum method is better than the row one
    % negative terms: the row method is better than the column one

    for ii = 2:M
        for jj = 1:(ii-1)
            d = eval(char(aaa_05{ii})) - eval(char(aaa_05{jj}));
            if (sum(imag(d)~=0) > 0)
                DM_05(ii,jj) = NaN;
            else
                nwse = sqrt(NeweyWest(d));
                DM_05(ii,jj) = mean(d)/nwse;
            end
            
            d = eval(char(aaa_1{ii})) - eval(char(aaa_1{jj}));
            if (sum(imag(d)~=0) > 0)
                DM_1(ii,jj) = NaN;
            else
                nwse = sqrt(NeweyWest(d));
                DM_1(ii,jj) = mean(d)/nwse;
            end
            
            d = eval(char(aaa_5{ii})) - eval(char(aaa_5{jj}));
            if (sum(imag(d)~=0) > 0)
                DM_5(ii,jj) = NaN;
            else
                nwse = sqrt(NeweyWest(d));
                DM_5(ii,jj) = mean(d)/nwse;
            end
        end
    end

    %% time varying     
    DM_v = NaN(M,1); % NaN will be removed while printing to tex
    DM_v_05 = diag(DM_v);
    DM_v_1 = diag(DM_v);
    DM_v_5 = diag(DM_v);
    % positive terms: the colum method is better than the row one
    % negative terms: the row method is better than the column one

    for ii = 2:M
        for jj = 1:(ii-1)
            d = eval(char(bbb_05{ii})) - eval(char(bbb_05{jj}));
            if (sum(imag(d)~=0) > 0)
                DM_v_05(ii,jj) = NaN;
            else
                nwse = sqrt(NeweyWest(d));
                DM_v_05(ii,jj) = mean(d)/nwse;
            end
            
            d = eval(char(bbb_1{ii})) - eval(char(bbb_1{jj}));
            if (sum(imag(d)~=0) > 0)
                DM_v_1(ii,jj) = NaN;
            else
                nwse = sqrt(NeweyWest(d));
                DM_v_1(ii,jj) = mean(d)/nwse;
            end
            
            d = eval(char(bbb_5{ii})) - eval(char(bbb_5{jj}));
            if (sum(imag(d)~=0) > 0)
                DM_v_5(ii,jj) = NaN;
            else
                nwse = sqrt(NeweyWest(d));
                DM_v_5(ii,jj) = mean(d)/nwse;
            end
        end
    end
    
    
%% %%%% SAVE    
name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'_DUPA.mat'];
if exist(name,'file')
   save(name,'-regexp','^mu','^Sigma','^draw',...%,'^mean','^median','^std',...
        '^accept','^mit',...%'^CV',...
        '^dens','^cdf','^C_score',...
        '^cdf_v_','^Cv_score',...
        '\w*_MLE','^THR',...
        '^DM',...
        '-append');
else
   save(name,'-regexp','^mu','^Sigma','^draw',...%,'^mean','^median','^std',...
        '^accept','^mit',...%'^CV',...
        '^dens','^cdf','^C_score',...
        '^cdf_v_','^Cv_score',...
        '\w*_MLE','^THR',...
        '^DM');    
end