%% %%%%%% ESTIMATION
%% Posterior
fprintf('*** Uncensored Posterior  ***\n');
kernel = @(xx) posterior_t_garch11_mex(xx, y(1:T), y_S, GamMat, hyper);
[draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
lnd = dmvgt(draw, mit, true, GamMat);
lnw = lnk - lnd;
lnw = lnw - max(lnw);
[ind, a] = fn_MH(lnw);
draw = draw(ind,:);
accept = a/(M+BurnIn);
draw = draw(BurnIn+1:BurnIn+M,:);    
%% Constant threshold
threshold = sort(y(1:T));
threshold = threshold(round(2*p_bar*T));
fprintf('*** Censored Posterior, threshold 10%% ***\n');
kernel = @(xx) C_posterior_t_garch11_2_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);
[draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
lnw_C = lnk_C - lnd_C;
lnw_C = lnw_C - max(lnw_C);
[ind, a] = fn_MH(lnw_C);
draw_C = draw_C(ind,:);
accept_C = a/(M+BurnIn);
draw_C = draw_C(BurnIn+1:BurnIn+M,:);

fprintf('*** Partially Censored Posterior, threshold 10%%, partition = 4 ***\n');
% mit_C: joint candidate for the joint censored posterior
partition = 4;
draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
[draw_PC, a_PC, lnw_PC] = sim_cond_mit_MH_outloop(mit_C, draw_short,...
partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
accept_PC = mean(a_PC); 
% fprintf('*** Partially Censored Posterior, threshold 10%%, partition = 3 ***\n');
% % mit_C: joint candidate for the joint censored posterior
% partition2 = 3;
% [draw_PC2, a_PC2, lnw_PC2] = sim_cond_mit_MH_outloop(mit_C, draw_short,...
% partition2, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
% accept_PC2 = mean(a_PC2); 


%% Time varying threshold

threshold_m = 0.1; %<---------- HiGhER?
quantile = tinv(threshold_m, mu_MLE(1));
fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold_m);
kernel = @(xx) C_posterior_t_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S, GamMat, hyper);
[draw_Cm, lnk_Cm] = fn_rmvgt_robust(M+BurnIn, mit_Cm, kernel, false);
lnd_Cm = dmvgt(draw_Cm, mit_Cm, true, GamMat);    
lnw_Cm = lnk_Cm - lnd_Cm;
lnw_Cm = lnw_Cm - max(lnw_Cm);
[ind, a] = fn_MH(lnw_Cm);
draw_Cm = draw_Cm(ind,:);
accept_Cm = a/(M+BurnIn);
draw_Cm = draw_Cm(BurnIn+1:BurnIn+M,:);

fprintf('*** Partially Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold_m);
% mit_C: joint candidate for the joint censored posterior    
draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
partition = 4; 
[draw_PCm, a_PCm, lnw_PCm] = sim_cond_mit_MH_outloop(mit_Cm, draw_short, ...
partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
accept_PCm = mean(a_PCm);

% fprintf('*** Partially Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold_m);
% % mit_C: joint candidate for the joint censored posterior    
% draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
% partition = 3; 
% [draw_PCm2, a_PCm2, lnw_PCm2] = sim_cond_mit_MH_outloop(mit_Cm, draw_short, ...
% partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
% accept_PCm2 = mean(a_PCm2); 



    
%% %%%% SAVE    
name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'_DUPA.mat'];
if exist(name,'file')
   save(name,'-regexp','^mu','^Sigma','^draw',...
        '^accept',...
        '\w*_MLE','^THR',...  
        '-append');
else
   save(name,'-regexp','^mu','^Sigma','^draw',...
        '^accept',...
        '\w*_MLE','^THR');    
end