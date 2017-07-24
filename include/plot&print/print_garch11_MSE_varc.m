clear all
addpath(genpath('include/'));

model = 'agarch11';

TT = [1000,2500];
  	
S = 50;
H = 100;
sigma1 = 1;
sigma2 = 2;
c = (sigma2 - sigma1)/sqrt(2*pi); % mean of eps
mu = 0;
mu2 = 0;
omega = 1;
alpha = 0.1;
beta = 0.8;
parameters = {'$\\mu$','$\\omega$','$\\alpha$','$\\beta$'};



%% T = 1000;
T = 1000;
T = 2500;
clear -regexp ^MSE ^mean ^draw

% VARC LOW (_5 variables are for 0.5% VaR)
name = ['results\',model,'\',model,'_1_2_T',num2str(T),'_H100_II10_PCP0_MC_(R2017a)_varc_low.mat'];
load(name,'param_true')
load(name,'-regexp','^MSE')
load(name,'-regexp','^mean')
load(name,'-regexp','^draw')

% draw2 = draw;
% draw_Cm2 = draw_C0;
% draw_PCm2 = draw_PC0;
% draw_Cah2 = draw_C;
% draw_PCah2 = draw_PC;
% 
% mean_draw2 = mean_draw;
% mean_draw_Cm2 = mean_draw_C0;
% mean_draw_PCm2 = mean_draw_PC0;
% mean_draw_Cah2 = mean_draw_C;
% mean_draw_PCah2 = mean_draw_PC;

MSE_05 = MSE_5;
MSE_05_post = MSE_5_post;
MSE_05_post_Cm = MSE_5_post_C0;
MSE_05_post_PCm = MSE_5_post_PC0;
MSE_05_post_Cah = MSE_5_post_C;
MSE_05_post_PCah = MSE_5_post_PC;

% VARC 
name = ['results\',model,'\',model,'_1_2_T',num2str(T),'_H100_II10_PCP0_MC_(R2017a)_varc.mat'];
load(name,'-regexp','^MSE')
load(name,'-regexp','^mean')
load(name,'-regexp','^draw')

draw_Cm = draw_C0;
draw_PCm = draw_PC0;
draw_Cah = draw_C;
draw_PCah = draw_PC;

mean_draw_Cm = mean_draw_C0;
mean_draw_PCm = mean_draw_PC0;
mean_draw_Cah = mean_draw_C;
mean_draw_PCah = mean_draw_PC;

MSE_5_post_Cm = MSE_5_post_C0;
MSE_5_post_PCm = MSE_5_post_PC0;
MSE_5_post_Cah = MSE_5_post_C;
MSE_5_post_PCah = MSE_5_post_PC;

MSE_1_post_Cm = MSE_1_post_C0;
MSE_1_post_PCm = MSE_1_post_PC0;
MSE_1_post_Cah = MSE_1_post_C;
MSE_1_post_PCah = MSE_1_post_PC;

clear 'draw_C0' 'draw_PC0' 'draw_C' 'draw_PC'
clear 'mean_draw_PC0' 'mean_draw_PC' 'mean_draw_C0' 'mean_draw_C'
clear 'MSE_1_post_PC0' 'MSE_1_post_PC' 'MSE_1_post_C0' 'MSE_1_post_C'
clear 'MSE_5_post_PC0' 'MSE_5_post_PC' 'MSE_5_post_C0' 'MSE_5_post_C'

mean(MSE_05_post_PCah(MSE_05_post_PCah<1))  % 0.0719
mean(MSE_05_post_PCm(MSE_05_post_PCm<1))    % 0.1826

mean(MSE_1_post_PCah(MSE_1_post_PCah<1))  %  0.0611
mean(MSE_1_post_PCm(MSE_1_post_PCm<1))    %  0.1664

mean(MSE_5_post_PCah(MSE_5_post_PCah<1))  % 0.0326
mean(MSE_5_post_PCm(MSE_5_post_PCm<1))    % 0.0531


MSE_1_post_PC0 = NaN;
MSE_1_post_PC = NaN;
MSE_1_post_C0 = NaN;
MSE_1_post_C = NaN;
MSE_5_post_PC0 = NaN;
MSE_5_post_PC = NaN;
MSE_5_post_C0 = NaN;
MSE_5_post_C = NaN;
MSE_05_post_PC0 = NaN;
MSE_05_post_PC = NaN;
MSE_05_post_C0 = NaN;
MSE_05_post_C = NaN;



%% Create table
fname = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_H',num2str(H),'_pcp_mc_varc.tex'];
FID = fopen(fname, 'w+');

fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');
fprintf(FID, '\\begin{sidewaystable} \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\begin{tabular}{cc cc cccc cccc} \n');

fprintf(FID, ['Value && True & Posterior & CP0  & PCP0 & CP10\\%%  & PCP10\\%% &', ...
'CP var ah & PCP var ah & CP var mle & PCP var mle \\\\ \\hline \n']); 
for T = TT
        fprintf(FID, '\\hline \n');
        fprintf(FID, ['\\multicolumn{12}{c}{$T =',num2str(T),'$}  \\\\ \n']);
        fprintf(FID, '\\hline \n');
        
        fprintf(FID,['VaR 0.5\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f &', ...
            '%6.4f & %6.4f & %6.4f & %6.4f \\\\ \n'], ...
            mean(MSE_05), mean(MSE_05_post), mean(MSE_05_post_C0), mean(MSE_05_post_PC0), ...
            mean(MSE_05_post_C), mean(MSE_05_post_PC), ...        
            mean(MSE_05_post_Cah(MSE_05_post_Cah<1)), mean(MSE_05_post_PCah(MSE_05_post_PCah<1)), ...     
            mean(MSE_05_post_Cm(MSE_05_post_Cm<1)), mean(MSE_05_post_PCm(MSE_05_post_PCm<1))); 
        fprintf(FID,['no. MSE>1 && %d & %d & %d & %d & %d & %d &', ...
            '%d & %d & %d & %d \\\\[1ex] \n'], ...
            sum(MSE_05>=1), sum(MSE_05_post>=1), sum(MSE_05_post_C0>=1), sum(MSE_05_post_PC0>=1), ...
            sum(MSE_05_post_C>=1), sum(MSE_05_post_PC>=1), ...        
            sum(MSE_05_post_Cah>=1), sum(MSE_05_post_PCah>=1), ...     
            sum(MSE_05_post_Cm>=1), sum(MSE_05_post_PCm>=1));     
        
        fprintf(FID,['VaR 1\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f &', ...
            '%6.4f & %6.4f & %6.4f & %6.4f \\\\ \n'], ...
            mean(MSE_1), mean(MSE_1_post), mean(MSE_1_post_C0), mean(MSE_1_post_PC0), ...
            mean(MSE_1_post_C), mean(MSE_1_post_PC), ...        
            mean(MSE_1_post_Cah(MSE_1_post_Cah>=1)), mean(MSE_1_post_PCah(MSE_1_post_PCah>=1)), ...     
            mean(MSE_1_post_Cm(MSE_1_post_Cm>=1)), mean(MSE_1_post_PCm(MSE_1_post_PCm>=1)));        
        fprintf(FID,['no. MSE>1 && %d & %d & %d & %d & %d & %d &', ...
            '%d & %d & %d & %d \\\\[1ex] \n'], ...
            sum(MSE_1>=1), sum(MSE_1_post>=1), sum(MSE_1_post_C0>=1), sum(MSE_1_post_PC0>=1), ...
            sum(MSE_1_post_C>=1), sum(MSE_1_post_PC>=1), ...        
            sum(MSE_1_post_Cah>=1), sum(MSE_1_post_PCah>=1), ...     
            sum(MSE_1_post_Cm>=1), sum(MSE_1_post_PCm>=1));                
        
        fprintf(FID,['VaR 5\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f &', ...
            '%6.4f & %6.4f & %6.4f & %6.4f \\\\ \n'], ...
            mean(MSE_5), mean(MSE_5_post), mean(MSE_5_post_C0), mean(MSE_5_post_PC0), ...
            mean(MSE_5_post_C), mean(MSE_5_post_PC), ...        
            mean(MSE_5_post_Cah(MSE_5_post_Cah<1)), mean(MSE_5_post_PCah(MSE_5_post_PCah<1)), ...     
            mean(MSE_5_post_Cm(MSE_5_post_Cm<1)), mean(MSE_5_post_PCm(MSE_5_post_PCm<1)));          
        fprintf(FID,['no. MSE>1 && %d & %d & %d & %d & %d & %d &', ...
            '%d & %d & %d & %d \\\\[1ex] \n'], ...
            sum(MSE_5>=1), sum(MSE_5_post>=1), sum(MSE_5_post_C0>=1), sum(MSE_5_post_PC0>=1), ...
            sum(MSE_5_post_C>=1), sum(MSE_5_post_PC>=1), ...        
            sum(MSE_5_post_Cah>=1), sum(MSE_5_post_PCah>=1), ...     
            sum(MSE_5_post_Cm>=1), sum(MSE_5_post_PCm>=1));           
end      
fprintf(FID, '\\hline \n');

fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{CP: Censored posterior.}}  \\\\ \n');
fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{PCP: Partially censored posterior.}} \\\\ \n');
fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{(P)CP0: Censoring with threshold 0.}} \\\\ \n'); 
fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{(P)CP10\\%%: Censoring with threshold 10\\%% sample quantile.}}  \\\\ \n');
fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{(P)CP var ad: Time Varying Censoring, ad hoc method.}} \\\\ \n'); 
fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{(P)CP var mle: Time Varying Censoring, MLE based method.}}  \\\\ \n');    
fprintf(FID, '\\end{tabular}\n ');

caption = ['\\caption{MSEs for VaR prediction for standard posterior, censored posterior and partially censored posterior',...
    ' (the latter two with two time-constant and two time-varying thresholds) for the ',...
    model,' zero mean split normal model with $\\sigma_{1} = ',num2str(sigma1),'$ and $\\sigma_{2} = ',num2str(sigma2),'$.',...
    ' For the censored and the partially censored posterior the focus is on the left tail.',...            
    ' Average MSEs (over 50 simulations) averaged over out-of-sample horizon of $H=',num2str(H),'$.',...
    ' (Mean) standard errors in parentheses, (Mean) MSEs in brackets.} \n'];        

fprintf(FID, caption);

label = ['\\label{tab:',model,'_pcp_var}  \n'];
fprintf(FID, label);

fprintf(FID, '\\end{sidewaystable}\n');

fprintf(FID, '}');
fclose(FID);        
  