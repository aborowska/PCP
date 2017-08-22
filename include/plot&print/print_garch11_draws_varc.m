clear all
addpath(genpath('include/'));

% model = 'arch1'; %'arch1';
% TT = [1000,2500];
% suffix = [];

model = 'garch11'; 
TT = [1000, 2500];
suffix = '_tunning';% [];
  	
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
if strcmp(model,'garch11')
    parameters = {'$\\mu$','$\\omega$','$\\alpha$','$\\beta$'};
    param_true = [0,1,0.1,0.8];
elseif strcmp(model,'agarch11')
    parameters = {'$\\mu_{1}$','$\\omega$','$\\mu_{2}$','$\\alpha$','$\\beta$'};  
    param_true = [0,1,0,0.1,0.8]; 
elseif strcmp(model,'arch1')
    parameters = {'$\\mu_{1}$','$\\omega$','$\\mu_{2}$','$\\alpha$'};  
    param_true = [0,1,0,0.1];     
end
% clear -regexp ^MSE ^mean ^draw


%% Create table
fname = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),...
    '_H',num2str(H),'_pcp_mc_varc_draws',suffix,'.tex'];
FID = fopen(fname, 'w+');

fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');
fprintf(FID, '\\begin{sidewaystable} \n');
fprintf(FID, '\\center \n');
fprintf(FID, '\\begin{tabular}{cc cc| cccc| cccc} \n');

fprintf(FID, ['Value && DGP & Posterior & CP0  & PCP0 & CP10\\%%  & PCP10\\%% &', ...
'CP var ah & PCP var ah & CP var mle & PCP var mle \\\\ \\hline \n']); 
for T = TT
    %   load variables
    %% VARC  
    varc = true;
    name = ['results\',model,'\',model,'_1_2_T',num2str(T),...
        '_H100_II10_PCP0_MC_(R2017a)_varc_low_es',suffix,'.mat'];
    load(name,'-regexp','^mean')
    load(name,'-regexp','^median')
    load(name,'-regexp','^std')

    %% TIME CONSTANT
    varc = false;
    name = ['results\',model,'\',model,'_1_2_T',num2str(T),...
        '_H100_II10_PCP0_MC_(R2017a)_low_es',suffix,'.mat'];
    load(name,'-regexp','^mean')
    load(name,'-regexp','^median')
    load(name,'-regexp','^std')
    
    %% Print table block
    fprintf(FID, '\\hline \n');
    fprintf(FID, ['\\multicolumn{12}{c}{$T =',num2str(T),'$}  \\\\ \n']);
    fprintf(FID, '\\hline \n');

    for ii = 1:length(parameters)
        param = parameters(ii);
        fprintf(FID,' \\rowcolor{LightCyan} \n' );
        fprintf(FID,[char(param),'&& %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f &', ...
            ' %6.4f & %6.4f & %6.4f & %6.4f  \\\\   \n'], ...
            param_true(ii), mean(mean_draw(:,ii)), ...
            mean(mean_draw_C0(:,ii)), mean(mean_draw_PC0(:,ii)), ...
            mean(mean_draw_C(:,ii)), mean(mean_draw_PC(:,ii)), ...
            mean(mean_draw_Cah(:,ii)), mean(mean_draw_PCah(:,ii)), ...
            mean(mean_draw_Cm(:,ii)), mean(mean_draw_PCm(:,ii)));  
        fprintf(FID,['median &&  & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] &', ...
            ' [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f]  \\\\ \n'], ...
            mean(median_draw(:,ii)), ...
            mean(median_draw_C0(:,ii)), mean(median_draw_PC0(:,ii)), ...
            mean(median_draw_C(:,ii)), mean(median_draw_PC(:,ii)), ...
            mean(median_draw_Cah(:,ii)), mean(median_draw_PCah(:,ii)), ...
            mean(median_draw_Cm(:,ii)), mean(median_draw_PCm(:,ii)));              
        fprintf(FID,['std &&   & (%6.4f) & (%6.4f) & (%6.4f) & (%6.4f) & (%6.4f) &', ...
            '(%6.4f) & (%6.4f) & (%6.4f) & (%6.4f)  \\\\ \n'], ...
            mean(std_draw(:,ii)), mean(std_draw_C0(:,ii)), mean(std_draw_PC0(:,ii)), ...
            mean(std_draw_C(:,ii)), mean(std_draw_PC(:,ii)), ... 
            mean(std_draw_Cah(:,ii)), mean(std_draw_PCah(:,ii)), ... 
            mean(std_draw_Cm(:,ii)), mean(std_draw_PCm(:,ii))); 
    end     
end      
fprintf(FID, '\\hline \n');

fprintf(FID, '%%\\multicolumn{12}{l}{\\footnotesize{CP: Censored posterior.}}  \\\\ \n');
fprintf(FID, '%%\\multicolumn{12}{l}{\\footnotesize{PCP: Partially censored posterior.}} \\\\ \n');
fprintf(FID, '%%\\multicolumn{12}{l}{\\footnotesize{(P)CP0: Censoring with threshold 0.}} \\\\ \n'); 
fprintf(FID, '%%\\multicolumn{12}{l}{\\footnotesize{(P)CP10\\%%: Censoring with threshold 10\\%% sample quantile.}}  \\\\ \n');
fprintf(FID, '%%\\multicolumn{12}{l}{\\footnotesize{(P)CP var ad: Time Varying Censoring, ad hoc method.}} \\\\ \n'); 
fprintf(FID, '%%\\multicolumn{12}{l}{\\footnotesize{(P)CP var mle: Time Varying Censoring, MLE based method.}}  \\\\ \n');    
fprintf(FID, '\\end{tabular}\n ');

caption = ['\\caption{Draws statistics (means, medians, standard deviations) for standard posterior, censored posterior and partially censored posterior',...
    ' (the latter two with two time-constant and two time-varying thresholds) for the ',...
    model,' zero mean split normal model with $\\sigma_{1} = ',num2str(sigma1),'$ and $\\sigma_{2} = ',num2str(sigma2),'$.',...
    ' For the censored and the partially censored posterior the focus is on the left tail.',...            
    ' Averages over 50 simulations of the simulation averages over 10,000 draws.} \n'];        

fprintf(FID, caption);

label = ['\\label{tab:',model,'_pcp_draws}  \n'];
fprintf(FID, label);

fprintf(FID, '\\end{sidewaystable}\n');

fprintf(FID, '}');
fclose(FID);        
  