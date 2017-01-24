function print_table_pcp_mc(model,parameters,sigma1,sigma2)
    fname = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_pcp_mc.tex'];
    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, '\\begin{tabular}{cc cccccc} \n');

    fprintf(FID, 'Value & & True/MC$^*$ & Posterior & CP0  & CP10\\%% & PCP0  & PCP10\\%% \\\\ \\hline \n'); 
    
    for T = [100,1000,10000]
        load(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_PCP0_MC.mat'])
        fprintf(FID, '\\hline \n');
        fprintf(FID, ['\\multicolumn{8}{c}{$T =',num2str(T),'$}  \\\\ \n']);
        fprintf(FID, '\\hline \n');
    
        fprintf(FID,'AR && -- & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\  \n', ...
            mean(accept), mean(accept_C0), mean(accept_C), mean(accept_PC0),mean(accept_PC));        

        for ii = 1:length(parameters)
             param = parameters(ii);
             fprintf(FID,[char(param),'&& %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n'], ...
                 param_true(ii), mean(draw(:,ii)), mean(draw_C0(:,ii)), mean(draw_C(:,ii)), mean(draw_PC0(:,ii)), mean(draw_PC(:,ii)));        
             fprintf(FID,'&&   & (%6.4f) & (%6.4f) & (%6.4f) & (%6.4f) & (%6.4f) \\\\ \n', ...
                 std(draw(:,ii)), std(draw_C0(:,ii)), std(draw_C(:,ii)), std(draw_PC0(:,ii)), std(draw_PC(:,ii)));        
        end
        fprintf(FID,'VaR 1\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
            mean(VaR_1), mean(VaR_1_post), mean(VaR_1_post_C0), mean(VaR_1_post_C), mean(VaR_1_post_PC0), mean(VaR_1_post_PC));        
        fprintf(FID,'  && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ \n', ...
            mean(MSE_1), mean(MSE_1_post), mean(MSE_1_post_C0), mean(MSE_1_post_C), mean(MSE_1_post_PC0), mean(MSE_1_post_PC));        
        fprintf(FID,'VaR 5\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n',  ...
            mean(VaR_5), mean(VaR_5_post), mean(VaR_5_post_C0), mean(VaR_5_post_C), mean(VaR_5_post_PC0), mean(VaR_5_post_PC));               
        fprintf(FID,' && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ \n', ...
            mean(MSE_5), mean(MSE_5_post), mean(MSE_5_post_C0), mean(MSE_5_post_C), mean(MSE_5_post_PC0), mean(MSE_5_post_PC));        
    end
    
    fprintf(FID, '\\hline \n');
    fprintf(FID, '\\multicolumn{6}{l}{\\footnotesize{$^*$ True value for parameters, MC mean value for VaRs.}}  \\\\ \n');
    fprintf(FID, '\\multicolumn{6}{l}{\\footnotesize{AR: acceptance rate for the independent MH.}}  \\\\ \n');
    fprintf(FID, '\\multicolumn{6}{l}{\\footnotesize{CP: Censored posterior, threshold 10\\%% sample quantile.}}  \\\\ \n');
    fprintf(FID, '\\multicolumn{6}{l}{\\footnotesize{PCP: Partially censored posterior,  threshold 10\\%% sample quantile.}} \\\\ \n');
    fprintf(FID, '\\multicolumn{6}{l}{\\footnotesize{(P)CP0: Censoring with threshold 0.}} \\\\ \n'); 
    fprintf(FID, '\\multicolumn{6}{l}{\\footnotesize{(P)CP10\\%%: Censoring with threshold 10\\%% sample quantile.}}  \\\\ \n');
    fprintf(FID, '\\end{tabular}\n ');
    
    caption = ['\\caption{Simulation results for standard posterior, censored posterior and partially censored posterior',...
        ' (the latter two with two threshold values) for the ',...
        model,' zero mean split normal model with $\\sigma_{1} = ',num2str(sigma1),'$ and $\\sigma_{2} = ',num2str(sigma2),'$.',...
        ' For the censored and the partially censored posterior the focus is on the left tail.',...
        ' Standard errors in parentheses, MSEs in brackets.} \n'];
    fprintf(FID, caption);

    label = ['\\label{tab:',model,'_pcp}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fclose(FID);
end