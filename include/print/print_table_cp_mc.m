function print_table_cp_mc(model,parameters,sigma1,sigma2)
    fname = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_mc.tex'];
    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, '\\begin{tabular}{cc cccc} \n');

    fprintf(FID, 'Value & & True/MC$^*$ & Posterior & CP10\\%% & CP0 \\\\ \\hline \n'); 
    
    for T = [100,1000,10000]
        load(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_MC.mat'])
        fprintf(FID, '\\hline \n');
        fprintf(FID, ['\\multicolumn{6}{c}{$T =',num2str(T),'$}  \\\\ \n']);
        fprintf(FID, '\\hline \n');
    
        fprintf(FID,'AR && -- & %6.4f & %6.4f & %6.4f  \\\\  \n', accept, accept_C,accept_C0);        

        for ii = 1:length(parameters)
             param = parameters(ii);
             fprintf(FID,[char(param),'&& %6.4f & %6.4f & %6.4f & %6.4f  \\\\ \n'], param_true(ii), mean(draw(:,ii)), mean(draw_C(:,ii)), mean(draw_C0(:,ii)));        
             fprintf(FID,'&&   & (%6.4f) & (%6.4f) & (%6.4f)  \\\\ \n', std(draw(:,ii)), std(draw_C(:,ii)), std(draw_C0(:,ii)));        
        end
        fprintf(FID,'VaR 1\\%% && %6.4f & %6.4f & %6.4f & %6.4f  \\\\ \n', mean(VaR_1),mean(VaR_1_post),mean(VaR_1_post_C),mean(VaR_1_post_C0));        
        fprintf(FID,'  && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f]  \\\\ \n', mean(MSE_1),mean(MSE_1_post),mean(MSE_1_post_C),mean(MSE_1_post_C0));        
        fprintf(FID,'VaR 5\\%% && %6.4f & %6.4f & %6.4f & %6.4f  \\\\ \n',  mean(VaR_5),mean(VaR_5_post),mean(VaR_5_post_C),mean(VaR_5_post_C0));               
        fprintf(FID,' && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f]  \\\\ \n', mean(MSE_5),mean(MSE_5_post),mean(MSE_5_post_C),mean(MSE_5_post_C0));        
    end
    
    fprintf(FID, '\\hline \n');
    fprintf(FID, '\\multicolumn{6}{l}{\\footnotesize{$^*$ True value for parameters, MC mean value for VaRs.}}  \\\\ \n');
    fprintf(FID, '\\multicolumn{6}{l}{\\footnotesize{AR: acceptance rate for the independent MH.}}  \\\\ \n');
    fprintf(FID, '\\multicolumn{6}{l}{\\footnotesize{CP10\\%%: Censored posterior, threshold 10\\%% sample quantile.}}  \\\\ \n');
    fprintf(FID, '\\multicolumn{6}{l}{\\footnotesize{CP0: Censored posterior, threshold 0.}} \n');
    
    fprintf(FID, '\\end{tabular}\n ');
    
    caption = ['\\caption{Simulation results for posterior without and with censoring (with different thresholds) for the ',...
        model,' zero mean split normal model with $\\sigma_{1} = ',num2str(sigma1),'$ and $\\sigma_{2} = ',num2str(sigma2),'$.',...
        ' For the censored posterior the focus is on the left tail.',...
        ' Standard errors in parantheses, MSEs in brackets.} \n'];
    fprintf(FID, caption);

    label = ['\\label{tab:',model,'}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fclose(FID);
end