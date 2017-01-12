function print_table_cp(model,parameters)
    fname = ['results/',model,'.tex'];
    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, '\\begin{tabular}{cc cccc} \n');

    fprintf(FID, 'Value & & True & Posterior & CP10\\%% & CP0 \\\\ \\hline \n'); 

    for T = [100,1000,1000]
        load(['results/',model,'_',num2str(T),'.mat']);
        fprintf(FID, '\\hline \n');
        fprintf(FID, ['\\multicolumn{6}{c}{$T =',num2str(T),'$}  \\\\ \n']);
        fprintf(FID, '\\hline \n');

        for ii = 1:length(parameters)
             param = parameters(ii);
             fprintf(FID,[char(param),'&& %6.4f & %6.4f & %6.4f & %6.4f  \\\\ \n'], param_true(ii), mean(draw(:,ii)), mean(draw_C(:,ii)), mean(draw_C0(:,ii)));        
             fprintf(FID,'&&   & (%6.4f) & (%6.4f) & (%6.4f)  \\\\ \n', std(draw(:,ii)), std(draw_C(:,ii)), std(draw_C0(:,ii)));        
        end
        fprintf(FID,'VaR 1\\%% && %6.4f & %6.4f & %6.4f & %6.4f  \\\\ \n', VaR_1,VaR_1_post,VaR_1_post_C,VaR_1_post_C0);        
        fprintf(FID,'VaR 5\\%% && %6.4f & %6.4f & %6.4f & %6.4f  \\\\ \n',  VaR_5,VaR_5_post,VaR_5_post_C,VaR_5_post_C0);               
    end
    
    fprintf(FID, '\\hline \n');
    
    fprintf(FID, '\\multicolumn{6}{l}{\\footnotesize{CP10\\%%: Censored posterior, threshold 10\\%% sample quantile.}}  \\\\ \n');
    fprintf(FID, '\\multicolumn{6}{l}{\\footnotesize{CP0: Censored posterior, threshold 0.}} \n');
    
    fprintf(FID, '\\end{tabular}\n ');
    
    caption = ['\\caption{Simulation results for posterior without and with censoring (with different thresholds) for the ',...
        model,' mean zero split normal model. For the censored posterior the focus is on the left tail. } \n'];
    fprintf(FID, caption);

    label = ['\\label{tab:',model,'}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fclose(FID);
end