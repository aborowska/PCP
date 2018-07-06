function print_table_estim_emp(data_name,parameters)

%     res_name = ['results/tgarch11/',data_name,'/PCP_emp_tgarch11_data_',data_name,'.mat'];
    res_name = ['results/tgarch11/',data_name,'/PCP_emp_tgarch11_data_',data_name,'_with_ah.mat'];
    load(res_name,'-regexp','^draw');
    
    fname = ['results/tgarch11/',data_name,'/',data_name,'_estim_results2.tex'];
    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
%     fprintf(FID, '\\begin{tabular}{cc ccccc} \n');
    fprintf(FID, '\\begin{tabular}{cc ccccccc} \n');

%    fprintf(FID, 'Parameter & & Posterior & CP & CP$_{var}$& PCP & PCP$_{var}$ \\\\ \\hline \n'); 
    fprintf(FID, ['Parameter & & Posterior & CP10\\%% & PCP10\\%% ', ...
        ' & CP$_{var,mf}$ & PCP$_{var,mf}$ & CP$_{var,mle}$ & PCP$_{var,mle}$ \\\\ \\hline \n']); 
    

    for ii = 1:length(parameters)
        param = parameters(ii);
        fprintf(FID,'\\rowcolor{LightCyan} \n');
%         fprintf(FID,[char(param),'&& %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n'], ...
%             mean(draw(:,ii)), mean(draw_C(:,ii)), mean(draw_Cm(:,ii)), ...
%             mean(draw_PC(:,ii)), mean(draw_PCm(:,ii)));        
%         fprintf(FID,'&&    (%6.4f) & (%6.4f) & (%6.4f)  & (%6.4f) & (%6.4f)  \\\\ \n', ...
%             std(draw(:,ii)), std(draw_C(:,ii)), std(draw_Cm(:,ii)), ...
%             std(draw_PC(:,ii)), std(draw_PCm(:,ii))); 
        fprintf(FID,[char(param),'&& %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n'], ...
            mean(draw(:,ii)), mean(draw_C(:,ii)), mean(draw_PC(:,ii)), ...
            mean(draw_Cah(:,ii)), mean(draw_PCah(:,ii)), ...        
            mean(draw_Cm(:,ii)), mean(draw_PCm(:,ii)));        
        fprintf(FID,'&&    (%6.4f) & (%6.4f) & (%6.4f)  & (%6.4f) & (%6.4f) & (%6.4f) & (%6.4f)  \\\\ \n', ...
            std(draw(:,ii)), std(draw_C(:,ii)), std(draw_PC(:,ii)), ...
            std(draw_Cah(:,ii)), std(draw_PCah(:,ii)), ...         
            std(draw_Cm(:,ii)), std(draw_PCm(:,ii)));         
    end
 
    
    fprintf(FID, '\\hline \n');
%     fprintf(FID, ['\\multicolumn{7}{l}{\\footnotesize{CP: Censored posterior.}}  \\\\ \n']);
%     fprintf(FID, ['\\multicolumn{7}{l}{\\footnotesize{PCP: Partially censored posterior.}} \\\\ \n']);
%     fprintf(FID, ['\\multicolumn{7}{l}{\\footnotesize{(P)CP: Censoring with constant threshold at 10\\%% sample quantile.}}  \\\\ \n']);
%     fprintf(FID, ['\\multicolumn{7}{l}{\\footnotesize{(P)CP$_{var}$: Censoring with time-varying threshold at 10\\%% MLE implied quantile.}} \\\\ \n']); 
    
    fprintf(FID, '\\end{tabular}\n ');
    
    data = data_name(1:(strfind(data_name,'_')-1));
      
    caption = ['\\caption{',data,': estimation results (posterior means and standard deviations) ',...
        ' for the GARCH(1,1)-$t$ model with  different estimation methods.',...
        ' For the censored methods the focus is on the left tail. }\n'];
    fprintf(FID, caption);

    label = ['\\label{tab:',data,'_estim}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fclose(FID);
end