function fname = print_table_DM(DM, model, p_bar, methods,v)

    fname = ['results/',model,'/',model,'_DM',v,'_',num2str(p_bar),'.tex'];

%     methods = {'True','Posterior','CP0','PCP0','CP10\%','PCP10\%'};
    met = length(methods);
    
    STAR = cell(met,met);
    [STAR{:}] = deal({'\phantom{$^{***}$}'});
    
    alpha = 0.10;
    SIG = logical(abs(DM) > norminv(1-alpha/2)); % 1.6449
    STAR(SIG) = {'$^{*}$\phantom{$^{**}$}'};

    alpha = 0.05;
    SIG = logical(abs(DM) > norminv(1-alpha/2)); % 1.9600
    STAR(SIG) = {'$^{**}$\phantom{$^{*}$}'};

    alpha = 0.01;
    SIG = logical(abs(DM) > norminv(1-alpha/2)); % 2.5758
    STAR(SIG) = {'$^{***}$'};
 

    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n {\\footnotesize \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, ['\\begin{tabular}{l | ', repelem('r',met), '} \n']);

        for ii = 1:met
            fprintf(FID, '& \\multicolumn{1}{c}{%s}', methods{ii});
        end
%         '& \\multicolumn{1}{c}{True}',...
%         ' & \\multicolumn{1}{c}{Posterior} & \\multicolumn{1}{c}{CP0}',...
%         ' & \\multicolumn{1}{c}{PCP0} & \\multicolumn{1}{c}{CP10\\%%}',...
%         ' & \\multicolumn{1}{c}{PCP10\\%%} 
    fprintf(FID,' \\\\ \\hline \n');
 
    fprintf(FID, ['\\multicolumn{',num2str(met+1),'}{c}{$\\tau = ',num2str(p_bar),'$} \\\\ \\hline \n']);
    
    
    for ii = 1:met
        fprintf(FID,'%s & ',methods{ii});
        for jj = 1:met
            if (jj == met)
                fprintf(FID,'%6.4f%s  ',DM(ii,jj),char(STAR{ii,jj})); 
            else
                fprintf(FID,'%6.4f%s & ',DM(ii,jj),char(STAR{ii,jj})); 
            end
        end
        fprintf(FID,' \\\\ \n');        
    end
     
    fprintf(FID, '\\hline \n');
    
%     fprintf(FID, ['\\multicolumn{',num2str(met+1),'}{l}{\\footnotesize{CP: Censored posterior.}}  \\\\ \n']);
%     fprintf(FID, ['\\multicolumn{',num2str(met+1),'}{l}{\\footnotesize{PCP: Partially censored posterior.}} \\\\ \n']);
%     fprintf(FID, ['\\multicolumn{',num2str(met+1),'}{l}{\\footnotesize{(P)CP: Censoring with constant threshold at 10\\%% sample quantile.}}  \\\\ \n']);
%     fprintf(FID, ['\\multicolumn{',num2str(met+1),'}{l}{\\footnotesize{(P)CP$_{var}$: Censoring with time-varying threshold at 10\\%% MLE implied quantile.}} \\\\ \n']); 
    
    fprintf(FID, '\\end{tabular}\n ');
    
%     if (met == 5)
%         ttt = 'and a time-varying one based on the 10\\%% MLE implied quantile';
%     elseif (met == 7)
%         ttt = 'and two time-varying ones, model-free and based on the 10\\%% MLE implied quantile';
%     end
    
    if isempty(v)
        vvv = ' VaR ';
    else
        vvv = ' ES';
    end
        
    caption = ['\\caption{',model,...
        ': Diebold-Mariano test statistics for pairwise method comparison ', ...
        ' of forecasting of the $\\tau$\\%% ', vvv,' based on ', ...
        ' the RMSE for $H=100$ out-of-sample periods, ', ...
        ' between  the regular posterior, censored posterior (CP) and ', ...
        ' partially censored posterior (PCP) with different thresholds, ', ...
        ' a time-constant one at the 10\\%% data percentile (CP10\\%%, PCP 10\\%%). }\n'];
    fprintf(FID, caption);

    label = ['\\label{tab:',model,'_DM_',v,'_',num2str(p_bar),'}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}}');
    fclose(FID);
    
    
    Remove_NaN(fname);

end