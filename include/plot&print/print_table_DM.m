function fname = print_table_DM(DM, model, T, S)
    fname = ['results/',model,'/',model,'_DM_T',num2str(T),'.tex'];
    
    methods = {'True','Posterior','CP0','PCP0','CP10\%','PCP10\%'};

    STAR = cell(6,6);
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
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, '\\begin{tabular}{l | rrr rrr} \n');

    fprintf(FID, ['VaR 1\\%% $\\setminus$ 5\\%% & \\multicolumn{1}{c}{True}',...
        ' & \\multicolumn{1}{c}{Posterior} & \\multicolumn{1}{c}{CP0}',...
        ' & \\multicolumn{1}{c}{PCP0} & \\multicolumn{1}{c}{CP10\\%%}',...
        ' & \\multicolumn{1}{c}{PCP10\\%%} \\\\ \\hline \n']);  
    for ii = 1:6
        fprintf(FID,'%s & ',methods{ii});
        for jj = 1:6
            if (jj == 6)
                fprintf(FID,'%6.4f%s ',DM(ii,jj),char(STAR{ii,jj})); 
            else
                fprintf(FID,'%6.4f%s & ',DM(ii,jj),char(STAR{ii,jj})); 
            end
        end
        fprintf(FID,' \\\\ \n');        
    end
     
    fprintf(FID, '\\hline \n');
    
    fprintf(FID, '\\multicolumn{7}{l}{\\footnotesize{CP: Censored posterior.}}  \\\\ \n');
    fprintf(FID, '\\multicolumn{7}{l}{\\footnotesize{PCP: Partially censored posterior.}} \\\\ \n');
    fprintf(FID, '\\multicolumn{7}{l}{\\footnotesize{(P)CP0: Censoring with threshold 0.}} \\\\ \n'); 
    fprintf(FID, '\\multicolumn{7}{l}{\\footnotesize{(P)CP10\\%%: Censoring with threshold 10\\%% sample quantile.}}  \\\\ \n');
    
    fprintf(FID, '\\end{tabular}\n ');
    
    caption = ['\\caption{Diebold-Mariano test statistics for ',...
        ' pairwise method comparision of forecasting performace based on ',...
        'loss diffential vectors of length $S = ',num2str(S),...
        '$, with the loss function ',...
        'defined as the RMSE over $H=100$ out-of-sample periods. ',...
        'In-sample period: $T = ',num2str(T),'$.} \n'];
    fprintf(FID, caption);

    label = ['\\label{tab:',model,'_DM_T_',num2str(T),'}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fclose(FID);
end