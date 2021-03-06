function print_table_pcp_mc(model,parameters,sigma1,sigma2,H, varc)
%     addpath(genpath('include/'));
%     model = 'ar1'; sigma1=1; sigma2=2; H=100; parameters ={'$\\mu$','$\\sigma$','$\\phi$'}; varc  = 1;
    rowcolor = @(ff) fprintf(ff, '\\rowcolor{LightCyan} \n');
 
    %% ESTIMATION
    if varc
        VV = '_varc';
    else
        VV = '';
    end
    fname = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_H',num2str(H),'_pcp_mc_est',VV,'.tex'];
    
    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');
    fprintf(FID, '{\\footnotesize \n');
    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, '\\begin{tabular}{cc cccccc} \n');

    if varc
        fprintf(FID, 'Value & & DGP & Posterior & CP$_{var,mf}$ & PCP$_{var,mf}$ & CP$_{var,mle}$ & PCP$_{var,mle }$\\\\ \\hline \n'); 
    else
        fprintf(FID, 'Value & & DGP & Posterior & CP0  & CP10\\%% & PCP0  & PCP10\\%% \\\\ \\hline \n'); 
    end
    
    for T = [100,1000,10000]
%         if (sigma2 == 2)
            name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_II10_PCP0_MC_(R2017a)',VV,'.mat'];
%         else
%             name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_II10_PCP0_MC_(R2015a).mat'];
%         end
        clear '-regexp' '^mean' '^std' '^VaR' '^ES' 'MSE' '^cdf''param_true'
        load(name,'-regexp','^mean','^std','^VaR','^ES','MSE','^cdf','param_true')
        MSE_ES;
        
        fprintf(FID, '\\hline \n');
        fprintf(FID, ['\\multicolumn{8}{c}{$T =',num2str(T),'$}  \\\\ \n']);
        fprintf(FID, '\\hline \n');
    
%         fprintf(FID,'AR && -- & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\  \n', ...
%             mean(accept), mean(accept_C0), mean(accept_C), mean(accept_PC0),mean(accept_PC));        
     
        for ii = 1:length(parameters)
            rowcolor(FID);        
            param = parameters(ii);
            if varc
                fprintf(FID,[char(param),'&& %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n'], ...
                    param_true(ii), mean(mean_draw(:,ii)), mean(mean_draw_Cah(:,ii)),...
                    mean(mean_draw_PCah(:,ii)), mean(mean_draw_Cm(:,ii)), mean(mean_draw_PCm(:,ii)));        
                fprintf(FID,'&&   & (%6.4f) & (%6.4f) & (%6.4f) & (%6.4f) & (%6.4f) \\\\ [1ex]\n', ...
                    mean(std_draw(:,ii)), mean(std_draw_Cah(:,ii)), mean(std_draw_PCah(:,ii)), ...
                    mean(std_draw_Cm(:,ii)), mean(std_draw_PCm(:,ii)));
            else
                fprintf(FID,[char(param),'&& %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n'], ...
                    param_true(ii), mean(mean_draw(:,ii)), mean(mean_draw_C0(:,ii)),...
                    mean(mean_draw_C(:,ii)), mean(mean_draw_PC0(:,ii)), mean(mean_draw_PC(:,ii)));        
                fprintf(FID,'&&   & (%6.4f) & (%6.4f) & (%6.4f) & (%6.4f) & (%6.4f) \\\\ [1ex]\n', ...
                    mean(std_draw(:,ii)), mean(std_draw_C0(:,ii)), mean(std_draw_C(:,ii)), ...
                    mean(std_draw_PC0(:,ii)), mean(std_draw_PC(:,ii)));
            end
        end    
    end
    
    fprintf(FID, '\\hline \n');
    fprintf(FID, '\\end{tabular}\n ');
    
    if sigma2 == 1
        tttt = 'Symmetric (correctly specified) ';
    else
        tttt = 'Asymmetric (misspecified) ';
    end
    if varc
        caption = ['\\caption{',tttt, model,' zero mean split normal model with $\\sigma_{1} = ', ...
                num2str(sigma1),'$ and $\\sigma_{2} = ',num2str(sigma2),'$: ',...
                ' simulation results for the regular posterior, censored posterior (CP) ', ...
                ' and partially censored posterior (PCP) with different time-varying thresholds, ', ...
                ' the mode-free one (CP$_{var,mf}$ and PCP$_{var,mf}$)',...
                ' and  the MLE-based one (CP$_{var,mle}$ and PCP$_{var,mle }). ' ,...
                ' For the censored methods the focus is on the left tail. ', ...
                ' All results averaged 50 MC replications, (mean) standard errors in parentheses.', ...
                '} \n'];        
    else
        caption = ['\\caption{',tttt, model,' zero mean split normal model with $\\sigma_{1} = ', ...
                num2str(sigma1),'$ and $\\sigma_{2} = ',num2str(sigma2),'$: ',...
                ' simulation results for the regular posterior, censored posterior (CP) ', ...
                ' and partially censored posterior (PCP) with different thresholds, ', ...
                ' at $0$ (CP0, PCP0) and at the 10\\%% data percentile (CP10\\%%, PCP 10\\%%). ' ,...
                ' For the censored methods the focus is on the left tail. ', ...
                ' All results averaged 100 MC replications, (mean) standard errors in parentheses.', ...
                '} \n'];
    end
    fprintf(FID, caption);

    label = ['\\label{tab:',model,'_s',num2str(sigma2),'_pcp_est',VV,'}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fprintf(FID, '}');
    fclose(FID);
    
    
    %% VaR & ES
    fname = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_H',num2str(H),'_pcp_mc_var_es',VV,'.tex'];    
    
    FID = fopen(fname, 'w+');
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');
    fprintf(FID, '{\\footnotesize \n');
    fprintf(FID, '\\begin{table} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, '\\begin{tabular}{cc cccccc} \n');

    if varc
        fprintf(FID, 'Value & & MC mean & Posterior & CP$_{var,mf}$ & PCP$_{var,mf}$ & CP$_{var,mle}$ & PCP$_{var,mle }$\\\\ \\hline \n'); 
    else
        fprintf(FID, 'Value & & MC mean & Posterior & CP0  & CP10\\%% & PCP0  & PCP10\\%% \\\\ \\hline \n'); 
    end
     
    for T = [100,1000,10000]
        name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_II10_PCP0_MC_(R2017a)',VV,'.mat'];
        clear '-regexp' '^mean' '^std' '^VaR' '^ES' 'MSE' '^cdf''param_true'
        load(name,'-regexp','^mean','^std','^VaR','^ES','MSE','^cdf','param_true')
        MSE_ES;
        
        fprintf(FID, '\\hline \n');
        fprintf(FID, ['\\multicolumn{8}{c}{$T =',num2str(T),'$}  \\\\ \n']);
        fprintf(FID, '\\hline \n');

        if varc
            rowcolor(FID);
            fprintf(FID,'VaR 0.5\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
                mean(mean(VaR_05)), mean(mean(VaR_05_post)), mean(mean(VaR_05_post_Cah)), ...
                mean(mean(VaR_05_post_PCah)), mean(mean(VaR_05_post_Cm)), mean(mean(VaR_05_post_PCm)));                     
            fprintf(FID,'  && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ \n', ...
                mean(MSE_05), mean(MSE_05_post), mean(MSE_05_post_Cah), mean(MSE_05_post_PCah), ...
                mean(MSE_05_post_Cm), mean(MSE_05_post_PCm));                         
             rowcolor(FID);
             fprintf(FID,'ES 0.5\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
                mean(mean(ES_05)), mean(mean(ES_05_post)), mean(mean(ES_05_post_Cah)), ...
                mean(mean(ES_05_post_PCah)), mean(mean(ES_05_post_Cm)), mean(mean(ES_05_post_PCm)));                      
            fprintf(FID,'  && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ [1ex]\n', ...
                mean(MSE_es_05), mean(MSE_es_05_post), mean(MSE_es_05_post_Cah), mean(MSE_es_05_post_PCah), ...
                mean(MSE_es_05_post_Cm), mean(MSE_es_05_post_PCm));                     

            rowcolor(FID);        
            fprintf(FID,'VaR 1\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
                mean(mean(VaR_1)), mean(mean(VaR_1_post)), mean(mean(VaR_1_post_Cah)), ...
                mean(mean(VaR_1_post_PCah)), mean(mean(VaR_1_post_Cm)), mean(mean(VaR_1_post_PCm)));                      
            fprintf(FID,'  && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ \n', ...
                mean(MSE_1), mean(MSE_1_post), mean(MSE_1_post_Cah), mean(MSE_1_post_PCah), ...
                mean(MSE_1_post_Cm), mean(MSE_1_post_PCm));        
            rowcolor(FID);
            fprintf(FID,'ES 1\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
                mean(mean(ES_1)), mean(mean(ES_1_post)), mean(mean(ES_1_post_Cah)), ...
                mean(mean(ES_1_post_PCah)), mean(mean(ES_1_post_Cm)), mean(mean(ES_1_post_PCm)));                      
            fprintf(FID,'  && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ [1ex]\n', ...
                mean(MSE_es_1), mean(MSE_es_1_post), mean(MSE_es_1_post_Cah), mean(MSE_es_1_post_PCah), ...
                mean(MSE_es_1_post_Cm), mean(MSE_es_1_post_PCm));        

            rowcolor(FID);         
            fprintf(FID,'VaR 5\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
                mean(mean(VaR_5)), mean(mean(VaR_5_post)), mean(mean(VaR_5_post_Cah)), ...
                mean(mean(VaR_5_post_PCah)), mean(mean(VaR_5_post_Cm)), mean(mean(VaR_5_post_PCm)));                                   
            fprintf(FID,' && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ \n', ...
                mean(MSE_5), mean(MSE_5_post), mean(MSE_5_post_Cah), mean(MSE_5_post_PCah), ...
                mean(MSE_5_post_Cm), mean(MSE_5_post_PCm));        
            rowcolor(FID);
            fprintf(FID,'ES 5\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
                mean(mean(ES_5)), mean(mean(ES_5_post)), mean(mean(ES_5_post_Cah)), ...
                mean(mean(ES_5_post_PCah)), mean(mean(ES_5_post_Cm)), mean(mean(ES_5_post_PCm)));                               
            fprintf(FID,' && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ \n', ...
                mean(MSE_es_5), mean(MSE_es_5_post), mean(MSE_es_5_post_Cah), mean(MSE_es_5_post_PCah), ...
                mean(MSE_es_5_post_Cm), mean(MSE_es_5_post_PCm)); 
        else
            rowcolor(FID);
            fprintf(FID,'VaR 0.5\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
                mean(mean(VaR_05)), mean(mean(VaR_05_post)), mean(mean(VaR_05_post_C0)), ...
                mean(mean(VaR_05_post_C)), mean(mean(VaR_05_post_PC0)), mean(mean(VaR_05_post_PC)));                     
            fprintf(FID,'  && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ \n', ...
                mean(MSE_05), mean(MSE_05_post), mean(MSE_05_post_C0), mean(MSE_05_post_C), ...
                mean(MSE_05_post_PC0), mean(MSE_05_post_PC));                         
             rowcolor(FID);
             fprintf(FID,'ES 0.5\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
                mean(mean(ES_05)), mean(mean(ES_05_post)), mean(mean(ES_05_post_C0)), ...
                mean(mean(ES_05_post_C)), mean(mean(ES_05_post_PC0)), mean(mean(ES_05_post_PC)));                      
            fprintf(FID,'  && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ [1ex]\n', ...
                mean(MSE_es_05), mean(MSE_es_05_post), mean(MSE_es_05_post_C0), mean(MSE_es_05_post_C), ...
                mean(MSE_es_05_post_PC0), mean(MSE_es_05_post_PC));                     

            rowcolor(FID);        
            fprintf(FID,'VaR 1\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
                mean(mean(VaR_1)), mean(mean(VaR_1_post)), mean(mean(VaR_1_post_C0)), ...
                mean(mean(VaR_1_post_C)), mean(mean(VaR_1_post_PC0)), mean(mean(VaR_1_post_PC)));                      
            fprintf(FID,'  && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ \n', ...
                mean(MSE_1), mean(MSE_1_post), mean(MSE_1_post_C0), mean(MSE_1_post_C), ...
                mean(MSE_1_post_PC0), mean(MSE_1_post_PC));        
            rowcolor(FID);
            fprintf(FID,'ES 1\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
                mean(mean(ES_1)), mean(mean(ES_1_post)), mean(mean(ES_1_post_C0)), ...
                mean(mean(ES_1_post_C)), mean(mean(ES_1_post_PC0)), mean(mean(ES_1_post_PC)));                      
            fprintf(FID,'  && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ [1ex]\n', ...
                mean(MSE_es_1), mean(MSE_es_1_post), mean(MSE_es_1_post_C0), mean(MSE_es_1_post_C), ...
                mean(MSE_es_1_post_PC0), mean(MSE_es_1_post_PC));        

            rowcolor(FID);         
            fprintf(FID,'VaR 5\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
                mean(mean(VaR_5)), mean(mean(VaR_5_post)), mean(mean(VaR_5_post_C0)), ...
                mean(mean(VaR_5_post_C)), mean(mean(VaR_5_post_PC0)), mean(mean(VaR_5_post_PC)));                                   
            fprintf(FID,' && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ \n', ...
                mean(MSE_5), mean(MSE_5_post), mean(MSE_5_post_C0), mean(MSE_5_post_C), ...
                mean(MSE_5_post_PC0), mean(MSE_5_post_PC));        
            rowcolor(FID);
            fprintf(FID,'ES 5\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', ...
                mean(mean(ES_5)), mean(mean(ES_5_post)), mean(mean(ES_5_post_C0)), ...
                mean(mean(ES_5_post_C)), mean(mean(ES_5_post_PC0)), mean(mean(ES_5_post_PC)));                               
            fprintf(FID,' && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] \\\\ \n', ...
                mean(MSE_es_5), mean(MSE_es_5_post), mean(MSE_es_5_post_C0), mean(MSE_es_5_post_C), ...
                mean(MSE_es_5_post_PC0), mean(MSE_es_5_post_PC)); 
        end
    end
    
    fprintf(FID, '\\hline \n');
    fprintf(FID, '\\end{tabular}\n ');
    
    if sigma2 == 1
        tttt = 'Symmetric (correctly specified) ';
    else
        tttt = 'Asymmetric (misspecified) ';
    end
    
    if varc
        caption = ['\\caption{',tttt, model,' zero mean split normal model with $\\sigma_{1} = ', ...
                num2str(sigma1),'$ and $\\sigma_{2} = ',num2str(sigma2),'$: ',...
                ' tail forecasting results for the regular posterior, censored posterior (CP) ', ...
                ' and partially censored posterior (PCP) with different time-varying thresholds, ', ...
                ' the mode-free one (CP$_{var,mf}$ and PCP$_{var,mf}$)',...
                ' and  the MLE-based one (CP$_{var,mle}$ and PCP$_{var,mle }). ' ,...
                ' For the censored methods the focus is on the left tail. ', ...
                ' Results averaged 100 MC replications and over the out-of-sample horizon of $H=100$, ', ...
                ' with (mean) MSEs in brackets.',...
                '} \n']; 
    else
        caption = ['\\caption{',tttt, model,' zero mean split normal model with $\\sigma_{1} = ', ...
                num2str(sigma1),'$ and $\\sigma_{2} = ',num2str(sigma2),'$: ',...
                ' tail forecasting results for the regular posterior, censored posterior (CP) ', ...
                ' and partially censored posterior (PCP) with different thresholds, ', ...
                ' at $0$ (CP0, PCP0) and at the 10\\%% data percentile (CP10\\%%, PCP 10\\%%). ' ,...
                ' For the censored methods the focus is on the left tail. ', ...
                ' Results averaged 100 MC replications and over the out-of-sample horizon of $H=100$, ', ...
                ' with (mean) MSEs in brackets.',...
                '} \n'];        
    end
    fprintf(FID, caption);

    label = ['\\label{tab:',model,'_s',num2str(sigma2),'_pcp_var_es',VV,'}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{table}\n');
    
    fprintf(FID, '}');
    fprintf(FID, '}');
    fclose(FID);    
end