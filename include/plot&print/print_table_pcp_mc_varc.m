function print_table_pcp_mc_varc(model,parameters,sigma1,sigma2,T,H, version)
    
    if nargin == 6
        version = '(R2017a)';
    end
       
    fname = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_H',num2str(H),'_pcp_mc_varc.tex'];
    FID = fopen(fname, 'w+');
    %% Create table
    fprintf(FID, '{ \\renewcommand{\\arraystretch}{1.2} \n');

    fprintf(FID, '\\begin{sidewaystable} \n');
    fprintf(FID, '\\center \n');
    fprintf(FID, '\\begin{tabular}{cc cc cccc cccc} \n');

    fprintf(FID, ['Value & & True/MC$^*$ & Posterior & CP0  & PCP0 & CP10\\%%  & PCP10\\%% &', ...
        'CP var mf & PCP var mf & CP var mle & PCP var mle \\\\ \\hline \n']); 
    
    for t = T
        % Load time varying threshold results
        load(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(t),'_H',num2str(H),'_II10_PCP0_MC_',version,'_varc.mat'])            
        
        %% Rename time varying threshold variables
        accept_Cah = accept_C;
        accept_PCah = accept_PC;
        mean_draw_Cah = mean_draw_C;
        mean_draw_PCah = mean_draw_PC;
        std_draw_Cah = std_draw_C;
        std_draw_PCah = std_draw_PC;        
        VaR_1_post_Cah = VaR_1_post_C;
        VaR_5_post_Cah = VaR_5_post_C;  
        VaR_1_post_PCah = VaR_1_post_PC;
        VaR_5_post_PCah = VaR_5_post_PC; 
        MSE_1_post_Cah = MSE_1_post_C;
        MSE_5_post_Cah = MSE_5_post_C;  
        MSE_1_post_PCah = MSE_1_post_PC;
        MSE_5_post_PCah = MSE_5_post_PC;       
        if ~exist('CV_Cm','var') 
            accept_Cm = accept_C0;
            accept_PCm = accept_PC0;
            mean_draw_Cm = mean_draw_C0;
            mean_draw_PCm = mean_draw_PC0;
            std_draw_Cm = std_draw_C0;
            std_draw_PCm = std_draw_PC0;
            VaR_1_post_Cm = VaR_1_post_C0;
            VaR_5_post_Cm = VaR_5_post_C0;  
            VaR_1_post_PCm = VaR_1_post_PC0;
            VaR_5_post_PCm = VaR_5_post_PC0; 
            MSE_1_post_Cm = MSE_1_post_C0;
            MSE_5_post_Cm = MSE_5_post_C0;  
            MSE_1_post_PCm = MSE_1_post_PC0;
            MSE_5_post_PCm = MSE_5_post_PC0;             
        end
        
        %% Load time constant threshold results

        load(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(t),'_H',num2str(H),'_II10_PCP0_MC_',version,'.mat'])            

        %% Fill-in table
        fprintf(FID, '\\hline \n');
        fprintf(FID, ['\\multicolumn{12}{c}{$T =',num2str(t),'$}  \\\\ \n']);
        fprintf(FID, '\\hline \n');
    
        fprintf(FID,['AR && -- & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f &', ...
            '%6.4f & %6.4f & %6.4f & %6.4f  \\\\  \n'], ...
            mean(accept), mean(accept_C0), mean(accept_PC0), mean(accept_C),mean(accept_PC), ...
            mean(accept_Cah), mean(accept_PCah), mean(accept_Cm),mean(accept_PCm));        

        for ii = 1:length(parameters)
            param = parameters(ii);
            fprintf(FID,[char(param),'&& %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f &', ...
                ' %6.4f & %6.4f & %6.4f & %6.4f  \\\\ \n'], ...
                param_true(ii), mean(mean_draw(:,ii)), mean(mean_draw_C0(:,ii)), mean(mean_draw_PC0(:,ii)), ...
                mean(mean_draw_C(:,ii)), mean(mean_draw_PC(:,ii)), ...
                mean(mean_draw_Cah(:,ii)), mean(mean_draw_PCah(:,ii)), ...
                mean(mean_draw_Cm(:,ii)), mean(mean_draw_PCm(:,ii)));        
            fprintf(FID,['&&   & (%6.4f) & (%6.4f) & (%6.4f) & (%6.4f) & (%6.4f) &', ...
                '(%6.4f) & (%6.4f) & (%6.4f) & (%6.4f)  \\\\ \n'], ...
                mean(std_draw(:,ii)), mean(std_draw_C0(:,ii)), mean(std_draw_PC0(:,ii)), ...
                mean(std_draw_C(:,ii)), mean(std_draw_PC(:,ii)), ... 
                mean(std_draw_Cah(:,ii)), mean(std_draw_PCah(:,ii)), ... 
                mean(std_draw_Cm(:,ii)), mean(std_draw_PCm(:,ii))); 
        end
          
        fprintf(FID,['VaR 1\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f &', ...
            '%6.4f & %6.4f & %6.4f & %6.4f \\\\ \n'], ...
            mean(mean(VaR_1)), mean(mean(VaR_1_post)), mean(mean(VaR_1_post_C0)), mean(mean(VaR_1_post_PC0)), ...
            mean(mean(VaR_1_post_C)), mean(mean(VaR_1_post_PC)), ...             
            mean(mean(VaR_1_post_Cah)), mean(mean(VaR_1_post_PCah)), ...           
            mean(mean(VaR_1_post_Cm)), mean(mean(VaR_1_post_PCm)));             

        fprintf(FID,['  && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] &', ...
            '[%6.4f] & [%6.4f] & [%6.4f] & [%6.4f]  \\\\ \n'], ...
            mean(MSE_1), mean(MSE_1_post), mean(MSE_1_post_C0), mean(MSE_1_post_PC0), ...
            mean(MSE_1_post_C), mean(MSE_1_post_PC), ...        
            mean(MSE_1_post_Cah), mean(MSE_1_post_PCah), ...     
            mean(MSE_1_post_Cm), mean(MSE_1_post_PCm));        

        fprintf(FID,['VaR 5\\%% && %6.4f & %6.4f & %6.4f & %6.4f & %6.4f & %6.4f &', ...
            '%6.4f & %6.4f & %6.4f & %6.4f \\\\ \n'], ...
            mean(mean(VaR_5)), mean(mean(VaR_5_post)), mean(mean(VaR_5_post_C0)), mean(mean(VaR_5_post_PC0)), ...
            mean(mean(VaR_5_post_C)), mean(mean(VaR_5_post_PC)), ...                          
            mean(mean(VaR_5_post_Cah)), mean(mean(VaR_5_post_PCah)), ...                       
            mean(mean(VaR_5_post_Cm)), mean(mean(VaR_5_post_PCm)));                          

        fprintf(FID,['  && [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] & [%6.4f] &', ...
            '[%6.4f] & [%6.4f] & [%6.4f] & [%6.4f]  \\\\ \n'], ...
            mean(MSE_5), mean(MSE_5_post), mean(MSE_5_post_C0), mean(MSE_5_post_PC0), ...
            mean(MSE_5_post_C), mean(MSE_5_post_PC), ...        
            mean(MSE_5_post_Cah), mean(MSE_5_post_PCah), ...     
            mean(MSE_5_post_Cm), mean(MSE_5_post_PCm));          
    end
    
    fprintf(FID, '\\hline \n');
    
    fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{$^*$ True value for parameters, MC mean value for the mean (over out-of-sample horizon) VaRs.}}  \\\\ \n');

    fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{AR: acceptance rate for the independent MH (M = 10,000, burn in of 1,000).}}  \\\\ \n');
    fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{CP: Censored posterior.}}  \\\\ \n');
    fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{PCP: Partially censored posterior.}} \\\\ \n');
    fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{(P)CP0: Censoring with threshold 0.}} \\\\ \n'); 
    fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{(P)CP10\\%%: Censoring with threshold 10\\%% sample quantile.}}  \\\\ \n');
    fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{(P)CP var ad: Time Varying Censoring, ad hoc method.}} \\\\ \n'); 
    fprintf(FID, '\\multicolumn{12}{l}{\\footnotesize{(P)CP var mle: Time Varying Censoring, MLE based method.}}  \\\\ \n');    
    fprintf(FID, '\\end{tabular}\n ');
    
    caption = ['\\caption{Simulation results for standard posterior, censored posterior and partially censored posterior',...
        ' (the latter two with two time-constant and two time-varying thresholds) for the ',...
        model,' zero mean split normal model with $\\sigma_{1} = ',num2str(sigma1),'$ and $\\sigma_{2} = ',num2str(sigma2),'$.',...
        ' For the censored and the partially censored posterior the focus is on the left tail.',...            
        ' All results averaged 100 MC replications.',...
        ' Results for the VaRs additionally averaged over out-of-sample horizon of $H=',num2str(H),'$.',...
        ' (Mean) standard errors in parentheses, (Mean) MSEs in brackets.} \n'];        

    fprintf(FID, caption);

    label = ['\\label{tab:',model,'_pcp_var}  \n'];
    fprintf(FID, label);
    
    fprintf(FID, '\\end{sidewaystable}\n');
    
    fprintf(FID, '}');
    fclose(FID);
end