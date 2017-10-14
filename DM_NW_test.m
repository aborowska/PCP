function DM = DM_NW_test(model, T, sigma1, sigma2, H, II)
%     model = 'ar1';
%     sigma1 = 1;
%     sigma2 = 2;
%     T = 1000;
%     H = 100;
%      model = 'garch11';

    model = 'tgarch11';
    data_name = 'MSFT';
    H = 2000;

    name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'.mat'];

    load(name, '-regexp','^C_score')
    aaa = who('-regexp','^C_score');
    M = length(aaa)/3;  % no of methods compared (div 3: 3 thresholds considered)

    load(name, '-regexp','^Cv_score')
    bbb = who('-regexp','^Cv_score');
    N = length(bbb)/3;
    
    aaa_05 = who('-regexp','C_score\w*_05$'); 
    aaa_1 = who('-regexp','C_score\w*_1$'); 
    aaa_5 = who('-regexp','C_score\w*_5$'); 

    bbb_05 = who('-regexp','Cv_score\w*_05$'); 
    bbb_1 = who('-regexp','Cv_score\w*_1$'); 
    bbb_5 = who('-regexp','Cv_score\w*_5$');    
    % Ordering: MC sampling from the true model,
    % (Standard) Posterior,
    % Censored Posterior and Partially Censored Posterior with threshold at 0,
    % Censored Posterior and Partially Censored Posterior with threshold at the 10th sample percentile


    

    %% DM test statistics based on Censored Score Rules differentials
% d = censored scoring rule of model 1 at time t - censored scoring rule of model 2 at time t.
% We test whether the sample mean of the d_t significantly differs from 0, 
% which is equivalent with a regression of d_t on a constant term:
% d_t  = mu + error_t  (*)
% where the OLS estimator of mu is obviously equal to the sample mean of the d_t. 
% Because the d_t may have serial correlation, it is useful to consider this in the form of regression model (*), 
% where the Newey-West standard error of the OLS estimate of mu can be used to get a standard error 
% (at the sample mean of the d_t) that takes into account the serial correlation in the d_t.
% So, the t-statistic (that for a large enough number of observations can be compared with critical values -1.96 and 1.96) is then given by
% t-statistic = (sample mean of the d_t)/(Newey-West standard error) = (OLS estimator of mu)/(Newey-West standard error).
    
    DM = NaN(M,1); % NaN will be removed while printing to tex
    DM_05 = diag(DM);
    DM_1 = diag(DM);
    DM_5 = diag(DM);
    % positive terms: the colum method is better than the row one
    % negative terms: the row method is better than the column one

    for ii = 2:M
        for jj = 1:(ii-1)
            d = eval(char(aaa_05{ii})) - eval(char(aaa_05{jj}));
            if (sum(imag(d)~=0) > 0)
                DM_05(ii,jj) = NaN;
            else
                nwse = sqrt(NeweyWest(d));
                DM_05(ii,jj) = mean(d)/nwse;
            end
            
            d = eval(char(aaa_1{ii})) - eval(char(aaa_1{jj}));
            if (sum(imag(d)~=0) > 0)
                DM_1(ii,jj) = NaN;
            else
                nwse = sqrt(NeweyWest(d));
                DM_1(ii,jj) = mean(d)/nwse;
            end
            
            d = eval(char(aaa_5{ii})) - eval(char(aaa_5{jj}));
            if (sum(imag(d)~=0) > 0)
                DM_5(ii,jj) = NaN;
            else
                nwse = sqrt(NeweyWest(d));
                DM_5(ii,jj) = mean(d)/nwse;
            end
        end
    end

    %% time varying     
    DM_v = NaN(N,1); % NaN will be removed while printing to tex
    DM_v_05 = diag(DM_v);
    DM_v_1 = diag(DM_v);
    DM_v_5 = diag(DM_v);
    % positive terms: the colum method is better than the row one
    % negative terms: the row method is better than the column one

    for ii = 2:N
        for jj = 1:(ii-1)
            d = eval(char(bbb_05{ii})) - eval(char(bbb_05{jj}));
            if (sum(imag(d)~=0) > 0)
                DM_v_05(ii,jj) = NaN;
            else
                nwse = sqrt(NeweyWest(d));
                DM_v_05(ii,jj) = mean(d)/nwse;
            end
            
            d = eval(char(bbb_1{ii})) - eval(char(bbb_1{jj}));
            if (sum(imag(d)~=0) > 0)
                DM_v_1(ii,jj) = NaN;
            else
                nwse = sqrt(NeweyWest(d));
                DM_v_1(ii,jj) = mean(d)/nwse;
            end
            
            d = eval(char(bbb_5{ii})) - eval(char(bbb_5{jj}));
            if (sum(imag(d)~=0) > 0)
                DM_v_5(ii,jj) = NaN;
            else
                nwse = sqrt(NeweyWest(d));
                DM_v_5(ii,jj) = mean(d)/nwse;
            end
        end
    end

    %% Combined
    if (N == M)
        DM_all_05 = DM_05; % NaN will be removed while printing to tex
        DM_all_1 = DM_1; % NaN will be removed while printing to tex
        DM_all_5 = DM_5; %
        
        for ii = 2:N
            for jj = 1:(ii-1)
                DM_all_05()
            end
        end
    end
    
    %% Print to tex file
    fname = print_table_DM(DM,model,T,S);
    Remove_NaN(fname);
end
