function DM = DM_NW_test(model, T, sigma1, sigma2, H, II)
%     model = 'ar1';
%     sigma1 = 1;
%     sigma2 = 2;
%     T = 1000;
%     H = 100;
%      model = 'garch11';

    model = 'tgarch11';
    data_name = 'IBM_T2000_crisis';
%     H = 1000;
    
    if strcmp(data_name, 'IBM_T2000_crisis')
        name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'_with_ah.mat'];
    else
        name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'.mat'];
    end
    
    load(name, '-regexp','^C_score')
    aaa = who('-regexp','^C_score');
    met = length(aaa)/3;  % no of methods compared (div 3: 3 thresholds considered)

    load(name, '-regexp','^Cv_score')
    bbb = who('-regexp','^Cv_score');
    N = length(bbb)/3;
    
    aaa_05 = who('-regexp','C_score\w*_05$'); 
    aaa_1 = who('-regexp','C_score\w*_1$'); 
    aaa_5 = who('-regexp','C_score\w*_5$'); 

    bbb_05 = who('-regexp','Cv_score\w*_05$'); 
    bbb_1 = who('-regexp','Cv_score\w*_1$'); 
    bbb_5 = who('-regexp','Cv_score\w*_5$');
    
    met_2 = (met-1)/2;
    ind = zeros(met,1);
    ind(1) = met;
    ind(2:2:met) = 1:(met_2);
    ind(3:2:met) = (met_2+1):(met-1);
    aaa_05 = aaa_05(ind);
    aaa_1 = aaa_1(ind);
    aaa_5 = aaa_5(ind); 
   
    bbb_05 = bbb_05(ind);
    bbb_1 = bbb_1(ind);
    bbb_5 = bbb_5(ind); 
%     [7, 1,4,2,5,3,6]
%     2,4,6 --> 1,2,3
%     3,5,7 --> 2,4,6

    % Ordering: MC sampling from the true model, <-- for SS only
    % (Standard) Posterior,
    % Censored Posterior and Partially Censored Posterior with threshold at 0,
    % Censored Posterior and Partially Censored Posterior with threshold at the 10th sample percentile
    % Censored Posterior and Partially Censored Posterior with model free time verying threshold  
    % Censored Posterior and Partially Censored Posterior with MLE-implied time verying threshold  
    

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
    
    DM = NaN(met,met); % NaN will be removed while printing to tex
    DM_05 = (DM);
    DM_1 = (DM);
    DM_5 = (DM);
    % positive terms: the row method is better than the columm one

    for ii = 2:met
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

    DM_05(1,:) = [];
    DM_1(1,:) = [];
    DM_5(1,:) = [];    
    %% time varying     
    DM_v = NaN(met,met); % NaN will be removed while printing to tex
    DM_v_05 = DM_v;
    DM_v_1 = DM_v;
    DM_v_5 = DM_v;
    % positive terms: the row method is better than the column one

    for ii = 2:met
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


    %% Print to tex file
%     fname = print_table_DM(DM,model,T,S);
%     methods = {'True','Posterior','CP0','PCP0','CP10\%','PCP10\%'};

    if (met == 5)
%         methods =  {'CP','CP$_{var}$','PCP','PCP$_{var}$','Post'};
        methods =  {'Posterior','CP10\%','PCP10\%',...
            'CP$_{var,mle}$','PCP$_{var,mle}$'};
    elseif (met == 7)
        methods =  {'Posterior','CP10\\%%','PCP10\\%%',...
            'CP$_{var,mf}$','PCP$_{var,mf}$',...
            'CP$_{var,mle}$','PCP$_{var,mle}$'};
    end
    print_table_DM_emp(DM_5, data_name, 5, methods,'');
    print_table_DM_emp(DM_1, data_name, 1, methods,'');
    print_table_DM_emp(DM_05, data_name, 0.5, methods,'');

    print_table_DM_emp(DM_v_5, data_name, 5, methods,'v');
    print_table_DM_emp(DM_v_1, data_name, 1, methods,'v');
    print_table_DM_emp(DM_v_05, data_name, 0.5, methods,'v');

end
