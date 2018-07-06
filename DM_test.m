function DM = DM_test(model, T, sigma1, sigma2, H, II)
%     addpath(genpath('include/'));
% clear all

%     model = 'ar1';
%     sigma1 = 1;
%     sigma2 = 2;
%     T = 10000;
%     H = 100;
%      model = 'garch11';

    model = 'ar1'; %agarch11';
    sigma1 = 1;
    sigma2 = 2;
    T = 1000;
    H = 100;
    varc = 0;
    
%     v_new = '(R2015a)';
    v_new = '(R2017a)';
    II = 10;
    if strcmp(model,'ar1')
        name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_II',num2str(II),'_PCP0_MC_',v_new,'.mat'];
    elseif strcmp(model,'garch11')
        name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_II',num2str(II),'_PCP0_MC_',v_new,'.mat'];
    else
        name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_II',num2str(II),'_PCP0_MC2_',v_new,'.mat'];    
    end
%     load(name, '-regexp','^q\w*')
%     load(name, '-regexp','^VaR\w*')
    load(name, '-regexp','^ES\w*')
    load(name, '-regexp','^MSE\w*')
    load(name, '-regexp','^cdf\w*')
    MSE_ES;
    
    % Ordering: MC sampling from the true model,
    % (Standard) Posterior,
    % Censored Posterior and Partially Censored Posterior with threshold at 0,
    % Censored Posterior and Partially Censored Posterior with threshold at the 10th sample percentile

%     S = size(MSE_1,1);
    aaa = who('-regexp','^MSE_es');
    met = length(aaa)/3;  % no of methods compared (div 3: 3 thresholds considered)

    bbb = who('-regexp','^MSE_[^es]');
    
    aaa_05 = who('-regexp','MSE_es_05\w*$'); 
    aaa_1 = who('-regexp','MSE_es_1\w*$'); 
    aaa_5 = who('-regexp','MSE_es_5\w*$'); 

    bbb_05 = who('-regexp','MSE_05\w*$'); 
    bbb_1 = who('-regexp','MSE_1\w*$'); 
    bbb_5 = who('-regexp','MSE_5\w*$');
  
    met_2 = (met-2)/2;
    ind = zeros(met,1);
    ind(1:2) = [1,2];
    ind(3:2:met) = 3:(met_2+2);
    ind(4:2:met) = (met_2+3):(met);
    
    aaa_05 = aaa_05(ind);
    aaa_1 = aaa_1(ind);
    aaa_5 = aaa_5(ind);    
    
    bbb_05 = bbb_05(ind);
    bbb_1 = bbb_1(ind);
    bbb_5 = bbb_5(ind);    
    

    %% DM test statistics based on loss differentials
    % FOR VAR
    %% DM test statistics based on loss differentials
    DM = NaN(met,met); % NaN will be removed while printing to tex
    DM_05 = (DM);
    DM_1 = (DM);
    DM_5 = (DM);
     
    % positive terms: the colum method is better than the row one
    % negative terms: the row method is better than the column one

    for ii = 2:met
        for jj = 1:(ii-1)
            % loss differnetial vector (loss = MSE over H out-of-sample periods)
            % elements of d are iid hence no Newey-West type correction needed
%             d = MSE_1(:,ii) - MSE_1(:,jj);
            d1 = eval(char(bbb_05{ii}));
            d2 = eval(char(bbb_05{jj}));
            d = sqrt(d1) - sqrt(d2);
            se = std(d)/sqrt(length(d));
            DM_05(ii,jj) = mean(d)/se;
 
%             d = MSE_5(:,jj) - MSE_5(:,ii);
%             DM(jj,ii) = mean(d)/(std(d)/sqrt(S));    
            d1 = eval(char(bbb_1{ii}));
            d2 = eval(char(bbb_1{jj}));
            d = sqrt(d1) - sqrt(d2);
            se = std(d)/sqrt(length(d));
            DM_1(ii,jj) = mean(d)/se;
            
            
            d1 = eval(char(bbb_5{ii}));
            d2 = eval(char(bbb_5{jj}));
            d = sqrt(d1) - sqrt(d2);
            se = std(d)/sqrt(length(d));
            DM_5(ii,jj) = mean(d)/se;
        end
    end    
    
    % FOR ES
    
    
    DM = NaN(met,met); % NaN will be removed while printing to tex
    DM_es_05 = (DM);
    DM_es_1 = (DM);
    DM_es_5 = (DM);
    
    % below the diagonal: for 99% VaR
    % above the diagonal: for 95% VaR
    % positive terms: the colum method is better than the row one
    % negative terms: the row method is better than the column one

    for ii = 2:met
        for jj = 1:(ii-1)
            % loss differnetial vector (loss = MSE over H out-of-sample periods)
            % elements of d are iid hence no Newey-West type correction needed
%             d = MSE_1(:,ii) - MSE_1(:,jj);
            d1 = eval(char(aaa_05{ii}));
            d2 = eval(char(aaa_05{jj}));
            d = sqrt(d1) - sqrt(d2);
            se = std(d)/sqrt(length(d));
            DM_es_05(ii,jj) = mean(d)/se;
 
%             d = MSE_5(:,jj) - MSE_5(:,ii);
%             DM(jj,ii) = mean(d)/(std(d)/sqrt(S));    
            d1 = eval(char(aaa_1{ii}));
            d2 = eval(char(aaa_1{jj}));
            d = sqrt(d1) - sqrt(d2);
            se = std(d)/sqrt(length(d));
            DM_es_1(ii,jj) = mean(d)/se;
            
            
            d1 = eval(char(aaa_5{ii}));
            d2 = eval(char(aaa_5{jj}));
            d = sqrt(d1) - sqrt(d2);
            se = std(d)/sqrt(length(d));
            DM_es_5(ii,jj) = mean(d)/se;
        end
    end
     
%     DM_05(1,:) = [];
%     DM_1(1,:) = [];
%     DM_5(1,:) = [];   
%     DM_es_05(1,:) = [];
%     DM_es_1(1,:) = [];
%     DM_es_5(1,:) = [];  
     
    
    if (met == 6)
%         methods =  {'CP','CP$_{var}$','PCP','PCP$_{var}$','Post'};
        methods =  {'True','Posterior','CP0','PCP0','CP10\%','PCP10\%'};
    end
%     print_table_DM(DM_5, model, 5, methods,'');
%     print_table_DM(DM_1, model, 1, methods,'');
%     print_table_DM(DM_05, model, 0.5, methods,'');
% 
%     
%     print_table_DM(DM_es_5, model, 5, methods,'es');
%     print_table_DM(DM_es_1, model, 1, methods,'es');
%     print_table_DM(DM_es_05, model, 0.5, methods,'es');
    
    DM = zeros(met, met, 3);
    DM(:,:,1) = DM_05;
    DM(:,:,2) = DM_1;
    DM(:,:,3) = DM_5;
    print_table_DM_comb(DM, T, model, [0.5, 1, 5], methods,'');

    DMes = zeros(met, met, 3);
    DMes(:,:,1) = DM_es_05;
    DMes(:,:,2) = DM_es_1;
    DMes(:,:,3) = DM_es_5;    
    print_table_DM_comb(DMes, T, model, [0.5, 1, 5], methods,'es');

    
    % Print to tex file
%     fname = print_table_DM(DM,model,T,S);
%     Remove_NaN(fname);
end
