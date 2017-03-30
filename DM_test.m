function DM = DM_test(model, T, sigma1, sigma2, H, II)
%     model = 'ar1';
%     sigma1 = 1;
%     sigma2 = 2;
%     T = 1000;
%     H = 100;
%      model = 'garch11';

    model = 'agarch11';
    sigma1 = 1;
    sigma2 = 2;
    T = 1000;
    H = 100;

    v_new = '(R2015a)';
    II = 10;
    if strcmp(model,'ar1')
        name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_II',num2str(II),'_PCP0_MC2_',v_new,'.mat'];
    elseif strcmp(model,'garch11')
        name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_II',num2str(II),'_PCP0_MC_',v_new,'.mat'];
    else
        name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_II',num2str(II),'_PCP0_MC2_',v_new,'.mat'];    
    end
    load(name, '-regexp','^q\w*')
    load(name, '-regexp','^VaR\w*')
    % Ordering: MC sampling from the true model,
    % (Standard) Posterior,
    % Censored Posterior and Partially Censored Posterior with threshold at 0,
    % Censored Posterior and Partially Censored Posterior with threshold at the 10th sample percentile

    S = size(q1,1);

    %% Squared Forecast Errors
    % SE_1 = zeros(S,H,6);
    % SE_1(:,:,1) = (VaR_1 - q1).^2;
    % SE_1(:,:,2) = (VaR_1_post - q1).^2;
    % SE_1(:,:,3) = (VaR_1_post_C0 - q1).^2;
    % SE_1(:,:,4) = (VaR_1_post_PC0 - q1).^2;
    % SE_1(:,:,5) = (VaR_1_post_C - q1).^2;
    % SE_1(:,:,6) = (VaR_1_post_PC - q1).^2;
    % 
    % SE_5 = zeros(S,H,6);
    % SE_5(:,:,1) = (VaR_5 - q5).^2;
    % SE_5(:,:,2) = (VaR_5_post - q5).^2;
    % SE_5(:,:,3) = (VaR_5_post_C0 - q5).^2;
    % SE_5(:,:,4) = (VaR_5_post_PC0 - q5).^2;
    % SE_5(:,:,5) = (VaR_5_post_C - q5).^2;
    % SE_5(:,:,6) = (VaR_5_post_PC - q5).^2;

    %% Mean Squared Errors 
    MSE_1 = zeros(S,6);
    MSE_1(:,1) = (mean((VaR_1 - q1).^2,2));
    MSE_1(:,2) = (mean((VaR_1_post - q1).^2,2));
    MSE_1(:,3) = (mean((VaR_1_post_C0 - q1).^2,2));
    MSE_1(:,4) = (mean((VaR_1_post_PC0 - q1).^2,2));
    MSE_1(:,5) = (mean((VaR_1_post_C - q1).^2,2));
    MSE_1(:,6) = (mean((VaR_1_post_PC - q1).^2,2));

    MSE_5 = zeros(S,6);
    MSE_5(:,1) = (mean((VaR_5 - q5).^2,2));
    MSE_5(:,2) = (mean((VaR_5_post - q5).^2,2));
    MSE_5(:,3) = (mean((VaR_5_post_C0 - q5).^2,2));
    MSE_5(:,4) = (mean((VaR_5_post_PC0 - q5).^2,2));
    MSE_5(:,5) = (mean((VaR_5_post_C - q5).^2,2));
    MSE_5(:,6) = (mean((VaR_5_post_PC - q5).^2,2));

%     % remove "bad" simulations
%     ind_imag_1 = (imag(sum(RMSE_1,2)) == 0);
%     ind_imag_5 = (imag(sum(RMSE_5,2)) == 0);
%     ind_imag = (ind_imag_1 & ind_imag_5);
%     RMSE_1 = RMSE_1(ind_imag,:);
%     RMSE_5 = RMSE_5(ind_imag,:);
%     S = sum(ind_imag);

    %% DM test statistics based on loss differentials
    DM = NaN(6,1); % NaN will be removed while printing to tex
    DM = diag(DM);
    % below the diagonal: for 99% VaR
    % above the diagonal: for 95% VaR
    % positive terms: the colum method is better than the row one
    % negative terms: the row method is better than the column one

    for ii = 2:6
        for jj = 1:(ii-1)
            % loss differnetial vector (loss = MSE over H out-of-sample periods)
            % elements of d are iid hence no Newey-West type correction needed
            d = MSE_1(:,ii) - MSE_1(:,jj);
            DM(ii,jj) = mean(d)/(std(d)/sqrt(S));

            d = MSE_5(:,jj) - MSE_5(:,ii);
            DM(jj,ii) = mean(d)/(std(d)/sqrt(S));    
        end
    end

    %% Print to tex file
    fname = print_table_DM(DM,model,T,S);
    Remove_NaN(fname);
end
