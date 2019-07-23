
for sigma2 = 1.5:0.5:3.5
    clearvars -except sigma2 
    close all

    addpath(genpath('include/'));

    % s = RandStream('mt19937ar','Seed',1);
    % RandStream.setGlobalStream(s); 

    model = 'iid';
    fprintf('**** Model: %s. *****\n',model)
    parameters = {'$\\mu$','$\\sigma$'};

    sigma1 = 1; 
    % sigma2 = 2;
    c = (sigma2 - sigma1)/(sqrt(2*pi));%1/sqrt(2*pi); %0.3989

    fprintf(' *Sigma2 : %4.2f.*\n',sigma2)

    S = 100; % number of MC replications

    VaR_1_post = zeros(S,1);
    VaR_1_post_C = zeros(S,1,4);
    VaR_1_post_C0 = zeros(S,1);

    VaR_5_post = zeros(S,1);
    VaR_5_post_C = zeros(S,1,4);
    VaR_5_post_C0 = zeros(S,1);

    VaR_05_post = zeros(S,1);
    VaR_05_post_C = zeros(S,1,4);
    VaR_05_post_C0 = zeros(S,1);

    accept = zeros(S,1);
    accept_C = zeros(S,4);
    accept_C0 = zeros(S,1);

    T = 1000; %time series length

    TAU = [0.05, 0.1, 0.15, 0.2];

    p_bar05 = 0.005;
    p_bar1 = 0.01;
    p_bar = 0.05;    

    % true VaRs
    q05 = norminv(p_bar05,c,sigma2);
    q1 = norminv(p_bar1,c,sigma2);
    q5 = norminv(p_bar,c,sigma2); 

    M = 10000; % number of draws 
    BurnIn = 1000;

    x_gam = (0:0.00001:50)'+0.00001;
    GamMat = gamma(x_gam);
 
    plot_on = false;
    save_on = true;

    options = optimset('Display','off');
    % w = warning('query','last');
    % id = w.identifier;
    id = 'optim:fminunc:SwitchingMethod';
    warning('off',id);

    for s = 1:S
        if (mod(s,10)==0)
            fprintf(['\n',model, ' simualtion no. %i\n'],s)
        end

        %% iid simulation, mean 0
        eps = randn(T,1);
        ind = (eps>0);
        eps(ind) = c + sigma1.*eps(ind);
        eps(~ind) = c + sigma2.*eps(~ind);
        % eps1 = c + sigma1.*eps(eps>0);
        % eps2 = c + sigma2.*eps(eps<0);
        % eps = [eps1;eps2];
        y = eps; % -mean(eps);
        y_true_sort = sort(y);
        Thresholds = y_true_sort(round(0.1*TAU*T));

        %% Uncensored Posterior
        % Misspecified model: normal with unknown mu and sigma
        % Metropolis-Hastings for the parameters

        % Uncensored likelihood
        kernel_init = @(xx) -loglik_iid(xx,y);
        [mu,~,~,~,~,Sigma] = fminunc(kernel_init,[0,1],options);
        Sigma = inv(T*Sigma);
        df = 5;
        draw = rmvt(mu,Sigma,df,M+BurnIn);
        kernel = @(ss) posterior_iid(ss,y);
        lnk = kernel(draw);

        lnd = dmvgt_mex(draw, mu, Sigma, df, 1, GamMat, double(1));
        lnw = lnk - lnd;
        lnw = lnw - max(lnw);
        [ind, a] = fn_MH(lnw);
        draw = draw(ind,:);
        accept(s) = a/(M+BurnIn);
        draw = draw(BurnIn+1:BurnIn+M,:);    

        y_post = draw(:,1) + draw(:,2).*randn(M,1);
        y_post = sort(y_post);
        VaR_05_post(s,1) = y_post(p_bar05*M); 
        VaR_1_post(s,1) = y_post(p_bar1*M); 
        VaR_5_post(s,1) = y_post(p_bar*M);

        %% Censored posterior: take values below the threshold
        % Misspecified model: N(mu,sigma)
        Draws_C = zeros(M,2,4);

        for tau = 1:4
            fprintf('tau = %4.2f\n',TAU(tau));
            threshold = Thresholds(tau);
            % 1. Threshold = 10% perscentile of the data sample
            kernel_init = @(xx) - C_posterior_iid(xx, y, threshold)/T;
            [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu,options);
            Sigma_C = inv(T*Sigma_C);
            draw_C = rmvt(mu_C,Sigma_C,df,M+BurnIn);
            kernel = @(ss) C_posterior_iid(ss, y, threshold);
            lnk_C = kernel(draw_C);
            lnd_C = dmvgt_mex(draw_C, mu_C, Sigma_C, df, 1, GamMat, double(1));
            lnw_C = lnk_C - lnd_C;
            lnw_C = lnw_C - max(lnw_C);
            [ind, a] = fn_MH(lnw_C);
            draw_C = draw_C(ind,:);
            accept_C(s,tau) = a/(M+BurnIn);
            draw_C = draw_C(BurnIn+1:BurnIn+M,:);

            Draws_C(:,:,tau) = draw_C;

            y_post_C = draw_C(:,1) + draw_C(:,2).*randn(M,1);
            y_post_C = sort(y_post_C);
            VaR_05_post_C(s,1,tau) = y_post_C(p_bar05*M); 
            VaR_1_post_C(s,1,tau) = y_post_C(p_bar1*M); 
            VaR_5_post_C(s,1,tau) = y_post_C(p_bar*M); 
        end

        % 2. Threshold = 0
        threshold0 = 0;
        kernel_init = @(xx) - C_posterior_iid(xx, y, threshold0)/T;
        [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu,options);
        Sigma_C0 = inv(T*Sigma_C0);
        draw_C0 = rmvt(mu_C0,Sigma_C0,df,M+BurnIn);
        kernel = @(ss) C_posterior_iid(ss, y, threshold0);
        lnk_C0 = kernel(draw_C0);
        lnd_C0 = dmvgt_mex(draw_C0, mu_C0, Sigma_C0, df, 1, GamMat, double(1));
        lnw_C0 = lnk_C0 - lnd_C0;
        lnw_C0 = lnw_C0 - max(lnw_C0);
        [ind, a] = fn_MH(lnw_C0);
        draw_C0 = draw_C0(ind,:);
        accept_C0(s) = a/(M+BurnIn);
        draw_C0 = draw_C0(BurnIn+1:BurnIn+M,:);

        y_post_C0 = draw_C0(:,1) + draw_C0(:,2).*randn(M,1);
        y_post_C0 = sort(y_post_C0);
        VaR_05_post_C0(s,1) = y_post_C0(p_bar05*M); 
        VaR_1_post_C0(s,1) = y_post_C0(p_bar1*M); 
        VaR_5_post_C0(s,1) = y_post_C0(p_bar*M); 
    end

    MSE_1_post = sum((VaR_1_post - q1).^2)/S;
    MSE_1_post_C = zeros(4,1);
    for tau = 1:4
        MSE_1_post_C(tau) = sum((squeeze(VaR_1_post_C(:,:,tau)) - q1).^2)/S;
    end
    MSE_1_post_C0 = sum((VaR_1_post_C0 - q1).^2)/S;

    MSE_5_post = sum((VaR_5_post - q5).^2)/S;
    MSE_5_post_C = zeros(4,1);
    for tau = 1:4
        MSE_5_post_C(tau) = sum((squeeze(VaR_5_post_C(:,:,tau)) - q5).^2)/S;
    end
    MSE_5_post_C0 = sum((VaR_5_post_C0 - q5).^2)/S;

    MSE_05_post = sum((VaR_05_post - q05).^2)/S;
    MSE_05_post_C = zeros(4,1);
    for tau = 1:4
        MSE_05_post_C(tau) = sum((squeeze(VaR_05_post_C(:,:,tau)) - q05).^2)/S;
    end
    MSE_05_post_C0 = sum((VaR_05_post_C0 - q05).^2)/S;


    param_true = [c,sigma2];
    if save_on
        save(['results/',model,'/',model,'_',num2str(sigma1),'_',...
            num2str(sigma2),'_T',num2str(T),'_MC_diff_thresholds.mat'],...
            'y','param_true','q1','q5','q05',...,
            'draw','Draws_C','draw_C0',...
            'accept','accept_C','accept_C0',...
            'VaR_1_post','VaR_1_post_C','VaR_1_post_C0',...
            'VaR_5_post','VaR_5_post_C','VaR_5_post_C0',...
            'VaR_05_post','VaR_05_post_C','VaR_05_post_C0',...
            'MSE_1_post','MSE_1_post_C','MSE_1_post_C0',...
            'MSE_5_post','MSE_5_post_C','MSE_5_post_C0',...
            'MSE_05_post','MSE_05_post_C','MSE_05_post_C0')

    end
end
%print_table_cp_mc(model,parameters,sigma1,sigma2)