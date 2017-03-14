function PCP_agarch11_MC_fun(T, sigma2, S, II)
% T = 1000; S = 1; II = 10; sigma2 = 2 ; % <------------------ !!! 

    % clear all
    close all

    addpath(genpath('include/'));

    s = RandStream('mt19937ar','Seed',1);
    RandStream.setGlobalStream(s); 

    model = 'agarch11';     partition = 3;
    fprintf('Model: %s.\n',model)
    parameters = {'$\\mu$','$\\gamma$','$\\omega$','$\\alpha$','$\\beta$'};

    sigma1 = 1;
    % sigma2 = 2;
    c = (sigma2 - sigma1)/sqrt(2*pi); % mean of eps
    kappa = 0.5*(sigma1^2 + sigma2^2 - ((sigma2-sigma1)^2)/pi); % var of eps
    sigma1_k = sigma1/sqrt(kappa);
    sigma2_k = sigma2/sqrt(kappa);
    
    gama = 0; % "typo" on purpose: not to confuse with the gamma function
    omega = 1;
    alpha = 0.1;
    beta = 0.8;
    % theta  = [mu, gama, omega, alpha, beta]
    mu_true = [0, 0, omega, alpha, beta];
    param_true = [c,sigma2,gama,omega,alpha,beta];
    mu_init = [0, -0.1, 1, 0.05, 0.85];
    % mu_init = [mean(y), -0.01, 0.01, 0.05, 0.85];

    % S = 20; % number of MC replications
    H = 100;

    % quantiles of interest
    p_bar1 = 0.01;
    p_bar = 0.05;
    % theoretical quantiles
    q1 = zeros(S,H);
    q5 = zeros(S,H);

    %% simulated quantiles: 
    % true model, posterior, censored posterior 10%, partially censored posterior 10%, 
    % censored posterior at 0, partially censored posterior at 0
    VaR_1 = zeros(S,H);
    VaR_1_post = zeros(S,H);
    VaR_1_post_C = zeros(S,H);
    VaR_1_post_PC = zeros(S,H);
    VaR_1_post_C0 = zeros(S,H);
    VaR_1_post_PC0 = zeros(S,H);

    VaR_5 = zeros(S,H);
    VaR_5_post = zeros(S,H);
    VaR_5_post_C = zeros(S,H);
    VaR_5_post_PC = zeros(S,H);
    VaR_5_post_C0 = zeros(S,H);
    VaR_5_post_PC0 = zeros(S,H);

    %% simulated parameters:
    % true model, posterior, censored posterior 10%, partially censored posterior 10%, 
    % censored posterior at 0, partially censored posterior at 0
    mean_draw = zeros(S,5);
    mean_draw_C = zeros(S,5);
    mean_draw_PC = zeros(S,5);
    mean_draw_C0 = zeros(S,5);
    mean_draw_PC0 = zeros(S,5);

    median_draw = zeros(S,5);
    median_draw_C = zeros(S,5);
    median_draw_PC = zeros(S,5);
    median_draw_C0 = zeros(S,5);
    median_draw_PC0 = zeros(S,5);

    std_draw = zeros(S,5);
    std_draw_C = zeros(S,5);
    std_draw_PC = zeros(S,5);
    std_draw_C0 = zeros(S,5);
    std_draw_PC0 = zeros(S,5);

    accept = zeros(S,1);
    accept_C = zeros(S,1);
    accept_PC = zeros(S,1);
    accept_C0 = zeros(S,1);
    accept_PC0 = zeros(S,1);

    %% MitISEM results: mits and CVs
    mit = cell(S,1);
    mit_C = cell(S,1);
    mit_C0 = cell(S,1);

    CV = cell(S,1);
    CV_C = cell(S,1);
    CV_C0 = cell(S,1);

    %%
    % T = 1000; % time series length
    fprintf('Time series length T = %d.\n',T)

    % Metropolis-Hastings for the parameters
    M = 10000; % number of draws 
    BurnIn = 1000;

    x_gam = (0:0.00001:50)'+0.00001;
    GamMat = gamma(x_gam);

    cont = MitISEM_Control;
    cont.mit.CV_max = 1; %2?
    cont.mit.iter_max = 10;
    cont.mit.Hmax = 6;
    cont.mit.dfnc = 5;
    df = 5; % default df for a mit

    if (nargin == 3)
        II = 10;
    end
    %% various display options
    cont.disp = true; %false;

    v_new = ver('symbolic');
    v_new = v_new.Release;
    if strcmp(v_new,'(R2014a)')
        fn_hist = @(xx) hist(xx,20);
    else
        fn_hist = @(xx) histogram(xx,20);
    end

    plot_on = false;
    save_on = true;

    options = optimset('Display','off');
    % w = warning('query','last');
    % id = w.identifier;
    id = 'optim:fminunc:SwitchingMethod';
    warning('off',id);

    %% MC Simulations
    tic
    % for s = 1:S
    s = 0;
    while s < S    
        %     if (mod(s,10)==0)
                fprintf(['\n',model, ' simulation no. %i\n'],s)
        %     end
        try
            results = PCP_agarch11_run(c, sigma1, sigma2, kappa, omega, alpha, beta, p_bar1, p_bar, T, H, M, BurnIn, mu_init, df, cont, options, partition, II, GamMat);
            s = s+1;

            y = results.y;
            q1(s,:) = results.q1;
            q5(s,:) = results.q5;   

            VaR_1(s,:) = results.VaR_1;
            VaR_5(s,:) = results.VaR_5;

            mit{s,1} = results.mit;
            CV{s,1} = results.CV;
            draw = results.draw;
            accept(s,1) = results.accept;
            mean_draw(s,:) = results.mean_draw;
            median_draw(s,:) = results.median_draw;
            std_draw(s,:) =  results.std_draw;
            VaR_1_post(s,:) = results.VaR_1_post; 
            VaR_5_post(s,:) = results.VaR_5_post;

            mit_C{s,1} = results.mit_C;
            CV_C{s,1} = results.CV_C;        
            draw_C = results.draw_C;
            accept_C(s,1) = results.accept_C;
            mean_draw_C(s,:) = results.mean_draw_C;
            median_draw_C(s,:) = results.median_draw_C;
            std_draw_C(s,:) =  results.std_draw_C;
            VaR_1_post_C(s,:) = results.VaR_1_post_C; 
            VaR_5_post_C(s,:) = results.VaR_5_post_C; 

            draw_PC = results.draw_PC;
            accept_PC(s,1) = results.accept_PC;
            mean_draw_PC(s,:) = results.mean_draw_PC;
            median_draw_PC(s,:) = results.median_draw_PC;
            std_draw_PC(s,:) =  results.std_draw_PC;
            VaR_1_post_PC(s,:) = results.VaR_1_post_PC; 
            VaR_5_post_PC(s,:) = results.VaR_5_post_PC; 

            mit_C0{s,1} = results.mit_C0;
            CV_C0{s,1} = results.CV_C0;        
            draw_C0 = results.draw_C0;
            accept_C0(s,1) = results.accept_C0;
            mean_draw_C0(s,:) = results.mean_draw_C0;
            median_draw_C0(s,:) = results.median_draw_C0;
            std_draw_C0(s,:) =  results.std_draw_C0;
            VaR_1_post_C0(s,:) = results.VaR_1_post_C0; 
            VaR_5_post_C0(s,:) = results.VaR_5_post_C0; 

            draw_PC0 = results.draw_PC0;
            accept_PC0(s,1) = results.accept_PC0;
            mean_draw_PC0(s,:) = results.mean_draw_PC0;
            median_draw_PC0(s,:) = results.median_draw_PC0;
            std_draw_PC0(s,:) =  results.std_draw_PC0;
            VaR_1_post_PC0(s,:) = results.VaR_1_post_PC0; 
            VaR_5_post_PC0(s,:) = results.VaR_5_post_PC0;                 
        end
    end

    % MSEs
    MSE_1 = mean((VaR_1 - q1).^2,2);
    MSE_1_post = mean((VaR_1_post - q1).^2,2);
    MSE_1_post_C = mean((VaR_1_post_C - q1).^2,2);
    MSE_1_post_PC = mean((VaR_1_post_PC - q1).^2,2);
    MSE_1_post_C0 = mean((VaR_1_post_C0 - q1).^2,2);
    MSE_1_post_PC0 = mean((VaR_1_post_PC0 - q1).^2,2);

    MSE_5 = mean((VaR_5 - q5).^2,2);
    MSE_5_post = mean((VaR_5_post - q5).^2,2);
    MSE_5_post_C = mean((VaR_5_post_C - q5).^2,2);
    MSE_5_post_PC = mean((VaR_5_post_PC - q5).^2,2);
    MSE_5_post_C0 = mean((VaR_5_post_C0 - q5).^2,2);
    MSE_5_post_PC0 = mean((VaR_5_post_PC0 - q5).^2,2);

    time_total = toc;

    if save_on
        name = ['results/',model,'/',model,'_',num2str(sigma1),'_',...
            num2str(sigma2),'_T',num2str(T),'_H',num2str(H),...
            '_II',num2str(II),'_PCP0_MC_',v_new,'.mat'];
        save(name,...
        'time_total',...
        'y','draw','draw_C','draw_PC','draw_C0','draw_PC0','param_true','q1','q5',...
        'mean_draw','mean_draw_C','mean_draw_PC','mean_draw_C0','mean_draw_PC0',...
        'median_draw','median_draw_C','median_draw_PC','median_draw_C0','median_draw_PC0',...
        'std_draw','std_draw_C','std_draw_PC','std_draw_C0','std_draw_PC0',...
        'accept','accept_C','accept_PC','accept_C0','accept_PC0',...
        'II','mit','CV','mit_C','CV_C','mit_C0','CV_C0',...
        'VaR_1','VaR_1_post','VaR_1_post_C','VaR_1_post_PC','VaR_1_post_C0','VaR_1_post_PC0',...
        'VaR_5','VaR_5_post','VaR_5_post_C','VaR_5_post_PC','VaR_5_post_C0','VaR_5_post_PC0',...
        'MSE_1','MSE_1_post','MSE_1_post_C','MSE_1_post_PC','MSE_1_post_C0','MSE_1_post_PC0',...
        'MSE_5','MSE_5_post','MSE_5_post_C','MSE_5_post_PC','MSE_5_post_C0','MSE_5_post_PC0')
    end
end