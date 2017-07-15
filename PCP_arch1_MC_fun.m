function PCP_arch1_MC_fun(T, sigma2, S, II)
% T = 2500; S = 1; II = 10; sigma2 = 2 ; % <------------------ !!! 

    % clear all
    close all

    addpath(genpath('include/'));

%     s = RandStream('mt19937ar','Seed',1);
%     RandStream.setGlobalStream(s); 

    model = 'arch1';    
    partition = 3;
    if ~exist('II','var')
        II = 10;
    end
    fprintf('Model: %s.\n',model)
    parameters = {'$\\mu$','$\\omega$','$\\mu2$','$\\alpha$'};
    fprintf('Time series length T = %d.\n',T);
    
    sigma1 = 1;
    % sigma2 = 2;
    c = (sigma2 - sigma1)/sqrt(2*pi); % mean of eps
    kappa = 0.5*(sigma1^2 + sigma2^2 - ((sigma2-sigma1)^2)/pi); % var of eps
    sigma1_k = sigma1/sqrt(kappa);
    sigma2_k = sigma2/sqrt(kappa);
    
    mu2 = 0; % gama = mu - mu2;
    omega = 1;
    alpha = 0.1;

    % theta  = [mu, omega, mu2, alpha]
    mu_true = [0, omega, 0, alpha];
    param_true = [c,sigma2,omega,mu2,alpha];
    mu_init = [0, 1, 0.1, 0.05];

    % S = 20; % number of MC replications
    H = 100;

    varc = true; % run the version with time varying threshold

    % quantiles of interest
    p_bar1 = 0.01;
    p_bar = 0.05;
    % theoretical quantiles
    q1 = zeros(S,H);
    q5 = zeros(S,H);

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

    %% various display options
    cont.disp = false; %true; %false;

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
    mean_draw = zeros(S,4);
    mean_draw_C = zeros(S,4);
    mean_draw_PC = zeros(S,4);
    mean_draw_C0 = zeros(S,4);
    mean_draw_PC0 = zeros(S,4);

    median_draw = zeros(S,4);
    median_draw_C = zeros(S,4);
    median_draw_PC = zeros(S,4);
    median_draw_C0 = zeros(S,4);
    median_draw_PC0 = zeros(S,4);

    std_draw = zeros(S,4);
    std_draw_C = zeros(S,4);
    std_draw_PC = zeros(S,4);
    std_draw_C0 = zeros(S,4);
    std_draw_PC0 = zeros(S,4);

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
    
    SDD = zeros(S,1);

    %% MC Simulations
    tic
    % for s = 1:S
    s = 0;
    sdd = 0;
    while s < S    
        sdd = sdd + 1;
        %     if (mod(s,10)==0)
                fprintf(['\n',model, ' simulation no. %i\n'],s)
        %     end
        try
            if varc            
                results = PCP_arch1_run_varc(sdd, c, sigma1, sigma2, kappa, omega, alpha, p_bar1, p_bar, T, H, M, BurnIn, mu_init, df, cont, options, partition, II, GamMat);
            else
                results = PCP_arch1_run(sdd, c, sigma1, sigma2, kappa, omega, alpha, p_bar1, p_bar, T, H, M, BurnIn, mu_init, df, cont, options, partition, II, GamMat);
            end
            s = s+1;
            SDD(s,:) = sdd;
 
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
            if ~varc
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
            else
                mit_C0{s,1} = results.mit_Cm;
                CV_C0{s,1} = results.CV_Cm;        
                draw_C0 = results.draw_Cm;
                accept_C0(s,1) = results.accept_Cm;
                mean_draw_C0(s,:) = results.mean_draw_Cm;
                median_draw_C0(s,:) = results.median_draw_Cm;
                std_draw_C0(s,:) =  results.std_draw_Cm;
                VaR_1_post_C0(s,:) = results.VaR_1_post_Cm; 
                VaR_5_post_C0(s,:) = results.VaR_5_post_Cm; 

                draw_PC0 = results.draw_PCm;
                accept_PC0(s,1) = results.accept_PCm;
                mean_draw_PC0(s,:) = results.mean_draw_PCm;
                median_draw_PC0(s,:) = results.median_draw_PCm;
                std_draw_PC0(s,:) =  results.std_draw_PCm;
                VaR_1_post_PC0(s,:) = results.VaR_1_post_PCm; 
                VaR_5_post_PC0(s,:) = results.VaR_5_post_PCm;                 
            end
        end
    end

    % MSEs
    MSE_1 = mean((VaR_1 - q1).^2,2);
    MSE_1_post = mean((VaR_1_post - q1).^2,2);
    MSE_1_post_C = mean((VaR_1_post_C - q1).^2,2);
    MSE_1_post_PC = mean((VaR_1_post_PC - q1).^2,2);
%     if ~varc
        MSE_1_post_C0 = mean((VaR_1_post_C0 - q1).^2,2);
        MSE_1_post_PC0 = mean((VaR_1_post_PC0 - q1).^2,2);
%     end
    
    MSE_5 = mean((VaR_5 - q5).^2,2);
    MSE_5_post = mean((VaR_5_post - q5).^2,2);
    MSE_5_post_C = mean((VaR_5_post_C - q5).^2,2);
    MSE_5_post_PC = mean((VaR_5_post_PC - q5).^2,2);
%     if ~varc
        MSE_5_post_C0 = mean((VaR_5_post_C0 - q5).^2,2);
        MSE_5_post_PC0 = mean((VaR_5_post_PC0 - q5).^2,2);
%     end
    time_total = toc;

    if save_on
        if varc        
            name = ['results/',model,'/',model,'_',num2str(sigma1),'_',...
                num2str(sigma2),...
                '_T',num2str(T),'_H',num2str(H),'_II',num2str(II)...
                '_PCP0_MC_',v_new,'_varc.mat'];
        else
            name = ['results/',model,'/',model,'_',num2str(sigma1),'_',...
            num2str(sigma2),'_T',num2str(T),'_H',num2str(H),...
            '_II',num2str(II),'_PCP0_MC_',v_new,'.mat'];            
        end

%         if varc
%             save(name,...
%             'time_total','SDD',...
%             'y','draw','draw_C','draw_PC','param_true','q1','q5',...
%             'mean_draw','mean_draw_C','mean_draw_PC',...
%             'std_draw','std_draw_C','std_draw_PC',...
%             'accept','accept_C','accept_PC',...
%             'II','mit','CV','mit_C','CV_C',...
%             'VaR_1','VaR_1_post','VaR_1_post_C','VaR_1_post_PC',...
%             'VaR_5','VaR_5_post','VaR_5_post_C','VaR_5_post_PC',...
%             'MSE_1','MSE_1_post','MSE_1_post_C','MSE_1_post_PC',...
%             'MSE_5','MSE_5_post','MSE_5_post_C','MSE_5_post_PC')
%         else
            save(name,...
            'time_total','SDD',...
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
%         end
    end
end