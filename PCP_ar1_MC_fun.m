function PCP_ar1_MC_fun(T, sigma2, II)
    close all

    addpath(genpath('include/'));

    s = RandStream('mt19937ar','Seed',1);
    RandStream.setGlobalStream(s); 

    model = 'ar1';
    fprintf('Model: %s.\n',model)
    parameters = {'$\\mu$','$\\sigma$','$\\phi$'};

    sigma1 = 1;
    % sigma2 = 2;    
    c = (sigma2 - sigma1)/sqrt(2*pi);
    rho = 0.8;
    param_true = [c, sigma2, rho];
    mu_init = [0,1,0.9];

    S = 100; % number of MC replications
    H = 100;
    
    threshold_c = true; % true = run the version with an additional theoretical threshold equal to c
    if (sigma1 == sigma2)
        threshold_c = false;
    end
    varc = false; % run the version with time varying threshold
 
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
    mean_draw = zeros(S,3);
    mean_draw_C = zeros(S,3);
    mean_draw_PC = zeros(S,3);
    mean_draw_C0 = zeros(S,3);
    mean_draw_PC0 = zeros(S,3);

    std_draw = zeros(S,3);
    std_draw_C = zeros(S,3);
    std_draw_PC = zeros(S,3);
    std_draw_C0 = zeros(S,3);
    std_draw_PC0 = zeros(S,3);

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

    if threshold_c 
        mit_Cc = cell(S,1);
        CV_Cc = cell(S,1);      
        accept_Cc = zeros(S,1);
        mean_draw_Cc = zeros(S,3);
        std_draw_Cc = zeros(S,3);
        VaR_1_post_Cc = zeros(S,H); 
        VaR_5_post_Cc = zeros(S,H);

        accept_PCc = zeros(S,1);
        mean_draw_PCc = zeros(S,3);
        std_draw_PCc =  zeros(S,3);
        VaR_1_post_PCc = zeros(S,H); 
        VaR_5_post_PCc = zeros(S,H); 
    end    
    
    %%
    % T = 10000; % time series length
    fprintf('Time series length T = %d.\n',T)

    % Metropolis-Hastings for the parameters
    M = 10000; % number of draws 
    BurnIn = 1000;

    x_gam = (0:0.00001:50)'+0.00001;
    GamMat = gamma(x_gam);

    df = 5; % default df for a mit
    cont = MitISEM_Control;
    cont.mit.iter_max = 10;
    cont.mit.Hmax = 6;
    cont.mit.dfnc = 5;
    
    if (nargin == 2)
        II = 10;
    end
    
    partition = 3;

    %% various display options
    cont.disp = false;

    plot_on = false;
    save_on = true;

    v_new = ver('symbolic');
    v_new = v_new.Release;
    if strcmp(v_new,'(R2014a)')
        fn_hist = @(xx) hist(xx,20);
    else
        fn_hist = @(xx) histogram(xx,50);
    end

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
            if threshold_c
                results = PCP_ar1_run_c(c, sigma1, sigma2, rho, p_bar1, p_bar, T, H, M, BurnIn, mu_init, df, cont, options, partition, II, GamMat);
            elseif varc
                results = PCP_ar1_run_varc(c, sigma1, sigma2, rho, p_bar1, p_bar, T, H, M, BurnIn, mu_init, df, cont, options, partition, II, GamMat);
            else
                results = PCP_ar1_run(c, sigma1, sigma2, rho, p_bar1, p_bar, T, H, M, BurnIn, mu_init, df, cont, options, partition, II, GamMat);
            end
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
            std_draw(s,:) =  results.std_draw;
            VaR_1_post(s,:) = results.VaR_1_post; 
            VaR_5_post(s,:) = results.VaR_5_post;

            mit_C{s,1} = results.mit_C;
            CV_C{s,1}  = results.CV_C;        
            draw_C = results.draw_C;
            accept_C(s,1) = results.accept_C;
            mean_draw_C(s,:) = results.mean_draw_C;
            std_draw_C(s,:) =  results.std_draw_C;
            VaR_1_post_C(s,:) = results.VaR_1_post_C; 
            VaR_5_post_C(s,:) = results.VaR_5_post_C; 

            draw_PC = results.draw_PC;
            accept_PC(s,1) = results.accept_PC;
            mean_draw_PC(s,:) = results.mean_draw_PC;
            std_draw_PC(s,:) =  results.std_draw_PC;
            VaR_1_post_PC(s,:) = results.VaR_1_post_PC; 
            VaR_5_post_PC(s,:) = results.VaR_5_post_PC; 

            mit_C0{s,1} = results.mit_C0;  
            CV_C0{s,1} = results.CV_C0;
            draw_C0 = results.draw_C0;
            accept_C0(s,1) = results.accept_C0;
            mean_draw_C0(s,:) = results.mean_draw_C0;
            std_draw_C0(s,:) =  results.std_draw_C0;
            VaR_1_post_C0(s,:) = results.VaR_1_post_C0; 
            VaR_5_post_C0(s,:) = results.VaR_5_post_C0; 

            draw_PC0 = results.draw_PC0;
            accept_PC0(s,1) = results.accept_PC0;
            mean_draw_PC0(s,:) = results.mean_draw_PC0;
            std_draw_PC0(s,:) =  results.std_draw_PC0;
            VaR_1_post_PC0(s,:) = results.VaR_1_post_PC0; 
            VaR_5_post_PC0(s,:) = results.VaR_5_post_PC0; 

            if threshold_c
                mit_Cc{s,1} = results.mit_Cc;
                CV_Cc{s,1} = results.CV_Cc;        
                draw_Cc = results.draw_Cc;
                accept_Cc(s,1) = results.accept_Cc;
                mean_draw_Cc(s,:) = results.mean_draw_Cc;
                std_draw_Cc(s,:) =  results.std_draw_Cc;
                VaR_1_post_Cc(s,:) = results.VaR_1_post_Cc; 
                VaR_5_post_Cc(s,:) = results.VaR_5_post_Cc; 

                draw_PCc = results.draw_PCc;
                accept_PCc(s,1) = results.accept_PCc;
                mean_draw_PCc(s,:) = results.mean_draw_PCc;
                std_draw_PCc(s,:) =  results.std_draw_PCc;
                VaR_1_post_PCc(s,:) = results.VaR_1_post_PCc; 
                VaR_5_post_PCc(s,:) = results.VaR_5_post_PCc;     
            end
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

    
    if threshold_c
        MSE_1_post_Cc = mean((VaR_1_post_Cc - q1).^2,2);
        MSE_1_post_PCc = mean((VaR_1_post_PCc - q1).^2,2);
        MSE_5_post_Cc = mean((VaR_5_post_Cc - q5).^2,2);
        MSE_5_post_PCc = mean((VaR_5_post_PCc - q5).^2,2);
    end
    
    time_total = toc;

    if save_on
        if threshold_c
            name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),...
                '_T',num2str(T),'_H',num2str(H),'_II',num2str(II)...
                '_PCP0_MC2_',v_new,'_c.mat'];
        elseif varc        
            name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),...
                '_T',num2str(T),'_H',num2str(H),'_II',num2str(II)...
                '_PCP0_MC2_',v_new,'_varc.mat'];
        else
            name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),...
                '_T',num2str(T),'_H',num2str(H),'_II',num2str(II)...
                '_PCP0_MC2_',v_new,'.mat'];
        end
        save(name,...
        'time_total',...
        'y','draw','draw_C','draw_PC','draw_C0','draw_PC0','param_true','q1','q5',...
        'mean_draw','mean_draw_C','mean_draw_PC','mean_draw_C0','mean_draw_PC0',...
        'std_draw','std_draw_C','std_draw_PC','std_draw_C0','std_draw_PC0',...
        'accept','accept_C','accept_PC','accept_C0','accept_PC0',...
        'II','mit','CV','mit_C','CV_C','mit_C0','CV_C0',...
        'VaR_1','VaR_1_post','VaR_1_post_C','VaR_1_post_PC','VaR_1_post_C0','VaR_1_post_PC0',...
        'VaR_5','VaR_5_post','VaR_5_post_C','VaR_5_post_PC','VaR_5_post_C0','VaR_5_post_PC0',...
        'MSE_1','MSE_1_post','MSE_1_post_C','MSE_1_post_PC','MSE_1_post_C0','MSE_1_post_PC0',...
        'MSE_5','MSE_5_post','MSE_5_post_C','MSE_5_post_PC','MSE_5_post_C0','MSE_5_post_PC0')
        if threshold_c
             save(name,...   
             'draw_Cc','draw_PCc','mean_draw_Cc','mean_draw_PCc','std_draw_Cc','std_draw_PCc',...
             'accept_Cc','accept_PCc','mit_Cc','CV_Cc',...
             'VaR_1_post_Cc','VaR_1_post_PCc','VaR_5_post_Cc','VaR_5_post_PCc',...
             'MSE_1_post_Cc','MSE_1_post_PC0',...
             'MSE_5_post_Cc','MSE_5_post_PCc','-append')
        end
    end

    % print_table_pcp_mc(model,parameters,sigma1,sigma2)
    % print_table_pcp_mc(model,parameters,sigma1,sigma2,100)

    if plot_on
        figure(21)
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
        hold on
        if strcmp(v_new,'(R2014a)')
            fn_hist([mean_draw(:,1), mean_draw_C(:,1), mean_draw_PC(:,1), mean_draw_C0(:,1), mean_draw_PC0(:,1)])    
        else
            fn_hist(mean_draw(:,1))
            fn_hist(mean_draw_C(:,1))
            fn_hist(mean_draw_PC(:,1))
            fn_hist(mean_draw_C0(:,1))
            fn_hist(mean_draw_PC0(:,1))
        end
        hold off
        xlabel('\mu','FontSize',11)
        plotTickLatex2D('FontSize',11); 
        leg = legend('Uncensored','CP 10\% ','PCP 10\% ','CP 0','PCP 0');
        set(leg,'Interpreter','latex','FontSize',11)

        figure(22)
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);    
        hold on
        if strcmp(v_new,'(R2014a)')
            fn_hist([mean_draw(:,2), mean_draw_C(:,2), mean_draw_PC(:,2), mean_draw_C0(:,2), mean_draw_PC0(:,2)])    
        else
            fn_hist(mean_draw(:,2))
            fn_hist(mean_draw_C(:,2))
            fn_hist(mean_draw_PC(:,2))
            fn_hist(mean_draw_C0(:,2))
            fn_hist(mean_draw_PC0(:,2))
        end
        hold off
        xlabel('\sigma','FontSize',11)
        plotTickLatex2D('FontSize',11); 
        leg = legend('Uncensored','CP 10\% ','PCP 10\% ','CP 0','PCP 0');
        set(leg,'Interpreter','latex','FontSize',11)

        figure(23)
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);    
        hold on    
        if strcmp(v_new,'(R2014a)')
            hist([mean_draw(:,3), mean_draw_C(:,3), mean_draw_PC(:,3), mean_draw_C0(:,3), mean_draw_PC0(:,3)],50)    
        else
    %         [pF, x] = ksdensity(mean_draw(:,3),'function','pdf');
    %         plot(x, pF/sum(pF))        
    %         [pF, x] = ksdensity(mean_draw_C(:,3),'function','pdf');
    %         plot(x, pF/sum(pF))        
    %         [pF, x] = ksdensity(mean_draw_PC(:,3),x,'function','pdf');
    %         plot(x, pF/sum(pF))        
    %         [pF, x] = ksdensity(mean_draw_C0(:,3),x,'function','pdf');
    %         plot(x, pF/sum(pF))        
    %         [pF, x] = ksdensity(mean_draw_PC0(:,3),x,'function','pdf');
    %         plot(x, pF/sum(pF))        
            fn_hist(mean_draw(:,3))
            fn_hist(mean_draw_C(:,3))
            fn_hist(mean_draw_PC(:,3))
            fn_hist(mean_draw_C0(:,3))
            fn_hist(mean_draw_PC0(:,3))
        end
        hold off
        xlabel('\phi','FontSize',11)
        plotTickLatex2D('FontSize',11); 
        leg = legend('Uncensored','CP 10\% ','PCP 10\% ','CP 0','PCP 0');
        set(leg,'Interpreter','latex','FontSize',11)    
    end


end