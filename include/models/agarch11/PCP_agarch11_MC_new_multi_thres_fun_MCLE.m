function PCP_agarch11_MC_new_multi_thres_fun_MCLE(T, sigma2, S, II)
% T = 1000; S = 20; II = 2; sigma2 = 2 ; % <------------------ !!! 
% T = 2000; S = 20; II = 2; sigma2 = 2 ; % <------------------ !!! 

    % clear all
    close all

    addpath(genpath('include/'));

    
 
    
%     s = RandStream('mt19937ar','Seed',1);
%     RandStream.setGlobalStream(s); 

    model = 'agarch11'; 
    fprintf('Model: %s.\n',model)
    parameters = {'$\\mu$','$\\omega$','$\\alpha$','$\\beta$'};
    fprintf('Time series length T = %d.\n',T)

    sigma1 = 1;
    c = (sigma2 - sigma1)/sqrt(2*pi); % mean of eps
    kappa = 0.5*(sigma1^2 + sigma2^2 - ((sigma2-sigma1)^2)/pi); % var of eps
    sigma1_k = sigma1/sqrt(kappa);
    sigma2_k = sigma2/sqrt(kappa);

    % gama = 0; % "typo" on purpose: not to confuse with the gamma function
    mu2 = 0; % gama = mu - mu2;    
    omega = 1;
    alpha = 0.1;
    beta = 0.8;
    % % theta  = [mu, gama, omega, alpha, beta]
    % mu_true = [0, 0, omega, alpha, beta];
    % param_true = [c,sigma2,gama,omega,alpha,beta];
    % mu_init = [0, -0.1, 1, 0.05, 0.85];
    % % mu_init = [mean(y), -0.01, 0.01, 0.05, 0.85];

    % theta  = [mu, omega, mu2, alpha, beta]
    mu_true = [0, omega, 0, alpha, beta];
    param_true = [c,sigma2,omega,mu2,alpha,beta];
    mu_init = [0, 1, 0.1, 0.05, 0.85];
    d = length(mu_init);
    
    % S = 20; % number of MC replications
    H = 1000;

    % quantiles of interest
    p_bar0 = 0.005;
    p_bar1 = 0.01;
    p_bar = 0.05; 
        
    x_gam = (0:0.00001:50)'+0.00001;
    GamMat = gamma(x_gam);

    cont = MitISEM_Control;
    cont.mit.CV_max = 1; %2?
    cont.mit.iter_max = 10;
    cont.mit.Hmax = 6;
    cont.mit.dfnc = 5;


    %% various display options
    cont.disp = false;

    save_on = true;

    options = optimset('Display','off');
    % w = warning('query','last');
    % id = w.identifier;
    id = 'optim:fminunc:SwitchingMethod';
    warning('off',id);
    
    %% MitISEM results: mits and CVs
    RES_MCLE = cell(S,1); 
   
    THRES = [ 0.1000    0.2000    0.3000    0.4000];
    
    %%  Sampling options
    % Metropolis-Hastings for the parameters
    M = 10000; % number of draws 
    BurnIn = 5000; %10000; %1000
    thinning = 10;
    BurnIn_PCP = BurnIn/5; %BurnIn/2; %BurnIn/10;
    df = 5; % default df for a mit   
    partition = 3;
    
    sampling_opt. M = M;
    sampling_opt.BurnIn = BurnIn;
    sampling_opt.BurnIn_PCP = BurnIn_PCP;
    sampling_opt.thinning = thinning;
    sampling_opt.df = df;
    sampling_opt.mu_init = mu_init;
    sampling_opt.partition = partition;
    sampling_opt.II = II;

%% REUSE
%     name_reuse = ['results/',model,'/',model,'_',num2str(sigma1),'_',...
%     num2str(sigma2),'_T',num2str(T),'_H',num2str(H),...
%     '_II',num2str(II),'_PCP_NEW_MULTI_THRES.mat'];    
    name_reuse = ['../../../../ownCloud2/ForDropbox/',model,'_',num2str(sigma1),'_',...
    num2str(sigma2),'_T',num2str(T),'_H',num2str(H),...
    '_II',num2str(II),'_PCP_NEW_MULTI_THRES.mat'];    
    load(name_reuse,'RES','SDD');
    
    %% MC Simulations
    tic
    for s = 1:S
        sdd = SDD(s); 
        %     if (mod(s,10)==0)
                fprintf(['\n',model, ' simulation no. %i\n'],s)
                fprintf('Time series length T = %d.\n',T)               
        %     end
        try
            RESULTS = RES{s,1};
            results = PCP_agarch11_new_multi_thres_run_MCLE(sdd,...
                c, sigma1, sigma2, kappa, omega, alpha, beta,...
                p_bar0, p_bar1, p_bar, T, H, ...
                sampling_opt, cont, options, GamMat,...
                THRES,...
                RESULTS);               
            
            RES_MCLE{s,1} = results;
        catch
            % nothing
        end
    end

    
    time_total = toc;

   
    name = ['results/',model,'/',model,'_',num2str(sigma1),'_',...
    num2str(sigma2),'_T',num2str(T),'_H',num2str(H),...
    '_II',num2str(II),'_PCP_NEW_MULTI_THRES_MCLE.mat'];
    
    save(name,'time_total','SDD','RES_MCLE');         

    
end