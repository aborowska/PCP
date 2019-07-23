function PCP_skt_agarch11_MC_new_multi_thres_fun_MCLE(T, lambda, nu, S, II)
% T = 1000; S = 20; II = 2; lambda = -0.3 ; nu = 5; % <------------------ !!! 
% T = 2000; S = 20; II = 2; lambda = -0.3 ; nu = 5; % <------------------ !!! 
% SIMULATION STUDY: DGP IS SKT-AGARCH BUT THE ESTIMATED MODEL IS NORMAL
% AGARCH HENCE THERE ARE 5 PARAMETERS TO BE ESTIMATED: MU, OMEGA, MU2,
% ALPHA, BETA


    % clear all
    close all

    addpath(genpath('include/'));

%     s = RandStream('mt19937ar','Seed',1);
%     RandStream.setGlobalStream(s); 

    model = 'skt_agarch11'; 
    fprintf('Model: %s.\n',model)
    parameters = {'$\\lambda$','$\\mu$','$\\omega$','$\\alpha$','$\\beta$'};
    fprintf('Time series length T = %d.\n',T)
    fprintf('Lambda  = %4.2f.\n',lambda)

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
    mu_true = [lambda,0, omega, 0, alpha, beta];
    param_true = [lambda,0,omega,mu2,alpha,beta];
    mu_init = [0, 1, 0.1, 0.05, 0.85];
    d = length(mu_init);
    
    % S = 20; % number of MC replications
    H = 100;

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
    RES = cell(S,1); 
    
    SDD = zeros(S,1); % RANDOM SEED NUMBER TO RECOVER GENERATED DATASETS
    
    THRES = [ 0.1000    0.2000    0.3000    0.4000];
%     THRES = [ 0.05000    0.1000    0.1500    0.2000];
    
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
    
    %% MC Simulations
    tic
%     s = 0;
%     sdd = 0;
%     while s < S  
    for s = 1:S  
        sdd = 20+s; %sdd + 1;
        %     if (mod(s,10)==0)
                fprintf(['\n',model, ' simulation no. %i\n'],s)
                fprintf('Time series length T = %d.\n',T)               
        %     end
        try
            results = PCP_skt_agarch11_new_multi_thres_run_MCLE(sdd,...
                lambda, nu, omega, alpha, beta,...
                p_bar0, p_bar1, p_bar, T, H, ...
                sampling_opt, cont, options, GamMat,...
                THRES);
         
%             s = s+1;
            SDD(s,:) = sdd;
            
            RES{s,1} = results;
        catch
            % nothing
        end
    end

    
    time_total = toc;

   
    name = ['results/',model,'/',model,'_',num2str(lambda),'_',num2str(nu),...
        '_T',num2str(T),'_H',num2str(H),'_II',num2str(II)...
        '_PCP_NEW_MULTI_THRES_MCLE_2.mat'];
    
    save(name,'time_total','SDD','RES','THRES');         

    
end