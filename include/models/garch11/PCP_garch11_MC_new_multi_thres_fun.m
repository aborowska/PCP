function PCP_garch11_MC_new_multi_thres_fun(T, sigma2, S, II)
% T = 1000; S = 20; II = 2; sigma2 = 2 ; % <------------------ !!! 
% T = 2000; S = 20; II = 2; sigma2 = 2 ; % <------------------ !!! 

    % clear all
    close all

    addpath(genpath('include/'));

%     s = RandStream('mt19937ar','Seed',1);
%     RandStream.setGlobalStream(s); 

    model = 'garch11'; 
    partition = 3;
%     model = 'garch11_v2';    partition = 2;
    fprintf('Model: %s.\n',model)
    parameters = {'$\\mu$','$\\omega$','$\\alpha$','$\\beta$'};
    fprintf('Time series length T = %d.\n',T)

    sigma1 = 1;
    % sigma2 = 2;
    c = (sigma2 - sigma1)/sqrt(2*pi); % mean of eps
    kappa = 0.5*(sigma1^2 + sigma2^2 - ((sigma2-sigma1)^2)/pi); % var of eps
    sigma1_k = sigma1/sqrt(kappa);
    sigma2_k = sigma2/sqrt(kappa);

    omega = 1;
    alpha = 0.1;
    beta = 0.8;
    mu_true = [0, omega, alpha, beta];
    param_true = [c,sigma2,omega,alpha,beta];
    mu_init = [0, 1, 0.05, 0.85];
    d = length(mu_init);
    
    % S = 20; % number of MC replications
    H = 100;

    % quantiles of interest
    p_bar0 = 0.005;
    p_bar1 = 0.01;
    p_bar = 0.05; 
        
    % Metropolis-Hastings for the parameters
    M = 10000; % number of draws 
    BurnIn = 5000; %10000; %1000

    x_gam = (0:0.00001:50)'+0.00001;
    GamMat = gamma(x_gam);

    cont = MitISEM_Control;
    cont.mit.CV_max = 1; %2?
    cont.mit.iter_max = 10;
    cont.mit.Hmax = 6;
    cont.mit.dfnc = 5;
    df = 5; % default df for a mit

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
    
    SDD = zeros(S,1);
    
    THRES = [ 0.1000    0.2000    0.3000    0.4000];
    %% MC Simulations
    tic
    s = 0;
    sdd = 40;
    while s < S  
        sdd = sdd + 1;
        %     if (mod(s,10)==0)
                fprintf(['\n',model, ' simulation no. %i\n'],s)
                fprintf('Time series length T = %d.\n',T)               
        %     end
        try
            results = PCP_garch11_new_multi_thres_run(sdd, c, sigma1, sigma2, kappa, omega,...
                alpha, beta, p_bar0, p_bar1, p_bar, T, H, M, BurnIn,...
                mu_init, df, cont, options, partition, II, GamMat,...
                THRES);
         
            s = s+1;
            SDD(s,:) = sdd;
            
            RES{s,1} = results;
        catch
            % nothing
        end
    end

    
    time_total = toc;

   
    name = ['results/',model,'/',model,'_',num2str(sigma1),'_',...
    num2str(sigma2),'_T',num2str(T),'_H',num2str(H),...
    '_II',num2str(II),'_PCP_NEW_MULTI_THRES_3.mat'];
    
    save(name,'time_total','SDD','RES');         

    
end