clear all
close all

addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

model = 'garch11';
fprintf('Model: %s.\n',model)
parameters = {'$\\mu$','$\\omega$','$\\alpha$','$\\beta$'};

sigma1 = 1;
sigma2 = 2;
c = (sigma2 - sigma1)/sqrt(2*pi); % mean of eps
kappa = 0.5*(sigma1^2 + sigma2^2 - ((sigma2-sigma1)^2)/pi); % var of eps
sigma1_k = sigma1/sqrt(kappa);
sigma2_k = sigma2/sqrt(kappa);

omega = 1;
alpha = 0.1;
beta = 0.8;
mu_true = [0, omega, alpha, beta];
param_true = [c,sigma2,omega,alpha,beta];
% partition = 3;
S = 100; % number of MC replications
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

%% simualted parameters:
% true model, posterior, censored posterior 10%, partially censored posterior 10%, 
% censored posterior at 0, partially censored posterior at 0
mean_draw = zeros(S,4);
mean_draw_C = zeros(S,4);
mean_draw_PC = zeros(S,4);
mean_draw_C0 = zeros(S,4);
mean_draw_PC0 = zeros(S,4);

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

%%
T = 1000; % time series length
fprintf('Time series length T = %d.\n',T)

% Metropolis-Hastings for the parameters
M = 10000; % number of draws 
BurnIn = 1000;

x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);

cont = MitISEM_Control;
cont.mit.iter_max = 10;
cont.mit.Hmax = 6;
cont.mit.dfnc = 5;
df = 5; % default df for a mit

%% various display options
cont.disp = true;%false;

v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014a)')
    fn_hist = @(xx) hist(xx,20);
else
    fn_hist = @(xx) histogram(xx,20);
end

plot_on = false;
save_on = false; %true;

options = optimset('Display','off');
% w = warning('query','last');
% id = w.identifier;
id = 'optim:fminunc:SwitchingMethod';
warning('off',id);

%% MC Simulations
tic
for s = 1:S
    %     if (mod(s,10)==0)
            fprintf(['\n',model, ' simulation no. %i\n'],s)
    %     end

    %% GARCH(1,1)
    eps = randn(T+H,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);
    eps = eps/sqrt(kappa);  
    y = zeros(T+H,1);
    h_true = zeros(T+H,1);

    for ii = 1:T+H
        if (ii == 1)
            h_true(ii,1) = omega;
        else
            h_true(ii,1) = omega*(1-alpha-beta) + alpha*(y(ii-1,1))^2 + beta*h_true(ii-1,1);
        end
        y(ii,1) = sqrt(h_true(ii,1))*eps(ii,1);
    end
    
    % true VaRs
    q1(s,:) = norminv(p_bar1, 0, h_true(T+1:T+H))';
    q5(s,:) = norminv(p_bar, 0, h_true(T+1:T+H))'; 
    
    % MC VaRs under the true model
    eps_sort = randn(M,H);
    ind = (eps_sort>0);
    eps_sort(ind) = c + sigma1.*eps_sort(ind);
    eps_sort(~ind) = c + sigma2.*eps_sort(~ind);
    eps_sort = eps_sort/sqrt(kappa); 

    y_sort = bsxfun(@times,eps_sort,sqrt(h_true(T+1:T+H,1))');
    y_sort = sort(y_sort);
    VaR_1(s,:) = y_sort(p_bar1*M,:); 
    VaR_5(s,:) = y_sort(p_bar*M,:); 

    %% Misspecified model: GARCH(1,1) normal 
    %% Uncensored Posterior
    % theta  = [mu, omega, alpha, beta]
    mu_init = [0, 0.1, 0.05, 0.85];
    y_S = var(y(1:T));
    mu_init(1,1) = mean(y(1:T));
    mu_init(1,2) = y_S;
   
    kernel_init = @(xx) -posterior_garch11_mex(xx, y(1:T), y_S)/T;
    kernel = @(xx) posterior_garch11_mex(xx, y(1:T), y_S);

    [mu,~,~,~,~,Sigma] = fminunc(kernel_init, mu_init);
    Sigma = inv(T*Sigma);
    % draw = rmvt(mu,Sigma,df,M+BurnIn);
    % lnk = kernel(draw);
    % lnd = dmvgt_mex(draw, mu, Sigma, df, 1, GamMat, double(1));
    [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
    [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
    lnd = dmvgt(draw, mit, true, GamMat); 

    lnw = lnk - lnd;
    lnw = lnw - max(lnw);
    [ind, a] = fn_MH(lnw);
    draw = draw(ind,:);
    lnw = lnw(ind);
    accept(s,1) = a/(M+BurnIn);
    draw = draw(BurnIn+1:BurnIn+M,:);    
    lnw = lnw(BurnIn+1:BurnIn+M,:);    
       
    h_post = volatility_garch11(draw,y,y_S,H);
    y_post = bsxfun(@times,randn(M,H),sqrt(h_post(T+1:T+H,1))');
    y_post = bsxfun(@plus,y_post,draw(:,1));
    y_post = sort(y_post);
    VaR_1_post(s,:) = y_post(p_bar1*M,:); 
    VaR_5_post(s,:) = y_post(p_bar*M,:); 
    mean_draw(s,:) = mean(draw);
    std_draw(s,:) = std(draw);
    
    %% Threshold = 10% perscentile of the data sample
    threshold = sort(y(1:T));
    threshold = threshold(2*p_bar*T);
    %% CENSORED
    % kernel_init = @(xx) - C_posterior_garch11(xx, y, threshold, y_S)/T;    
    % kernel = @(xx) C_posterior_garch11(xx, y, threshold, y_S);
    kernel_init = @(xx) - C_posterior_garch11_mex(xx, y(1:T,1), threshold, y_S)/T;    
    kernel = @(xx) C_posterior_garch11_mex(xx, y(1:T,1), threshold, y_S);

    try
        [mit_C, CV_C] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat); 
        if CV_C(end)>2
            [mit_C, CV_C] = MitISEM_new(mit_C, kernel, mu_init, cont, GamMat);   
        end        
        [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
        lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
    catch
        try
            [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
            Sigma_C = inv(T*Sigma_C);
            draw_C = rmvt(mu_C,Sigma_C,df,M+BurnIn);
            mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
            [mit_C, CV_C] = MitISEM_new(mit_C, kernel, mu_init, cont, GamMat);   
            if CV_C(end)>2
                [mit_C, CV_C] = MitISEM_new(mit_C, kernel, mu_init, cont, GamMat);   
            end
            [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
            lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
        catch
            mu_C = fminunc(kernel_init,mu_init,options);
            mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma,1,length(mu_C)^2),'df', df, 'p', 1);
            [mit_C, CV_C] = MitISEM_new(mit_C, kernel, mu_init, cont, GamMat);   
            if CV_C(end)>2
                [mit_C, CV_C] = MitISEM_new(mit_C, kernel, mu_init, cont, GamMat);   
            end
            [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
            lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
        end   
    end
    
    lnw_C = lnk_C - lnd_C;
    lnw_C = lnw_C - max(lnw_C);
    [ind, a] = fn_MH(lnw_C);
    draw_C = draw_C(ind,:);
    accept_C(s,1) = a/(M+BurnIn);
    draw_C = draw_C(BurnIn+1:BurnIn+M,:);

    h_post_C = volatility_garch11(draw_C,y,y_S,H);
    y_post_C = bsxfun(@times,randn(M,H),sqrt(h_post_C(T+1:T+H,1))');
    y_post_C = bsxfun(@plus,y_post_C,draw_C(:,1));
    y_post_C = sort(y_post_C);
    VaR_1_post_C(s,:) = y_post_C(p_bar1*M,:); 
    VaR_5_post_C(s,:) = y_post_C(p_bar*M,:); 
    mean_draw_C(s,:) = mean(draw_C);
    std_draw_C(s,:) = std(draw_C);
    
    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor mu and sigma
    % mit_C: joint candidate for the joint censored posterior    
    partition = 3; 

    II = 100;
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
    M_short = M/II;
    [draw_PC, a_PC] = sim_cond_mit_MH(mit_C, draw_short, partition, M_short, BurnIn, kernel, GamMat);
    accept_PC(s,1) = mean(a_PC);
    mean_draw_PC(s,:) = mean(draw_PC);
    std_draw_PC(s,:) = std(draw_PC);
    
    h_post_PC = volatility_garch11(draw_PC,y,y_S,H);
    y_post_PC = bsxfun(@times,randn(M,H),sqrt(h_post_PC(T+1:T+H,1))');
    y_post_PC = bsxfun(@plus,y_post_PC,draw_PC(:,1));
    y_post_PC = sort(y_post_PC);
    VaR_1_post_PC(s,:) = y_post_PC(p_bar1*M,:); 
    VaR_5_post_PC(s,:) = y_post_PC(p_bar*M,:); 

    %% Threshold = 0
    threshold0 = 0;
    %% CENSORED
    kernel_init = @(xx) - C_posterior_garch11_mex(xx, y(1:T,1), threshold0, y_S)/T;    
    kernel = @(xx) C_posterior_garch11_mex(xx, y(1:T,1), threshold0, y_S);
    try
        [mit_C0, CV_C0] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
        [draw_C0, lnk_C0] = fn_rmvgt_robust(M+BurnIn, mit_C0, kernel, false);
        lnd_C0 = dmvgt(draw_C0, mit_C0, true, GamMat);           
    catch
        try
            [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu_init,options);
            Sigma_C0 = inv(T*Sigma_C0);
            mit_C0 = struct('mu',mu_C0,'Sigma',reshape(Sigma_C0,1,length(mu_C0)^2),'df', df, 'p', 1);
            draw_C0 = rmvt(mu_C0,Sigma_C0,df,M+BurnIn);
            [mit_C0, CV_C0] = MitISEM_new(mit_C0, kernel, mu_init, cont, GamMat);   
            [draw_C0, lnk_C0] = fn_rmvgt_robust(M+BurnIn, mit_C0, kernel, false);
            lnd_C0 = dmvgt(draw_C0, mit_C0, true, GamMat);            
        catch
            mu_C0 = fminunc(kernel_init,mu_init,options);
            [~, ind_aux] = max(mit_C.p);
            Sigma_aux = mit_C.Sigma(ind_aux,:);
            mit_C0 = struct('mu',mu_C0,'Sigma',Sigma_aux,'df', df, 'p', 1);
            [mit_C0, CV_C0] = MitISEM_new(mit_C0, kernel, mu_init, cont, GamMat);   
            if CV_C0(end)>2
                [mit_C0, CV_C0] = MitISEM_new(mit_C0, kernel, mu_init, cont, GamMat);   
            end
            [draw_C0, lnk_C0] = fn_rmvgt_robust(M+BurnIn, mit_C0, kernel, false);
            lnd_C0 = dmvgt(draw_C0, mit_C0, true, GamMat);                
        end
    end
    lnw_C0 = lnk_C0 - lnd_C0;
    lnw_C0 = lnw_C0 - max(lnw_C0);
    [ind, a] = fn_MH(lnw_C0);
    draw_C0 = draw_C0(ind,:);
    accept_C0(s,1) = a/(M+BurnIn);
    draw_C0 = draw_C0(BurnIn+1:BurnIn+M,:);
    mean_draw_C0(s,:) = mean(draw_C0);
    std_draw_C0(s,:) = std(draw_C0);

    h_post_C0 = volatility_garch11(draw_C0,y,y_S,H);
    y_post_C0 = bsxfun(@times,randn(M,H),sqrt(h_post_C0(T+1:T+H,1))');
    y_post_C0 = bsxfun(@plus,y_post_C0,draw_C0(:,1));
    y_post_C0 = sort(y_post_C0);
    VaR_1_post_C0 = y_post_C0(p_bar1*M); 
    VaR_5_post_C0 = y_post_C0(p_bar*M); 
    
    %% PARTIAL CENSORING: keep alpha and beta uncensored, then censor mu and sigma
    % mit_C0: joint cnadidate for the joint censored posterior
    [draw_PC0, a_PC0] = sim_cond_mit_MH(mit_C0, draw_short, partition, M_short, BurnIn, kernel, GamMat);
    accept_PC0(s,1) = mean(a_PC0);
    mean_draw_PC0(s,:) = mean(draw_PC0);
    std_draw_PC0(s,:) = std(draw_PC0);    

    h_post_PC0 = volatility_garch11(draw_PC0,y,y_S,H);
    y_post_PC0 = bsxfun(@times,randn(M,H),sqrt(h_post_PC0(T+1:T+H,1))');
    y_post_PC0 = bsxfun(@plus,y_post_PC0,draw_PC0(:,1));
    y_post_PC0 = sort(y_post_PC0);
    VaR_1_post_PC0(s,:) = y_post_PC0(p_bar1*M,:); 
    VaR_5_post_PC0(s,:) = y_post_PC0(p_bar*M,:); 
end

% MSEs
MSE_1 = mean((VaR_1 - q1).^2);
MSE_1_post = mean((VaR_1_post - q1).^2);
MSE_1_post_C = mean((VaR_1_post_C - q1).^2);
MSE_1_post_PC = mean((VaR_1_post_PC - q1).^2);
MSE_1_post_C0 = mean((VaR_1_post_C0 - q1).^2);
MSE_1_post_PC0 = mean((VaR_1_post_PC0 - q1).^2);

MSE_5 = mean((VaR_5 - q5).^2);
MSE_5_post = mean((VaR_5_post - q5).^2);
MSE_5_post_C = mean((VaR_5_post_C - q5).^2);
MSE_5_post_PC = mean((VaR_5_post_PC - q5).^2);
MSE_5_post_C0 = mean((VaR_5_post_C0 - q5).^2);
MSE_5_post_PC0 = mean((VaR_5_post_PC0 - q5).^2);

time_total = toc;

if save_on
    name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_PCP0_MC2_',v_new,'.mat'];
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
end