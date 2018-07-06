close all 
clear all
addpath(genpath('include/'));


T = 1000; 
II = 10;
II_Q = 1;
sigma2 = 2 ;
fprintf('Time series length T = %d.\n',T)
   
model = 'garch11'; 
%     model = 'garch11_v2';    partition = 2;
fprintf('Model: %s.\n',model)
parameters = {'$\\mu$','$\\omega$','$\\alpha$','$\\beta$'};
partition = 3;

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
    
 
H = 1;

% quantiles of interest
p_bar0 = 0.005;
p_bar1 = 0.01;
p_bar = 0.05;

% Metropolis-Hastings for the parameters
M = 10000; % number of draws 
BurnIn = 1000;
thinning = 10;
BurnIn_PCP = BurnIn/10;
     
x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);

df = 5; % default df for a mit
cont = MitISEM_Control;
cont.mit.iter_max = 10;
cont.mit.Hmax = 6;
cont.mit.dfnc = 5;
cont.mit.CV_max = 1.9;

options = optimset('Display','off');
% w = warning('query','last');
% id = w.identifier;
id = 'optim:fminunc:SwitchingMethod';
warning('off',id);
    
S = 50;
VaR_post_PC = zeros(S,3);
VaR_post_PCQ = zeros(S,3);
VaR_post_PCQ_eps = zeros(S,3);
ES_post_PC = zeros(S,3);
ES_post_PCQ = zeros(S,3);
ES_post_PCQ_eps = zeros(S,3);   


%% GARCH(1,1)
sdd = 0;
s = RandStream('mt19937ar','Seed',sdd);
RandStream.setGlobalStream(s); 
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
y_S = var(y(1:T));

% Threshold = 10% perscentile of the data sample
threshold = sort(y(1:T));
threshold = threshold(2*p_bar*T);

% true VaRs ;
q05 = norminv(p_bar0, c, sigma2_k*sqrt(h_true(T+1:T+H)))';
q1 = norminv(p_bar1, c, sigma2_k*sqrt(h_true(T+1:T+H)))';
q5 = norminv(p_bar, c, sigma2_k*sqrt(h_true(T+1:T+H)))'; 

% true ESs
cdf05 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar0))/(p_bar0);
cdf1 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar1))/(p_bar1);
cdf5 = c - sigma2_k*sqrt(h_true(T+1:T+H))'*normpdf(norminv(1-p_bar))/(p_bar);    
    
parfor ss = 1:S
    s = RandStream('mt19937ar','Seed',ss);
    RandStream.setGlobalStream(s);   
    
    %% Misspecified model: GARCH(1,1) normal 
    %% Uncensored Posterior
    fprintf('*** Uncensored Posterior ***\n');
 
    kernel_init = @(xx) -posterior_garch11_mex(xx, y(1:T), y_S)/T;
    kernel = @(xx) posterior_garch11_mex(xx, y(1:T), y_S);
    try
        [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
        [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
        lnd = dmvgt(draw, mit, true, GamMat); 
    catch
        [mu,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init,options);
        Sigma = inv(T*Sigma);
        draw = rmvt(mu,Sigma,df,M+BurnIn);
        mit = struct('mu',mu,'Sigma',reshape(Sigma,1,length(mu)^2),'df', df, 'p', 1);
        [mit, CV] = MitISEM_new(mit, kernel, mu_init, cont, GamMat);            
        [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
        lnd = dmvgt(draw, mit, true, GamMat); 
    end

    lnw = lnk - lnd;
    lnw = lnw - max(lnw);
    [ind, a] = fn_MH(lnw);
    draw = draw(ind,:);
    accept = a/(M+BurnIn);
    draw = draw(BurnIn+1:BurnIn+M,:);    

 
    %% CENSORED
    fprintf('*** Censored Posterior, threshold 10%% ***\n');
    kernel_init = @(xx) - C_posterior_garch11_mex(xx, y(1:T,1), threshold, y_S)/T;    
    kernel = @(xx) C_posterior_garch11_mex(xx, y(1:T,1), threshold, y_S);
    
    CV_C = cont.mit.CV_old;
    while (CV_C(end) >= 2)
        try
            [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
            Sigma_C = inv(T*Sigma_C);
            draw_C = rmvt(mu_C,Sigma_C,df,M+BurnIn);
            mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
            [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            if CV_C(end)>2
                [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            end
            [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
            lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
        catch
            [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
            mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
            [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            if CV_C(end)>2
                [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            end
            [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
            lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
        end 
    end
    lnw_C = lnk_C - lnd_C;
    lnw_C = lnw_C - max(lnw_C);
    [ind, a] = fn_MH(lnw_C);
    draw_C = draw_C(ind,:);
    accept_C = a/(M+BurnIn);
    draw_C = draw_C(BurnIn+1:BurnIn+M,:);
    

    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior, threshold 10%% ***\n');
    % mit_C: joint candidate for the joint censored posterior
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
    [draw_PC, a_PC, lnw_PC] = sim_cond_mit_MH_outloop(mit_C, draw_short,...
        partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
    accept_PC = mean(a_PC); 
   
    eps_PC = randn(M,H);
    h_post_PC = volatility_garch11(draw_PC,y,y_S,H);
    y_post_PC = eps_PC.*sqrt(h_post_PC);
    y_post_PC = bsxfun(@plus,y_post_PC,draw_PC(:,1));
    [y_post_PC, ind_PC] = sort(y_post_PC);

    VaR_1_post_PC = y_post_PC(round(p_bar1*M),:); 
    VaR_5_post_PC = y_post_PC(round(p_bar*M),:); 
    VaR_05_post_PC = y_post_PC(round(p_bar0*M),:); 
    VaR_post_PC(ss,:) = [VaR_05_post_PC, VaR_1_post_PC, VaR_5_post_PC];

    ES_1_post_PC = mean(y_post_PC(1:round(p_bar1*M),:)); 
    ES_5_post_PC = mean(y_post_PC(1:round(p_bar*M),:)); 
    ES_05_post_PC = mean(y_post_PC(1:round(p_bar0*M),:));
    ES_post_PC(ss,:) = [ES_05_post_PC, ES_1_post_PC, ES_5_post_PC];
  

    %% QERMIT WITH PARTIAL CENSORING
    fprintf('*** QERMit for PCP ***\n');

    eps_PC_hl = eps_PC(ind_PC,:);
    draw_PC_hl = [eps_PC(ind_PC,:),draw_PC(ind_PC,:)];
    hl_boundary = y_post_PC(round(10*p_bar1*M)); % 10% lowest
    eps_PC_hl = eps_PC_hl(1:round(10*p_bar1*M),:);
    draw_PC_hl = draw_PC_hl(1:round(10*p_bar1*M),:);

    mu_PC_hl = mean(draw_PC_hl);
    Sigma_PC_hl = cov(draw_PC_hl);
    mit0_PC_hl = struct('mu',mu_PC_hl,'Sigma',reshape(Sigma_PC_hl,1,length(mu_PC_hl)^2),...
        'df', df, 'p', 1);
    mit_PC_hl = fn_optimt(draw_PC_hl, mit0_PC_hl, ones(size(draw_PC_hl,1),1), cont, GamMat);


    mu_eps_hl = mean(eps_PC_hl);
    Sigma_eps_hl = cov(eps_PC_hl);
    mit0_eps_hl = struct('mu',mu_eps_hl,'Sigma',reshape(Sigma_eps_hl,1,length(mu_eps_hl)^2),...
        'df', df, 'p', 1);
    mit_eps_hl = fn_optimt(eps_PC_hl, mit0_eps_hl, ones(size(eps_PC_hl,1),1), cont, GamMat);

    if false
        params = [{'$\\epsilon$'},parameters];
set(0,'defaultTextInterpreter','latex')
        subplot(3,2,1)
        hold on
        scatter(draw_PC(:,2-1),draw_PC(:,3-1))
        scatter(draw_PC_hl(:,2),draw_PC_hl(:,3))       
        xlabel(params{2})
        ylabel(params{3});
        hold off

        subplot(3,2,2)
        hold on
        scatter(draw_PC(:,2-1),draw_PC(:,4-1))
        scatter(draw_PC_hl(:,2),draw_PC_hl(:,4))       
        xlabel(params{2})
        ylabel(params{4});
        hold off

        subplot(3,2,3)
        hold on
        scatter(draw_PC(:,2-1),draw_PC(:,5-1))
        scatter(draw_PC_hl(:,2),draw_PC_hl(:,5))       
        xlabel(params{2})
        ylabel(params{5});
        hold off     
        

        subplot(3,2,4)
        hold on
        scatter(draw_PC(:,3-1),draw_PC(:,4-1))
        scatter(draw_PC_hl(:,3),draw_PC_hl(:,4))       
        xlabel(params{3})
        ylabel(params{4});
        hold off

        subplot(3,2,5)
        hold on
        scatter(draw_PC(:,3-1),draw_PC(:,5-1))
        scatter(draw_PC_hl(:,3),draw_PC_hl(:,5))       
        xlabel(params{3})
        ylabel(params{5});
        hold off

        subplot(3,2,6)
        hold on
        scatter(draw_PC(:,4-1),draw_PC(:,5-1))
        scatter(draw_PC_hl(:,4),draw_PC_hl(:,5))       
        xlabel(params{4})
        ylabel(params{5});
        hold off             
%         legend('All direct PCP','HL direct PCP')
    suptitle('Blue: all direct PCP, orange: HL direct PCP')
    end    
    
    kernel = @(xx) -0.5*(log(2*pi) + xx(:,1).^2);
    [draw_PCQ, LND_cond, LND_marg] = sim_cond_mit_outloop(mit_PC_hl, ...
        [zeros(M/II_Q,1), draw_PC((1:II_Q:M)',:)], 2, II_Q, GamMat);
    eps_PCQ = draw_PCQ(:,1);
    draw_PCQ = draw_PCQ(:,2:end);
    LNK_cond = kernel(eps_PCQ);
    LNW_cond = LNK_cond - LND_cond;

      
    h_post_PCQ = volatility_garch11(draw_PCQ,y,y_S,H);
    y_post_PCQ = eps_PCQ.*sqrt(h_post_PCQ);
    y_post_PCQ = bsxfun(@plus,y_post_PCQ,draw_PCQ(:,1));
    [y_post_PCQ, ind_PCQ] = sort(y_post_PCQ);
 
  
    w_PCQ = exp(LNW_cond)/M;
    w_PCQ = w_PCQ(ind_PCQ,:);
    cum_w = cumsum(w_PCQ);
    ind_var_05 = min(find(cum_w >= p_bar0))-1; 
    ind_var_1 = min(find(cum_w >= p_bar1))-1; 
    ind_var_5 = min(find(cum_w >= p_bar))-1;     

    VaR_1_post_PCQ = y_post_PCQ(ind_var_1,:); 
    VaR_5_post_PCQ = y_post_PCQ(ind_var_5,:); 
    VaR_05_post_PCQ = y_post_PCQ(ind_var_05,:); 
    VaR_post_PCQ(ss,:) = [VaR_05_post_PCQ, VaR_1_post_PCQ, VaR_5_post_PCQ];

    ES_1_post_PCQ = sum((w_PCQ(1:ind_var_1)/sum(w_PCQ(1:ind_var_1))).*y_post_PCQ(1:ind_var_1,:)); 
    ES_5_post_PCQ = sum((w_PCQ(1:ind_var_5)/sum(w_PCQ(1:ind_var_5))).*y_post_PCQ(1:ind_var_5,:)); 
    ES_05_post_PCQ = sum((w_PCQ(1:ind_var_05)/sum(w_PCQ(1:ind_var_05))).*y_post_PCQ(1:ind_var_05,:)); 
    ES_post_PCQ(ss,:) = [ES_05_post_PCQ, ES_1_post_PCQ, ES_5_post_PCQ];

     %% JUST EPS from HL independent from THETA
    eps_PCQ_eps = rmvgt2(M, mit_eps_hl.mu, mit_eps_hl.Sigma, mit_eps_hl.df, mit_eps_hl.p);
    draw_PCQ_eps = draw_PC;

    LNK_eps = kernel(eps_PCQ_eps);
    LND_eps = dmvgt(eps_PCQ_eps, mit_eps_hl, true, GamMat);
    LNW_eps = LNK_eps - LND_eps;        
    w_eps = exp(LNW_eps)/M;
   
      
    h_post_PCQ_eps = volatility_garch11(draw_PCQ_eps,y,y_S,H);
    y_post_PCQ_eps = eps_PCQ_eps.*sqrt(h_post_PCQ_eps);
    y_post_PCQ_eps = bsxfun(@plus,y_post_PCQ_eps,draw_PCQ_eps(:,1));
    [y_post_PCQ_eps, ind_PCQ_eps] = sort(y_post_PCQ_eps);
 
    w_eps = w_eps(ind_PCQ_eps,:);
    cum_w_eps = cumsum(w_eps);

    ind_var_05_eps = min(find(cum_w_eps >= p_bar0))-1; 
    ind_var_1_eps = min(find(cum_w_eps >= p_bar1))-1; 
    ind_var_5_eps = min(find(cum_w_eps >= p_bar))-1;     

    VaR_1_post_PCQ_eps = y_post_PCQ_eps(ind_var_1_eps,:); 
    VaR_5_post_PCQ_eps = y_post_PCQ_eps(ind_var_5_eps,:); 
    VaR_05_post_PCQ_eps = y_post_PCQ_eps(ind_var_05_eps,:); 
    VaR_post_PCQ_eps(ss,:) = [VaR_05_post_PCQ_eps, VaR_1_post_PCQ_eps, VaR_5_post_PCQ_eps];

    ES_1_post_PCQ_eps = sum((w_eps(1:ind_var_1_eps)/sum(w_eps(1:ind_var_1_eps))).*y_post_PCQ_eps(1:ind_var_1_eps,:)); 
    ES_5_post_PCQ_eps = sum((w_eps(1:ind_var_5_eps)/sum(w_eps(1:ind_var_5_eps))).*y_post_PCQ_eps(1:ind_var_5_eps,:)); 
    ES_05_post_PCQ_eps = sum((w_eps(1:ind_var_05_eps)/sum(w_eps(1:ind_var_05_eps))).*y_post_PCQ_eps(1:ind_var_05_eps,:)); 
    ES_post_PCQ_eps(ss,:) = [ES_05_post_PCQ_eps, ES_1_post_PCQ_eps, ES_5_post_PCQ_eps];   
end
 
save('PCPQ_garch11.mat','-regexp','^y','^q','^cdf','^VaR_post_PC','^ES_post_PC')

P_BARS = [p_bar0,p_bar1,p_bar]';

RES_VaR0 = [P_BARS,mean(VaR_post_PC)',mean(VaR_post_PCQ)',mean(VaR_post_PCQ_eps)'];
RES_ES0 = [P_BARS,mean(ES_post_PC)',mean(ES_post_PCQ)',mean(ES_post_PCQ_eps)'];

RES_VaR = [P_BARS,std(VaR_post_PC)',std(VaR_post_PCQ)',std(VaR_post_PCQ_eps)'];
RES_ES = [P_BARS,std(ES_post_PC)',std(ES_post_PCQ)',std(ES_post_PCQ_eps)'];

RES = [RES_VaR0;RES_VaR;RES_ES0;RES_ES];

RES = RES([1,4,7,10,2,5,8,11,3,6,9,12],:);