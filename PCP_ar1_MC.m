clear all
close all

addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

model = 'ar1';
fprintf('Model: %s.\n',model)
parameters = {'$\\mu$','$\\sigma$','$\\phi$'};

sigma1 = 1;
sigma2 = 2;
c = (sigma2 - sigma1)/sqrt(2*pi);

partition = 3;
S = 1; % number of MC replications
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

%%
T = 1000; % time series length
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

%% various display options
cont.disp = false;

v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014a)')
    fn_hist = @(xx) hist(xx,20);
else
    fn_hist = @(xx) histogram(xx);
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
for s = 1:S
%     if (mod(s,10)==0)
        fprintf(['\n',model, ' simulation no. %i\n'],s)
%     end

    %% simple AR(1)
    eps = randn(T+H,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);

    rho = 0.8;
    param_true = [c,sigma2,rho];
    y = zeros(T+H,1);

    y(1,1) = eps(1,1);
    for ii = 2:T+H
        y(ii,1) = rho*y(ii-1,1) + eps(ii,1);
    end
    % true VaRs
    q1(s,:) = norminv(p_bar1,c+rho*y(T:(T+H-1),1),sigma2)';
    q5(s,:) = norminv(p_bar,c+rho*y(T:(T+H-1),1),sigma2)'; 
     
    % MC VaRs under the true model
    eps_sort = randn(M,H);
    ind = (eps_sort>0);
    eps_sort(ind) = c + sigma1.*eps_sort(ind);
    eps_sort(~ind) = c + sigma2.*eps_sort(~ind);

    y_sort = bsxfun(@plus,eps_sort,rho*y(T:(T+H-1),1)');
    y_sort = sort(y_sort);
    VaR_1(s,:)  = y_sort(p_bar1*M,:); 
    VaR_5(s,:)  = y_sort(p_bar*M,:); 

    %% Misspecified model: AR1 normal with unknown mu and sigma
    %% UNCENSORED posterior
%     kernel_init = @(xx) -loglik_ar1(xx,y);
%     mu_init = [0,1,0.9];
   
    mu_init = [0,1,0.9];
    kernel_init = @(xx) -posterior_ar1_mex(xx,y(1:T))/T;
    kernel = @(xx) posterior_ar1_mex(xx,y(1:T));
    
    try
        [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
    catch
        [mu,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init,options);
        Sigma = inv(T*Sigma);
        mit = struct('mu',mu,'Sigma',reshape(Sigma,1,length(mu)^2),'df', df, 'p', 1);
        [mit, CV] = MitISEM_new(mit, kernel, mu_init, cont, GamMat);
    end
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

    y_post = bsxfun(@times,randn(M,H),draw(:,2));
    y_post = bsxfun(@plus,y_post,draw(:,1));
    y_post = y_post + draw(:,3)*y(T:(T+H-1),1)';
    y_post = sort(y_post);
    VaR_1_post(s,:)  = y_post(p_bar1*M,:); 
    VaR_5_post(s,:)  = y_post(p_bar*M,:); 
    mean_draw(s,:) = mean(draw);
    std_draw(s,:) = std(draw);
    
    %% Threshold = 10% perscentile of the data sample
    threshold = sort(y(1:T));
    threshold = threshold(2*p_bar*T);
    %% CENSORED
    kernel_init = @(xx) - C_posterior_ar1_mex(xx, y(1:T), threshold)/T;    
    kernel = @(xx) C_posterior_ar1_mex(xx, y(1:T), threshold);
    try
        [mit_C, CV_C] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);   
    catch
        [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
        Sigma_C = inv(T*Sigma_C);
        mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
        CV = cont.mit.CV_old;
    end
    [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
    lnd_C = dmvgt(draw_C, mit_C, true, GamMat);     
    lnw_C = lnk_C - lnd_C;
    lnw_C = lnw_C - max(lnw_C);
    [ind, a] = fn_MH(lnw_C);
    draw_C = draw_C(ind,:);
    accept_C(s,1) = a/(M+BurnIn);
    draw_C = draw_C(BurnIn+1:BurnIn+M,:);
    mean_draw_C(s,:) = mean(draw_C);
    std_draw_C(s,:) = std(draw_C);
    
%     y_post_C = draw_C(:,1) + draw_C(:,3).*y(T,1) + draw_C(:,2).*randn(M,1);
%     y_post_C = sort(y_post_C);
%     VaR_1_post_C(s,1) = y_post_C(p_bar1*M); 
%     VaR_5_post_C(s,1) = y_post_C(p_bar*M);   
    
    y_post_C = bsxfun(@times,randn(M,H),draw_C(:,2));
    y_post_C = bsxfun(@plus,y_post_C,draw_C(:,1));
    y_post_C = y_post_C + draw_C(:,3)*y(T:(T+H-1),1)';
    y_post_C = sort(y_post_C);
    VaR_1_post_C(s,:)  = y_post_C(p_bar1*M,:); 
    VaR_5_post_C(s,:)  = y_post_C(p_bar*M,:); 
    
    %% PARTIAL CENSORING: keep rho uncensored, then censor mu and sigma
    % mit_C: joint candidate for the joint censored posterior
    % Short version
    II = 100;
    % draw_short = draw(1:II,:);
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos 
    M_short = M/II;
    [draw_PC, a_PC] = sim_cond_mit_MH(mit_C, draw_short, partition, M_short, BurnIn, kernel, GamMat);
    

% M = M_org
% BurnIn = BurnIn_org;
% draw = draw_org;
    II2 = 10;
    draw_short2 = draw((1:II2:M)',:); % thinning - to get hight quality rhos 
    thinning = 10;
    % H_org = H;
    % M_org = M;
    % M = II2;
    % draw_org = draw;
    % draw = draw_short2;
    % mit_org = mit;
    % mit = mit_C;    
    % BurnIn_org = BurnIn;
    % BurnIn = 2*BurnIn;
    [draw_PC2, a_PC2] = sim_cond_mit_MH_outloop(mit_C, draw_short2, partition, II2, 2*BurnIn, kernel, GamMat, thinning, cont.disp);
    % H = H_org;
    accept_PC(s,1) = mean(a_PC);
    mean_draw_PC(s,:) = mean(draw_PC);
    std_draw_PC(s,:) = std(draw_PC);

if plot_on
    ff = figure(22);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
     
    subplot(1,3,1)
    hold on
    fn_hist(draw(:,1))
    fn_hist(draw_PC(:,1))
    fn_hist(draw_PC2(:,1))
    hold off
    xlabel('\mu','FontSize',11)
    plotTickLatex2D('FontSize',11);
    
    subplot(1,3,2)
    hold on
    fn_hist(draw(:,2))
    fn_hist(draw_PC(:,2))
    fn_hist(draw_PC2(:,2))
    hold off
    xlabel('\sigma','FontSize',11)
    plotTickLatex2D('FontSize',11);
    
    subplot(1,3,3)
    hold on
    fn_hist(draw(:,3))
    fn_hist(draw_PC(:,3))
    fn_hist(draw_PC2(:,3))
    hold off
    xlabel('\rho','FontSize',11)
    plotTickLatex2D('FontSize',11); 
    
    leg = legend('Uncensored','PCP 10\% II=100','PCP 10\% II=10');
    set(leg,'Interpreter','latex','FontSize',11)
 
    if save_on
        name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP_IIcomp.eps'];
        set(gcf,'PaperPositionMode','auto');
        print_fail = 1;
        while print_fail 
            try                 
                print(ff,name,'-depsc','-r0')
                print_fail = 0;
            catch
                print_fail = 1;
            end
        end
    end   
end
    
%     y_post_PC = draw_PC(:,1) + draw_PC(:,3).*y(T,1) + draw_PC(:,2).*randn(M,1);
%     y_post_PC = sort(y_post_PC);
%     VaR_1_post_PC(s,1) = y_post_PC(p_bar1*M); 
%     VaR_5_post_PC(s,1) = y_post_PC(p_bar*M); 

    y_post_PC = bsxfun(@times,randn(M,H),draw_PC(:,2));
    y_post_PC = bsxfun(@plus,y_post_PC,draw_PC(:,1));
    y_post_PC = y_post_PC + draw_PC(:,3)*y(T:(T+H-1),1)';
    y_post_PC = sort(y_post_PC);
    VaR_1_post_PC(s,:)  = y_post_PC(p_bar1*M,:); 
    VaR_5_post_PC(s,:)  = y_post_PC(p_bar*M,:); 
    
    %% Threshold = 0
    threshold0 = 0;
    %% CENSORED
    kernel_init = @(xx) - C_posterior_ar1_mex(xx, y(1:T), threshold0)/T; 
    kernel = @(xx) C_posterior_ar1_mex(xx, y(1:T), threshold0);

    try
        [mit_C0, CV_C0] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
    catch
        [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu_init,options);
        Sigma_C0 = inv(T*Sigma_C0);
        mit_C0 = struct('mu',mu_C0,'Sigma',reshape(Sigma_C0,1,length(mu_C0)^2),'df', df, 'p', 1);
        CV_0 = cont.mit.CV_old;
    end
    [draw_C0, lnk_C0] = fn_rmvgt_robust(M+BurnIn, mit_C0, kernel, false);
    lnd_C0 = dmvgt(draw_C0, mit_C0, true, GamMat);    
    lnw_C0 = lnk_C0 - lnd_C0;
    lnw_C0 = lnw_C0 - max(lnw_C0);
    [ind, a] = fn_MH(lnw_C0);
    draw_C0 = draw_C0(ind,:);
    accept_C0(s,1) = a/(M+BurnIn);
    draw_C0 = draw_C0(BurnIn+1:BurnIn+M,:);
    mean_draw_C0(s,:) = mean(draw_C0);
    std_draw_C0(s,:) = std(draw_C0);
    
%     y_post_C0 = draw_C0(:,1) + draw_C0(:,3).*y(T,1) + draw_C0(:,2).*randn(M,1);
%     y_post_C0 = sort(y_post_C0);
%     VaR_1_post_C0(s,1) = y_post_C0(p_bar1*M); 
%     VaR_5_post_C0(s,1) = y_post_C0(p_bar*M); 
%     
    y_post_C0 = bsxfun(@times,randn(M,H),draw_C0(:,2));
    y_post_C0 = bsxfun(@plus,y_post_C0,draw_C0(:,1));
    y_post_C0 = y_post_C0 + draw_C0(:,3)*y(T:(T+H-1),1)';
    y_post_C0 = sort(y_post_C0);
    VaR_1_post_C0(s,:) = y_post_C0(p_bar1*M,:); 
    VaR_5_post_C0(s,:) = y_post_C0(p_bar*M,:); 

    %% PARTIAL CENSORING: keep rho uncensored, then censor mu and sigma
    % mit_C0: joint candidate for the joint censored posterior
    
    % Short version
%     II = 100;
%     draw_short = draw(1:II,:);
%     M_short = M/II;

    [draw_PC0, a_PC0] = sim_cond_mit_MH(mit_C0, draw_short, partition, M_short, BurnIn, kernel, GamMat);
    accept_PC0(s,1) = mean(a_PC0);
    mean_draw_PC0(s,:) = mean(draw_PC0);
    std_draw_PC0(s,:) = std(draw_PC0);    

%     y_post_PC0 = draw_PC0(:,1) + draw_PC0(:,3).*y(T,1) + draw_PC0(:,2).*randn(M,1);
%     y_post_PC0 = sort(y_post_PC0);
%     VaR_1_post_PC0(s,1) = y_post_PC0(p_bar1*M); 
%     VaR_5_post_PC0(s,1) = y_post_PC0(p_bar*M); 

    y_post_PC0 = bsxfun(@times,randn(M,H),draw_PC0(:,2));
    y_post_PC0 = bsxfun(@plus,y_post_PC0,draw_PC0(:,1));
    y_post_PC0 = y_post_PC0 + draw_PC0(:,3)*y(T:(T+H-1),1)';
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
    name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_PCP0_MC_',v_new,'.mat'];
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

% print_table_pcp_mc(model,parameters,sigma1,sigma2)
% print_table_pcp_mc(model,parameters,sigma1,sigma2,H)
% print_table_pcp_mc(model,parameters,sigma1,1,H)

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


