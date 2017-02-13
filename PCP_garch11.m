clear all
close all

addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

model = 'garch11';
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
% theta  = [mu, omega, alpha, beta]
mu_true = [0, omega, alpha, beta];
param_true = [c,sigma2,omega,alpha,beta];
mu_init = [0, 1, 0.05, 0.85];


T = 1000; % time series length
p_bar1 = 0.01;
p_bar = 0.05;
% Metropolis-Hastings for the parameters
M = 10000; % number of draws 
BurnIn = 1000;

x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);

cont = MitISEM_Control;
cont.mit.iter_max = 10;

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


%% GARCH(1,1)
eps = randn(T,1);
ind = (eps>0);
eps(ind) = c + sigma1.*eps(ind);
eps(~ind) = c + sigma2.*eps(~ind);
eps = eps/sqrt(kappa);  

omega = 1;
alpha = 0.1;
beta = 0.8;
mu_true = [0, omega, alpha, beta];

param_true = [c,sigma2,omega,alpha,beta];
y = zeros(T,1);
h_true = zeros(T,1);

for ii = 1:T
    if (ii == 1)
        h_true(ii,1) = omega;
    else
        h_true(ii,1) = omega*(1-alpha-beta) + alpha*(y(ii-1,1))^2 + beta*h_true(ii-1,1);
    end
    y(ii,1) = sqrt(h_true(ii,1))*eps(ii,1);
end

if plot_on
    figure(1)
    hold on
    plot(eps,'Color',[0.64 0.64 0.74])
    plot(y,'k')
    plot(h_true,'r')
    hold off
    leg = legend('$\varepsilon$','$y$','$h$');
    set(leg,'Interpreter','latex','FontSize',11)
    plotTickLatex2D('FontSize',12);
end

% MC VaRs under the true model
eps_sort = randn(M,1);
ind = (eps_sort>0);
eps_sort(ind) = c + sigma1.*eps_sort(ind);
eps_sort(~ind) = c + sigma2.*eps_sort(~ind);
eps_sort = eps_sort/sqrt(kappa); 

y_sort = sqrt(h_true(T,1))*eps_sort;
y_sort = sort(y_sort);
VaR_1 = y_sort(p_bar1*M); 
VaR_5 = y_sort(p_bar*M); 

%% Misspecified model: GARCH(1,1) normal 
%% Uncensored Posterior
y_S = var(y);
kernel_init = @(xx) -posterior_garch11_mex(xx, y, y_S)/T;
kernel = @(xx) posterior_garch11_mex(xx, y, y_S);
 
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
lnw = lnw(ind);
accept = a/(M+BurnIn);
draw = draw(BurnIn+1:BurnIn+M,:);    
lnw = lnw(BurnIn+1:BurnIn+M,:);    

h_post = volatility_garch11(draw,y,y_S);
y_post = draw(:,1) + sqrt(h_post).*randn(M,1);
y_post = sort(y_post);
VaR_1_post = y_post(p_bar1*M); 
VaR_5_post = y_post(p_bar*M); 

%% Threshold = 10% perscentile of the data sample
threshold = sort(y);
threshold = threshold(2*p_bar*T);
%% CENSORED
% kernel_init = @(xx) - C_posterior_garch11(xx, y, threshold, y_S)/T;    
% kernel = @(xx) C_posterior_garch11(xx, y, threshold, y_S);
kernel_init = @(xx) - C_posterior_garch11_mex(xx, y, threshold, y_S)/T;    
kernel = @(xx) C_posterior_garch11_mex(xx, y, threshold, y_S);

try
    [mit_C, CV_C] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat); 
    [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
    lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
catch
    try
        [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
        Sigma_C = inv(T*Sigma_C);
        draw_C = rmvt(mu_C,Sigma_C,df,M+BurnIn);
        mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
        [mit_C, CV_C] = MitISEM_new(mit_C, kernel, mu_init, cont, GamMat);   
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
accept_C = a/(M+BurnIn);
draw_C = draw_C(BurnIn+1:BurnIn+M,:);

h_post_C = volatility_garch11(draw_C, y, y_S);
y_post_C = draw_C(:,1) + sqrt(h_post_C).*randn(M,1);
y_post_C = sort(y_post_C);
VaR_1_post_C = y_post_C(p_bar1*M); 
VaR_5_post_C = y_post_C(p_bar*M); 


if save_on
    save(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'.mat'],...
    'y','draw','draw_C','param_true',...
    'accept','accept_C',...
    'mit','CV','mit_C','CV_C',...
    'VaR_1','VaR_1_post','VaR_1_post_C',...
    'VaR_5','VaR_5_post','VaR_5_post_C')
end


if plot_on
    ff = figure(1);
    if strcmp(v_new,'(R2014a)')
        set(gcf,'units','normalized','outerposition',[0.1 0.05 0.45 0.75]);
        
        ax1 = axes('Position',[0.05 0.60 0.36 0.38],'Visible','on');        
        hold on
        fn_hist(draw(:,1))
        fn_hist(draw_C(:,1))
        xlabel('\mu','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([c c], YL,'Color','r','LineWidth',3); 
        hold off
  
        ax2 = axes('Position',[0.46 0.60 0.36 0.38],'Visible','on');        
        hold on
        fn_hist(draw(:,2))
        fn_hist(draw_C(:,2))
        xlabel('\omega','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([omega omega], YL,'Color','r','LineWidth',3); 
        hold off

        
        ax3 = axes('Position',[0.05 0.10 0.36 0.38],'Visible','on');        
        hold on
        fn_hist(draw(:,3))
        fn_hist(draw_C(:,3))
        xlabel('\alpha','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([alpha alpha], YL,'Color','r','LineWidth',3); 
        hold off

        
        ax4 = axes('Position',[0.46 0.10 0.36 0.38],'Visible','on');        
        hold on
        fn_hist(draw(:,4))
        fn_hist(draw_C(:,4))
        xlabel('\beta','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([beta beta], YL,'Color','r','LineWidth',3); 
        hold off

        leg = legend('Uncensored','Thr. = 10\% ','True');
        set(leg,'Interpreter','latex','FontSize',11,'position',[0.85 0.42 0.14 0.2])
        
    else
        set(gcf,'units','normalized','outerposition',[0.1 0.05 0.65 0.95]);
        
        ax1 = axes('Position',[0.05 0.60 0.36 0.38],'Visible','on');        
        fn_hist([draw(:,1),draw_C(:,1)])    
        xlabel('\mu','FontSize',10)
        plotTickLatex2D('FontSize',10);  
        YL = get(gca,'YLim');
        line([c c], YL,'Color','r','LineWidth',3); 
  
        ax2 = axes('Position',[0.46 0.60 0.36 0.38],'Visible','on');        
        fn_hist([draw(:,2),draw_C(:,2)])    
        xlabel('\omega','FontSize',10)
        plotTickLatex2D('FontSize',10);  
        YL = get(gca,'YLim');
        line([omega omega], YL,'Color','r','LineWidth',3); 

        ax3 = axes('Position',[0.05 0.10 0.36 0.38],'Visible','on');        
        fn_hist([draw(:,3),draw_C(:,3)])    
        xlabel('\alpha','FontSize',10)
        plotTickLatex2D('FontSize',10);  
        YL = get(gca,'YLim');
        line([alpha alpha], YL,'Color','r','LineWidth',3); 

        ax4 = axes('Position',[0.46 0.10 0.36 0.38],'Visible','on');        
        fn_hist([draw(:,4),draw_C(:,4)])    
        xlabel('\beta','FontSize',10)
        plotTickLatex2D('FontSize',10);  
        YL = get(gca,'YLim');
        line([beta beta], YL,'Color','r','LineWidth',3); 

        leg = legend('Uncensored','Thr. = 10\% ','True');
        set(leg,'Interpreter','latex','FontSize',10,'position',[0.85 0.42 0.14 0.2])
    end
    
    if save_on
        name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'.eps'];
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



%% PARTIALLY CENSORED
partition = 3; 

II = 100; 
draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
M_short = M/II;

% draw = draw_short;
% M = M/II;
[draw_PC, a_PC] = sim_cond_mit_MH(mit_C, draw_short, partition, M_short, BurnIn, kernel, GamMat);

accept_PC = mean(a_PC);
h_post_PC = volatility_garch11(draw_PC, y, y_S);
y_post_PC = draw_PC(:,1) + sqrt(h_post_PC).*randn(M,1);
y_post_PC = sort(y_post_PC);
VaR_1_post_PC = y_post_PC(p_bar1*M); 
VaR_5_post_PC = y_post_PC(p_bar*M); 

  

if save_on
    save(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'.mat'],...
    'y','draw','draw_C','param_true',...
    'accept','accept_C',...
    'mit','CV','mit_C','CV_C',...
    'VaR_1','VaR_1_post','VaR_1_post_C',...
    'VaR_5','VaR_5_post','VaR_5_post_C','-append')
end


if plot_on
      set(gcf,'units','normalized','outerposition',[0.1 0.05 0.45 0.75]);
        
        ax1 = axes('Position',[0.05 0.60 0.36 0.38],'Visible','on');        
        hold on
        fn_hist(draw(:,1))
        fn_hist(draw_C(:,1))
        fn_hist(draw_PC(:,1))
        xlabel('\mu','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([c c], YL,'Color','r','LineWidth',3); 
        hold off
  
        ax2 = axes('Position',[0.46 0.60 0.36 0.38],'Visible','on');        
        hold on
        fn_hist(draw(:,2))
        fn_hist(draw_C(:,2))
        fn_hist(draw_PC(:,2))
        xlabel('\omega','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([omega omega], YL,'Color','r','LineWidth',3); 
        hold off

        
        ax3 = axes('Position',[0.05 0.10 0.36 0.38],'Visible','on');        
        hold on
        fn_hist(draw(:,3))
        fn_hist(draw_C(:,3))
        fn_hist(draw_PC(:,3))        
        xlabel('\alpha','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([alpha alpha], YL,'Color','r','LineWidth',3); 
        hold off

        
        ax4 = axes('Position',[0.46 0.10 0.36 0.38],'Visible','on');        
        hold on
        fn_hist(draw(:,4))
        fn_hist(draw_C(:,4))
        fn_hist(draw_PC(:,4))        
        xlabel('\beta','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([beta beta], YL,'Color','r','LineWidth',3); 
        hold off

        leg = legend('Uncensored','CP10\%','PCP10\%','True');
        set(leg,'Interpreter','latex','FontSize',11,'position',[0.85 0.42 0.14 0.2])
        
end
%% ADDITIONAL PARAMS
mu_add_init = [0, 1, 1];
% mu_add = zeros(20,3);
% Sigma_add = zeros(20,3*3);
% 
% mu_add_mex = zeros(20,3);
% Sigma_add_mex = zeros(20,3*3);
% 
% for ii = 1:20
% %     tic
% %     kernel_init = @(aa) - C_addparam_posterior_garch11(aa, draw(10*ii,:), y, threshold, y_S)/T;
% %     [mu_add(ii,:),~,~,~,~,S_add] = fminunc(kernel_init, mu_add_init);
% %     Sigma_add(ii,:) = reshape(inv(T*S_add),1,9);
% %     toc
% 
% %     tic
%     kernel_init = @(aa) - C_addparam_posterior_garch11_mex(aa, draw(10*ii,:), y, threshold, y_S)/T;   
%     [mu_add_mex(ii,:),~,~,~,~,S_add] = fminunc(kernel_init, mu_add_init);
%     Sigma_add_mex(ii,:) = reshape(inv(T*S_add),1,9);
% %     toc
% end

mu_add = zeros(M,3);
Sigma_add = zeros(M,3*3);

for ii = 1:M
    if (mod(ii,100)==0)
        fprintf('Add param iter = %d\n',ii)
    end
    kernel_init = @(aa) - C_addparam_posterior_garch11_mex(aa, draw(ii,:), y, threshold, y_S)/T;   
    [mu_add(ii,:),~,~,~,~,S_add] = fminunc(kernel_init, mu_add_init,options);
    Sigma_add(ii,:) = reshape(inv(T*S_add),1,9);
end

if plot_on
    figure(2)
    ff = figure(2);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.4]);
       
    subplot(1,3,1)
    fn_hist(mu_add(:,1))    
    xlabel('Additional a','Interpreter','latex','FontSize',11)
    plotTickLatex2D('FontSize',11);       
    
    subplot(1,3,2)
    fn_hist(mu_add(:,2))    
    xlabel('Additional b','Interpreter','latex','FontSize',11)
    plotTickLatex2D('FontSize',11);  

    subplot(1,3,3)    
    fn_hist(mu_add(:,3))    
    xlabel('Additional c','Interpreter','latex','FontSize',11)
    plotTickLatex2D('FontSize',11);  
    if save_on
        name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_additional_params.eps'];
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

kernel = @(aa, xx) C_addparam_posterior_garch11(aa, xx, y, y_S);



[mu,~,~,~,~,Sigma] = fminunc(kernel_init, mu_init);



%% Threshold = 0
threshold0 = 0;
%% CENSORED
kernel_init = @(xx) - C_posterior_garch11_mex(xx, y, threshold0, y_S)/T;    
kernel = @(xx) C_posterior_garch11_mex(xx, y, threshold0, y_S);
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
accept_C0 = a/(M+BurnIn);
draw_C0 = draw_C0(BurnIn+1:BurnIn+M,:);

h_post_C0 = volatility_garch11(draw_C0,y,y_S);
y_post_C0 = draw_C0(:,1) + sqrt(h_post_C0).*randn(M,1);
VaR_1_post_C0) = y_post_C0(p_bar1*M); 
VaR_5_post_C0 = y_post_C0(p_bar*M); 




if save_on
    save(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_PCP0.mat'],...
    'y','draw','draw_C','draw_PC','draw_C0','draw_PC0','param_true','q1','q5',...
    'accept','accept_C','accept_PC','accept_C0','accept_PC0',...
    'II','mit','CV','mit_C','CV_C','mit_C0','CV_C0',...
    'VaR_1','VaR_1_post','VaR_1_post_C','VaR_1_post_PC','VaR_1_post_C0','VaR_1_post_PC0',...
    'VaR_5','VaR_5_post','VaR_5_post_C','VaR_5_post_PC','VaR_5_post_C0','VaR_5_post_PC0',...
    'MSE_1','MSE_1_post','MSE_1_post_C','MSE_1_post_PC','MSE_1_post_C0','MSE_1_post_PC0',...
    'MSE_5','MSE_5_post','MSE_5_post_C','MSE_5_post_PC','MSE_5_post_C0','MSE_5_post_PC0')
end