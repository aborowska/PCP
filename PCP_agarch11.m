clear all
close all
%%%%% https://www.kevinsheppard.com/images/9/95/MFE_Toolbox_Documentation.pdf
% % addpath(genpath('../MFEToolbox/'));
% % mex composite_likelihood.c
% % mex armaxerrors.c
% % mex agarch_core.c
% % mex egarch_core.c
% % mex igarch_core.c
% % mex tarch_core.c

% startingvals =[0.0591 0.1040 0.8393];

% % parameters = agarch(y-mean(y),1,1)
% % [PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS] = agarch(y-mean(y),1,1)
% omega alpha (-gamma) beta
%     0.0575  0.1025  0.0639  0.8418
% omega2 = omega/(1-alpha-beta) 
% omega2 = PARAMETERS(1)/(1-PARAMETERS(2)-PARAMETERS(4))
%   1.0319

% mu:
%     -0.0739    1.0318    0.1027    0.8414
%  omega2 = omega*(1-alpha-beta) 


addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

model = 'agarch11'; 
partition = 3;
% parameters = {'$\\mu$','$\\gamma$','$\\omega$','$\\alpha$','$\\beta$'};
parameters = {'$\\mu$','$\\omega$','$\\mu2$','$\\alpha$','$\\beta$'};

sigma1 = 1;
sigma2 = 2;
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


%% GARCH(1,1)
eps = randn(T,1);
ind = (eps>0);
eps(ind) = c + sigma1.*eps(ind);
eps(~ind) = c + sigma2.*eps(~ind);
eps = eps/sqrt(kappa);  
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
kernel_init = @(xx) -posterior_agarch11_mex(xx, y, y_S)/T;
kernel = @(xx) posterior_agarch11_mex(xx, y, y_S);
[mu,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init,options);
Sigma = inv(T*Sigma);
% sigma2 = 1 ==>  -0.0188   -0.0739    1.0318    0.1027    0.8414
% sigma2 = 2 ==>  -0.0121   -0.1388    1.0135    0.0887    0.8364

kernel_init = @(xx) - posterior_agarch11_mex(transform_param_agarch(xx, 'back'), y, y_S)/T;    
[mu2,~,~,~,~,Sigma2] = fminunc(kernel_init,transform_param_agarch(mu_init, 'opt'),options);
mu2 = transform_param_agarch(mu2, 'back');


kernel_init = @(xx) -posterior_agarch11_mex(xx, y, 0)/T;
kernel = @(xx) posterior_agarch11_mex(xx, y, 0);
[mu,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init,options);
Sigma = inv(T*Sigma);
% sigma2 = 1 ==> -0.0186   -0.0738    1.0159    0.1015    0.8418
% sigma2 = 2 ==> -0.0118   -0.1384    1.0032    0.0878    0.8366


try
    [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
    [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
    lnd = dmvgt(draw, mit, true, GamMat); 
catch
    [mu,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init,options);
    %    -0.0188   -0.0739    1.0318    0.1027    0.8414
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

h_post = volatility_agarch11(draw,y,y_S);
y_post = draw(:,1) + sqrt(h_post).*randn(M,1);
y_post = sort(y_post);
VaR_1_post = y_post(p_bar1*M); 
VaR_5_post = y_post(p_bar*M); 

%% Threshold = 10% perscentile of the data sample
threshold = sort(y);
threshold = threshold(2*p_bar*T);
%% CENSORED
kernel_init = @(xx) - C_posterior_agarch11_mex(xx, y, threshold, y_S)/T;    
kernel = @(xx) C_posterior_agarch11_mex(xx, y, threshold, y_S);
[mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
Sigma_C = inv(T*Sigma_C);
% sigma2 = 1 ==> -0.0779   -0.1001    0.9484    0.1321    0.6173
% sigma2 = 2 ==>  0.1767    0.1498    2.2550    0.2103    0.5550


kernel_init = @(xx) - C_posterior_agarch11_mex(xx, y, threshold, 0)/T;    
kernel = @(xx) C_posterior_agarch11_mex(xx, y, threshold, 0);
[mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
Sigma_C = inv(T*Sigma_C);
% unconditional variance: 
%h1 = omega + (alpha.*gama.^2)./(1-alpha-beta); 
h1 = mu_init(1,3) + (mu_init(1,4).*mu_init(1,2).^2)./(1-mu_init(1,4)-mu_init(1,5)); % 1.005
% sigma2 = 1 ==> -0.0887   -0.1077    0.9261    0.1295    0.6214
% sigma2 = 2 ==>  0.1655    0.1305    2.1720    0.2104    0.5371


% WITH TRANSFORMATION

kernel_init = @(xx) - C_posterior_agarch11_mex(transform_param_agarch(xx, 'back'), y, threshold, y_S)/T;    
[mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,transform_param_agarch(mu_init, 'opt'),options);
mu_C = transform_param_agarch(mu_C, 'back');
% sigma2 = 1 ==> -0.0778   -0.1001    0.9484    0.1321    0.6173
% sigma2 = 2 ==>  0.1767    0.1498    2.2550    0.2103    0.5550


kernel_init = @(xx) - C_posterior_agarch11_mex(transform_param_agarch(xx, 'back'), y, threshold, 0)/T;    
[mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,transform_param_agarch(mu_init, 'opt'),options);
mu_C = transform_param_agarch(mu_C, 'back');
% sigma2 = 1 ==> -0.0886   -0.1077    0.9261    0.1295    0.6214
% sigma2 = 2 ==>  0.1655    0.1305    2.1720    0.2104    0.5371


try
    mu_init_C = mu_init;
    mu_init_C(1,1) = 0;
    mu_init_C(1,2) = 0;
    [mit_C, CV_C] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat); 
    [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
    lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
catch
    try
        [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
        Sigma_C = inv(T*Sigma_C);
        draw_C = rmvt(mu_C,Sigma_C,df,M+BurnIn);
        mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
%         [mit_C, CV_C] = MitISEM_new(mit_C, kernel, mu_init, cont, GamMat);   
        cont.mit. CV_max = 1.9;
        [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
        [draw_C, lnk_C] = fn_rmvgt_robust(M+BurnIn, mit_C, kernel, false);
        lnd_C = dmvgt(draw_C, mit_C, true, GamMat);    
    catch
        mu_C = fminunc(kernel_init,mu_init,options);
        mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma,1,length(mu_C)^2),'df', df, 'p', 1);
%         [mit_C, CV_C] = MitISEM_new(mit_C, kernel, mu_init, cont, GamMat);   
        cont.mit. CV_max = 1.9;
        [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);
        if CV_C(end)>2
%             [mit_C, CV_C] = MitISEM_new(mit_C, kernel, mu_init, cont, GamMat);   
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

h_post_C = volatility_agarch11(draw_C, y, y_S);
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

%% PARTIALLY CENSORED
% partition = 3; 

II = 10; 
draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
M_short = M/II;

% draw = draw_short;
% M = M/II;
[draw_PC, a_PC] = sim_cond_mit_MH(mit_C, draw_short, partition, M_short, BurnIn, kernel, GamMat);

accept_PC = mean(a_PC);
h_post_PC = volatility_agarch11(draw_PC, y, y_S);
y_post_PC = draw_PC(:,1) + sqrt(h_post_PC).*randn(M,1);
y_post_PC = sort(y_post_PC);
VaR_1_post_PC = y_post_PC(p_bar1*M); 
VaR_5_post_PC = y_post_PC(p_bar*M); 
  

if save_on
    save(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'.mat'],...
    'draw_PC','II','accept_PC','VaR_1_post_PC','VaR_5_post_PC','-append')
end

if plot_on
    ff = figure(10);
    if strcmp(v_new,'(R2014a)')
        set(gcf,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
      
        subplot(2,3,1)
%         ax1 = axes('Position',[0.05 0.60 0.36 0.38],'Visible','on');        
        hold on
        fn_hist(draw(:,1))
        fn_hist(draw_C(:,1))
%         fn_hist(draw_PC(:,1))     
        xlabel('\mu','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([c c], YL,'Color','r','LineWidth',3); 
        hold off

        subplot(2,3,2)
%         ax2 = axes('Position',[0.46 0.60 0.36 0.38],'Visible','on');        
        hold on
        fn_hist(draw(:,2))
        fn_hist(draw_C(:,2))
%         fn_hist(draw_PC(:,2))        
        xlabel('\mu_2','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([mu2  mu2], YL,'Color','r','LineWidth',3); 
        hold off

 
        subplot(2,3,3)
%         ax3 = axes('Position',[0.05 0.10 0.36 0.38],'Visible','on');        
        hold on
        fn_hist(draw(:,3))
        fn_hist(draw_C(:,3))
%         fn_hist(draw_PC(:,3))        
        xlabel('\omega','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([omega omega], YL,'Color','r','LineWidth',3); 
        hold off

        subplot(2,3,4)
%         ax4 = axes('Position',[0.46 0.10 0.36 0.38],'Visible','on');                
        hold on
        fn_hist(draw(:,4))
        fn_hist(draw_C(:,4))
%         fn_hist(draw_PC(:,4))        
        xlabel('\alpha','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([alpha alpha], YL,'Color','r','LineWidth',3); 
        hold off  
              
        subplot(2,3,5)
        hold on
        fn_hist(draw(:,5))
        fn_hist(draw_C(:,5))
%         fn_hist(draw_PC(:,5))        
        xlabel('\beta','FontSize',11)
        plotTickLatex2D('FontSize',11);  
        YL = get(gca,'YLim');
        line([beta beta], YL,'Color','r','LineWidth',3); 
        hold off

        leg = legend('Uncens.','Cens. 10\%','Part. Cens. 10\%','True');
        set(leg,'Interpreter','latex','FontSize',11,'position',[0.7 0.22 0.14 0.2])
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