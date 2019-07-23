clear all
close all

addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

model = 'ar1';
parameters = {'$\\mu$','$\\sigma$','$\\phi$'};

sigma1 = 1;
sigma2 = 2;
c = (sigma2 - sigma1)/sqrt(2*pi);

T = 1000; % time series length
p_bar1 = 0.01;
p_bar = 0.05;
% Metropolis-Hastings for the parameters
M = 10000; % number of draws 
BurnIn = 1000;

x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);

v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014a)')
    fn_hist = @(xx) hist(xx,20);
else
    fn_hist = @(xx) histogram(xx,20);
end

plot_on = false;
save_on = true;

%% simple AR(1)
eps = randn(T,1);
ind = (eps>0);
eps(ind) = c + sigma1.*eps(ind);
eps(~ind) = c + sigma2.*eps(~ind);

rho = 0.8;
param_true = [c,sigma2,rho];
y = zeros(T,1);

y(1,1) = eps(1,1);
for ii = 2:T
    y(ii,1) = rho*y(ii-1,1) + eps(ii,1);
end

% MC VaRs under the true model
eps_sort = randn(M,1);
ind = (eps_sort>0);
eps_sort(ind) = c + sigma1.*eps_sort(ind);
eps_sort(~ind) = c + sigma2.*eps_sort(~ind);

y_sort = rho*y(T,1) + eps_sort;
y_sort = sort(y_sort);
VaR_1 = y_sort(p_bar1*M); 
VaR_5 = y_sort(p_bar*M); 

%% Misspecified model: AR1 normal with unknown mu and sigma
%% Uncensored Posterior

% kernel_init = @(xx) -loglik_ar1(xx,y);
% [mu,~,~,~,~,Sigma] = fminunc(kernel_init,[0,1,0.9]);
% Sigma = inv(T*Sigma);
% df = 5;
% draw = rmvt(mu,Sigma,df,M+BurnIn);

mu_init = [0,1,0.9];
kernel_init = @(xx) -posterior_ar1_mex(xx,y)/T;
% kernel = @(xx) posterior_ar1(xx,y);
kernel = @(xx) posterior_ar1_mex(xx,y);
% tic
% lnk = kernel(draw);
% toc    
[mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
% tic
[draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
% toc
% Elapsed time is 1.178022 seconds.
% Elapsed time is 0.124007 seconds.

% lnd = dmvgt_mex(draw, mu, Sigma, df, 1, GamMat, double(1));
lnd = dmvgt(draw, mit, true, GamMat);   
lnw = lnk - lnd;
lnw = lnw - max(lnw);
[ind, a] = fn_MH(lnw);
draw = draw(ind,:);
lnw = lnw(ind);
accept = a/(M+BurnIn);
draw = draw(BurnIn+1:BurnIn+M,:);    
lnw = lnw(BurnIn+1:BurnIn+M,:);    
 
y_post = draw(:,1) + draw(:,3).*y(T,1) + draw(:,2).*randn(M,1);
y_post = sort(y_post);
VaR_1_post = y_post(p_bar1*M); 
VaR_5_post = y_post(p_bar*M); 

%% PARTIAL CENSORING: keep rho uncensored, then censor mu i sigma
partition = 3;

%% THRESHOLD = 10% perscentile of the data sample
threshold = sort(y);
threshold = threshold(2*p_bar*T);
% conditional candidate 
% joint cnadidate for the joint censored posterior
% kernel_init = @(xx) - C_posterior_ar1(xx, y, threshold);
% kernel = @(xx) C_posterior_ar1(xx, y, threshold);
kernel_init = @(xx) - C_posterior_ar1_mex(xx, y, threshold);
kernel = @(xx) C_posterior_ar1_mex(xx, y, threshold);
mu_init = [0, 1, 0.5];
cont = MitISEM_Control;
cont.mit.iter_max = 10;
[mit_C, CV_C] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);

if save_on
    save(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP.mat'],...
    'y','draw','accept','param_true',...
    'mit','CV','mit_C','CV_C','cont',...
    'VaR_1','VaR_1_post','VaR_5','VaR_5_post');
end

% load(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP.mat'])
cont1 = cont;
cont1.mit.iter_max = 0;
[mit_C1, CV_C1] = MitISEM_new(kernel_init, kernel, mu_init, cont1, GamMat);

if save_on
    save(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP.mat'],...
    'mit_C1','CV_C1','cont1','-append')
end

%% True rho
[draw_partial_true, accept_partial_true] = sim_cond_mit_MH(mit_C, [0,0,0.8], partition, M, BurnIn, kernel, GamMat);
[draw_partial_true1, accept_partial_true1] = sim_cond_mit_MH(mit_C1, [0,0,0.8], partition, M, BurnIn, kernel, GamMat);

% Plot partial true with 1 and 2 comp mixtures
if plot_on
    ff = figure(11);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.4]);
    ax1 = axes('Position',[0.06 0.15 0.32 0.8],'Visible','on');
    axes(ax1)
    hold on
    if strcmp(v_new,'(R2014a)')
        hist([draw(:,1), draw_partial_true(:,1),draw_partial_true1(:,1)],20)
    else
        fn_hist(draw(:,1))
        fn_hist(draw_partial_true(:,1))
        fn_hist(draw_partial_true1(:,1))
    end
    YL = get(gca,'YLim');
    line([c c], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\mu','FontSize',12)

    ax2 = axes('Position',[0.44 0.15 0.32 0.8],'Visible','on');
    axes(ax2)
    hold on
    if strcmp(v_new,'(R2014a)')
        hist([draw(:,2), draw_partial_true(:,2),draw_partial_true1(:,2)],20)
    else
        fn_hist(draw(:,2))
        fn_hist(draw_partial_true(:,2))
        fn_hist(draw_partial_true1(:,2))
    end
    
    YL = get(gca,'YLim');
    line([sigma2 sigma2], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\sigma','FontSize',12)
    leg = legend('Uncensored','Partial true, 2 comp.','Partial true, 1 comp.','True');%,)
    set(leg,'Interpreter','latex','FontSize',11,'position',[0.78 0.42 0.18 0.2])

    if save_on
        name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP_true.eps'];
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

if save_on
    save(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP.mat'],...
    'draw_partial_true','draw_partial_true1','accept_partial_true','accept_partial_true1','-append')
end

%% Short version
II = 100;
% draw_org = draw;
% M_org = M;
% draw_short = draw(1:II,:);
draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
M_short = M/II;

% draw = draw_short;
% M = M/II;
[draw_partial_short, accept_partial_short] = sim_cond_mit_MH(mit_C, draw_short, partition, M_short, BurnIn, kernel, GamMat);
% draw_partial_short = draw_partial;
% accept_partial_short = accept;

if save_on
    save(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP.mat'],...
    'draw_partial_short','accept_partial_short','-append')
end

% Plot partial true and short
if plot_on
    ff = figure(12);
    if strcmp(v_new,'(R2014a)')
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.7 0.55]);
        ax1 = axes('Position',[0.05 0.15 0.23 0.8],'Visible','on');
    else
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.4]);
        ax1 = axes('Position',[0.04 0.15 0.24 0.8],'Visible','on');        
    end
    axes(ax1)
    hold on
    if strcmp(v_new,'(R2014a)')
        fn_hist([draw(:,1), draw_partial_true(:,1),draw_partial_short(:,1)])
    else
        fn_hist(draw(:,1))    
        fn_hist(draw_partial_true(:,1)) 
        fn_hist(draw_partial_short(:,1)) 
%         hist([draw(:,1), draw_C(:,1),draw_C0(:,1)],20)
    end
    YL = get(gca,'YLim');
    line([c c], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\mu','FontSize',12)

    if strcmp(v_new,'(R2014a)')
        ax2 = axes('Position',[0.33 0.15 0.23 0.8],'Visible','on');
    else    
        ax2 = axes('Position',[0.32 0.15 0.24 0.8],'Visible','on');
    end
    axes(ax2)
    hold on
    if strcmp(v_new,'(R2014a)')
        fn_hist([draw(:,2), draw_partial_true(:,2),draw_partial_short(:,2)])     
    else
        fn_hist(draw(:,2))    
        fn_hist(draw_partial_true(:,2))  
        fn_hist(draw_partial_short(:,2))    
    end
    YL = get(gca,'YLim');
    line([sigma2 sigma2], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\sigma','FontSize',12)

    if strcmp(v_new,'(R2014a)')   
        ax3 = axes('Position',[0.61 0.15 0.23 0.8],'Visible','on');
    else
        ax3 = axes('Position',[0.61 0.15 0.24 0.8],'Visible','on');
    end
    axes(ax3)
    hold on
    if strcmp(v_new,'(R2014a)')   
        fn_hist([draw(:,3),[],draw_partial_short(:,3)])
    else        
        fn_hist(draw(:,3))    
        fn_hist(draw_partial_true(1,3))  
        fn_hist(draw_partial_short(:,3))    
    end
    YL = get(gca,'YLim');
    line([rho rho], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\rho','FontSize',12)
    
    leg = legend('Uncensored','Partial true, 2 comp.','Partial 100 rhos.','True');%,)
    if strcmp(v_new,'(R2014a)')   
        set(leg,'Interpreter','latex','FontSize',10,'position',[0.85 0.42 0.14 0.2])
    else
        set(leg,'Interpreter','latex','FontSize',11,'position',[0.85 0.42 0.14 0.2])
    end
    
    if save_on
        name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP_short.eps'];
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

%% Highest weights rhos
[~,ind_max] = sort(lnw);
draw_max = draw(ind_max,:);
draw_max = draw_max(M-II+1:M,:);

[draw_partial_max, accept_partial_max] = sim_cond_mit_MH(mit_C, draw_max, partition, M_short, BurnIn, kernel, GamMat);

% Plot partial short and max
if plot_on
    ff = figure(13);
    if strcmp(v_new,'(R2014a)')
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.7 0.55]);
        ax1 = axes('Position',[0.05 0.15 0.23 0.8],'Visible','on');
    else
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.4]);
        ax1 = axes('Position',[0.04 0.15 0.24 0.8],'Visible','on');        
    end
    axes(ax1)
    hold on
    if strcmp(v_new,'(R2014a)')
        fn_hist([draw(:,1),draw_partial_short(:,1),draw_partial_max(:,1)])
    else
        fn_hist(draw(:,1))    
%         fn_hist(draw_partial_true(:,1)) 
        fn_hist(draw_partial_short(:,1)) 
        fn_hist(draw_partial_max(:,1))            
    end
    YL = get(gca,'YLim');
    line([c c], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\mu','FontSize',12)

    if strcmp(v_new,'(R2014a)')
        ax2 = axes('Position',[0.33 0.15 0.23 0.8],'Visible','on');
    else    
        ax2 = axes('Position',[0.32 0.15 0.24 0.8],'Visible','on');
    end
    axes(ax2)
    hold on
    if strcmp(v_new,'(R2014a)')
        fn_hist([draw(:,2),draw_partial_short(:,2),draw_partial_max(:,2)])     
    else
        fn_hist(draw(:,2))    
%         fn_hist(draw_partial_true(:,2))  
        fn_hist(draw_partial_short(:,2))   
        fn_hist(draw_partial_max(:,2))                    
    end
    YL = get(gca,'YLim');
    line([sigma2 sigma2], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\sigma','FontSize',12)

    if strcmp(v_new,'(R2014a)')   
        ax3 = axes('Position',[0.61 0.15 0.23 0.8],'Visible','on');
    else
        ax3 = axes('Position',[0.61 0.15 0.24 0.8],'Visible','on');
    end
    axes(ax3)
    hold on
    if strcmp(v_new,'(R2014a)')   
        fn_hist([draw(:,3),draw_partial_short(:,3),draw_partial_max(:,3)])
    else        
        fn_hist(draw(:,3))    
%         fn_hist(draw_partial_true(1,3))  
        fn_hist(draw_partial_short(:,3))    
        fn_hist(draw_partial_max(:,3))            
    end
    YL = get(gca,'YLim');
    line([rho rho], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\rho','FontSize',12)
    
    leg = legend('Uncensored','Partial 100 rhos','Partial 100 max','True');%,)
    if strcmp(v_new,'(R2014a)')   
        set(leg,'Interpreter','latex','FontSize',10,'position',[0.85 0.42 0.14 0.2])
    else
        set(leg,'Interpreter','latex','FontSize',11,'position',[0.85 0.42 0.14 0.2])
    end
    
    if save_on
        name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP_short.eps'];
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

 %% THRESHOLD = 0
threshold0 = 0;
% conditional candidate 
% joint candidate for the joint censored posterior
% kernel_init = @(xx) - C_posterior_ar1(xx, y, threshold0);
% kernel = @(xx) C_posterior_ar1(xx, y, threshold0);
kernel_init = @(xx) - C_posterior_ar1_mex(xx, y, threshold0);
kernel = @(xx) C_posterior_ar1_mex(xx, y, threshold0);
[mit_C0, CV_C0] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);

% True rho
[draw_partial0_true, accept_partial0_true] = sim_cond_mit_MH(mit_C0, [0,0,0.8], partition, M, BurnIn, kernel, GamMat);
 
% Short version
[draw_partial0_short, accept_partial0_short] = sim_cond_mit_MH(mit_C0, draw_short, partition, M_short, BurnIn, kernel, GamMat);

% Highest weights rhos
[draw_partial0_max, accept_partial0_max] = sim_cond_mit_MH(mit_C0, draw_max, partition, M_short, BurnIn, kernel, GamMat);

 if save_on
    save(['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP.mat'],...
    'mit_C0','CV_C0',...
    'draw_partial0_true','accept_partial0_true',...
    'draw_partial0_short','accept_partial0_short',...
    'draw_partial0_max','accept_partial0_max',...
    '-append')
 end 
 
% Plot partial true both thresholds 
if plot_on
    ff = figure(100);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.4]);
    ax1 = axes('Position',[0.06 0.15 0.32 0.8],'Visible','on');
    axes(ax1)
    hold on
    if strcmp(v_new,'(R2014a)')
        hist([draw(:,1), draw_partial_true(:,1),draw_partial0_true(:,1)],20)
    else
        fn_hist(draw(:,1))
        fn_hist(draw_partial_true(:,1))
        fn_hist(draw_partial0_true(:,1))
    end
    YL = get(gca,'YLim');
    line([c c], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\mu','FontSize',12)

    ax2 = axes('Position',[0.44 0.15 0.32 0.8],'Visible','on');
    axes(ax2)
    hold on
    if strcmp(v_new,'(R2014a)')
        hist([draw(:,2), draw_partial_true(:,2),draw_partial0_true(:,2)],20)
    else
        fn_hist(draw(:,2))
        fn_hist(draw_partial_true(:,2))
        fn_hist(draw_partial0_true(:,2))
    end
    
    YL = get(gca,'YLim');
    line([sigma2 sigma2], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\sigma','FontSize',12)
    leg = legend('Uncensored','Partial true, thr. 10\%.','Partial true, thr. 0','True');%,)
    set(leg,'Interpreter','latex','FontSize',11,'position',[0.78 0.42 0.18 0.2])

    if save_on
        name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP_thresholds_true.eps'];
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

% Plot partial short both thresholds 
if plot_on
    ff = figure(110);
    if strcmp(v_new,'(R2014a)')
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.7 0.55]);
        ax1 = axes('Position',[0.05 0.15 0.23 0.8],'Visible','on');
    else
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.4]);
        ax1 = axes('Position',[0.04 0.15 0.24 0.8],'Visible','on');        
    end
    axes(ax1)
    hold on
    if strcmp(v_new,'(R2014a)')
        fn_hist([draw(:,1), draw_partial_short(:,1), draw_partial0_short(:,1)])
    else
        fn_hist(draw(:,1))    
        fn_hist(draw_partial_short(:,1)) 
        fn_hist(draw_partial0_short(:,1)) 
    end
    YL = get(gca,'YLim');
    line([c c], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\mu','FontSize',12)

    if strcmp(v_new,'(R2014a)')
        ax2 = axes('Position',[0.33 0.15 0.23 0.8],'Visible','on');
    else    
        ax2 = axes('Position',[0.32 0.15 0.24 0.8],'Visible','on');
    end
    axes(ax2)
    hold on
    if strcmp(v_new,'(R2014a)')
        fn_hist([draw(:,2), draw_partial_short(:,2),draw_partial0_short(:,2)])     
    else
        fn_hist(draw(:,2))    
        fn_hist(draw_partial_short(:,2))  
        fn_hist(draw_partial0_short(:,2))  
    end
    YL = get(gca,'YLim');
    line([sigma2 sigma2], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\sigma','FontSize',12)

    if strcmp(v_new,'(R2014a)')   
        ax3 = axes('Position',[0.60 0.15 0.23 0.8],'Visible','on');
    else
        ax3 = axes('Position',[0.60 0.15 0.24 0.8],'Visible','on');
    end
    axes(ax3)
    hold on
    if strcmp(v_new,'(R2014a)')   
        fn_hist([draw(:,3),draw_partial_short(:,3),draw_partial0_short(:,3)])
    else        
        fn_hist(draw(:,3))    
        fn_hist(draw_partial_short(1,3))  
        fn_hist(draw_partial0_short(:,3))    
    end
    YL = get(gca,'YLim');
    line([rho rho], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\rho','FontSize',12)
    
    leg = legend('Uncensored','Partial short thr.=10\%','Partial short thr.=0','True');
    if strcmp(v_new,'(R2014a)')   
        set(leg,'Interpreter','latex','FontSize',10,'position',[0.84 0.42 0.14 0.2])
    else
        set(leg,'Interpreter','latex','FontSize',10,'position',[0.84 0.42 0.14 0.2])
    end
    
    if save_on
        name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_PCP_thresholds.eps'];
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