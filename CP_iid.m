clear all
close all

addpath(genpath('include/'));

% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s); 

model = 'iid';
parameters = {'$\\mu$','$\\sigma$'};

sigma1 = 1; 
sigma2 = 1; %2;
c = (sigma2 - sigma1)/(sqrt(2*pi));

T = 1000; % time series length
p_bar1 = 0.01;
p_bar = 0.05;
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

plot_on = true;
save_on = false;

%% split normal, mode 0
% if plot_on
%     x = linspace(-5,5,100);
%     fn_norm_pdf = @(xx,s) exp(-xx.^2/(2*s^2))/sqrt(2*pi*s^2);
%     A = sqrt(2/pi)/(sigma1+sigma2);
%     split1 = fn_norm_pdf(x(1:50),sigma2);
%     split2= fn_norm_pdf(x(51:100),sigma1);
%     ff = figure(1)
%     plot(x(1:50),split1,'linewidth',2)
%     hold on
%     plot(x(51:100),split2,'linewidth',2)
%     scatter(-c,0,'MarkerFaceColor','r','MarkerEdgeColor','r')
%     text(-0.75,0.015,'E[X]','FontSize',12,'Interpreter','latex')
%     hold off  
% 	plotTickLatex2D('FontSize',12);
%     
%     if save_on
%         name = [figures_path,model,'_hor_direct_H', num2str(H),'.png'];
%         name = ['figures/iid/Idd_SplitNormal_poster.eps'];
%         set(gcf,'PaperPositionMode','auto');
%         print_fail = 1;
%         while print_fail 
%             try 
%                 print(name,'-dpng','-r0')
%                 print(name,'-depsc','-r0')                    
%                 print(ff,name,'-depsc','-r0')
%                 print_fail = 0;
%             catch
%                 print_fail = 1;
%             end
%         end
%     end
% 
% end

%% iid simulation, mean 0
eps = randn(T,1);
ind = (eps>0);
eps(ind) = c + sigma1.*eps(ind);
eps(~ind) = c + sigma2.*eps(~ind);
% eps1 = c + sigma1.*eps(eps>0);
% eps2 = c + sigma2.*eps(eps<0);
% eps = [eps1;eps2];
y = eps; % -mean(eps);
median(y) %0.3716
mean(y)  % -0.0134

% true VaRs
q1 = norminv(p_bar1,c,sigma2);
q5 = norminv(p_bar,c,sigma2); 

eps_sort = randn(M,1);
ind = (eps_sort>0);
eps_sort(ind) = c + sigma1.*eps_sort(ind);
eps_sort(~ind) = c + sigma2.*eps_sort(~ind);

y_sort = sort(eps_sort);
VaR_1 = y_sort(p_bar1*M); 
VaR_5 = y_sort(p_bar*M);
threshold = y_sort(2*p_bar*M);

%% Uncensored Posterior
% Misspecified model: normal with unknown mu and sigma
% Metropolis-Hastings for the parameters

% Uncensored likelihood
kernel_init = @(xx) -loglik_iid(xx,y);
[mu,~,~,~,~,Sigma] = fminunc(kernel_init,[0,1]);
Sigma = inv(T*Sigma);
df = 5;
draw = rmvt(mu,Sigma,df,M+BurnIn);
kernel = @(ss) posterior_iid(ss,y);
lnk = kernel(draw);

lnd = dmvgt_mex(draw, mu, Sigma, df, 1, GamMat, double(1));
lnw = lnk - lnd;
lnw = lnw - max(lnw);
[ind, a] = fn_MH(lnw);
draw = draw(ind,:);
accept = a/(M+BurnIn);
draw = draw(BurnIn+1:BurnIn+M,:);    

% if plot_on
%     subplot(1,2,1)
%     fn_hist(draw(:,1))    
%     xlabel('\mu')
%     subplot(1,2,2)
%     fn_hist(draw(:,2))    
%     xlabel('\sigma')
% end

y_post = draw(:,1) + draw(:,2).*randn(M,1);
y_post = sort(y_post);
VaR_1_post = y_post(p_bar1*M); 
VaR_5_post = y_post(p_bar*M); 

%% Censored posterior: take values below the threshold
% Misspecified model: N(mu,sigma)

% 1. Threshold = 10% perscentile of the data sample
kernel_init = @(xx) - C_posterior_iid(xx, y, threshold)/T;
[mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu);
Sigma_C = inv(T*Sigma_C);
draw_C = rmvt(mu_C,Sigma_C,df,M+BurnIn);
kernel = @(ss) C_posterior_iid(ss, y, threshold);
lnk_C = kernel(draw_C);
lnd_C = dmvgt_mex(draw_C, mu_C, Sigma_C, df, 1, GamMat, double(1));
lnw_C = lnk_C - lnd_C;
lnw_C = lnw_C - max(lnw_C);
[ind, a] = fn_MH(lnw_C);
draw_C = draw_C(ind,:);
accept_C = a/(M+BurnIn);
draw_C = draw_C(BurnIn+1:BurnIn+M,:);

y_post_C = draw_C(:,1) + draw_C(:,2).*randn(M,1);
y_post_C = sort(y_post_C);
VaR_1_post_C = y_post_C(p_bar1*M); 
VaR_5_post_C = y_post_C(p_bar*M); 

% if plot_on
%     subplot(1,2,1)
%     fn_hist(draw_C(:,1))    
%     xlabel('\mu')
%     subplot(1,2,2)
%     fn_hist(draw_C(:,2))    
%     xlabel('\sigma')
% end


% 2. Threshold = 0             
threshold = 0;
kernel_init = @(xx) - C_posterior_iid(xx, y, threshold)/T;
[mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu);
Sigma_C0 = inv(T*Sigma_C0);
draw_C0 = rmvt(mu_C0,Sigma_C0,df,M+BurnIn);
kernel = @(ss) C_posterior_iid(ss, y, threshold);
lnk_C0 = kernel(draw_C0);
lnd_C0 = dmvgt_mex(draw_C0, mu_C0, Sigma_C0, df, 1, GamMat, double(1));
lnw_C0 = lnk_C0 - lnd_C0;
lnw_C0 = lnw_C0 - max(lnw_C0);
[ind, a] = fn_MH(lnw_C0);
draw_C0 = draw_C0(ind,:);
accept_C0 = a/(M+BurnIn);
draw_C0 = draw_C0(BurnIn+1:BurnIn+M,:);

y_post_C0 = draw_C0(:,1) + draw_C0(:,2).*randn(M,1);
y_post_C0 = sort(y_post_C0);
VaR_1_post_C0 = y_post_C0(p_bar1*M); 
VaR_5_post_C0 = y_post_C0(p_bar*M); 

% if plot_on
%     subplot(1,2,1)
%     hold on
%     fn_hist([draw(:,1),draw_C(:,1),draw_C0(:,1)])    
%     YL = get(gca,'YLim');
%     line([c c], YL,'Color','r','LineWidth',3); 
%     hold off
%     xlabel('\mu')
%     
%     subplot(1,2,2)
%     hold on
%     fn_hist([draw(:,2),draw_C(:,2),draw_C0(:,2)])    
%     YL = get(gca,'YLim');    
%     line([sigma2 sigma2], YL,'Color','r','LineWidth',3); 
%     hold off
%     xlabel('\sigma')
% end

%% save and plot
param_true = [c,sigma2];
if save_on
    save(['results/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'.mat'],...
    'y','draw','draw_C','draw_C0','param_true',...
    'accept','accept_C','accept_C0',...
    'VaR_1','VaR_1_post','VaR_1_post_C','VaR_1_post_C0',...
    'VaR_5','VaR_5_post','VaR_5_post_C','VaR_5_post_C0')
end

if plot_on
    ff = figure(11);
%     set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.4]);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.6]);
%     ax1 = axes('Position',[0.06 0.15 0.32 0.8],'Visible','on');
    ax1 = axes('Position',[0.03 0.15 0.36 0.8],'Visible','on');
    axes(ax1)
    [f,xi] = ksdensity(draw(:,1));
    [f_C,xi_C] = ksdensity(draw_C(:,1));
    [f_C0,xi_C0] = ksdensity(draw_C0(:,1));
    hold on
    plot(xi,f,'k','linewidth',2)
    plot(xi_C,f_C,'color',[         0    0.4470    0.7410],'linewidth',2)
    plot(xi_C0,f_C0,'color',[    0.3010    0.7450    0.9330],'linewidth',2)
    YL = get(gca,'YLim');
    line([c c], YL,'Color','r','LineWidth',2); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\mu','FontSize',12)

    ax2 = axes('Position',[0.43 0.15 0.36 0.8],'Visible','on');
    axes(ax2)
    hold on
    [f,xi] = ksdensity(draw(:,2));
    [f_C,xi_C] = ksdensity(draw_C(:,2));
    [f_C0,xi_C0] = ksdensity(draw_C0(:,2));
    hold on
    plot(xi,f,'k','linewidth',2)
    plot(xi_C,f_C,'color',[         0    0.4470    0.7410],'linewidth',2)
    plot(xi_C0,f_C0,'color',[    0.3010    0.7450    0.9330],'linewidth',2)
    YL = get(gca,'YLim');
    line([sigma2 sigma2], YL,'Color','r','LineWidth',3); 
    hold off
    plotTickLatex2D('FontSize',12);
    xlabel('\sigma','FontSize',12)
    leg = legend('Uncensored','Thr. = 10\%','Thr. = 0','True');%,)
    set(leg,'Interpreter','latex','FontSize',12,'position',[0.8 0.42 0.18 0.16])

    if save_on
        name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'_dens.eps'];
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