addpath(genpath('include/'));
clear all
close all


model = 'iid';
% model = 'ar1';
T = 10000;
sigma1 = 1;
% sigma2 = 2;
sigma2 = 1;
c = (sigma2 - sigma1)/(sqrt(2*pi));%1/sqrt(2*pi); %0.3989
 
p_bar05 = 0.005;
p_bar1 = 0.01;
p_bar = 0.05;    

% true VaRs
q05 = norminv(p_bar05,c,sigma2);
q1 = norminv(p_bar1,c,sigma2);
q5 = norminv(p_bar,c,sigma2); 


if strcmp(model,'ar1')
    if (sigma2 == 2)
        name = ['results/',model,'/',model,'_1_',num2str(sigma2),'_T',num2str(T),...
            '_H100_II10_PCP0_MC_(R2017a)_varc.mat'];
        load(name, '-regexp','^draw','^param')
        name = ['results/',model,'/',model,'_1_',num2str(sigma2),'_T',num2str(T),...
            '_H100_II10_PCP0_MC_(R2017a).mat'];
        load(name, '-regexp','^draw','^param')
    else
        name = ['results/',model,'/',model,'_1_',num2str(sigma2),'_T',num2str(T),...
            '_H100_II10_PCP0_MC_(R2015b)_varc.mat'];
        load(name, '-regexp','^draw','^param')
        name = ['results/',model,'/',model,'_1_',num2str(sigma2),'_T',num2str(T),...
            '_H100_II10_PCP0_MC_(R2015b).mat'];
        load(name, '-regexp','^draw','^param')
    end
elseif strcmp(model,'iid')
    name = ['results/',model,'/',model,'_1_',num2str(sigma2),'_T',num2str(T),...
        '_MC.mat'];
    load(name, '-regexp','^draw','^param')
%     name = ['results/',model,'/',model,'_1_',num2str(sigma2),'_T',num2str(T),...
%         '_MC.mat'];
%     load(name, '-regexp','^draw','^param')    
end
    
lw = 1;

ColPalette = [     0    0.4470    0.7410    % deep blue
              0.3010    0.7450    0.9330    % light blue
              0.0000    0.5000    0.0000    % deep green        
              0.4660    0.6740    0.1880];  % light green

Qs = norminv([p_bar05, p_bar1, p_bar]);     

VaR_05 = draw(:,1) + Qs(1)*draw(:,2);
VaR_1 = draw(:,1) + Qs(2)*draw(:,2);
VaR_5 = draw(:,1) + Qs(3)*draw(:,2);

VaR_C_05 = draw_C(:,1) + Qs(1)*draw_C(:,2);
VaR_C_1 = draw_C(:,1) + Qs(2)*draw_C(:,2);
VaR_C_5 = draw_C(:,1) + Qs(3)*draw_C(:,2);

VaR_C0_05 = draw_C0(:,1) + Qs(1)*draw_C0(:,2);
VaR_C0_1 = draw_C0(:,1) + Qs(2)*draw_C0(:,2);
VaR_C0_5 = draw_C0(:,1) + Qs(3)*draw_C0(:,2);

%% IID
ff = figure(11);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.48 0.42]);    

% set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.4]);
% ax1 = axes('Position',[0.03 0.15 0.36 0.8],'Visible','on');
% axes(ax1)

subplot(1,3,1)
hold on
% histogram(VaR_05,100)
% histogram(VaR_C_05,100)
% histogram(VaR_C0_05,100)
[f,xi] = ksdensity(VaR_05);%,'bandwidth',bw);%
[f_C,xi_C] = ksdensity(VaR_C_05,'bandwidth',0.4);
[f_C0,xi_C0] = ksdensity(VaR_C0_05);%,'bandwidth',bw);
hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_C0,f_C0,'color',ColPalette(2,:),'linewidth',lw)
P = get(gca,'Position');
set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
% if (T == 100)
%     if (sigma2 == 1)
%         xlim([-2 4])
%     else
%         xlim([-2 7])
%     end
% elseif (T == 10000)
%     if (sigma2 == 1)
%         xlim([-0.25 0.25])
%     else
%         xlim([-0.25 1])
%     end
% end
YL = get(gca,'YLim');
line([q05 q05], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
xlabel('VaR 99.5\%','FontSize',12,'Interpreter','latex')


subplot(1,3,2)
hold on
% histogram(VaR_05,100)
% histogram(VaR_C_05,100)
% histogram(VaR_C0_05,100)
[f,xi] = ksdensity(VaR_1);%,'bandwidth',bw);%
[f_C,xi_C] = ksdensity(VaR_C_1,'bandwidth',0.4);
[f_C0,xi_C0] = ksdensity(VaR_C0_1);%,'bandwidth',bw);
hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_C0,f_C0,'color',ColPalette(2,:),'linewidth',lw)
P = get(gca,'Position');
set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
% if (T == 100)
%     if (sigma2 == 1)
%         xlim([-2 4])
%     else
%         xlim([-2 7])
%     end
% elseif (T == 10000)
%     if (sigma2 == 1)
%         xlim([-0.25 0.25])
%     else
%         xlim([-0.25 1])
%     end
% end
YL = get(gca,'YLim');
line([q1 q1], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
xlabel('VaR 99\%','FontSize',12,'Interpreter','latex')



subplot(1,3,3)
hold on
% histogram(VaR_05,100)
% histogram(VaR_C_05,100)
% histogram(VaR_C0_05,100)
[f,xi] = ksdensity(VaR_5);%,'bandwidth',bw);%
[f_C,xi_C] = ksdensity(VaR_C_5,'bandwidth',0.4);
[f_C0,xi_C0] = ksdensity(VaR_C0_5);%,'bandwidth',bw);
hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_C0,f_C0,'color',ColPalette(2,:),'linewidth',lw)
P = get(gca,'Position');
set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
% if (T == 100)
%     if (sigma2 == 1)
%         xlim([-2 4])
%     else
%         xlim([-2 7])
%     end
% elseif (T == 10000)
%     if (sigma2 == 1)
%         xlim([-0.25 0.25])
%     else
%         xlim([-0.25 1])
%     end
% end
YL = get(gca,'YLim');
line([q5 q5], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
xlabel('VaR 95\%','FontSize',12,'Interpreter','latex')


leg = legend('Uncensored','CP 10\%','CP 0','True'); 
set(leg,'Interpreter','latex','FontSize',12, 'location',...  
        'southoutside','Orientation','horizontal');

if save_on
    name = ['figures/',model,'/',model,'_1_',num2str(sigma2),'_T',num2str(T),'_dens_VaR.eps'];
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
