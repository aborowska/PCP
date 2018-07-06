addpath(genpath('include/'));
clear all
close all


%model = 'garch';
model = 'tgarch';
% T = 1000;
ver = version('-release');

data_name = 'IBM_T2000_crisis';
    
lw = 1;

ColPalette = [     0    0.4470    0.7410    % deep blue
              0.3010    0.7450    0.9330    % light blue
              0.0000    0.5000    0.0000    % deep green        
              0.4660    0.6740    0.1880];  % light green
%%
ff = figure(11);
if strcmp(ver,'2015b')
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.8]);    
else
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.35]);
    % set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.4]);
end
% ax1 = axes('Position',[0.03 0.15 0.36 0.8],'Visible','on');
% axes(ax1)
subplot(2,3,1)
[f,xi] = ksdensity(draw(:,1));%,'bandwidth',bw);%
[f_C,xi_C] = ksdensity(draw_C(:,1)); %,'bandwidth',0.4);
[f_Cm,xi_Cm] = ksdensity(draw_Cm(:,1));%,'bandwidth',bw);
[f_PC,xi_PC] = ksdensity(draw_PC(:,1)); %,'bandwidth',0.4);
[f_PCm,xi_PCm] = ksdensity(draw_PCm(:,1));%,'bandwidth',bw);

hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)

plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)

% P = get(gca,'Position');
% set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
xlim([-10 120])
% YL = get(gca,'YLim');
% line([param_true(1) param_true(1)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
xlabel('\nu','FontSize',12)

% ax2 = axes('Position',[0.43 0.15 0.36 0.8],'Visible','on');
% axes(ax2)
subplot(2,3,2)

[f,xi] = ksdensity(draw(:,2));%,'bandwidth',bw);%
[f_C,xi_C] = ksdensity(draw_C(:,2)); %,'bandwidth',0.4);
[f_Cm,xi_Cm] = ksdensity(draw_Cm(:,2));%,'bandwidth',bw);
[f_PC,xi_PC] = ksdensity(draw_PC(:,2)); %,'bandwidth',0.4);
[f_PCm,xi_PCm] = ksdensity(draw_PCm(:,2));%,'bandwidth',bw);

hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)

plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)

xlim([-0.5 2.5])

% P = get(gca,'Position');
% set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
% YL = get(gca,'YLim');
% line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
xlabel('\mu','FontSize',12)

subplot(2,3,3)

[f,xi] = ksdensity(draw(:,3));%,'bandwidth',bw);%
[f_C,xi_C] = ksdensity(draw_C(:,3)); %,'bandwidth',0.4);
[f_Cm,xi_Cm] = ksdensity(draw_Cm(:,3));%,'bandwidth',bw);
[f_PC,xi_PC] = ksdensity(draw_PC(:,3)); %,'bandwidth',0.4);
[f_PCm,xi_PCm] = ksdensity(draw_PCm(:,3));%,'bandwidth',bw);

hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)

plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)

xlim([0 1000])

% P = get(gca,'Position');
% set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
% YL = get(gca,'YLim');
% line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
xlabel('\omega','FontSize',12)

subplot(2,3,4)
bw = 0.01;
[f,xi,b] = ksdensity(draw(:,4),'bandwidth',bw);%
[f_C,xi_C,b_C] = ksdensity(draw_C(:,4),'bandwidth',2*bw);
[f_Cm,xi_Cm,b_Cm] = ksdensity(draw_Cm(:,4),'bandwidth',2*bw);
[f_PC,xi_PC,b_PC] = ksdensity(draw_PC(:,4),'bandwidth',bw);
[f_PCm,xi_PCm,b_PCm] = ksdensity(draw_PCm(:,4),'bandwidth',bw);

hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)

plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)

% xlim([-0.5 2.5])
% P = get(gca,'Position');
% set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
% YL = get(gca,'YLim');
% line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
xlabel('\alpha','FontSize',12)



subplot(2,3,5)
bw = 0.01;
[f,xi,b] = ksdensity(draw(:,5),'bandwidth',bw);%
[f_C,xi_C,b_C] = ksdensity(draw_C(:,5),'bandwidth',2*bw);
[f_Cm,xi_Cm,b_Cm] = ksdensity(draw_Cm(:,5),'bandwidth',2*bw);
[f_PC,xi_PC,b_PC] = ksdensity(draw_PC(:,5),'bandwidth',bw);
[f_PCm,xi_PCm,b_PCm] = ksdensity(draw_PCm(:,5),'bandwidth',bw);

hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)

plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)

% xlim([-0.5 2.5])
% P = get(gca,'Position');
% set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
% YL = get(gca,'YLim');
% line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
xlabel('\beta','FontSize',12)



leg = legend('Uncensored','CP','CP$_{var}$','PCP','PCP$_{var}$'); 
set(leg,'Interpreter','latex','FontSize',12, 'location', ...%,'position',[0.63 0.15 0.23 0.3])
        'eastoutside');%,'Orientation','horizontal');

if save_on
    name = ['figures/',model,'/',data_name,'/',data_name,'_ksdens.eps'];
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
