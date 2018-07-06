addpath(genpath('include/'));
clear all
close all


% model = 'iid';
model = 'ar1';
T = 100;
% sigma2 = 2;
sigma2 = 1;

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
%% IID
ff = figure(11);
if strcmp(version('-release'),'2015b')
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.48 0.42]);    
else
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.35]);
end
% set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.4]);
% ax1 = axes('Position',[0.03 0.15 0.36 0.8],'Visible','on');
% axes(ax1)
subplot(1,2,1)
[f,xi] = ksdensity(draw(:,1));%,'bandwidth',bw);%
[f_C,xi_C] = ksdensity(draw_C(:,1)); %,'bandwidth',0.4);
[f_C0,xi_C0] = ksdensity(draw_C0(:,1));%,'bandwidth',bw);
hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_C0,f_C0,'color',ColPalette(2,:),'linewidth',lw)
P = get(gca,'Position');
set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
if (T == 100)
    if (sigma2 == 1)
        xlim([-2 4])
    else
        xlim([-2 7])
    end
elseif (T == 10000)
    if (sigma2 == 1)
        xlim([-0.25 0.25])
    else
        xlim([-0.25 1])
    end
end
YL = get(gca,'YLim');
line([param_true(1) param_true(1)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
xlabel('\mu','FontSize',12)

% ax2 = axes('Position',[0.43 0.15 0.36 0.8],'Visible','on');
% axes(ax2)
subplot(1,2,2)
hold on
[f,xi] = ksdensity(draw(:,2));%,'bandwidth',bw);%
[f_C,xi_C] = ksdensity(draw_C(:,2)); %,'bandwidth',0.4);%
[f_C0,xi_C0] = ksdensity(draw_C0(:,2));%,'bandwidth',bw);%
hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:), 'linewidth',lw)
plot(xi_C0,f_C0,'color',ColPalette(2,:) ,'linewidth',lw)
P = get(gca,'Position');
set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
if (sigma2 == 1)
    if (T == 1000)
        xlim([0.65 1.35])
    elseif (T == 100)
        xlim([0.25 4])
    else
        xlim([0.85 1.15])
    end
end
YL = get(gca,'YLim');
line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
xlabel('\sigma','FontSize',12)

leg = legend('Uncensored','CP 10\%','CP 0','True'); 
set(leg,'Interpreter','latex','FontSize',12, 'location',...  
        'northoutside','Orientation','horizontal');

if save_on
    name = ['figures/',model,'/',model,'_1_',num2str(sigma2),'_T',num2str(T),'_dens.eps'];
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


%% AR1
ff = figure(11);
% % set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.55]);
% % ax1 = axes('Position',[0.03 0.15 0.26 0.8],'Visible','on');
if strcmp(version('-release'),'2015b')
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.72 0.42]);    
else
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.35]);
end
% set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.6]);
% ax1 = axes('Position',[0.03 0.58 0.45 0.38],'Visible','on');
% axes(ax1)
subplot(1,3,1)
[f,xi] = ksdensity(draw(:,1),'function','pdf');
[f_C,xi_C] = ksdensity(draw_C(:,1),'function','pdf');
[f_C0,xi_C0] = ksdensity(draw_C0(:,1),'function','pdf');
[f_PC,xi_PC] = ksdensity(draw_PC(:,1),'function','pdf');
[f_PC0,xi_PC0] = ksdensity(draw_PC0(:,1),'function','pdf');
hold on
plot(xi	,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_C0,f_C0,'color',ColPalette(2,:) ,'linewidth',lw)
plot(xi_PC,f_PC,'color', ColPalette(3,:) ,'linewidth',lw)
plot(xi_PC0,f_PC0,'color',ColPalette(4,:) ,'linewidth',lw)
P = get(gca,'Position');
set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
if (T == 100)
    if (sigma2 == 1)
        xlim([-2.5 2.5])
    else
        xlim([-5 7.5])
    end
elseif ((T == 10000) && (sigma2  == 1))
    xlim([-0.25 0.25])    
else
    if (sigma2 == 1)
        xlim([-0.75 0.75])        
    else
        xlim([-1 2])
    end
end
YL = get(gca,'YLim');
line([param_true(1) param_true(1)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
xlabel('\mu','FontSize',12);
% pos = get(XL,'pos'); % Read position [x y z]
% set(XL,'pos',pos+[0 0.8 0]) % Move label up 

% ax2 = axes('Position',[0.32 0.15 0.26 0.8],'Visible','on');
%ax2 = axes('Position',[0.53 0.58 0.45 0.38],'Visible','on');
%axes(ax2)
subplot(1,3,2)
[f,xi] = ksdensity(draw(:,2),'function','pdf'); %,'bandwidth',0.04);
[f_C,xi_C] = ksdensity(draw_C(:,2),'function','pdf');%,'bandwidth',0.1);
[f_C0,xi_C0] = ksdensity(draw_C0(:,2),'function','pdf');%,'bandwidth',0.5);
[f_PC,xi_PC] = ksdensity(draw_PC(:,2),'function','pdf');%,'bandwidth',0.5);
[f_PC0,xi_PC0] = ksdensity(draw_PC0(:,2),'function','pdf');%,'bandwidth',0.1);
hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_C0,f_C0,'color',ColPalette(2,:) ,'linewidth',lw)
plot(xi_PC,f_PC,'color', ColPalette(3,:) ,'linewidth',lw)
plot(xi_PC0,f_PC0,'color',ColPalette(4,:) ,'linewidth',lw)
P = get(gca,'Position');
set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
if (T == 100)
    if (sigma2 == 1)
        xlim([0 3])
    else
        xlim([0 6])
    end
end

YL = get(gca,'YLim');
line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
XL = xlabel('\sigma','FontSize',12);
% pos = get(XL,'pos'); % Read position [x y z]
% set(XL,'pos',pos+[0 0.8 0]) % Move label up 

% ax3 = axes('Position',[0.61 0.15 0.26 0.8],'Visible','on');
% ax3 = axes('Position',[0.03 0.1 0.45 0.38],'Visible','on');
% axes(ax3)
subplot(1,3,3)
[f,xi] = ksdensity(draw(:,3));
if ((T == 100) && (sigma2 == 1))
    [f_C,xi_C,bw_C] = ksdensity(draw_C(:,3),'function','pdf','bandwidth',0.06);%);
else
    [f_C,xi_C] = ksdensity(draw_C(:,3),'function','pdf');%,'bandwidth',0.04);%);
end
[f_C0,xi_C0] = ksdensity(draw_C0(:,3),'function','pdf');%,'bandwidth',0.01);
[f_PC,xi_PC] = ksdensity(draw_PC(:,3),'function','pdf');%,'bandwidth',0.01);%);
[f_PC0,xi_PC0] = ksdensity(draw_PC0(:,3),'function','pdf');%,'bandwidth',0.01);%);
hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_C0,f_C0,'color',ColPalette(2,:) ,'linewidth',lw)
plot(xi_PC,f_PC,'color', ColPalette(3,:) ,'linewidth',lw)
plot(xi_PC0,f_PC0,'color',ColPalette(4,:) ,'linewidth',lw)
P = get(gca,'Position');
set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])         
if (T == 100)
    xlim([0 1.3])
end
YL = get(gca,'YLim');
line([param_true(3) param_true(3)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
XL = xlabel('\rho','FontSize',12);
% pos = get(XL,'pos'); % Read position [x y z]
% set(XL,'pos',pos+[0 2 0]) % Move label up 


leg = legend('Uncensored','CP 10\%','CP 0','PCP 10\%','PCP 0','True');%,)
% set(leg,'Interpreter','latex','FontSize',12,'position',[0.890 0.42 0.10 0.16])
set(leg,'Interpreter','latex','FontSize',12, 'location', ...%,'position',[0.63 0.15 0.23 0.3])
        'northoutside','Orientation','horizontal');

if save_on
    name = ['figures/',model,'/',model,'_1_',num2str(sigma2),'_T',num2str(T),'_dens.eps'];
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


%% AR(1) time varying
ff = figure(12);
% % set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.55]);
% % ax1 = axes('Position',[0.03 0.15 0.26 0.8],'Visible','on');

set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.35]);
% set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.6]);
% ax1 = axes('Position',[0.03 0.58 0.45 0.38],'Visible','on');
% axes(ax1)
subplot(1,3,1)
[f,xi] = ksdensity(draw(:,1));
[f_C,xi_C] = ksdensity(draw_Cm(:,1));
[f_C0,xi_C0] = ksdensity(draw_Cah(:,1));
[f_PC,xi_PC] = ksdensity(draw_PCm(:,1));
[f_PC0,xi_PC0] = ksdensity(draw_PCah(:,1));
hold on
plot(xi	,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_C0,f_C0,'color',ColPalette(2,:) ,'linewidth',lw)
plot(xi_PC,f_PC,'color', ColPalette(3,:) ,'linewidth',lw)
plot(xi_PC0,f_PC0,'color',ColPalette(4,:) ,'linewidth',lw)
xlim([-1 2])
YL = get(gca,'YLim');
line([param_true(1) param_true(1)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
XL = xlabel('\mu','FontSize',12);
% pos = get(XL,'pos'); % Read position [x y z]
% set(XL,'pos',pos+[0 0.8 0]) % Move label up 
P = get(gca,'Position');
set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])

% ax2 = axes('Position',[0.32 0.15 0.26 0.8],'Visible','on');
%ax2 = axes('Position',[0.53 0.58 0.45 0.38],'Visible','on');
%axes(ax2)
subplot(1,3,2)
[f,xi] = ksdensity(draw(:,2));
[f_C,xi_C] = ksdensity(draw_Cm(:,2));
[f_C0,xi_C0] = ksdensity(draw_Cah(:,2));
[f_PC,xi_PC] = ksdensity(draw_PCm(:,2));
[f_PC0,xi_PC0] = ksdensity(draw_PCah(:,2));
hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_C0,f_C0,'color',ColPalette(2,:) ,'linewidth',lw)
plot(xi_PC,f_PC,'color', ColPalette(3,:) ,'linewidth',lw)
plot(xi_PC0,f_PC0,'color',ColPalette(4,:) ,'linewidth',lw)
YL = get(gca,'YLim');
line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
XL = xlabel('\sigma','FontSize',12);
% pos = get(XL,'pos'); % Read position [x y z]
% set(XL,'pos',pos+[0 0.8 0]) % Move label up 
P = get(gca,'Position');
set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])

% ax3 = axes('Position',[0.61 0.15 0.26 0.8],'Visible','on');
% ax3 = axes('Position',[0.03 0.1 0.45 0.38],'Visible','on');
% axes(ax3)
subplot(1,3,3)
[f,xi] = ksdensity(draw(:,3));
[f_C,xi_C] = ksdensity(draw_Cm(:,3),'bandwidth',0.02);%);
[f_C0,xi_C0] = ksdensity(draw_Cah(:,3));
[f_PC,xi_PC] = ksdensity(draw_PCm(:,3),'bandwidth',0.01);%);
[f_PC0,xi_PC0] = ksdensity(draw_PCah(:,3),'bandwidth',0.01);%);
hold on
plot(xi,f,'k','linewidth',lw)
plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
plot(xi_C0,f_C0,'color',ColPalette(2,:) ,'linewidth',lw)
plot(xi_PC,f_PC,'color', ColPalette(3,:) ,'linewidth',lw)
plot(xi_PC0,f_PC0,'color',ColPalette(4,:) ,'linewidth',lw)
YL = get(gca,'YLim');
line([param_true(3) param_true(3)], YL,'Color','r','LineWidth',lw); 
hold off
% plotTickLatex2D('FontSize',12);
set(gca,'TickLabelInterpreter','latex');
XL = xlabel('\rho','FontSize',12);
% pos = get(XL,'pos'); % Read position [x y z]
% set(XL,'pos',pos+[0 2 0]) % Move label up 
P = get(gca,'Position');
set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])



leg = legend('Uncensored','CP var mle','CP var mf','PCP var mle','PCP var mf','True');%,)
% set(leg,'Interpreter','latex','FontSize',12,'position',[0.890 0.42 0.10 0.16])
set(leg,'Interpreter','latex','FontSize',12, 'location', ...%,'position',[0.63 0.15 0.23 0.3])
        'northoutside','Orientation','horizontal');

if save_on
    name = ['figures/',model,'/',model,'_1_',num2str(sigma2),'_T',num2str(T),'_var_dens.eps'];
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