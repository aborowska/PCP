bw = 0.008;
bw = 0.05;
bw = 0.1;

ff = figure(11);
% set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.6]);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.4]);
ax1 = axes('Position',[0.03 0.15 0.36 0.8],'Visible','on');
axes(ax1)
[f,xi,bw] = ksdensity(draw(:,1));%,'bandwidth',bw);%
[f_C,xi_C,bw] = ksdensity(draw_C(:,1));%,'bandwidth',bw);
[f_C0,xi_C0,bw] = ksdensity(draw_C0(:,1));%,'bandwidth',bw);
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
[f,xi] = ksdensity(draw(:,2));%,'bandwidth',bw);%
[f_C,xi_C] = ksdensity(draw_C(:,2));%,'bandwidth',bw);%
[f_C0,xi_C0] = ksdensity(draw_C0(:,2));%,'bandwidth',bw);%
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


%% AR1
ff = figure(11);
% % set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.55]);
% % ax1 = axes('Position',[0.03 0.15 0.26 0.8],'Visible','on');

% set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.9]);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.6]);
ax1 = axes('Position',[0.03 0.58 0.45 0.38],'Visible','on');
axes(ax1)
[f,xi] = ksdensity(draw(:,1));
[f_C,xi_C] = ksdensity(draw_C(:,1));
[f_C0,xi_C0] = ksdensity(draw_C0(:,1));
[f_PC,xi_PC] = ksdensity(draw_PC(:,1));
[f_PC0,xi_PC0] = ksdensity(draw_PC0(:,1));
hold on
plot(xi,f,'k','linewidth',2)
plot(xi_C,f_C,'color',[         0    0.4470    0.7410],'linewidth',2)
plot(xi_C0,f_C0,'color',[    0.3010    0.7450    0.9330],'linewidth',2)
plot(xi_PC,f_PC,'color', [       0.4660    0.6740    0.1880],'linewidth',2)
plot(xi_PC0,f_PC0,'color',[     0 1 0],'linewidth',2)
YL = get(gca,'YLim');
line([c c], YL,'Color','r','LineWidth',2); 
hold off
plotTickLatex2D('FontSize',12);
XL = xlabel('\mu','FontSize',12);
pos = get(XL,'pos'); % Read position [x y z]
set(XL,'pos',pos+[0 0.8 0]) % Move label up 


% ax2 = axes('Position',[0.32 0.15 0.26 0.8],'Visible','on');
ax2 = axes('Position',[0.53 0.58 0.45 0.38],'Visible','on');
axes(ax2)
[f,xi] = ksdensity(draw(:,2));
[f_C,xi_C] = ksdensity(draw_C(:,2));
[f_C0,xi_C0] = ksdensity(draw_C0(:,2));
[f_PC,xi_PC] = ksdensity(draw_PC(:,2));
[f_PC0,xi_PC0] = ksdensity(draw_PC0(:,2));
hold on
plot(xi,f,'k','linewidth',2)
plot(xi_C,f_C,'color',[         0    0.4470    0.7410],'linewidth',2)
plot(xi_C0,f_C0,'color',[    0.3010    0.7450    0.9330],'linewidth',2)
plot(xi_PC,f_PC,'color', [       0.4660    0.6740    0.1880],'linewidth',2)
plot(xi_PC0,f_PC0,'color',[     0 1 0],'linewidth',2)
YL = get(gca,'YLim');
line([sigma2 sigma2], YL,'Color','r','LineWidth',2); 
hold off
plotTickLatex2D('FontSize',12);
XL = xlabel('\sigma','FontSize',12);
pos = get(XL,'pos'); % Read position [x y z]
set(XL,'pos',pos+[0 0.8 0]) % Move label up 

% ax3 = axes('Position',[0.61 0.15 0.26 0.8],'Visible','on');
ax3 = axes('Position',[0.03 0.1 0.45 0.38],'Visible','on');
axes(ax3)
[f,xi] = ksdensity(draw(:,3));
[f_C,xi_C] = ksdensity(draw_C(:,3),'bandwidth',0.02);%);
[f_C0,xi_C0] = ksdensity(draw_C0(:,3));
[f_PC,xi_PC] = ksdensity(draw_PC(:,3),'bandwidth',0.01);%);
[f_PC0,xi_PC0] = ksdensity(draw_PC0(:,3),'bandwidth',0.01);%);
hold on
plot(xi,f,'k','linewidth',2)
plot(xi_C,f_C,'color',[         0    0.4470    0.7410],'linewidth',2)
plot(xi_C0,f_C0,'color',[    0.3010    0.7450    0.9330],'linewidth',2)
plot(xi_PC,f_PC,'color', [       0.4660    0.6740    0.1880],'linewidth',2)
plot(xi_PC0,f_PC0,'color',[     0 1 0],'linewidth',2)
YL = get(gca,'YLim');
line([rho rho], YL,'Color','r','LineWidth',2); 
hold off
plotTickLatex2D('FontSize',12);
XL = xlabel('\rho','FontSize',12);
pos = get(XL,'pos'); % Read position [x y z]
set(XL,'pos',pos+[0 2 0]) % Move label up 



leg = legend('Uncensored','CP 10\%','CP 0','PCP 10\%','PCP 0','True');%,)
% set(leg,'Interpreter','latex','FontSize',12,'position',[0.890 0.42 0.10 0.16])
set(leg,'Interpreter','latex','FontSize',12,'position',[0.63 0.15 0.23 0.3])

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