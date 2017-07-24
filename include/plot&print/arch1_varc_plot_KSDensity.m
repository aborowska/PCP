%% mu
ff = figure(1);
set(gcf,'units','normalized','outerposition',[0.2 0.2 0.3 0.4]);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

[f,xi] = ksdensity(draw(:,1));
[f_Cah,xi_Cah] = ksdensity(draw_Cah(:,1));
[f_Cm,xi_Cm] = ksdensity(draw_Cm(:,1));
[f_PCah,xi_PCah] = ksdensity(draw_PCah(:,1));
[f_PCm,xi_PCm] = ksdensity(draw_PCm(:,1));
hold on
plot(xi,f,'k','linewidth',1)
plot(xi_Cah,f_Cah,'color',[0    0.4470    0.7410],'linewidth',1)
plot(xi_Cm,f_Cm,'color',[0.3010    0.7450    0.9330],'linewidth',1)
plot(xi_PCah,f_PCah,'color', [0.4660    0.6740    0.1880],'linewidth',1)
plot(xi_PCm,f_PCm,'color',[0 1 0],'linewidth',1)
YL = get(gca,'YLim');
line([c c], YL,'Color','r','LineWidth',1); 
hold off
% plotTickLatex2D('FontSize',10);
XL = xlabel('\mu','FontSize',10);
% pos = get(XL,'pos'); % Read position [x y z]
% set(XL,'pos',pos+[0 0.8 0]) % Move label up 

leg = legend('Uncensored','CP var ah','CP var mle','PC var ah','PC var mle','True');%,)
% set(leg,'Interpreter','latex','FontSize',10,'position',[0.890 0.42 0.10 0.16])
set(leg,'Interpreter','latex','FontSize',10,'location','eastoutside')


name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_dens_mu_varc.eps'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-depsc','-r0')



%%
ff = figure(2);
set(gcf,'units','normalized','outerposition',[0.2 0.2 0.3 0.4]);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

[f,xi] = ksdensity(draw(:,2));
[f_Cah,xi_Cah] = ksdensity(draw_Cah(:,2));
[f_Cm,xi_Cm] = ksdensity(draw_Cm(:,2));
[f_PCah,xi_PCah] = ksdensity(draw_PCah(:,2));
[f_PCm,xi_PCm] = ksdensity(draw_PCm(:,2));
hold on
plot(xi,f,'k','linewidth',1)
plot(xi_Cah,f_Cah,'color',[         0    0.4470    0.7410],'linewidth',1)
plot(xi_Cm,f_Cm,'color',[    0.3010    0.7450    0.9330],'linewidth',1)
plot(xi_PCah,f_PCah,'color', [       0.4660    0.6740    0.1880],'linewidth',1)
plot(xi_PCm,f_PCm,'color',[     0 1 0],'linewidth',1)
YL = get(gca,'YLim');
line([omega omega], YL,'Color','r','linewidth',1); 
hold off
% plotTickLatex2D('FontSize',10);
XL = xlabel('\omega','FontSize',10);
% pos = get(XL,'pos'); % Read position [x y z]
% set(XL,'pos',pos+[0 0.8 0]) % Move label up 

leg = legend('Uncensored','CP var ah','CP var mle','PC var ah','PC var mle','True');%,)
% set(leg,'Interpreter','latex','FontSize',10,'position',[0.890 0.42 0.10 0.16])
set(leg,'Interpreter','latex','FontSize',10,'location','eastoutside')

name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_dens_omega_varc.eps'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-depsc','-r0')

%%  mu2
ff = figure(3);
set(gcf,'units','normalized','outerposition',[0.2 0.2 0.3 0.4]);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

[f,xi] = ksdensity(draw(:,3));
[f_Cah,xi_Cah] = ksdensity(draw_Cah(:,3),'bandwidth',0.02);%);
[f_Cm,xi_Cm] = ksdensity(draw_Cm(:,3));
[f_PCah,xi_PCah] = ksdensity(draw_PCah(:,3),'bandwidth',0.01);%);
[f_PCm,xi_PCm] = ksdensity(draw_PCm(:,3),'bandwidth',0.01);%);
hold on
plot(xi,f,'k','linewidth',1)
plot(xi_Cah,f_Cah,'color',[         0    0.4470    0.7410],'linewidth',1)
plot(xi_Cm,f_Cm,'color',[    0.3010    0.7450    0.9330],'linewidth',1)
plot(xi_PCah,f_PCah,'color', [       0.4660    0.6740    0.1880],'linewidth',1)
plot(xi_PCm,f_PCm,'color',[     0 1 0],'linewidth',1)
YL = get(gca,'YLim');
line([mu2 mu2], YL,'Color','r','linewidth',1); 
hold off
% plotTickLatex2D('FontSize',10);
XL = xlabel('\mu_{2}','FontSize',10);
% pos = get(XL,'pos'); % Read position [x y z]
% set(XL,'pos',pos+[0 2 0]) % Move label up 

leg = legend('Uncensored','CP var ah','CP var mle','PC var ah','PC var mle','True');%,)
% set(leg,'Interpreter','latex','FontSize',10,'position',[0.890 0.42 0.10 0.16])
set(leg,'Interpreter','latex','FontSize',10,'location','eastoutside')

name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_dens_mu2_varc.eps'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-depsc','-r0')

%%  alpha
ff = figure(4);
set(gcf,'units','normalized','outerposition',[0.2 0.2 0.3 0.4]);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

[f,xi] = ksdensity(draw(:,4));
[f_Cah,xi_Cah] = ksdensity(draw_Cah(:,4),'bandwidth',0.02);%);
[f_Cm,xi_Cm] = ksdensity(draw_Cm(:,4));
[f_PCah,xi_PCah] = ksdensity(draw_PCah(:,4),'bandwidth',0.01);%);
[f_PCm,xi_PCm] = ksdensity(draw_PCm(:,4),'bandwidth',0.01);%);
hold on
plot(xi,f,'k','linewidth',1)
plot(xi_Cah,f_Cah,'color',[         0    0.4470    0.7410],'linewidth',1)
plot(xi_Cm,f_Cm,'color',[    0.3010    0.7450    0.9330],'linewidth',1)
plot(xi_PCah,f_PCah,'color', [       0.4660    0.6740    0.1880],'linewidth',1)
plot(xi_PCm,f_PCm,'color',[     0 1 0],'linewidth',1)
YL = get(gca,'YLim');
line([alpha alpha], YL,'Color','r','linewidth',1); 
hold off
% plotTickLatex2D('FontSize',10);
XL = xlabel('\alpha','FontSize',10);
% pos = get(XL,'pos'); % Read position [x y z]
% set(XL,'pos',pos+[0 2 0]) % Move label up 

leg = legend('Uncensored','CP var ah','CP var mle','PC var ah','PC var mle','True');%,)
% set(leg,'Interpreter','latex','FontSize',10,'position',[0.890 0.42 0.10 0.16])
set(leg,'Interpreter','latex','FontSize',10,'location','eastoutside')

name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_dens_alpha_varc.eps'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-depsc','-r0')
