% bw = 0.008;
% bw = 0.05;
% bw = 0.1;

%% ARCH(1)

model = 'arch1';
T = 2500;
% % set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.55]);
% % ax1 = axes('Position',[0.03 0.15 0.26 0.8],'Visible','on');

name = ['C:\Users\aba228\Dropbox\New Projects\Censored Posterior\PCP\results\',model,'\',model,'_1_2_T',num2str(T),'_H100_II10_PCP0_MC_(R2017a)_varc.mat'];
load(name,'-regexp','^draw')
load(name,'-regexp','^mean')

draw_Cm = draw_C0;
draw_PCm = draw_PC0;
draw_Cah = draw_C;
draw_PCah = draw_PC;

mean_draw_Cm = mean_draw_C0;
mean_draw_PCm = mean_draw_PC0;
mean_draw_Cah = mean_draw_C;
mean_draw_PCah = mean_draw_PC;

name = ['C:\Users\aba228\Dropbox\New Projects\Censored Posterior\PCP\results\',model,'\',model,'_1_2_T',num2str(T),'_H100_II10_PCP0_MC_(R2015a)_add.mat'];
load(name,'-regexp','^draw')
load(name,'-regexp','^mean')
load(name,'param_true')

sigma1 = 1;
sigma2 = 2;
mu = 0;    
mu2 = 0; % gama = mu - mu2;
omega = 1;
alpha = 0.1;
c = (sigma2 - sigma1)/sqrt(2*pi); % mean of eps


 
%% mu
ff = figure(1)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

[f,xi] = ksdensity(draw(:,1));
[f_C,xi_C] = ksdensity(draw_C(:,1));
[f_C0,xi_C0] = ksdensity(draw_C0(:,1));
[f_PC,xi_PC] = ksdensity(draw_PC(:,1));
[f_PC0,xi_PC0] = ksdensity(draw_PC0(:,1));
hold on
plot(xi,f,'k','linewidth',1)
plot(xi_C,f_C,'color',[0    0.4470    0.7410],'linewidth',1)
plot(xi_C0,f_C0,'color',[0.3010    0.7450    0.9330],'linewidth',1)
plot(xi_PC,f_PC,'color', [0.4660    0.6740    0.1880],'linewidth',1)
plot(xi_PC0,f_PC0,'color',[0 1 0],'linewidth',1)
YL = get(gca,'YLim');
line([c c], YL,'Color','r','LineWidth',1); 
hold off
% plotTickLatex2D('FontSize',10);
XL = xlabel('\mu','FontSize',10);
% pos = get(XL,'pos'); % Read position [x y z]
% set(XL,'pos',pos+[0 0.8 0]) % Move label up 

leg = legend('Uncensored','CP 10\%','CP 0','PCP 10\%','PCP 0','True');%,)
% set(leg,'Interpreter','latex','FontSize',10,'position',[0.890 0.42 0.10 0.16])
set(leg,'Interpreter','latex','FontSize',10,'location','southeastoutside')


%%
ff = figure(2)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

[f,xi] = ksdensity(draw(:,2));
[f_C,xi_C] = ksdensity(draw_C(:,2));
[f_C0,xi_C0] = ksdensity(draw_C0(:,2));
[f_PC,xi_PC] = ksdensity(draw_PC(:,2));
[f_PC0,xi_PC0] = ksdensity(draw_PC0(:,2));
hold on
plot(xi,f,'k','linewidth',1)
plot(xi_C,f_C,'color',[         0    0.4470    0.7410],'linewidth',1)
plot(xi_C0,f_C0,'color',[    0.3010    0.7450    0.9330],'linewidth',1)
plot(xi_PC,f_PC,'color', [       0.4660    0.6740    0.1880],'linewidth',1)
plot(xi_PC0,f_PC0,'color',[     0 1 0],'linewidth',1)
YL = get(gca,'YLim');
line([omega omega], YL,'Color','r','linewidth',1); 
hold off
% plotTickLatex2D('FontSize',10);
XL = xlabel('\omega','FontSize',10);
pos = get(XL,'pos'); % Read position [x y z]
set(XL,'pos',pos+[0 0.8 0]) % Move label up 

leg = legend('Uncensored','CP 10\%','CP 0','PCP 10\%','PCP 0','True');%,)
% set(leg,'Interpreter','latex','FontSize',10,'position',[0.890 0.42 0.10 0.16])
set(leg,'Interpreter','latex','FontSize',10,'position',[0.63 0.15 0.23 0.3])

%%  mu2
ff = figure(3)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

[f,xi] = ksdensity(draw(:,3));
[f_C,xi_C] = ksdensity(draw_C(:,3),'bandwidth',0.02);%);
[f_C0,xi_C0] = ksdensity(draw_C0(:,3));
[f_PC,xi_PC] = ksdensity(draw_PC(:,3),'bandwidth',0.01);%);
[f_PC0,xi_PC0] = ksdensity(draw_PC0(:,3),'bandwidth',0.01);%);
hold on
plot(xi,f,'k','linewidth',1)
plot(xi_C,f_C,'color',[         0    0.4470    0.7410],'linewidth',1)
plot(xi_C0,f_C0,'color',[    0.3010    0.7450    0.9330],'linewidth',1)
plot(xi_PC,f_PC,'color', [       0.4660    0.6740    0.1880],'linewidth',1)
plot(xi_PC0,f_PC0,'color',[     0 1 0],'linewidth',1)
YL = get(gca,'YLim');
line([mu2 mu2], YL,'Color','r','linewidth',1); 
hold off
% plotTickLatex2D('FontSize',10);
XL = xlabel('\mu_{2}','FontSize',10);
pos = get(XL,'pos'); % Read position [x y z]
set(XL,'pos',pos+[0 2 0]) % Move label up 

leg = legend('Uncensored','CP 10\%','CP 0','PCP 10\%','PCP 0','True');%,)
% set(leg,'Interpreter','latex','FontSize',10,'position',[0.890 0.42 0.10 0.16])
set(leg,'Interpreter','latex','FontSize',10,'position',[0.63 0.15 0.23 0.3])


%%  mu2
ff = figure(4)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

[f,xi] = ksdensity(draw(:,3));
[f_C,xi_C] = ksdensity(draw_C(:,3),'bandwidth',0.02);%);
[f_C0,xi_C0] = ksdensity(draw_C0(:,3));
[f_PC,xi_PC] = ksdensity(draw_PC(:,3),'bandwidth',0.01);%);
[f_PC0,xi_PC0] = ksdensity(draw_PC0(:,3),'bandwidth',0.01);%);
hold on
plot(xi,f,'k','linewidth',1)
plot(xi_C,f_C,'color',[         0    0.4470    0.7410],'linewidth',1)
plot(xi_C0,f_C0,'color',[    0.3010    0.7450    0.9330],'linewidth',1)
plot(xi_PC,f_PC,'color', [       0.4660    0.6740    0.1880],'linewidth',1)
plot(xi_PC0,f_PC0,'color',[     0 1 0],'linewidth',1)
YL = get(gca,'YLim');
line([alpha alpha], YL,'Color','r','linewidth',1); 
hold off
% plotTickLatex2D('FontSize',10);
XL = xlabel('\alpha}','FontSize',10);
pos = get(XL,'pos'); % Read position [x y z]
set(XL,'pos',pos+[0 2 0]) % Move label up 

leg = legend('Uncensored','CP 10\%','CP 0','PCP 10\%','PCP 0','True');%,)
% set(leg,'Interpreter','latex','FontSize',10,'position',[0.890 0.42 0.10 0.16])
set(leg,'Interpreter','latex','FontSize',10,'position',[0.63 0.15 0.23 0.3])






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