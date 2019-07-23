addpath(genpath('include/'));


nu = 5;
xxx = -5:0.01:5;
p_bar05 = 0.005;
p_bar1 = 0.01;
p_bar = 0.05;


Lambda = [-0.5,-0.3,-0.1,0,0.1,0.3,0.5];

for ii = 1:length(Lambda)
    
    lambda = Lambda(ii);
    q05 = skewtinv(p_bar05,nu,lambda); % -3.5685    -3.75313
    q1 = skewtinv(p_bar1,nu,lambda); % -2.9420  -3.0798
    q5 = skewtinv(p_bar,nu,lambda); 

    ff = figure(ii);
    set(gcf,'defaulttextinterpreter','latex');

    pdfxxx = skewtpdf(xxx,nu,lambda);
    hold on
    plot(xxx,pdfxxx);    
    scatter([q05,q05],[0,skewtpdf(q05,nu,lambda)],20,[0.6350, 0.0780, 0.1840],'filled')
    scatter([q1,q1],[0,skewtpdf(q1,nu,lambda)],20,[0.8500, 0.3250, 0.0980],'filled')
    scatter([q5,q5],[0,skewtpdf(q5,nu,lambda)],20,[0.9290, 0.6940, 0.1250],'filled')

    line([q05,q05],[0,skewtpdf(q05,nu,lambda)],'color',[0.6350, 0.0780, 0.1840])
    line([q1,q1],[0,skewtpdf(q1,nu,lambda)],'color',[0.8500, 0.3250, 0.0980])
    line([q5,q5],[0,skewtpdf(q5,nu,lambda)],'color',[0.9290, 0.6940, 0.1250])
    legend({'pdf','0.5\% percentile','1\% percentile','5\% percentile'},...
        'Interpreter','latex','fontsize',11);
    set(gca,'TickLabelInterpreter','latex')
            
    set(gcf,'PaperPositionMode','auto');
    name = ['figures/skewed_t_quant_lambda',num2str(lambda),'.eps'];
    print(ff,name,'-depsc','-r0')
    name = ['figures/skewed_t_quant_lambda',num2str(lambda),'.png'];
    print(ff,name,'-dpng','-r0')                  
end





for ii = 1:length(Lambda)
    
    lambda = Lambda(ii);
    q05 = skewtinv(p_bar05,nu,lambda); % -3.5685    -3.7531
    q1 = skewtinv(p_bar1,nu,lambda); % -2.9420  -3.0798
    q5 = skewtinv(p_bar,nu,lambda); 

    ff = figure(ii);
    set(gcf,'defaulttextinterpreter','latex');

    pdfxxx = skewtpdf(xxx,nu,lambda);
    hold on
    plot(xxx,pdfxxx);    
    scatter([q05,q05],[0,skewtpdf(q05,nu,lambda)],20,[0.6350, 0.0780, 0.1840],'filled')
    scatter([q1,q1],[0,skewtpdf(q1,nu,lambda)],20,[0.8500, 0.3250, 0.0980],'filled')
    scatter([q5,q5],[0,skewtpdf(q5,nu,lambda)],20,[0.9290, 0.6940, 0.1250],'filled')

    line([q05,q05],[0,skewtpdf(q05,nu,lambda)],'color',[0.6350, 0.0780, 0.1840])
    line([q1,q1],[0,skewtpdf(q1,nu,lambda)],'color',[0.8500, 0.3250, 0.0980])
    line([q5,q5],[0,skewtpdf(q5,nu,lambda)],'color',[0.9290, 0.6940, 0.1250])
    legend({'pdf','0.5th percentile','1th percentile','5th percentile'},...
        'Interpreter','latex','fontsize',11);
    set(gca,'TickLabelInterpreter','latex')
            
    set(gcf,'PaperPositionMode','auto');
    name = ['figures/skewed_t_quant_lambda_0',num2str(10*lambda),'.eps'];
    print(ff,name,'-depsc','-r0')
    name = ['figures/skewed_t_quant_lambda',num2str(lambda),'.png'];
    print(ff,name,'-dpng','-r0')                  
end