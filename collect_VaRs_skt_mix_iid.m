% clear all;
model = 'skt_iid';%'mix_iid'; %'mix_iid'; %'skt_iid';

switch model
    case 'skt_iid'
        nu = 5;
        param = nu;
%         Lambda = -0.6 : 0.1: 0.2;
        Lambda = -0.6 : 0.1: 0.6;
    case 'mix_iid'
        rho = 0.75;
        param = rho;
%         Lambda = 0.1 : 0.1: 0.9;
        Lambda = 0.05:0.05:0.45;
end
R = length(Lambda);

 TT = [100,1000,10000];


for t = 1:3
    close all
    T = TT(t);
    
    VaRs_05 = zeros(50,4,R);
    VaRs_1 = zeros(50,4,R);
    VaRs_5 = zeros(50,4,R);
    Draws = zeros(2,3,R);

   
    
    for ii = 1:R
        lambda = Lambda(ii);

        load(['results/',model,'/',model,'_',num2str(lambda),'_',num2str(param),'_T',num2str(T),'_MC.mat'],...
            '-regexp','^VaR','^q','^draw')
         
        %% 99.5%
        VaRs_05(:,1,ii) = q05;
        VaRs_1(:,1,ii) = q1;
        VaRs_5(:,1,ii) = q5;
        
        VaRs_05(:,2,ii) = VaR_05_post(1:50);
        VaRs_05(:,3,ii) = VaR_05_post_C0(1:50);
        VaRs_05(:,4,ii) = VaR_05_post_C(1:50);
        
        VaRs_1(:,2,ii) = VaR_1_post(1:50);
        VaRs_1(:,3,ii) = VaR_1_post_C0(1:50);
        VaRs_1(:,4,ii) = VaR_1_post_C(1:50);
        
        VaRs_5(:,2,ii) = VaR_5_post(1:50);
        VaRs_5(:,3,ii) = VaR_5_post_C0(1:50);
        VaRs_5(:,4,ii) = VaR_5_post_C(1:50);
             
        Draws(:,1,ii) = mean(draw);
        Draws(:,2,ii) = mean(draw_C0);
        Draws(:,3,ii) = mean(draw_C);
    end
 

    VaRmeans_05 = squeeze(mean(VaRs_05));
    VaRmeans_1 = squeeze(mean(VaRs_1));
    VaRmeans_5 = squeeze(mean(VaRs_5));
    VaRstd_05 = squeeze(std(VaRs_05));
    VaRstd_1 = squeeze(std(VaRs_1));
    VaRstd_5 = squeeze(std(VaRs_5));

    Draw_mu = squeeze(Draws(1,:,:));
    Draw_sigma = squeeze(Draws(2,:,:));

    %% Parameter means 
    ff = figure(10);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.35]);
    set(gcf,'defaulttextinterpreter','latex');
    subplot(1,2,1)
    bar(Lambda,Draw_mu')
    xlabel('DGP $\lambda$')
    title('Estimated mean $\mu$')
    legend({'Posterior','CP0','CP10\%'},'Interpreter','latex','fontsize',10,'location','NorthEast');
    set(gca,'TickLabelInterpreter','latex')
    subplot(1,2,2)
    bar(Lambda,Draw_sigma')
    xlabel('DGP $\lambda$')
    title('Estimated mean $\sigma$')
    legend({'Posterior','CP0','CP10\%'},'Interpreter','latex','fontsize',10,'location','NorthEast');
    set(gca,'TickLabelInterpreter','latex')

    set(gcf,'PaperPositionMode','auto');
    name = ['figures/',model,'/',model,'_T',num2str(T),'means_bar.eps'];
    print(ff,name,'-depsc','-r0')
    name = ['figures/',model,'/',model,'_T',num2str(T),'means_bar.png'];
    print(ff,name,'-dpng','-r0')                

%% 99.5% VaR true predicted
    ff = figure(7);
    set(gcf,'defaulttextinterpreter','latex');
    hold on
    scatter(VaRmeans_05(1,:),VaRmeans_05(1,:),30,[0 0 0],'filled')
    scatter(VaRmeans_05(1,:),VaRmeans_05(2,:),'filled')
    scatter(VaRmeans_05(1,:),VaRmeans_05(3,:),'filled')
    scatter(VaRmeans_05(1,:),VaRmeans_05(4,:),'filled')
    for jj=1:length(Lambda)
        line([VaRmeans_05(1,jj),VaRmeans_05(1,jj)],[VaRmeans_05(2,jj)-VaRstd_05(2,jj),VaRmeans_05(2,jj)+VaRstd_05(2,jj)],'color',[0, 0.4, 0.7],'Linewidth',2)
        line([VaRmeans_05(1,jj),VaRmeans_05(1,jj)],[VaRmeans_05(3,jj)-VaRstd_05(3,jj),VaRmeans_05(3,jj)+VaRstd_05(3,jj)],'color',[0.85, 0.32, 0.1],'Linewidth',2)
        line([VaRmeans_05(1,jj),VaRmeans_05(1,jj)],[VaRmeans_05(4,jj)-VaRstd_05(4,jj),VaRmeans_05(4,jj)+VaRstd_05(4,jj)],'color',[0.9, 0.6, 0.1],'Linewidth',2)
    end

    xlabel('True')
    ylabel('Predicted')
    ll=legend({'True','Posterior','CP0','CP10\%'},'Interpreter','latex','fontsize',11,'location','SouthEast');
    set(gca,'TickLabelInterpreter','latex')
    % scatter(VaRmeans_05(1,:),VaRmeans_05(1,:),30,[0 0 0],'filled')
    set(gcf,'PaperPositionMode','auto');
    name = ['figures/',model,'/',model,'_T',num2str(T),'_VaRs_05_true_predicted.eps'];
    print(ff,name,'-depsc','-r0')
    name = ['figures/',model,'/',model,'_T',num2str(T),'_VaRs_05_true_predicted.png'];
    print(ff,name,'-dpng','-r0')      
    
%% 99% VaR true predicted
    ff = figure(8);
    set(gcf,'defaulttextinterpreter','latex');
    hold on
    scatter(VaRmeans_1(1,:),VaRmeans_1(1,:),30,[0 0 0],'filled')
    scatter(VaRmeans_1(1,:),VaRmeans_1(2,:),'filled')
    scatter(VaRmeans_1(1,:),VaRmeans_1(3,:),'filled')
    scatter(VaRmeans_1(1,:),VaRmeans_1(4,:),'filled')
    for jj=1:length(Lambda)
        line([VaRmeans_1(1,jj),VaRmeans_1(1,jj)],[VaRmeans_1(2,jj)-VaRstd_1(2,jj),VaRmeans_1(2,jj)+VaRstd_1(2,jj)],'color',[0, 0.4, 0.7],'Linewidth',2)
        line([VaRmeans_1(1,jj),VaRmeans_1(1,jj)],[VaRmeans_1(3,jj)-VaRstd_1(3,jj),VaRmeans_1(3,jj)+VaRstd_1(3,jj)],'color',[0.85, 0.32, 0.1],'Linewidth',2)
        line([VaRmeans_1(1,jj),VaRmeans_1(1,jj)],[VaRmeans_1(4,jj)-VaRstd_1(4,jj),VaRmeans_1(4,jj)+VaRstd_1(4,jj)],'color',[0.9, 0.6, 0.1],'Linewidth',2)
    end

    xlabel('True')
    ylabel('Predicted')
    ll=legend({'True','Posterior','CP0','CP10\%'},'Interpreter','latex','fontsize',11,'location','SouthEast');
    set(gca,'TickLabelInterpreter','latex')
    % scatter(VaRmeans_1(1,:),VaRmeans_1(1,:),30,[0 0 0],'filled')
    set(gcf,'PaperPositionMode','auto');
    name = ['figures/',model,'/',model,'_T',num2str(T),'_VaRs_1_true_predicted.eps'];
    print(ff,name,'-depsc','-r0')
    name = ['figures/',model,'/',model,'_T',num2str(T),'_VaRs_1_true_predicted.png'];
    print(ff,name,'-dpng','-r0')       
    
%% 95% VaR true predicted
    ff = figure(9);
    set(gcf,'defaulttextinterpreter','latex');
    hold on
    scatter(VaRmeans_5(1,:),VaRmeans_5(1,:),30,[0 0 0],'filled')
    scatter(VaRmeans_5(1,:),VaRmeans_5(2,:),'filled')
    scatter(VaRmeans_5(1,:),VaRmeans_5(3,:),'filled')
    scatter(VaRmeans_5(1,:),VaRmeans_5(4,:),'filled')
    for jj=1:length(Lambda)
        line([VaRmeans_5(1,jj),VaRmeans_5(1,jj)],[VaRmeans_5(2,jj)-VaRstd_5(2,jj),VaRmeans_5(2,jj)+VaRstd_5(2,jj)],'color',[0, 0.4, 0.7],'Linewidth',2)
        line([VaRmeans_5(1,jj),VaRmeans_5(1,jj)],[VaRmeans_5(3,jj)-VaRstd_5(3,jj),VaRmeans_5(3,jj)+VaRstd_5(3,jj)],'color',[0.85, 0.32, 0.1],'Linewidth',2)
        line([VaRmeans_5(1,jj),VaRmeans_5(1,jj)],[VaRmeans_5(4,jj)-VaRstd_5(4,jj),VaRmeans_5(4,jj)+VaRstd_5(4,jj)],'color',[0.9, 0.6, 0.1],'Linewidth',2)
    end

    xlabel('True')
    ylabel('Predicted')
    ll=legend({'True','Posterior','CP0','CP10\%'},'Interpreter','latex','fontsize',11,'location','SouthEast');
    set(gca,'TickLabelInterpreter','latex')
    % scatter(VaRmeans_5(1,:),VaRmeans_5(1,:),30,[0 0 0],'filled')
    set(gcf,'PaperPositionMode','auto');
    name = ['figures/',model,'/',model,'_T',num2str(T),'_VaRs_5_true_predicted.eps'];
    print(ff,name,'-depsc','-r0')
    name = ['figures/',model,'/',model,'_T',num2str(T),'_VaRs_5_true_predicted.png'];
    print(ff,name,'-dpng','-r0')      
        

end
%%

 for t = 1:3
    T = TT(t);
    
    for ii = 1:R
        lambda = round(-0.6 + (ii-1)*0.1,1);

        load(['results/',model,'/',model,'_',num2str(lambda),'_',num2str(nu),'_T',num2str(T),'_MC.mat'],...
            '-regexp','^MSE')               
        
        temp = [mean(MSE_05_post);
                mean(MSE_05_post_C0);
                mean(MSE_05_post_C)];
        meanMSEs_05(:,ii,t) = temp;

        temp = [mean(MSE_1_post);
            mean(MSE_1_post_C0);
            mean(MSE_1_post_C)];
        meanMSEs_1(:,ii,t) = temp;

        temp = [mean(MSE_5_post);
            mean(MSE_5_post_C0);
            mean(MSE_5_post_C)];
        meanMSEs_5(:,ii,t) = temp;


        temp = [std(MSE_05_post);
            std(MSE_05_post_C0);
            std(MSE_05_post_C)];
        stdMSEs_05(:,ii,t) = temp;

        temp = [std(MSE_1_post);
            std(MSE_1_post_C0);
            std(MSE_1_post_C)];
        stdMSEs_1(:,ii,t) = temp;

        temp = [std(MSE_5_post);
            std(MSE_5_post_C0);
            std(MSE_5_post_C)];
        stdMSEs_5(:,ii,t) = temp;
    end
 end


 