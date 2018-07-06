addpath(genpath('include/'));
clear all
close all


model = 'garch11';
% model = 'agarch11';
T = 1000; 
% T = 2500;

file_path = ['results/',model,'/'];
save_on = false;
if strcmp(model,'garch11')
    name = [file_path,model,'_1_2_T',num2str(T),'_H100_II10_PCP0_MC_(R2017a)_low_es_tunning.mat'];
    load(name,'-regexp','^draw','param_true');
    name = [file_path,model,'_1_2_T',num2str(T),'_H100_II10_PCP0_MC_(R2017a)_varc_low_es_tunning.mat'];
    load(name,'-regexp','^draw');
else
    name = [file_path,model,'_1_2_T',num2str(T),'_all_draws.mat'];
    load(name,'-regexp','^draw');
end
ver = version('-release');
    
lw = 1;

ColPalette = [     0    0.4470    0.7410    % deep blue
              0.3010    0.7450    0.9330    % light blue
              0.0000    0.5000    0.0000    % deep green        
              0.4660    0.6740    0.1880];  % light green
          
%% GARCH
if strcmp(model,'garch11')

    ff = figure(11);
    if strcmp(ver,'2015b')
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.8]);    
    else
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.35]);
        % set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.4]);
    end

    % ax2 = axes('Position',[0.43 0.15 0.36 0.8],'Visible','on');
    % axes(ax2)
    subplot(2,3,1)
    bw = 0.0245;
    [f,xi] = ksdensity(draw(:,1));%,'bandwidth',bw);%
    [f_C,xi_C] = ksdensity(draw_C(:,1)); %,'bandwidth',0.4);
    [f_Cm,xi_Cm,bw_Cm] = ksdensity(draw_Cm(:,1),'bandwidth',2*bw);
    [f_PC,xi_PC] = ksdensity(draw_PC(:,1)); %,'bandwidth',0.4);
    [f_PCm,xi_PCm,bw_PCm] = ksdensity(draw_PCm(:,1),'bandwidth',2*bw);

    hold on
    plot(xi,f,'k','linewidth',lw)
    plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
    plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)
    plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
    plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)
    hold off
    if (T == 1000)
        xlim([-0.5 1])
    else
        xlim([-0.25 0.75])
    end
    % P = get(gca,'Position');
    % set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
    % YL = get(gca,'YLim');
    % line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
    set(gca,'TickLabelInterpreter','latex');
    xlabel('\mu','FontSize',12)

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
    hold off
    if (T == 1000)
        xlim([0 20])
    else
        xlim([0.5 3])    
    end
    % P = get(gca,'Position');
    % set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
    % YL = get(gca,'YLim');
    % line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
    set(gca,'TickLabelInterpreter','latex');
    xlabel('\omega','FontSize',12)

    subplot(2,3,4)
    bw = 0.01;
    [f,xi,b] = ksdensity(draw(:,3));%,'bandwidth',bw);%
    [f_C,xi_C,b_C] = ksdensity(draw_C(:,3),'bandwidth',2*bw);
    [f_Cm,xi_Cm,b_Cm] = ksdensity(draw_Cm(:,3),'bandwidth',2*bw);
    [f_PC,xi_PC,b_PC] = ksdensity(draw_PC(:,3),'bandwidth',bw);
    [f_PCm,xi_PCm,b_PCm] = ksdensity(draw_PCm(:,3),'bandwidth',bw);
    hold on
    plot(xi,f,'k','linewidth',lw)
    plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
    plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)
    plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
    plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)
    hold off
    xlim([0 0.6])
    % P = get(gca,'Position');
    % set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
    % YL = get(gca,'YLim');
    % line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
    set(gca,'TickLabelInterpreter','latex');
    xlabel('\alpha','FontSize',12)


    subplot(2,3,5)
    bw = 0.01;
    [f,xi,b] = ksdensity(draw(:,4),'bandwidth',1.5*bw);%
    [f_C,xi_C,b_C] = ksdensity(draw_C(:,4),'bandwidth',2.5*bw);
    [f_Cm,xi_Cm,b_Cm] = ksdensity(draw_Cm(:,4),'bandwidth',3*bw);
    [f_PC,xi_PC,b_PC] = ksdensity(draw_PC(:,4),'bandwidth',1.5*bw);
    [f_PCm,xi_PCm,b_PCm] = ksdensity(draw_PCm(:,4),'bandwidth',1.5*bw);
    hold on
    plot(xi,f,'k','linewidth',lw)
    plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
    plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)
    plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
    plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)
    hold off
    xlim([0 1])
    % P = get(gca,'Position');
    % set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
    % YL = get(gca,'YLim');
    % line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
    set(gca,'TickLabelInterpreter','latex');
    xlabel('\beta','FontSize',12)


    leg = legend('Uncensored','CP','CP$_{var}$','PCP','PCP$_{var}$'); 
    set(leg,'Interpreter','latex','FontSize',12, 'location', ...%,'position',[0.63 0.15 0.23 0.3])
            'eastoutside');%,'Orientation','horizontal');

    if save_on
        name = ['figures/',model,'/',model,'_T',num2str(T),'_ksdens.eps'];
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

%% AGARCH
else 
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
    if (T == 1000)
        bw = 0.005;  
        [f,xi,bw] = ksdensity(draw(:,1));%,'bandwidth',bw);%
        [f_C,xi_C,bw_C] = ksdensity(draw_C(:,1),'bandwidth',0.06);
        [f_Cm,xi_Cm,bw_Cm] = ksdensity(draw_Cm(:,1),'bandwidth',0.07);
        [f_PC,xi_PC,bw_PC] = ksdensity(draw_PC(:,1),'bandwidth',0.04);
        [f_PCm,xi_PCm,bw_PCm] = ksdensity(draw_PCm(:,1),'bandwidth',0.04);
    else
        bw = 0.005;  
        [f,xi,bw] = ksdensity(draw(:,1));%,'bandwidth',bw);%
        [f_C,xi_C,bw_C] = ksdensity(draw_C(:,1),'bandwidth',2*0.021);
        [f_Cm,xi_Cm,bw_Cm] = ksdensity(draw_Cm(:,1),'bandwidth',2*0.017);
        [f_PC,xi_PC,bw_PC] = ksdensity(draw_PC(:,1),'bandwidth',2*0.019);
        [f_PCm,xi_PCm,bw_PCm] = ksdensity(draw_PCm(:,1),'bandwidth',2*0.016);
    end
    hold on
    plot(xi,f,'k','linewidth',lw)
    plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
    plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)
    plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
    plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)
    hold off
    % P = get(gca,'Position');
    % set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
    if (T == 1000)
        xlim([-0.5 1.50])
    else
        xlim([-0.10 0.8])
    end
    % YL = get(gca,'YLim');
    % line([param_true(1) param_true(1)], YL,'Color','r','LineWidth',lw); 
    set(gca,'TickLabelInterpreter','latex');
    xlabel('$\mu_{1}$','FontSize',12,'Interpreter','latex');


  subplot(2,3,2)
    if (T == 1000)
        bw = 0.05;
        [f,xi,bw] = ksdensity(draw(:,2),'bandwidth',1.1*0.0213);%
        [f_C,xi_C,bw_C] = ksdensity(draw_C(:,2),'bandwidth',1.5*0.1986);
        [f_Cm,xi_Cm,bw_Cm] = ksdensity(draw_Cm(:,2),'bandwidth',2*0.4564);
        [f_PC,xi_PC,bw_PC] = ksdensity(draw_PC(:,2),'bandwidth',1.5*0.1265);
        [f_PCm,xi_PCm,bw_PCm] = ksdensity(draw_PCm(:,2),'bandwidth',1.5*0.1492);
    else
        [f,xi,bw] = ksdensity(draw(:,2));%,'bandwidth',1.1*0.0213);%
        [f_C,xi_C,bw_C] = ksdensity(draw_C(:,2),'bandwidth',2.5*0.11);
        [f_Cm,xi_Cm,bw_Cm] = ksdensity(draw_Cm(:,2),'bandwidth',2.5*0.11);
        [f_PC,xi_PC,bw_PC] = ksdensity(draw_PC(:,2),'bandwidth',2.5*0.07);
        [f_PCm,xi_PCm,bw_PCm] = ksdensity(draw_PCm(:,2),'bandwidth',2.5*0.06);        
    end
    hold on
    plot(xi,f,'k','linewidth',lw)
    plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
    plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)
    plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
    plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)
    hold off
    if (T == 1000)
        xlim([0 20])
    else
        xlim([0 7])    
    end
    % P = get(gca,'Position');
    % set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
    % YL = get(gca,'YLim');
    % line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
    set(gca,'TickLabelInterpreter','latex');
    xlabel('\omega','FontSize',12)
    
    
    % ax2 = axes('Position',[0.43 0.15 0.36 0.8],'Visible','on');
    % axes(ax2)
    subplot(2,3,4)
    bw = 0.0245;
    if (T == 1000)
        [f,xi] = ksdensity(draw(:,3),'bandwidth',2.5*bw);%
        [f_C,xi_C,bw_C] = ksdensity(draw_C(:,3),'bandwidth',7*bw);
        [f_Cm,xi_Cm,bw_Cm] = ksdensity(draw_Cm(:,3),'bandwidth',6*bw);
        [f_PC,xi_PC,bw_PC] = ksdensity(draw_PC(:,3),'bandwidth',2.5*bw);
        [f_PCm,xi_PCm,bw_PCm] = ksdensity(draw_PCm(:,3),'bandwidth',2.5*bw);
    else
        [f,xi,bw] = ksdensity(draw(:,3),'bandwidth',1.5*0.017);%
        [f_C,xi_C,bw_C] = ksdensity(draw_C(:,3),'bandwidth',2*0.03);
        [f_Cm,xi_Cm,bw_Cm] = ksdensity(draw_Cm(:,3),'bandwidth',2*0.03);
        [f_PC,xi_PC,bw_PC] = ksdensity(draw_PC(:,3),'bandwidth',1*0.017);
        [f_PCm,xi_PCm,bw_PCm] = ksdensity(draw_PCm(:,3),'bandwidth',1*0.017);      
    end
    hold on
    plot(xi,f,'k','linewidth',lw)
    plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
    plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)
    plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
    plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)
    hold off
    if (T == 1000)
        xlim([-1.5 2.0])
    else
%         xlim([-0.25 0.75])
    end
    % P = get(gca,'Position');
    % set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
    % YL = get(gca,'YLim');
    % line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
    set(gca,'TickLabelInterpreter','latex');
    xlabel('$\mu_{2}$','FontSize',12,'Interpreter','latex')

    
  
    subplot(2,3,5)
    if (T == 1000)
        %     bw = 0.01;
        [f,xi,bw] = ksdensity(draw(:,4));%,'bandwidth',bw);%
        [f_C,xi_C,bw_C] = ksdensity(draw_C(:,4),'bandwidth',3*0.018);
        [f_Cm,xi_Cm,bw_Cm] = ksdensity(draw_Cm(:,4),'bandwidth',3*0.0258);
        [f_PC,xi_PC,bw_PC] = ksdensity(draw_PC(:,4));%,'bandwidth',bw);
        [f_PCm,xi_PCm,bw_PCm] = ksdensity(draw_PCm(:,4));%,'bandwidth',bw);
    else
        %     bw = 0.01;
        [f,xi,bw] = ksdensity(draw(:,4));%,'bandwidth',bw);%
        [f_C,xi_C,bw_C] = ksdensity(draw_C(:,4),'bandwidth',3*0.0127);
        [f_Cm,xi_Cm,bw_Cm] = ksdensity(draw_Cm(:,4),'bandwidth',2*0.0118);
        [f_PC,xi_PC,bw_PC] = ksdensity(draw_PC(:,4));%,'bandwidth',bw);
        [f_PCm,xi_PCm,bw_PCm] = ksdensity(draw_PCm(:,4));%,'bandwidth',bw);
    end
    hold on
    plot(xi,f,'k','linewidth',lw)
    plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
    plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)
    plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
    plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)
    hold off
    xlim([0 1.0])
    % P = get(gca,'Position');
    % set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
    % YL = get(gca,'YLim');
    % line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
    set(gca,'TickLabelInterpreter','latex');
    xlabel('\alpha','FontSize',12)



    subplot(2,3,6)
    if (T == 1000)
        bw = 0.01;
        [f,xi,bw] = ksdensity(draw(:,5));%,'bandwidth',1.5*bw);%
        [f_C,xi_C,bw_C] = ksdensity(draw_C(:,5),'bandwidth',1.5*0.0342);
        [f_Cm,xi_Cm,bw_Cm] = ksdensity(draw_Cm(:,5),'bandwidth',2.5*0.028);
        [f_PC,xi_PC,bw_PC] = ksdensity(draw_PC(:,5),'bandwidth',2*0.0095);
        [f_PCm,xi_PCm,bw_PCm] = ksdensity(draw_PCm(:,5),'bandwidth',2*0.0095);
    else
        bw = 0.01;
        [f,xi,bw] = ksdensity(draw(:,5));%,'bandwidth',1.5*bw);%
        [f_C,xi_C,bw_C] = ksdensity(draw_C(:,5),'bandwidth',2*0.0174);
        [f_Cm,xi_Cm,bw_Cm] = ksdensity(draw_Cm(:,5),'bandwidth',2*0.0175);
        [f_PC,xi_PC,bw_PC] = ksdensity(draw_PC(:,5));%,'bandwidth',2*0.0095);
        [f_PCm,xi_PCm,bw_PCm] = ksdensity(draw_PCm(:,5));%,'bandwidth',2*0.0095);
    end
    hold on
    plot(xi,f,'k','linewidth',lw)
    plot(xi_C,f_C,'color',ColPalette(1,:) ,'linewidth',lw)
    plot(xi_Cm,f_Cm,'color',ColPalette(2,:),'linewidth',lw)
    plot(xi_PC,f_PC,'color',ColPalette(3,:) ,'linewidth',lw)
    plot(xi_PCm,f_PCm,'color',ColPalette(4,:),'linewidth',lw)
    hold off
    xlim([0 1])
    % P = get(gca,'Position');
    % set(gca,'Position',[P(1) P(2)+0.15 P(3) P(4)-0.15])
    % YL = get(gca,'YLim');
    % line([param_true(2) param_true(2)], YL,'Color','r','LineWidth',lw); 
    set(gca,'TickLabelInterpreter','latex');
    xlabel('\beta','FontSize',12)


    leg = legend('Uncensored','CP','CP$_{var}$','PCP','PCP$_{var}$'); 
    set(leg,'Interpreter','latex','FontSize',12, 'location', ...%,'position',[0.63 0.15 0.23 0.3])
            'eastoutside');%,'Orientation','horizontal');

    if save_on
        name = ['figures/',model,'/',model,'_T',num2str(T),'_ksdens.eps'];
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