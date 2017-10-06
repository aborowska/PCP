ff = figure(101);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.65 0.5]);
% set(gcf,'units','normalized','outerposition',[0.1 0.1 0.66 0.4]);
set(gcf,'defaulttextinterpreter','latex');


ii = 1;
subplot(2,3,[1,4])
hold on
histogram(draw(:,ii))
histogram(draw_PC(:,ii))
histogram(draw_PCm(:,ii))
%     histogram(draw_PC2(:,ii))
hold off
plotTickLatex2D('FontSize',11);
title(['Draws of ',params{ii}])
ll = legend('Posterior','PCP10\%','PCPm');
set(ll,'interpreter','latex','FontSize',11)

ii = 2;
subplot(2,3,[2,5])
hold on
histogram(draw(:,ii))
histogram(draw_PC(:,ii))
histogram(draw_PCm(:,ii))
%     histogram(draw_PC2(:,ii))
hold off
plotTickLatex2D('FontSize',11);
title(['Draws of ',params{ii}])

ii = 3;
subplot(2,3,3)
hold on
histogram(draw(:,ii))
histogram(draw_PCm(:,ii),'FaceColor',  [0.9290    0.6940    0.1250])

% histogram(draw_PC(:,ii))
%     histogram(draw_PC2(:,ii))
hold off
plotTickLatex2D('FontSize',11);
title(['Draws of ',params{ii}])


subplot(2,3,6)
hold on
% histogram(draw(:,ii))
histogram(draw_PC(:,ii),'FaceColor',  [0.8500    0.3250    0.0980])
%     histogram(draw_PC2(:,ii))
hold off
plotTickLatex2D('FontSize',11);
title(['Draws of ',params{ii}])


suptitle('Parameters from $\theta_2$')

if save_on 
    name_plot = [figures_path,'/Hist_theta2_',data_name,'_.eps'];
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
            print(ff,name_plot,'-depsc','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end

