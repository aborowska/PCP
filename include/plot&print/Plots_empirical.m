
if plot_on
    figure(1)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');

    hold on
    histogram(draw(:,1))
    histogram(draw_PC(:,1))
    histogram(draw_PCm(:,1))
    hold off 
    leg = legend('Posterior','PC 10%','PC var MLE');   
    set(leg,'Interpreter','latex','FontSize',10)
    title('Draws of \nu in GARCH-t')
end




if plot_on   
    ff = figure(21);
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    subplot(1,3,1)
        plot(C_score_post_05 - C_score_PC_05)
        title('q=0.005\%','interpreter','latex')
    subplot(1,3,2)
        plot(C_score_post_1 - C_score_PC_1)
        title('q=0.01\%','interpreter','latex')
    subplot(1,3,3)
        plot(C_score_post_5 - C_score_PC_5)
        title('q=0.05\%','interpreter','latex')
    suptitle(['[',data_name,'] Differences in censored (at q) SR: posterior vs. PCP (part.4)'])
    

    if save_on 
        name_plot = [figures_path,'/C_score_diff_PCP(p4)_',data_name,'_.eps'];
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

    ff = figure(22)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    subplot(1,3,1)
        plot(C_score_post_05 - C_score_PC_05)
        title('q=0.005\%','interpreter','latex')
    subplot(1,3,2)
        plot(C_score_post_1 - C_score_PC_1)
        title('q=0.01\%','interpreter','latex')
    subplot(1,3,3)
        plot(C_score_post_5 - C_score_PC_5)
        title('q=0.05\%','interpreter','latex')
    suptitle(['[',data_name,'] Differences in censored (at q) SR: posterior vs. PCP (part.3)'])
    
    if save_on 
        name_plot = [figures_path,'/C_score_diff_PCP(p3)_',data_name,'_.eps'];
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
end


ff = figure(88);
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    COL = get(gca,'colororder');

    hold on
    plot(THR_MLE')
    set(gca, 'ColorOrder', COL(1:3,:))
    plot(repmat(THR_emp,1,H)','--')
    hold off
    leg = legend('MLE thr 5\%','MLE thr 1\%', 'MLE thr 0.5\%',...
        'EMP thr 5\%','EMP thr 1\%', 'EMP thr 0.5\%');
    set(leg,'Interpreter','latex','FontSize',10,'location','southwest');
    title([data_name,' - evaluation thresholds'],'interpreter','latex')
    if save_on 
        name_plot = [figures_path,'/Thr_MLE_emp_',data_name,'_.eps'];
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
