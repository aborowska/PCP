%% CSR: const 5\%
ff = figure(232);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');

subplot(2,3,1)
plot(C_score_post_5 - C_score_PC_5)
set(gca,'TickLabelInterpreter','latex');
title('Post - PCP')

subplot(2,3,4)
plot(C_score_post_5 - C_score_PCm_5)
set(gca,'TickLabelInterpreter','latex');
title('Post - PCP$_{var}$')

subplot(2,3,2)
plot(C_score_C_5 - C_score_PC_5)
set(gca,'TickLabelInterpreter','latex');
title('CP - PCP')

subplot(2,3,5)
plot(C_score_C_5 - C_score_PCm_5)
set(gca,'TickLabelInterpreter','latex');
title('CP - PCP$_{var}$')

subplot(2,3,3)
plot(C_score_Cm_5 - C_score_PC_5)
set(gca,'TickLabelInterpreter','latex');
title('CP$_{var}$ - PCP')

subplot(2,3,6)
plot(C_score_Cm_5 - C_score_PCm_5)
set(gca,'TickLabelInterpreter','latex');
title('CP$_{var}$  - PCP$_{var}$')

suptitle('CSR: const 5\%')


if save_on
    name_fig = [figures_path,data_name,'_C_diffs_5.eps'];
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
            print(ff,name_fig,'-depsc','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end

%% CSR: const 1\%
ff = figure(233);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');

subplot(2,3,1)
plot(C_score_post_1 - C_score_PC_1)
set(gca,'TickLabelInterpreter','latex');
title('Post - PCP')

subplot(2,3,4)
plot(C_score_post_1 - C_score_PCm_1)
set(gca,'TickLabelInterpreter','latex');
title('Post - PCP$_{var}$')

subplot(2,3,2)
plot(C_score_C_1 - C_score_PC_1)
set(gca,'TickLabelInterpreter','latex');
title('CP - PCP')

subplot(2,3,5)
plot(C_score_C_1 - C_score_PCm_1)
set(gca,'TickLabelInterpreter','latex');
title('CP - PCP$_{var}$')

subplot(2,3,3)
plot(C_score_Cm_1 - C_score_PC_1)
set(gca,'TickLabelInterpreter','latex');
title('CP$_{var}$ - PCP')

subplot(2,3,6)
plot(C_score_Cm_1 - C_score_PCm_1)
set(gca,'TickLabelInterpreter','latex');
title('CP$_{var}$  - PCP$_{var}$')

suptitle('CSR: const 1\%')


if save_on
    name_fig = [figures_path,data_name,'_C_diffs_1.eps'];
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
            print(ff,name_fig,'-depsc','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end

%% CSR: const 0.5\%
ff = figure(234);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');

subplot(2,3,1)
plot(C_score_post_05 - C_score_PC_05)
set(gca,'TickLabelInterpreter','latex');
% title('Post - PC; CSR: const 0.5\%')
title('Post - PCP')

subplot(2,3,4)
plot(C_score_post_05 - C_score_PCm_05)
set(gca,'TickLabelInterpreter','latex');
% title('Post - PCm; CSR: const 0.5\%')
title('Post - PCP$_{var}$')

subplot(2,3,2)
plot(C_score_C_05 - C_score_PC_05)
set(gca,'TickLabelInterpreter','latex');
% title('C - PC; CSR: const 0.5\%')
title('CP - PCP')

subplot(2,3,5)
plot(C_score_C_05 - C_score_PCm_05)
set(gca,'TickLabelInterpreter','latex');
% title('C - PCm; CSR: const 0.5\%')
title('CP - PCP$_{var}$')

subplot(2,3,3)
plot(C_score_Cm_05 - C_score_PC_05)
set(gca,'TickLabelInterpreter','latex');
% title('Cm - PC; CSR: const 0.5\%')
title('CP$_{var}$ - PCP')

subplot(2,3,6)
plot(C_score_Cm_05 - C_score_PCm_05)
set(gca,'TickLabelInterpreter','latex');
% title('Cm - PCm; CSR: const 0.5\%')
title('CP$_{var}$  - PCP$_{var}$')

suptitle('CSR: const 0.5\%')


if save_on
    name_fig = [figures_path,data_name,'_C_diffs_05.eps'];
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
            print(ff,name_fig,'-depsc','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end

 
%% %%%%%%%%%%%%%%%%%%%%%%

 %% CSR: var 5\%
ff = figure(332);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');

subplot(2,3,1)
plot(Cv_score_post_5 - Cv_score_PC_5)
set(gca,'TickLabelInterpreter','latex');
title('Post - PCP')

subplot(2,3,4)
plot(Cv_score_post_5 - Cv_score_PCm_5)
set(gca,'TickLabelInterpreter','latex');
title('Post - PCP$_{var}$')

subplot(2,3,2)
plot(Cv_score_C_5 - Cv_score_PC_5)
set(gca,'TickLabelInterpreter','latex');
title('CP - PCP')

subplot(2,3,5)
plot(Cv_score_C_5 - Cv_score_PCm_5)
set(gca,'TickLabelInterpreter','latex');
title('CP - PCP$_{var}$')

subplot(2,3,3)
plot(Cv_score_Cm_5 - Cv_score_PC_5)
set(gca,'TickLabelInterpreter','latex');
title('CP$_{var}$ - PCP')

subplot(2,3,6)
plot(Cv_score_Cm_5 - Cv_score_PCm_5)
set(gca,'TickLabelInterpreter','latex');
title('CP$_{var}$  - PCP$_{var}$')

suptitle('CSR: var 5\%')


if save_on
    name_fig = [figures_path,data_name,'_Cv_diffs_5.eps'];
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
            print(ff,name_fig,'-depsc','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end

%% CSR: var 1\%
ff = figure(333);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');

subplot(2,3,1)
plot(Cv_score_post_1 - Cv_score_PC_1)
set(gca,'TickLabelInterpreter','latex');
title('Post - PCP')

subplot(2,3,4)
plot(Cv_score_post_1 - Cv_score_PCm_1)
set(gca,'TickLabelInterpreter','latex');
title('Post - PCP$_{var}$')

subplot(2,3,2)
plot(Cv_score_C_1 - Cv_score_PC_1)
set(gca,'TickLabelInterpreter','latex');
title('CP - PCP')

subplot(2,3,5)
plot(Cv_score_C_1 - Cv_score_PCm_1)
set(gca,'TickLabelInterpreter','latex');
title('CP - PCP$_{var}$')

subplot(2,3,3)
plot(Cv_score_Cm_1 - Cv_score_PC_1)
set(gca,'TickLabelInterpreter','latex');
title('CP$_{var}$ - PCP')

subplot(2,3,6)
plot(Cv_score_Cm_1 - Cv_score_PCm_1)
set(gca,'TickLabelInterpreter','latex');
title('CP$_{var}$  - PCP$_{var}$')

suptitle('CSR: var 1\%')


if save_on
    name_fig = [figures_path,data_name,'_Cv_diffs_1.eps'];
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
            print(ff,name_fig,'-depsc','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end

%% CSR: var 0.5\%
ff = figure(334);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');

subplot(2,3,1)
plot(Cv_score_post_05 - Cv_score_PC_05)
set(gca,'TickLabelInterpreter','latex');
title('Post - PCP')

subplot(2,3,4)
plot(Cv_score_post_05 - Cv_score_PCm_05)
set(gca,'TickLabelInterpreter','latex');
title('Post - PCP$_{var}$')

subplot(2,3,2)
plot(Cv_score_C_05 - Cv_score_PC_05)
set(gca,'TickLabelInterpreter','latex');
title('CP - PCP')

subplot(2,3,5)
plot(Cv_score_C_05 - Cv_score_PCm_05)
set(gca,'TickLabelInterpreter','latex');
title('CP - PCP$_{var}$')

subplot(2,3,3)
plot(Cv_score_Cm_05 - Cv_score_PC_05)
set(gca,'TickLabelInterpreter','latex');
title('CP$_{var}$ - PCP')

subplot(2,3,6)
plot(Cv_score_Cm_05 - Cv_score_PCm_05)
set(gca,'TickLabelInterpreter','latex');
title('CP$_{var}$  - PCP$_{var}$')

suptitle('CSR: var 0.5\%')


if save_on
    name_fig = [figures_path,data_name,'_Cv_diffs_05.eps'];
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
            print(ff,name_fig,'-depsc','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end

 

 



%%
figure(456)
subplot(2,3,1)
plot(C_score_post_5 - C_score_PC_5)
title('Post - PCP; CSR: const 5%')

subplot(2,3,4)
plot(Cv_score_post_5 - Cv_score_PC_5)
title('Post - PCP; CSR: var 5%')

subplot(2,3,2)
plot(C_score_C_5 - C_score_PC_5)
title('CP - PCP; CSR: const 5%')

subplot(2,3,5)
plot(Cv_score_C_5 - Cv_score_PC_5)
title('CP - PCP; CSR: var 5%')

subplot(2,3,3)
plot(C_score_C_5 - C_score_PC_5)
title('CP - PCP; CSR: const 5%')

subplot(2,3,6)
plot(Cv_score_C_5 - Cv_score_PC_5)
title('CP - PCP; CSR: var 5%')



%% Slides/Poster
ff = figure(567)
set(gcf,'defaulttextinterpreter','latex');

subplot(2,3,1)
plot(1:H,C_score_post_05 - C_score_PC_05,'linewidth',2)
xlim([1,H])
title('Post - PC; CSR: const 0.5\%','interpreter','latex')
plotTickLatex2D('FontSize',11);

subplot(2,3,4)
plot(1:H,Cv_score_post_05 - Cv_score_PC_05,'linewidth',2)
xlim([1,H])
title('Post - PC; CSR: var 0.5\%','interpreter','latex')
plotTickLatex2D('FontSize',11);

subplot(2,3,2)
plot(1:H,C_score_C_05 - C_score_PC_05,'linewidth',2)
xlim([1,H])
title('C - PC; CSR: const 0.5\%','interpreter','latex')
plotTickLatex2D('FontSize',11);

subplot(2,3,5)
plot(1:H,Cv_score_C_05 - Cv_score_PC_05,'linewidth',2)
xlim([1,H])
title('C - PC; CSR: var 0.5\%','interpreter','latex')
plotTickLatex2D('FontSize',11);

subplot(2,3,3)
plot(1:H,C_score_C_05 - C_score_PC_05,'linewidth',2)
xlim([1,H])
title('C - PC; CSR: const 0.5\%','interpreter','latex')
plotTickLatex2D('FontSize',11);

subplot(2,3,6)
plot(1:H,Cv_score_C_05 - Cv_score_PC_05,'linewidth',2)
xlim([1,H])
title('C - PC; CSR: var 0.5\%','interpreter','latex')
plotTickLatex2D('FontSize',11);


if save_on
%         name = [figures_path,model,'_hor_direct_H', num2str(H),'.png'];
    name_fig = [figures_path,data_name,'_C_diffs.eps'];
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try 
%                 print(name,'-dpng','-r0')
%                 print(name,'-depsc','-r0')                    
            print(ff,name_fig,'-depsc','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end