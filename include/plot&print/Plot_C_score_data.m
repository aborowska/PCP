function data_stats = Plot_C_score_data(C_score_diff,y,H, time,save_on,figures_path)
% garch_old: time = [1998,2008]
% gas_crisis: time=[2005, 2016.5]
    close all

    T = length(y);
    in_sample = 1:(T-H);
    out_sample = (T-H+1):T;

    y_out = y(out_sample);
    y_in = y(in_sample);

    xx = linspace(time(1,1),time(1,2),length(y));
    xx_out = xx(out_sample);
    time_out = [xx_out(1),xx_out(end)];

    ff = figure(1);
    set(gcf,'units','normalized')%,'outerposition',[0.1 0.1 0.49 0.3]);
    % set(gcf,'units','normalized','outerposition',[0.1 0.1 0.66 0.4]);
    set(gcf,'defaulttextinterpreter','latex');

    [hAx,hLine1,hLine2] = plotyy(xx_out,y_out',xx_out,C_score_diff);
    hLine1.Color = [0.3,0.3,0.3];
    hLine2.Color = [1, 0, 0]; 
    xlim(hAx(1),time_out)
    xlim(hAx(2),time_out)
    plotTickLatex2D('FontSize',11);

    title('Difference in CSR and out-of-sample data')
    
    if save_on
%         name = [figures_path,model,'_hor_direct_H', num2str(H),'.png'];
        name_fig = [figures_path,'C_score_data.eps'];
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
end