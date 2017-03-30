clear all 
close all

save_on = false;

methods = {'Post.','CP0','PCP0','CP10%','PCP10%'};
model = 'agarch11';
v_new = '(R2015a)';
II = 10;
T = 1000;
sigma1 = 1;
sigma2 = 2;
H = 100;

name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_II',num2str(II),'_PCP0_MC_',v_new,'.mat'];
load(name, '-regexp','^mean\w*')
load(name, '-regexp','^std\w*')
load(name, '-regexp','^q\w*')
load(name, '-regexp','^VaR\w*')

%% Boxplots: MEAN DRAWS
ff = figure(100);
set(gcf,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
subplot(2,3,1)
boxplot([mean_draw(:,1),mean_draw_C0(:,1),mean_draw_PC0(:,1),mean_draw_C(:,1),mean_draw_PC(:,1)],...
    'labels' ,methods)
xlabel('\mu')

subplot(2,3,2)
boxplot([mean_draw(:,2),mean_draw_C0(:,2),mean_draw_PC0(:,2),mean_draw_C(:,2),mean_draw_PC(:,2)],...
    'labels' ,methods)
% xlabel('\gamma')
xlabel('\omega')

subplot(2,3,3)
boxplot([mean_draw(:,3),mean_draw_C0(:,3),mean_draw_PC0(:,3),mean_draw_C(:,3),mean_draw_PC(:,3)],...
    'labels' ,methods)
% xlabel('\omega')
xlabel('\mu2')

subplot(2,3,4)
boxplot([mean_draw(:,4),mean_draw_C0(:,4),mean_draw_PC0(:,4),mean_draw_C(:,4),mean_draw_PC(:,4)],...
    'labels' ,methods)
xlabel('\alpha')

subplot(2,3,5)
boxplot([mean_draw(:,5),mean_draw_C0(:,5),mean_draw_PC0(:,5),mean_draw_C(:,5),mean_draw_PC(:,5)],...
    'labels' ,methods)
xlabel('\beta')


subplot(2,3,6)
boxplot([mean_draw(:,5)+mean_draw(:,4),mean_draw_C0(:,5)+mean_draw_C0(:,4),mean_draw_PC0(:,5)+mean_draw_PC0(:,4),...
    mean_draw_C(:,5)+ mean_draw_C(:,4),mean_draw_PC(:,5)+mean_draw_PC(:,4)],...
    'labels' ,methods)
xlabel('\alpha + \beta')
suptitle([model,': Mean draws'])

if save_on 
    name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),...
        '_T',num2str(T),'_H',num2str(H),'_II',num2str(II),'_box_mean.eps'];
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

%% Boxplots: MEDIAN DRAWS
ff = figure(104);
set(gcf,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
subplot(2,3,1)
boxplot([median_draw(:,1),median_draw_C0(:,1),median_draw_PC0(:,1),median_draw_C(:,1),median_draw_PC(:,1)],...
    'labels' ,methods)
xlabel('\mu')

subplot(2,3,2)
boxplot([median_draw(:,2),median_draw_C0(:,2),median_draw_PC0(:,2),median_draw_C(:,2),median_draw_PC(:,2)],...
    'labels' ,methods)
% xlabel('\gamma')
xlabel('\omega')

subplot(2,3,3)
boxplot([median_draw(:,3),median_draw_C0(:,3),median_draw_PC0(:,3),median_draw_C(:,3),median_draw_PC(:,3)],...
    'labels' ,methods)
% xlabel('\omega')
xlabel('\mu2')

subplot(2,3,4)
boxplot([median_draw(:,4),median_draw_C0(:,4),median_draw_PC0(:,4),median_draw_C(:,4),median_draw_PC(:,4)],...
    'labels' ,methods)
xlabel('\alpha')

subplot(2,3,5)
boxplot([median_draw(:,5),median_draw_C0(:,5),median_draw_PC0(:,5),median_draw_C(:,5),median_draw_PC(:,5)],...
    'labels' ,methods)
xlabel('\beta')

suptitle([model,': Median draws'])

if save_on 
    name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),...
        '_T',num2str(T),'_H',num2str(H),'_II',num2str(II),'_box_median.eps'];
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



%% Boxplots: STD DRAWS
ff = figure(102);
set(gcf,'units','normalized','outerposition',[0.1 0.05 0.5 0.6]);
subplot(2,3,1)
boxplot([std_draw(:,1),std_draw_C0(:,1),std_draw_PC0(:,1),std_draw_C(:,1),std_draw_PC(:,1)],...
    'labels' ,methods)
xlabel('\mu')

subplot(2,3,2)
boxplot([std_draw(:,2),std_draw_C0(:,2),std_draw_PC0(:,2),std_draw_C(:,2),std_draw_PC(:,2)],...
    'labels' ,methods)
% xlabel('\gamma')
xlabel('\omega')

subplot(2,3,3)
boxplot([std_draw(:,3),std_draw_C0(:,3),std_draw_PC0(:,3),std_draw_C(:,3),std_draw_PC(:,3)],...
    'labels' ,methods)
% xlabel('\omega')
xlabel('\mu2')

subplot(2,3,4)
boxplot([std_draw(:,4),std_draw_C0(:,4),std_draw_PC0(:,4),std_draw_C(:,4),std_draw_PC(:,4)],...
    'labels' ,methods)
xlabel('\alpha')
suptitle('Std draw')

subplot(2,3,5)
boxplot([std_draw(:,4),std_draw_C0(:,4),std_draw_PC0(:,4),std_draw_C(:,4),std_draw_PC(:,4)],...
    'labels' ,methods)
xlabel('\beta')
suptitle('Std draw')

if save_on 
    name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),...
        '_T',num2str(T),'_H',num2str(H),'_II',num2str(II),'_box_std.eps'];
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



%% MEANS

mean_draw1 = [mean_draw(:,1),mean_draw_C0(:,1),mean_draw_PC0(:,1),mean_draw_C(:,1),mean_draw_PC(:,1)];
mean_draw2 = [mean_draw(:,2),mean_draw_C0(:,2),mean_draw_PC0(:,2),mean_draw_C(:,2),mean_draw_PC(:,2)];
mean_draw3 = [mean_draw(:,3),mean_draw_C0(:,3),mean_draw_PC0(:,3),mean_draw_C(:,3),mean_draw_PC(:,3)];
mean_draw4 = [mean_draw(:,4),mean_draw_C0(:,4),mean_draw_PC0(:,4),mean_draw_C(:,4),mean_draw_PC(:,4)];

ff = figure(31);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

subplot(2,2,1)
plot(mean_draw1)
xlabel('\mu')

subplot(2,2,2)
plot(mean_draw2)
xlabel('\omega')

subplot(2,2,3)
plot(mean_draw3)
xlabel('\alpha')

subplot(2,2,4)
plot(mean_draw4)
xlabel('\beta')

suptitle('Mean draw')


%% MSEs
methods = {'MC','Post.','CP0','PCP0','CP10%','PCP10%'};

MSE1 = [MSE_1,MSE_1_post,MSE_1_post_C0,MSE_1_post_PC0,MSE_1_post_C,MSE_1_post_PC];
MSE5 = [MSE_5,MSE_5_post,MSE_5_post_C0,MSE_5_post_PC0,MSE_5_post_C,MSE_5_post_PC];

ff = figure(98);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

subplot(2,2,1)
boxplot(MSE1,'labels' ,methods)
title('MSE 1')


subplot(2,2,2)
boxplot(MSE5,'labels' ,methods)
title('MSE 5')

subplot(2,2,3)
plot(MSE1)

subplot(2,2,4)
plot(MSE5)
leg = legend({'MSE MC','MSE post','MSE post_C_0','MSE post_P_C_0','MSE post_C','MSE post_P_C'});


%% CV_last
CV_last = cellfun(@(xx) xx(end),CV);
CV_last_C = cellfun(@(xx) xx(end),CV_C);
CV_last_C0 = cellfun(@(xx) xx(end),CV_C0);

ff = figure(56);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
hold all
plot(CV_last)
plot(CV_last_C0)
plot(CV_last_C)
hold off
xlabel('CVs')
leg = legend({'Post', 'C0', 'C'});