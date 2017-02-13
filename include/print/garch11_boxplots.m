methods = {'Post.','CP0','PCP0','CP10%','PCP10%'};
model = 'garch11';
v_new = '(R2015a)';
%     II = 100;
name = ['results/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_T',num2str(T),'_H',num2str(H),'_II',num2str(II),'_PCP0_MC2_',v_new,'.mat'];
load(name, '-regexp','^mean\w*')
load(name, '-regexp','^std\w*')
load(name, '-regexp','^q\w*')
load(name, '-regexp','^VaR\w*')

%% Boxplots: MEAN DRAWS
figure(100)
set(gcf,'units','normalized','outerposition',[0.1 0.05 0.5 0.6]);
subplot(2,2,1)
boxplot([mean_draw(:,1),mean_draw_C0(:,1),mean_draw_PC0(:,1),mean_draw_C(:,1),mean_draw_PC(:,1)],...
    'labels' ,methods)
xlabel('\mu')

subplot(2,2,2)
boxplot([mean_draw(:,2),mean_draw_C0(:,2),mean_draw_PC0(:,2),mean_draw_C(:,2),mean_draw_PC(:,2)],...
    'labels' ,methods)
xlabel('\omega')

subplot(2,2,3)
boxplot([mean_draw(:,3),mean_draw_C0(:,3),mean_draw_PC0(:,3),mean_draw_C(:,3),mean_draw_PC(:,3)],...
    'labels' ,methods)
xlabel('\alpha')

subplot(2,2,4)
boxplot([mean_draw(:,4),mean_draw_C0(:,4),mean_draw_PC0(:,4),mean_draw_C(:,4),mean_draw_PC(:,4)],...
    'labels' ,methods)
xlabel('\beta')
suptitle('Mean draw')


%% Boxplots: STD DRAWS
figure(102)
set(gcf,'units','normalized','outerposition',[0.1 0.05 0.5 0.6]);
subplot(2,2,1)
boxplot([std_draw(:,1),std_draw_C0(:,1),std_draw_PC0(:,1),std_draw_C(:,1),std_draw_PC(:,1)],...
    'labels' ,methods)
xlabel('\mu')

subplot(2,2,2)
boxplot([std_draw(:,2),std_draw_C0(:,2),std_draw_PC0(:,2),std_draw_C(:,2),std_draw_PC(:,2)],...
    'labels' ,methods)
xlabel('\omega')

subplot(2,2,3)
boxplot([std_draw(:,3),std_draw_C0(:,3),std_draw_PC0(:,3),std_draw_C(:,3),std_draw_PC(:,3)],...
    'labels' ,methods)
xlabel('\alpha')

subplot(2,2,4)
boxplot([std_draw(:,4),std_draw_C0(:,4),std_draw_PC0(:,4),std_draw_C(:,4),std_draw_PC(:,4)],...
    'labels' ,methods)
xlabel('\beta')
suptitle('Std draw')