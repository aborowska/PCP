%% h_true vs h
figure(2)
set(gcf,'units','normalized','outerposition',[0.0 0.00 0.5 1]);

subplot(2,2,1)
ii = 1;
    hold on
    plot(h_true(T+1:T+H),'r')
    plot(h_post(ii,:),'b')    
%     plot(h_post_C(ii,:),'g')    
%     plot(h_post_PC(ii,:),'k')
    hold off
    params = ['draw: ',sprintf('%5.3f ',draw(ii,:))];
    t = text(0.1,1,params,'units','normalized');
    t.Color='b';    
%     params = ['draw_PC: ',sprintf('%5.3f ',draw_PC(ii,:))];
%     t = text(0.1,0.9,params,'units','normalized');
%     t.Color='g';
%     params = ['draw_PC: ',sprintf('%5.3f ',draw_PC(ii,:))];
%     t = text(0.1,0.8,params,'units','normalized');
%     t.Color='k';
    
subplot(2,2,2)
ii = 2501;
    hold on
    plot(h_true(T+1:T+H),'r')
    plot(h_post(ii,:),'b')    
%     plot(h_post_C(ii,:),'g')    
%     plot(h_post_PC(ii,:),'k')
    hold off
    params = ['draw: ',sprintf('%5.3f ',draw(ii,:))];
    t = text(0.1,1,params,'units','normalized');
    t.Color='b';    
%     params = ['draw_PC: ',sprintf('%5.3f ',draw_PC(ii,:))];
%     t = text(0.1,0.9,params,'units','normalized');
%     t.Color='g';
%     params = ['draw_PC: ',sprintf('%5.3f ',draw_PC(ii,:))];
%     t = text(0.1,0.8,params,'units','normalized');
%     t.Color='k';

subplot(2,2,3)
ii = 5001;
    hold on
    plot(h_true(T+1:T+H),'r')
    plot(h_post(ii,:),'b')    
%     plot(h_post_C(ii,:),'g')    
%     plot(h_post_PC(ii,:),'k')
    hold off
    params = ['draw: ',sprintf('%5.3f ',draw(ii,:))];
    t = text(0.1,1,params,'units','normalized');
    t.Color='b';    
%     params = ['draw_PC: ',sprintf('%5.3f ',draw_PC(ii,:))];
%     t = text(0.1,0.9,params,'units','normalized');
%     t.Color='g';
%     params = ['draw_PC: ',sprintf('%5.3f ',draw_PC(ii,:))];
%     t = text(0.1,0.8,params,'units','normalized');
%     t.Color='k';

subplot(2,2,4)
ii = 7501;
    hold on
    plot(h_true(T+1:T+H),'r')
    plot(h_post(ii,:),'b')    
%     plot(h_post_C(ii,:),'g')    
%     plot(h_post_PC(ii,:),'k')
    hold off
    params = ['draw: ',sprintf('%5.3f ',draw(ii,:))];
    t = text(0.1,1,params,'units','normalized');
    t.Color='b';    
%     params = ['draw_PC: ',sprintf('%5.3f ',draw_PC(ii,:))];
%     t = text(0.1,0.9,params,'units','normalized');
%     t.Color='g';
%     params = ['draw_PC: ',sprintf('%5.3f ',draw_PC(ii,:))];
%     t = text(0.1,0.8,params,'units','normalized');
%     t.Color='k';
    
leg = legend({'h true','h post','h post_C','h post_P_C'},...
    'Position',[0.4505 0.4599 0.1091 0.0924],'units', 'normalized');
suptitle([model, '\sigma_2 = ',num2str(sigma2)])

%% q vs. VaR  single
ff = figure(1);
set(gcf,'units','normalized','outerposition',[0.0 0.00 1 1]);

ii = 1;
    hold on
    plot(q5(ii,:),'r')
    plot(VaR_5(ii,:),'m')
    plot(VaR_5_post(ii,:),'b')
    plot(VaR_5_post_C(ii,:),'g')
    plot(VaR_5_post_PC(ii,:),'k')
    
    plot(q1(ii,:),'r:')
    plot(VaR_1(ii,:),'m:')
    plot(VaR_1_post(ii,:),'b:')
    plot(VaR_1_post_C(ii,:),'g:')
%     plot(VaR_1_post_PC(ii,:),'k:')    
    hold off
leg = legend({'q (- 5%, : 1%)','VaR true','VaR post','VaR post_C','VaR post_P_C'},...
    'Position',[0.8054 0.1318 0.1829 0.1815],'units', 'normalized');
suptitle([model, ' \sigma_2 = ',num2str(sigma2)])

name = ['figures/',model,'/',model,'_',num2str(sigma1),'_',num2str(sigma2),'_',num2str(T),'q_VaR_single.eps'];
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


%% q vs. VaR 
figure(10)
set(gcf,'units','normalized','outerposition',[0.1 0.05 0.9 0.9]);

subplot(2,2,1)
ii = 3;%1;
    hold on
    plot(q5(ii,:),'r')
    plot(VaR_5(ii,:),'m')
    plot(VaR_5_post(ii,:),'b')
    plot(VaR_5_post_C(ii,:),'g')
    plot(VaR_5_post_PC(ii,:),'k')
    
    plot(q1(ii,:),'r:')
    plot(VaR_1(ii,:),'m:')
    plot(VaR_1_post(ii,:),'b:')
    plot(VaR_1_post_C(ii,:),'g:')
    plot(VaR_1_post_PC(ii,:),'K:')    
    hold off

subplot(2,2,2)
ii = 32;%10;
    hold on
    plot(q5(ii,:),'r')
    plot(VaR_5(ii,:),'m')
    plot(VaR_5_post(ii,:),'b')
    plot(VaR_5_post_C(ii,:),'g')
    plot(VaR_5_post_PC(ii,:),'k')
    
    plot(q1(ii,:),'r:')
    plot(VaR_1(ii,:),'m:')
    plot(VaR_1_post(ii,:),'b:')
    plot(VaR_1_post_C(ii,:),'g:')
    plot(VaR_1_post_PC(ii,:),'K:')    
    hold off

subplot(2,2,3)
 ii = 17;%30;
    hold on
    plot(q5(ii,:),'r')
    plot(VaR_5(ii,:),'m')
    plot(VaR_5_post(ii,:),'b')
    plot(VaR_5_post_C(ii,:),'g')
    plot(VaR_5_post_PC(ii,:),'k')
    
    plot(q1(ii,:),'r:')
    plot(VaR_1(ii,:),'m:')
    plot(VaR_1_post(ii,:),'b:')
    plot(VaR_1_post_C(ii,:),'g:')
    plot(VaR_1_post_PC(ii,:),'K:')    
    hold off

subplot(2,2,4)
ii = 14;%50;
    hold on
    plot(q5(ii,:),'r')
    plot(VaR_5(ii,:),'m')
    plot(VaR_5_post(ii,:),'b')
    plot(VaR_5_post_C(ii,:),'g')
    plot(VaR_5_post_PC(ii,:),'k')
    
    plot(q1(ii,:),'r:')
    plot(VaR_1(ii,:),'m:')
    plot(VaR_1_post(ii,:),'b:')
    plot(VaR_1_post_C(ii,:),'g:')
    plot(VaR_1_post_PC(ii,:),'K:')    
    hold off

leg = legend({'q (- 5%, : 1%)','VaR true','VaR post','VaR post_C','VaR post_P_C'},...
    'Position',[0.8054 0.1318 0.1829 0.1815],'units', 'normalized');
suptitle([model, ' \sigma_2 = ',num2str(sigma2)])
    