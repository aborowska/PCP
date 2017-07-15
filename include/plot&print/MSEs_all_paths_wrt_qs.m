figure(123)

subplot(2,2,1)
title('P')
hold on
plot(h_true(T+1:T+H),'r')
plot(q1,'r')
plot(q5,'r')
plot(y_post((1:10:250)',:)','b:')
plot(y_post((100)',:)','k')
plot(y_post((250)',:)','k')

subplot(2,2,2)
title('P add')
hold on
plot(h_true(T+1:T+H),'r')
plot(q1,'r')
plot(q5,'r')
plot(y_post_add((1:10:250)',:)','b:')
plot(y_post_add((100)',:)','k')
plot(y_post_add((250)',:)','k')

subplot(2,2,3)
title('CP')
hold on
plot(h_true(T+1:T+H),'r')
plot(q1,'r')
plot(q5,'r')
plot(y_post_C((1:10:250)',:)','b:')
plot(y_post_C((100)',:)','k')
plot(y_post_C((250)',:)','k')


subplot(2,2,4)
title('PCP')
hold on
plot(h_true(T+1:T+H),'r')
plot(q1,'r')
plot(q5,'r')
plot(y_post_PC((1:10:250)',:)','b:')
plot(y_post_PC((100)',:)','k')
plot(y_post_PC((250)',:)','k')


%%


MSE_1_post = mean(bsxfun(@minus,y_post,q1).^2,2);
MSE_1_post_C = mean(bsxfun(@minus,y_post_C,q1).^2,2);
MSE_1_post_PC = mean(bsxfun(@minus,y_post_PC,q1).^2,2);
MSE_1_post_C0 = mean(bsxfun(@minus,y_post_C0,q1).^2,2);
MSE_1_post_PC0 = mean(bsxfun(@minus,y_post_PC0,q1).^2,2);
% MSE_1_post_add = mean(bsxfun(@minus,y_post_add,q1).^2,2);

MSE_5_post = mean(bsxfun(@minus,y_post,q5).^2,2);
MSE_5_post_C = mean(bsxfun(@minus,y_post_C,q5).^2,2);
MSE_5_post_PC = mean(bsxfun(@minus,y_post_PC,q5).^2,2);
MSE_5_post_C0 = mean(bsxfun(@minus,y_post_C0,q5).^2,2);
MSE_5_post_PC0 = mean(bsxfun(@minus,y_post_PC0,q5).^2,2);
% MSE_5_post_add = mean(bsxfun(@minus,y_post_add,q5).^2,2);


figure(345)
subplot(1,2,1)
hold all
plot(MSE_1_post)
plot(MSE_1_post_C)
plot(MSE_1_post_PC)
plot(MSE_1_post_C0)
plot(MSE_1_post_PC0)
% plot(MSE_1_post_add)
hold off


subplot(1,2,2)
hold all
plot(MSE_5_post)
plot(MSE_5_post_C)
plot(MSE_5_post_PC)
plot(MSE_5_post_C0)
plot(MSE_5_post_PC0)
% plot(MSE_5_post_add)
hold off

legend('Post','C','PC','C0','PC0')
% legend('Post','C','PC','C0','PC0','P add')
suptitle([model,', sigma2= ',num2str(sigma2),', ',num2str(100*p_bar1),'%, ',num2str(100*p_bar),'%'])


%%
figure(5555)
hold on
plot(q1,'r')
plot(y_post(46,:),'b')
plot(y_post(100,:),'g')
plot(y_post(30,:),'g')


%%

 h_post2 = volatility_arch1(draw,y,y_S,H);   
 y_post2 = randn(M,H).*sqrt(h_post2);
 y_post2 = bsxfun(@plus,y_post,draw(:,1));
    y_post = sort(y_post);
    VaR_1_post = y_post(p_bar1*M,:); 
    VaR_5_post = y_post(p_bar*M,:); 
 