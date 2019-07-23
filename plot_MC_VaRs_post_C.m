lambda = 0.4;
T = 1000;
load(['results/skt_iid/skt_iid_',num2str(lambda),'_5_T',num2str(T),'_MC.mat'])


subplot(3,1,1)
plot(VaR_05_post)
hold on
plot(VaR_05_post_C)
plot(VaR_05_post_C0)
plot(VaR_05_post*0+q05,'k')
legend('post','C10%','C0','true')
title('99.5% VaR')


subplot(3,1,2)
plot(VaR_1_post)
hold on
plot(VaR_1_post_C)
plot(VaR_1_post_C0)
plot(VaR_1_post*0+q1,'k')
legend('post','C10%','C0','true')
title('99% VaR')


subplot(3,1,3)
plot(VaR_5_post)
hold on
plot(VaR_5_post_C)
plot(VaR_5_post_C0)
plot(VaR_5_post*0+q5,'k')
legend('post','C10%','C0','true')
title('95% VaR')

suptitle(['Sk-t iid T=',num2str(T),', lambda=',...
    num2str(lambda),' (right skewness)'])