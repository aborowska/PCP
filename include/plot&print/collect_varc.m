clear all
models = {'ar1','arch1','garch11','agarch11'};
Ts = [1000,2500];
H = 100;
II = 10;
v_new = '(R2015a)';
sigma1 = 1;
sigma2 = 2;

MSEs = zeros(16,4);
ii = 0;
for model = models
    for T = Ts
        ii = ii+1;
        name = ['results/',char(model),'/',char(model),'_',num2str(sigma1),'_',num2str(sigma2),...
            '_T',num2str(T),'_H',num2str(H),'_II',num2str(II)...
            '_PCP0_MC_',v_new,'_varc.mat'];
        load(name,'MSE_1','MSE_1_post','MSE_1_post_C','MSE_1_post_PC',...
            'MSE_5','MSE_5_post','MSE_5_post_C','MSE_5_post_PC')

        MSEs(2*ii-1,1) = mean(MSE_1);
        MSEs(2*ii-1,2) = mean(MSE_1_post);
        MSEs(2*ii-1,3) = mean(MSE_1_post_C);
        MSEs(2*ii-1,4) = mean(MSE_1_post_PC);

        MSEs(2*ii,1) = mean(MSE_5);
        MSEs(2*ii,2) = mean(MSE_5_post);
        MSEs(2*ii,3) = mean(MSE_5_post_C);
        MSEs(2*ii,4) = mean(MSE_5_post_PC);

        clear MSE_1 MSE_1_post MSE_1_post_C MSE_1_post_PC MSE_5 MSE_5_post MSE_5_post_C MSE_5_post_PC;
    end
end