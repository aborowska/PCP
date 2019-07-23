% clear all;
model = 'skt_iid';%'mix_iid'; %'mix_iid'; %'skt_iid';

switch model
    case 'skt_iid'
        nu = 5;
        param = nu;
%         Lambda = -0.6 : 0.1: 0.2;
        Lambda = -0.6 : 0.1: 0.6;
    case 'mix_iid'
        rho = 0.75;
        param = rho;
%         Lambda = 0.1 : 0.1: 0.9;
        Lambda = 0.05:0.05:0.45;
end
R = length(Lambda);


meanMSEs_05 = zeros(3,R,3);
meanMSEs_1 = zeros(3,R,3);
meanMSEs_5 = zeros(3,R,3);

stdMSEs_05 = zeros(3,R,3);
stdMSEs_1 = zeros(3,R,3);
stdMSEs_5 = zeros(3,R,3);

TT = [100,1000,10000];

DM_05 = NaN(3,R,3);
DM_1 = NaN(3,R,3);
DM_5 = NaN(3,R,3);

for t = 1:3
    T = TT(t);
    
    for ii = 1:R
        lambda = Lambda(ii);

        load(['results/',model,'/',model,'_',num2str(lambda),'_',num2str(param),'_T',num2str(T),'_MC.mat'],...
            '-regexp','^MSE')
         
        %% 99.5%
        d = sqrt(MSE_05_post) - sqrt(MSE_05_post_C0);
        se = std(d)/sqrt(length(d));
        DM_05(1,ii,t) = mean(d)/se; % Post - CP0 > 0 GOOD
         
        d = sqrt(MSE_05_post) - sqrt(MSE_05_post_C);
        se = std(d)/sqrt(length(d));
        DM_05(2,ii,t) = mean(d)/se; % Post - CP > 0 GOOD

        d = sqrt(MSE_05_post_C0) - sqrt(MSE_05_post_C);
        se = std(d)/sqrt(length(d));
        DM_05(3,ii,t) = mean(d)/se; % CP0 - CP > 0 more censoring better
        
        %% 99%
        d = sqrt(MSE_1_post) - sqrt(MSE_1_post_C0);
        se = std(d)/sqrt(length(d));
        DM_1(1,ii,t) = mean(d)/se; % Post - CP0 > 0 GOOD
         
        d = sqrt(MSE_1_post) - sqrt(MSE_1_post_C);
        se = std(d)/sqrt(length(d));
        DM_1(2,ii,t) = mean(d)/se; % Post - CP > 0 GOOD        
         
        d = sqrt(MSE_1_post_C0) - sqrt(MSE_1_post_C);
        se = std(d)/sqrt(length(d));
        DM_1(3,ii,t) = mean(d)/se; % CP0 - CP > 0 more censoring better
        
        %% 95%
        d = sqrt(MSE_5_post) - sqrt(MSE_5_post_C0);
        se = std(d)/sqrt(length(d));
        DM_5(1,ii,t) = mean(d)/se; % Post - CP0 > 0 GOOD
         
        d = sqrt(MSE_5_post) - sqrt(MSE_5_post_C);
        se = std(d)/sqrt(length(d));
        DM_5(2,ii,t) = mean(d)/se; % Post - CP > 0 GOOD       
         
        d = sqrt(MSE_5_post_C0) - sqrt(MSE_5_post_C);
        se = std(d)/sqrt(length(d));
        DM_5(3,ii,t) = mean(d)/se; % CP0 - CP > 0 more censoring better        
    end
end



 for t = 1:3
    T = TT(t);
    
    for ii = 1:R
        lambda = round(-0.6 + (ii-1)*0.1,1);

        load(['results/',model,'/',model,'_',num2str(lambda),'_',num2str(nu),'_T',num2str(T),'_MC.mat'],...
            '-regexp','^MSE')               
        
        temp = [mean(MSE_05_post);
                mean(MSE_05_post_C0);
                mean(MSE_05_post_C)];
        meanMSEs_05(:,ii,t) = temp;

        temp = [mean(MSE_1_post);
            mean(MSE_1_post_C0);
            mean(MSE_1_post_C)];
        meanMSEs_1(:,ii,t) = temp;

        temp = [mean(MSE_5_post);
            mean(MSE_5_post_C0);
            mean(MSE_5_post_C)];
        meanMSEs_5(:,ii,t) = temp;


        temp = [std(MSE_05_post);
            std(MSE_05_post_C0);
            std(MSE_05_post_C)];
        stdMSEs_05(:,ii,t) = temp;

        temp = [std(MSE_1_post);
            std(MSE_1_post_C0);
            std(MSE_1_post_C)];
        stdMSEs_1(:,ii,t) = temp;

        temp = [std(MSE_5_post);
            std(MSE_5_post_C0);
            std(MSE_5_post_C)];
        stdMSEs_5(:,ii,t) = temp;
    end
 end


 