%     addpath(genpath('include/'));
% clear all

% Lambda = [-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5];
Lambda = [-0.2,0.2];
R = length(Lambda);
  
met = 5;
    
DM_05 = NaN(met,met,R);
DM_1 =  NaN(met,met,R);
DM_5 =  NaN(met,met,R);
 
model = 'skt_ar1';

nu = 5; 
T = 10000; 
II = 10; 
H = 100;
 
for rr = 1:R
    lambda = Lambda(rr); 
    
    clear MSE* bbb* d se 

    try
        name = ['results/',model,'/',model,'_',num2str(lambda),'_',num2str(nu),'_T',num2str(T),'_H',num2str(H),'_II',num2str(II),'_PCP0_MC.mat'];    
        load(name, '-regexp','^MSE\w*')

        bbb = who('-regexp','^MSE_[^es]');

        bbb_05 = who('-regexp','MSE_05\w*$'); 
        bbb_1 = who('-regexp','MSE_1\w*$'); 
        bbb_5 = who('-regexp','MSE_5\w*$');

        ind = [1,3,5,2,4];

        bbb_05 = bbb_05(ind);
        bbb_1 = bbb_1(ind);
        bbb_5 = bbb_5(ind);    

        for ii = 2:met
            for jj = 1:(ii-1)
                % loss differnetial vector (loss = MSE over H out-of-sample periods)
                % elements of d are iid hence no Newey-West type correction needed
    %             d = MSE_1(:,ii) - MSE_1(:,jj);
                d1 = eval(char(bbb_05{ii}));
                d2 = eval(char(bbb_05{jj}));
                d = sqrt(d1) - sqrt(d2);
                se = std(d)/sqrt(length(d));
                DM_05(ii,jj,rr) = mean(d)/se;

    %             d = MSE_5(:,jj) - MSE_5(:,ii);
    %             DM(jj,ii) = mean(d)/(std(d)/sqrt(S));    
                d1 = eval(char(bbb_1{ii}));
                d2 = eval(char(bbb_1{jj}));
                d = sqrt(d1) - sqrt(d2);
                se = std(d)/sqrt(length(d));
                DM_1(ii,jj,rr) = mean(d)/se;


                d1 = eval(char(bbb_5{ii}));
                d2 = eval(char(bbb_5{jj}));
                d = sqrt(d1) - sqrt(d2);
                se = std(d)/sqrt(length(d));
                DM_5(ii,jj,rr) = mean(d)/se;
            end
        end 
        
    catch
        
    end
end    
       