% clear

%% 
% T = 1000; 
% S = 60;
met = 19;
methods = {'post','CP0','PCP0',...
    'CP10','PCP10','CP20','PCP20','CP30','PCP30','CP40','PCP40',...
    'CPv10','PCPv10','CPv20','PCPv20','CPv30','PCPv30','CPv40','PCPv40'}';
% 
% 
% load(['..\..\..\..\ownCloud2\ForDropbox\garch11_1_2_T',num2str(T),'_H100_II2_PCP_NEW_MULTI_THRES.mat'])
% 
% RES1 = RES;
% SDD1 = SDD;
% time_total1 = time_total;
% 
% load(['..\..\..\..\ownCloud2\ForDropbox\garch11_1_2_T',num2str(T),'_H100_II2_PCP_NEW_MULTI_THRES_2.mat'])
% RES2 = RES;
% SDD2 = SDD;
% time_total2 = time_total;
% 
% load(['..\..\..\..\ownCloud2\ForDropbox\garch11_1_2_T',num2str(T),'_H100_II2_PCP_NEW_MULTI_THRES_3.mat'])
% RES3 = RES;
% SDD3 = SDD;
% time_total3 = time_total;
% 
% RES = [RES1;RES2;RES3];


%% Check bad draws (not satisfying prior)
S = 20;
D = 5;

BAD_dr = zeros(S,met);
mean_draws = -inf*ones(S,met,D);

for ii = 1:S

    mean_draws(ii,1,:) = mean(RES{ii,1}.post.draw_post);
    % all(draws05(1,:,1) == RES{ii,1}.post.draw_post(1,:))
    try
        mean_draws(ii,2,:) = mean(RES{ii,1}.thr_0.draw_C);
        mean_draws(ii,3,:) = mean(RES{ii,1}.thr_0.draw_PC); 
    catch
        mean_draws(ii,2,:) = NaN;
        mean_draws(ii,3,:) = NaN;
        BAD_dr(ii,2:3) = 1;
    end
    try
        mean_draws(ii,4,:) = mean(RES{ii,1}.thr_10.draw_C);
        mean_draws(ii,5,:) = mean(RES{ii,1}.thr_10.draw_PC); 
    catch
        mean_draws(ii,4,:) = NaN;
        mean_draws(ii,5,:) = NaN;         
        BAD_dr(ii,4:5) = 1;
    end
    try
        mean_draws(ii,6,:) = mean(RES{ii,1}.thr_20.draw_C);
        mean_draws(ii,7,:) = mean(RES{ii,1}.thr_20.draw_PC); 
    catch
        mean_draws(ii,6,:) = NaN;
        mean_draws(ii,7,:) = NaN; 
        BAD_dr(ii,6:7) = 1;        
    end
    try
        mean_draws(ii,8,:) = mean(RES{ii,1}.thr_30.draw_C);
        mean_draws(ii,9,:) = mean(RES{ii,1}.thr_30.draw_PC); 
    catch
        mean_draws(ii,8,:) = NaN;
        mean_draws(ii,9,:) = NaN;         
        BAD_dr(ii,8:9) = 1;        
    end
    try
        mean_draws(ii,10,:) = mean(RES{ii,1}.thr_40.draw_C);
        mean_draws(ii,11,:) = mean(RES{ii,1}.thr_40.draw_PC); 
    catch
        mean_draws(ii,10,:) = NaN;
        mean_draws(ii,11,:) = NaN;    
        BAD_dr(ii,10:11) = 1;        
    end
    try
        mean_draws(ii,12,:) = mean(RES{ii,1}.thr_v10.draw_C);
        mean_draws(ii,13,:) = mean(RES{ii,1}.thr_v10.draw_PC); 
    catch
        mean_draws(ii,12,:) = NaN;
        mean_draws(ii,13,:) = NaN;   
        BAD_dr(ii,12:13) = 1;
    end
    try
        mean_draws(ii,14,:) = mean(RES{ii,1}.thr_v20.draw_C);
        mean_draws(ii,15,:) = mean(RES{ii,1}.thr_v20.draw_PC); 
    catch
        mean_draws(ii,14,:) = NaN;
        mean_draws(ii,15,:) = NaN;  
        BAD_dr(ii,14:15) = 1;        
    end
    try
        mean_draws(ii,16,:) = mean(RES{ii,1}.thr_v30.draw_C);
        mean_draws(ii,17,:) = mean(RES{ii,1}.thr_v30.draw_PC); 
    catch
        mean_draws(ii,16,:) = NaN;
        mean_draws(ii,17,:) = NaN;
        BAD_dr(ii,16:17) = 1;        
    end
    try
        mean_draws(ii,18,:) = mean(RES{ii,1}.thr_v40.draw_C);
        mean_draws(ii,19,:) = mean(RES{ii,1}.thr_v40.draw_PC); 
    catch
        mean_draws(ii,18,:) = NaN;
        mean_draws(ii,19,:) = NaN;    
        BAD_dr(ii,18:19) = 1;        
    end

end
 




%% Collect VaR values and compute the MSE matrix 
VaR_MSEs = zeros(S,met,3);
VaR_MTEs = zeros(S,met,3); % Total error to see whether undr/overestimation
tmp = NaN(100,3);
BAD_var = zeros(S,met);

for ii = 1:S
    Q_theor = RES{ii,1}.q_theor;

    VaRs = zeros(met,100,3);
    VaRs(1,:,:) = RES{ii,1}.post.VaR_post';
    % all(VaRs05(1,:,1) == RES{ii,1}.post.VaR_post(1,:))
    try
        VaRs(2,:,:) = RES{ii,1}.thr_0.VaR_C';
        VaRs(3,:,:) = RES{ii,1}.thr_0.VaR_PC'; 
    catch
        VaRs(2,:,:) = tmp;
        VaRs(3,:,:) = tmp;
        BAD_var(ii,2:3) = 1;
    end
    try
        VaRs(4,:,:) = RES{ii,1}.thr_10.VaR_C';
        VaRs(5,:,:) = RES{ii,1}.thr_10.VaR_PC'; 
    catch
        VaRs(4,:,:) = tmp;
        VaRs(5,:,:) = tmp;         
        BAD_var(ii,4:5) = 1;
    end
    try
        VaRs(6,:,:) = RES{ii,1}.thr_20.VaR_C';
        VaRs(7,:,:) = RES{ii,1}.thr_20.VaR_PC'; 
    catch
        VaRs(6,:,:) = tmp;
        VaRs(7,:,:) = tmp; 
        BAD_var(ii,6:7) = 1;        
    end
    try
        VaRs(8,:,:) = RES{ii,1}.thr_30.VaR_C';
        VaRs(9,:,:) = RES{ii,1}.thr_30.VaR_PC'; 
    catch
        VaRs(8,:,:) = tmp;
        VaRs(9,:,:) = tmp;         
        BAD_var(ii,8:9) = 1;        
    end
    try
        VaRs(10,:,:) = RES{ii,1}.thr_40.VaR_C';
        VaRs(11,:,:) = RES{ii,1}.thr_40.VaR_PC'; 
    catch
        VaRs(10,:,:) = tmp;
        VaRs(11,:,:) = tmp;    
        BAD_var(ii,10:11) = 1;        
    end
    try
        VaRs(12,:,:) = RES{ii,1}.thr_v10.VaR_C';
        VaRs(13,:,:) = RES{ii,1}.thr_v10.VaR_PC'; 
    catch
        VaRs(12,:,:) = tmp;
        VaRs(13,:,:) = tmp;   
        BAD_var(ii,12:13) = 1;
    end
    try
        VaRs(14,:,:) = RES{ii,1}.thr_v20.VaR_C';
        VaRs(15,:,:) = RES{ii,1}.thr_v20.VaR_PC'; 
    catch
        VaRs(14,:,:) = tmp;
        VaRs(15,:,:) = tmp;  
        BAD_var(ii,14:15) = 1;        
    end
    try
        VaRs(16,:,:) = RES{ii,1}.thr_v30.VaR_C';
        VaRs(17,:,:) = RES{ii,1}.thr_v30.VaR_PC'; 
    catch
        VaRs(16,:,:) = tmp;
        VaRs(17,:,:) = tmp;
        BAD_var(ii,16:17) = 1;        
    end
    try
        VaRs(18,:,:) = RES{ii,1}.thr_v40.VaR_C';
        VaRs(19,:,:) = RES{ii,1}.thr_v40.VaR_PC'; 
    catch
        VaRs(18,:,:) = tmp;
        VaRs(19,:,:) = tmp;    
        BAD_var(ii,18:19) = 1;        
    end
    
    
    
    MSEs = zeros(met,3);
    MTEs = zeros(met,3);
    for mm = 1:met
        for pp = 1:3
            MSEs(mm,pp) = mean((VaRs(mm,:,pp) - Q_theor(pp,:)).^2,'omitnan');
            MTEs(mm,pp) = mean((VaRs(mm,:,pp) - Q_theor(pp,:)),'omitnan');
        end
    end


    VaR_MSEs(ii,:,:) = MSEs;
    VaR_MTEs(ii,:,:) = MTEs;
end



%% Compute the DM stats

    
DM_05 = NaN(met,met);
DM_1 =  NaN(met,met);
DM_5 =  NaN(met,met);

fn_loss = @(xx) sqrt(xx);
% fn_loss = @(xx) xx;

for ii = 2:met
    for jj = 1:(ii-1)
        % loss differnetial vector (loss = MSE over H out-of-sample periods)
        % elements of d are iid hence no Newey-West type correction needed
%             d = MSE_1(:,ii) - MSE_1(:,jj);
        d1 = VaR_MSEs(:,ii,1); %eval(char(bbb_05{ii}));
        d2 = VaR_MSEs(:,jj,1); %eval(char(bbb_05{jj}));
        d = fn_loss(d1) - fn_loss(d2);
        d = d(~isnan(d));
        se = std(d)/sqrt(length(d));
        DM_05(ii,jj) = mean(d)/se;

%             d = MSE_1(:,ii) - MSE_1(:,jj);
        d1 = VaR_MSEs(:,ii,2); %eval(char(bbb_05{ii}));
        d2 = VaR_MSEs(:,jj,2); %eval(char(bbb_05{jj}));
        d = fn_loss(d1) - fn_loss(d2);
        d = d(~isnan(d));
        se = std(d)/sqrt(length(d));
        DM_1(ii,jj) = mean(d)/se;        
%             d = MSE_1(:,ii) - MSE_1(:,jj);

        d1 = VaR_MSEs(:,ii,3); %eval(char(bbb_05{ii}));
        d2 = VaR_MSEs(:,jj,3); %eval(char(bbb_05{jj}));
        d = fn_loss(d1) - fn_loss(d2);
        d = d(~isnan(d));
        se = std(d)/sqrt(length(d));
        DM_5(ii,jj) = mean(d)/se;

    end
end 



%% Do the same for ES
collect_DM_ES_MC_garch


%% Save
name_s = ['results/agarch11/draws_agarch11_1_2_T',num2str(T),'_H100_II2_S',num2str(S),'.mat'];
save(name_s, '-regexp', '^DM','VaR_MSEs','ES_MSEs');
