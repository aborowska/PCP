% S = 60 - 2; %1;
% 
% % RES(25,:) = [];
% met = 19;

ES_MSEs = zeros(S,met,3);
tmp = NaN(H,3);
BAD_es = zeros(S,met);

% methods = {'post','CP0','PCP0',...
%     'CP10','PCP10','CP20','PCP20','CP30','PCP30','CP40','PCP40',...
%     'CPv10','PCPv10','CPv20','PCPv20','CPv30','PCPv30','CPv40','PCPv40'}';

for ii = 1:S
    CDF_theor = RES{ii,1}.cdf_theor;

    ESs = zeros(met,H,3);
    ESs(1,:,:) = RES{ii,1}.post.ES_post';
    % all(ESs05(1,:,1) == RES{ii,1}.post.ES_post(1,:))
    try
        ESs(2,:,:) = RES{ii,1}.thr_0.ES_C';
        ESs(3,:,:) = RES{ii,1}.thr_0.ES_PC'; 
    catch
        ESs(2,:,:) = tmp;
        ESs(3,:,:) = tmp;
        BAD_es(ii,2:3) = 1;
    end
    try
        ESs(4,:,:) = RES{ii,1}.thr_10.ES_C';
        ESs(5,:,:) = RES{ii,1}.thr_10.ES_PC'; 
    catch
        ESs(4,:,:) = tmp;
        ESs(5,:,:) = tmp;         
        BAD_es(ii,4:5) = 1;
    end
    try
        ESs(6,:,:) = RES{ii,1}.thr_20.ES_C';
        ESs(7,:,:) = RES{ii,1}.thr_20.ES_PC'; 
    catch
        ESs(6,:,:) = tmp;
        ESs(7,:,:) = tmp; 
        BAD_es(ii,6:7) = 1;        
    end
    try
        ESs(8,:,:) = RES{ii,1}.thr_30.ES_C';
        ESs(9,:,:) = RES{ii,1}.thr_30.ES_PC'; 
    catch
        ESs(8,:,:) = tmp;
        ESs(9,:,:) = tmp;         
        BAD_es(ii,8:9) = 1;        
    end
    try
        ESs(10,:,:) = RES{ii,1}.thr_40.ES_C';
        ESs(11,:,:) = RES{ii,1}.thr_40.ES_PC'; 
    catch
        ESs(10,:,:) = tmp;
        ESs(11,:,:) = tmp;    
        BAD_es(ii,10:11) = 1;        
    end
    try
        ESs(12,:,:) = RES{ii,1}.thr_v10.ES_C';
        ESs(13,:,:) = RES{ii,1}.thr_v10.ES_PC'; 
    catch
        ESs(12,:,:) = tmp;
        ESs(13,:,:) = tmp;   
        BAD_es(ii,12:13) = 1;
    end
    try
        ESs(14,:,:) = RES{ii,1}.thr_v20.ES_C';
        ESs(15,:,:) = RES{ii,1}.thr_v20.ES_PC'; 
    catch
        ESs(14,:,:) = tmp;
        ESs(15,:,:) = tmp;  
        BAD_es(ii,14:15) = 1;        
    end
    try
        ESs(16,:,:) = RES{ii,1}.thr_v30.ES_C';
        ESs(17,:,:) = RES{ii,1}.thr_v30.ES_PC'; 
    catch
        ESs(16,:,:) = tmp;
        ESs(17,:,:) = tmp;
        BAD_es(ii,16:17) = 1;        
    end
    try
        ESs(18,:,:) = RES{ii,1}.thr_v40.ES_C';
        ESs(19,:,:) = RES{ii,1}.thr_v40.ES_PC'; 
    catch
        ESs(18,:,:) = tmp;
        ESs(19,:,:) = tmp;    
        BAD_es(ii,18:19) = 1;        
    end
    
    if (met > 19)

        try
            ESs(20,:,:) = RES{ii,1}.thr_vcp10.ES_C';
            ESs(21,:,:) = RES{ii,1}.thr_vcp10.ES_PC'; 
        catch
            ESs(20,:,:) = tmp;
            ESs(21,:,:) = tmp;   
            BAD_es(ii,20:21) = 1;
        end
        try
            ESs(22,:,:) = RES{ii,1}.thr_vcp20.ES_C';
            ESs(23,:,:) = RES{ii,1}.thr_vcp20.ES_PC'; 
        catch
            ESs(22,:,:) = tmp;
            ESs(23,:,:) = tmp;  
            BAD_es(ii,22:23) = 1;        
        end
        try
            ESs(24,:,:) = RES{ii,1}.thr_vcp30.ES_C';
            ESs(25,:,:) = RES{ii,1}.thr_vcp30.ES_PC'; 
        catch
            ESs(24,:,:) = tmp;
            ESs(25,:,:) = tmp;
            BAD_es(ii,24:25) = 1;        
        end
        try
            ESs(26,:,:) = RES{ii,1}.thr_vcp40.ES_C';
            ESs(27,:,:) = RES{ii,1}.thr_vcp40.ES_PC'; 
        catch
            ESs(26,:,:) = tmp;
            ESs(27,:,:) = tmp;    
            BAD_es(ii,26:27) = 1;        
        end       
    end
    
    MSEs = zeros(met,3);
    for mm = 1:met
        for pp = 1:3
            MSEs(mm,pp) = mean((ESs(mm,:,pp) - CDF_theor(pp,:)).^2,'omitnan');
        end
    end


    ES_MSEs(ii,:,:) = MSEs;
end



%% DM statistics    
DMes_05 = NaN(met,met);
DMes_1 =  NaN(met,met);
DMes_5 =  NaN(met,met);

fn_loss = @(xx) sqrt(xx);
% fn_loss = @(xx) xx;

for ii = 2:met
    for jj = 1:(ii-1)
        % loss differnetial vector (loss = MSE over H out-of-sample periods)
        % elements of d are iid hence no Newey-West type correction needed
%             d = MSE_1(:,ii) - MSE_1(:,jj);
        d1 = ES_MSEs(:,ii,1); %eval(char(bbb_05{ii}));
        d2 = ES_MSEs(:,jj,1); %eval(char(bbb_05{jj}));
        d = fn_loss(d1) - fn_loss(d2);
        d = d(~isnan(d));
        se = std(d)/sqrt(length(d));
        DMes_05(ii,jj) = mean(d)/se;

%             d = MSE_1(:,ii) - MSE_1(:,jj);
        d1 = ES_MSEs(:,ii,2); %eval(char(bbb_05{ii}));
        d2 = ES_MSEs(:,jj,2); %eval(char(bbb_05{jj}));
        d = fn_loss(d1) - fn_loss(d2);
        d = d(~isnan(d));
        se = std(d)/sqrt(length(d));
        DMes_1(ii,jj) = mean(d)/se;        
%             d = MSE_1(:,ii) - MSE_1(:,jj);

        d1 = ES_MSEs(:,ii,3); %eval(char(bbb_05{ii}));
        d2 = ES_MSEs(:,jj,3); %eval(char(bbb_05{jj}));
        d = fn_loss(d1) - fn_loss(d2);
        d = d(~isnan(d));
        se = std(d)/sqrt(length(d));
        DMes_5(ii,jj) = mean(d)/se;

    end
end 

 