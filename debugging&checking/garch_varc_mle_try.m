% P(y_t>=threshold|y_1,...,y_(t-1)) = 1 - P(y_t<threshold|y_1,...,y_(t-1))

% threshold_t 
% P(y_t>=threshold_t|y_1,...,y_(t-1)) = 1 - P(y_t<threshold_t|y_1,...,y_(t-1))
  
thr_const = sort(y(1:T));
thr_const = thr_const(0.5*T);
threshold = 0.5; %<---------- HiGhER?
quantile = norminv(threshold); % quantile is the fixed quantile of the standard normal distribution 

cond = zeros(T,1);
h_MLE = zeros(T,1);
THR = zeros(T,1);
for jj=1:T
    if (jj==1)
        h_MLE(jj) = y_S;            
    else
        mu = y(jj-1) - mu_MLE(1);
        h_MLE(jj) = mu_MLE(2)*(1.0-mu_MLE(3)-mu_MLE(4)) + mu_MLE(3)*(mu*mu) + mu_MLE(4)*h_MLE(jj-1);                  
    end       
    THR(jj) = (y(jj) - mu_MLE(1))/sqrt(h_MLE(jj));
    if (THR(jj) >= quantile) % OTHER WAY ROUND AS IN MEX:COND IS WHEN TO CENSOR
        cond(jj) = 1;
    end
    THR(jj) = mu_MLE(1) + sqrt(h_MLE(jj))*quantile; % time varying threshold for a nonstandardised variable
                                                    % to put to the cdf function
end

THR_cond = THR; 
THR_cond(~cond) = NaN;
y_cond = y;
y_cond(~cond) = NaN;

figure(1)
set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
hold on
plot(y(1:T),'k')
plot(y_cond(1:T),'m')
plot(cond,'.r')
plot(THR*0 + thr_const,'r')
plot(THR,'g')
plot(THR_cond,'b')
hold off
title('garch(1,1), sigma2=2, MLE based time varying vs constant threshold (0.5 quantile)')
legend('y','y censored','censoring cond','time constant thres.',...
    'time var. thres. theor.','time var. thres. applied')

figure(2)
set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
hold on
plot(200:400,y(200:400),'k')
plot(200:400,y_cond(200:400),'m')
plot(200:400,cond(200:400),'.r')
plot(200:400,THR(200:400)*0 + thr_const,'r')
plot(200:400,THR(200:400),'g')
plot(200:400,THR_cond(200:400),'b')
hold off
title('zoom in: garch(1,1), sigma2=2, MLE based time varying vs constant threshold (0.5 quantile)')
legend('y','y censored','censoring condition','time constant thres.',...
    'time var. thres. theor.','time var. thres. applied')

MUS_CM = zeros(100,4);
for ii = 1:100
    kernel_init = @(xx) - C_posterior_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, norminv(ii/100), y_S)/T;    
    MUS_CM(ii,:) = fminunc(kernel_init,mu_init,options);
end

figure(444)
set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
param = {'\mu','\omega','\alpha','\beta'};
for ii=1:4
    subplot(2,2,ii)
    plot(0.01:0.01:1,MUS_CM(:,ii))
    xlabel('censoring quantile')
    title(param{ii})
end
suptitle('MLEs for the (regular) MLE based time varying threshold censored likelihoods as a function of quantile for censoring')


%%

threshold = 0.1; %<---------- HiGhER?
quantile = norminv(threshold);
fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold);
kernel_init = @(xx) - C_posterior_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S)/T;    
kernel = @(xx) C_posterior_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S);
[mu_Cm01,~,~,~,~,Sigma_Cm01] = fminunc(kernel_init,mu_init,options);
Sigma_Cm01 = inv(T*Sigma_Cm01);

%     cont.mit.CV_tol = 0.3; 
cont.mit.CV_max = 1.5; %1.9;
CV_Cm01 = cont.mit.CV_old;
while (CV_Cm01(end) >= 2)
    try
        [mu_Cm01,~,~,~,~,Sigma_Cm01] = fminunc(kernel_init,mu_init,options);
        Sigma_Cm01 = inv(T*Sigma_Cm01);
        draw_Cm01 = rmvt(mu_Cm01,Sigma_Cm01,df,M+BurnIn);
        mit_Cm01 = struct('mu',mu_Cm01,'Sigma',reshape(Sigma_Cm01,1,length(mu_Cm01)^2),'df', df, 'p', 1);
        [mit_Cm01, CV_Cm01] = MitISEM_new2(mit_Cm01, kernel, mu_init, cont, GamMat);   
        if CV_C(end)>2
            [mit_Cm01, CV_Cm01] = MitISEM_new2(mit_Cm01, kernel, mu_init, cont, GamMat);   
        end
        [draw_Cm01, lnk_Cm01] = fn_rmvgt_robust(M+BurnIn, mit_Cm01, kernel, false);
        lnd_Cm01 = dmvgt(draw_Cm01, mit_Cm01, true, GamMat);    
    catch
        mu_Cm01 = fminunc(kernel_init,mu_init,options);
        mit_Cm01 = struct('mu',mu_Cm01,'Sigma',reshape(Sigma,1,length(mu_Cm01)^2),'df', df, 'p', 1);
        [mit_Cm01, CV_Cm01] = MitISEM_new2(mit_Cm01, kernel, mu_init, cont, GamMat);   
        if CV_Cm01(end)>2
            [mit_Cm01, CV_Cm01] = MitISEM_new2(mit_Cm01, kernel, mu_init, cont, GamMat);   
        end
        [draw_Cm01, lnk_Cm01] = fn_rmvgt_robust(M+BurnIn, mit_Cm01, kernel, false);
        lnd_Cm01 = dmvgt(draw_Cm01, mit_Cm01, true, GamMat);    
    end 
end
lnw_Cm01 = lnk_Cm01 - lnd_Cm01;
lnw_Cm01 = lnw_Cm01 - max(lnw_Cm01);
[ind, a] = fn_MH(lnw_Cm01);
draw_Cm01 = draw_Cm01(ind,:);
accept_Cm01 = a/(M+BurnIn);
draw_Cm01 = draw_Cm01(BurnIn+1:BurnIn+M,:);

h_post_Cm01 = volatility_garch11(draw_Cm01,y,y_S,H);
y_post_Cm01 = randn(M,H).*sqrt(h_post_Cm01);
y_post_Cm01 = bsxfun(@plus,y_post_Cm01,draw_Cm01(:,1));
y_post_Cm01 = sort(y_post_Cm01);
VaR_1_post_Cm01 = y_post_Cm01(p_bar1*M,:); 
VaR_5_post_Cm01 = y_post_Cm01(p_bar*M,:); 
mean_draw_Cm01 = mean(draw_Cm01);
median_draw_Cm01 = median(draw_Cm01);
std_draw_Cm01 = std(draw_Cm01);

%%%%%%%%%%%%%%%%%%%
threshold = 0.5; %<---------- HiGhER?
quantile = norminv(threshold);
fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold);
kernel_init = @(xx) - C_posterior_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S)/T;    
kernel = @(xx) C_posterior_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S);
[mu_Cm05,~,~,~,~,Sigma_Cm05] = fminunc(kernel_init,mu_init,options);
Sigma_Cm05 = inv(T*Sigma_Cm05);

%     cont.mit.CV_tol = 0.3; 
cont.mit.CV_max = 1.5; %1.9;
CV_Cm05 = cont.mit.CV_old;
while (CV_Cm05(end) >= 2)
    try
        [mu_Cm05,~,~,~,~,Sigma_Cm05] = fminunc(kernel_init,mu_init,options);
        Sigma_Cm05 = inv(T*Sigma_Cm05);
        draw_Cm05 = rmvt(mu_Cm05,Sigma_Cm05,df,M+BurnIn);
        mit_Cm05 = struct('mu',mu_Cm05,'Sigma',reshape(Sigma_Cm05,1,length(mu_Cm05)^2),'df', df, 'p', 1);
        [mit_Cm05, CV_Cm05] = MitISEM_new2(mit_Cm05, kernel, mu_init, cont, GamMat);   
        if CV_C(end)>2
            [mit_Cm05, CV_Cm05] = MitISEM_new2(mit_Cm05, kernel, mu_init, cont, GamMat);   
        end
        [draw_Cm05, lnk_Cm05] = fn_rmvgt_robust(M+BurnIn, mit_Cm05, kernel, false);
        lnd_Cm05 = dmvgt(draw_Cm05, mit_Cm05, true, GamMat);    
    catch
        mu_Cm05 = fminunc(kernel_init,mu_init,options);
        mit_Cm05 = struct('mu',mu_Cm05,'Sigma',reshape(Sigma,1,length(mu_Cm05)^2),'df', df, 'p', 1);
        [mit_Cm05, CV_Cm05] = MitISEM_new2(mit_Cm05, kernel, mu_init, cont, GamMat);   
        if CV_Cm05(end)>2
            [mit_Cm05, CV_Cm05] = MitISEM_new2(mit_Cm05, kernel, mu_init, cont, GamMat);   
        end
        [draw_Cm05, lnk_Cm05] = fn_rmvgt_robust(M+BurnIn, mit_Cm05, kernel, false);
        lnd_Cm05 = dmvgt(draw_Cm05, mit_Cm05, true, GamMat);    
    end 
end
lnw_Cm05 = lnk_Cm05 - lnd_Cm05;
lnw_Cm05 = lnw_Cm05 - max(lnw_Cm05);
[ind, a] = fn_MH(lnw_Cm05);
draw_Cm05 = draw_Cm05(ind,:);
accept_Cm05 = a/(M+BurnIn);
draw_Cm05 = draw_Cm05(BurnIn+1:BurnIn+M,:);

h_post_Cm05 = volatility_garch11(draw_Cm05,y,y_S,H);
y_post_Cm05 = randn(M,H).*sqrt(h_post_Cm05);
y_post_Cm05 = bsxfun(@plus,y_post_Cm05,draw_Cm05(:,1));
y_post_Cm05 = sort(y_post_Cm05);
VaR_1_post_Cm05 = y_post_Cm05(p_bar1*M,:); 
VaR_5_post_Cm05 = y_post_Cm05(p_bar*M,:); 
mean_draw_Cm05 = mean(draw_Cm05);
median_draw_Cm05 = median(draw_Cm05);
std_draw_Cm05 = std(draw_Cm05);
