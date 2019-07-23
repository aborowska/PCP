% mu_C = [9.5292    0.0673    1.9216    0.0583    0.9400];
% kernel = @(xx) C_posterior_t_garch11(xx, y(1:T,1), threshold, y_S, hyper);
% kernel2 = @(xx) C_posterior_t_garch11_2(xx, y(1:T,1), threshold, y_S, hyper);
% kernel_mex = @(xx) C_posterior_t_garch11_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);
% kernel_mex2 = @(xx) C_posterior_t_garch11_2_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);
% 
% tic; kernel(mu_C); toc;
% tic; kernel2(mu_C); toc;
% tic; kernel_mex(mu_C); toc;
% tic; kernel_mex2(mu_C); toc;
% % Elapsed time is 0.230514 seconds. MATLAB with cdf every time
% % Elapsed time is 0.041971 seconds. MATLAB with cdf vectorised
% % Elapsed time is 0.240488 seconds. MEX with cdf every time
% % Elapsed time is 0.002031 seconds. mex with cdf vectorised
clear all

save_on = false;
figures_path='';
data = 44; % MSFT 44; % GSPC 22;
H = 2275 +248;
[y, T, y_plot, data_name, time] = Load_data_empirical(data, H);   
time = [2002,2016];       
% data_stats = Plot_data(y_plot,H,time,save_on,figures_path);

p_bar0 = 0.005;
p_bar1 = 0.01;
p_bar = 0.05; 
P_bars = [p_bar0, p_bar1, p_bar];

mu_init = [8,0,1,0.1,0.8];
y_S = var(y(1:T));
hyper = 0.01; 
options = optimset('Display','off');

kernel_init = @(xx) - posterior_t_garch11(xx, y(1:T), y_S, hyper)/T;
[mu_MLE1,~,exitflag,output,~,Sigma1] = fminunc(kernel_init,mu_init,options);
[mu_MLE,~,exitflag,output,~,Sigma] = fminunc(kernel_init,mu_MLE1,options);


threshold = sort(y(1:T));
threshold = threshold(round(2*p_bar*T));
fprintf('*** Censored Posterior, threshold 10%% ***\n');
%     kernel_init = @(xx) - C_posterior_t_garch11_mex(xx, y(1:T,1), threshold, y_S, GamMat, hyper)/T;    
kernel_init = @(xx) - C_posterior_t_garch11_2(xx, y(1:T,1), threshold, y_S, hyper)/T;    
[mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
[mu_C2,~,~,~,~,Sigma_C2] = fminunc(kernel_init,mu_C,options);


threshold_m = 0.1; %<---------- HiGhER?
quantile = tinv(threshold_m, mu_MLE(1));
fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold_m);
kernel_init = @(xx) - C_posterior_t_garch11_varc_mle(xx, y(1:T,1), mu_MLE, quantile, y_S, hyper)/T;    
kernel = @(xx) C_posterior_t_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S, GamMat, hyper);

[mu_Cm2,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_init,options);
[mu_Cm,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_Cm2,options);
  
save([data_name,'_short_try.mat'],'-regexp','^mu','^Sigma');
if false
%% CENSORED MLE PARAMETERS    
    threshold_m = 0.1; %<---------- HiGhER?
    quantile = tinv(threshold_m, mu_MLE(1));
    fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',threshold_m);
    kernel_init = @(xx) - C_posterior_t_garch11_varc_mle(xx, y(1:T,1), mu_MLE, quantile, y_S, hyper)/T;    
    kernel = @(xx) C_posterior_t_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S, GamMat, hyper);
  
    kernel_init(mu_MLE)
    draw_PCm = sim_cond_inv_trans(draw, grid, kernel);
 

    mean_draw_PCm = mean(draw_PCm);
    median_draw_PCm = median(draw_PCm);
    std_draw_PCm = std(draw_PCm);
    
    
    % compute the implied volatility for the last in-sample period
    hT_PCm = volatility_t_garch11(draw_PCm,y(1:T),y_S,0); 
    


   kernel1 = @(xx) C_posterior_t_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S, GamMat, hyper);
   kernel2 = @(xx) C_posterior_t_garch11_varc_mle(xx, y(1:T,1), mu_MLE, quantile, y_S, hyper);

   kernel1(mu_MLE)
   kernel2(mu_MLE)
   

threshold = sort(y(1:T));   
threshold = threshold(round(2*p_bar*T));
mu_MLE = [ 5.3233    0.0164   66.2934    0.0559    0.9440];

kernel = @(xx) C_posterior_t_garch11(xx, y(1:T,1), threshold, y_S, hyper);
kernel2 = @(xx) C_posterior_t_garch11_2(xx, y(1:T,1), threshold, y_S, hyper);
kernel_mex = @(xx) C_posterior_t_garch11_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);
kernel_mex2 = @(xx) C_posterior_t_garch11_2_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);
 
   kernel(mu_MLE)
   kernel2(mu_MLE)
   kernel_mex(mu_MLE)
   kernel_mex2(mu_MLE)
end