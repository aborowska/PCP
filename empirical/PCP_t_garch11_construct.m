%% Uncensored Posterior
kernel_init = @(xx) - posterior_t_garch11(xx, y(1:T), y_S, hyper)/T;
kernel = @(xx) posterior_t_garch11_mex(xx, y(1:T), y_S, GamMat, hyper);
[mu_MLE1,~,exitflag,output,~,Sigma] = fminunc(kernel_init,mu_init,options);
[mu_MLE,~,exitflag,output,~,Sigma] = fminunc(kernel_init,mu_MLE1,options);

Sigma = inv(T*Sigma);  
r = fn_testSigma(reshape(Sigma,1,5^2)); % if r==1 then there is a problem with Sigma_C
if r
    Sigma_start = Sigma;        
else
    Sigma_start = csvread('MSFT_Sigma_mle.csv',1,1);
end
cont.mit.CV_max = 1.9;

mit = struct('mu',mu_MLE,'Sigma',reshape(Sigma_start,1,length(mu_MLE)^2),'df', df, 'p', 1);
[mit, CV] = MitISEM_new2(mit, kernel, mu_init, cont, GamMat);            

%% Censored Posterio constant threshold
kernel_init = @(xx) - C_posterior_t_garch11_2(xx, y(1:T,1), threshold, y_S, hyper)/T;    
[mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
[mu_C2,~,~,~,~,Sigma_C2] = fminunc(kernel_init,mu_C,options);
Sigma_C = inv(T*Sigma_C);
    
kernel = @(xx) C_posterior_t_garch11_2_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);

r = fn_testSigma(reshape(Sigma_C,1,5^2)); % if r==1 then there is a problem with Sigma_C
if ~r
    Sigma_start = Sigma_C;   
else
    Sigma_start = nearestSPD(Sigma_C);      
end

cont.mit.CV_max = 1.9;

df_C = 5;
mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_start,1,length(mu_C)^2),'df', df_C, 'p', 1);
[mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
 
  
%% Censored Posterio time-varying threshold
threshold_m = 0.1; %<---------- HiGhER?
quantile = tinv(threshold_m, mu_MLE(1));
kernel_init = @(xx) - C_posterior_t_garch11_varc_mle(xx, y(1:T,1), mu_MLE, quantile, y_S, hyper)/T;    
[mu_Cm2,~,~,~,~,Sigma_Cm2] = fminunc(kernel_init,mu_init,options);
[mu_Cm,~,~,~,~,Sigma_Cm] = fminunc(kernel_init,mu_Cm2,options);
Sigma_Cm = inv(T*Sigma_Cm);
    
kernel = @(xx) C_posterior_t_garch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, quantile, y_S, GamMat, hyper);
r = fn_testSigma(reshape(Sigma_Cm,1,5^2)); % if r==1 then there is a problem with Sigma_C
if ~r
    Sigma_start = Sigma_Cm;   
else
    Sigma_start = nearestSPD(Sigma_Cm);      
end
% cont.mit.CV_max = 2.1; 
% CV_Cm = cont.mit.CV_old;                    
mit_Cm = struct('mu',mu_Cm,'Sigma',reshape(Sigma_start,1,length(mu_Cm)^2),'df', df, 'p', 1);
[mit_Cm, CV_Cm] = MitISEM_new2(mit_Cm, kernel, mu_init, cont, GamMat);   
         

    
%% %%%% SAVE    
name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'_DUPA.mat'];
if exist(name,'file')
   save(name,'-regexp','^mit',...
        '-append');
else
   save(name,'-regexp','^mit','^CV');    
end