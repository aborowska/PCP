% function results = PCP_skt_gas_empirical(sdd, data, p_bar0, p_bar1, p_bar, T, H, M, BurnIn, ...
%     mu_init, df, cont, options, partition, II, GamMat)
% addpath(genpath('include/'));
% 
% clear all
% sdd = 1;
% plot_on = false;
% save_on = false; %true;
addpath(genpath('include/'));

clear% all

save_on = false;
plot_on = true;

sdd = 1;        
s = RandStream('mt19937ar','Seed',sdd);
RandStream.setGlobalStream(s); 
 
data = 103; %104; %44; % MSFT 4; % GSPC 2;
H = 1500; %1000; %2275 +248;

model = 'skt_gas';
figures_path = ['figures/',model,'/'];
parameters = {'$\\lambda$','$\\nu$','$\\mu$','$\\omega$','$A$','$B$'};
params = {'$\lambda$','$\nu$','$\mu$','$\omega$','$A$','$B$'};

% quantiles of interest
p_bar0 = 0.005;
p_bar1 = 0.01;
p_bar = 0.05; 
P_bars = [p_bar0, p_bar1, p_bar];

%theta[i] <= lambda
%  theta[i+N] <= nu
%  theta[i+2*N] <= mu
%  theta[i+3*N] <= omega
%  theta[i+4*N] <= alpha
%  theta[i+5*N] <= beta 
mu_init = [0, 8,0,1,0.1,0.8];
D = length(mu_init);
x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);

M = 10000; % number of draws 
BurnIn = 1000;

cont = MitISEM_Control;
cont.mit.CV_max = 1; %2?
cont.mit.iter_max = 10;
cont.mit.Hmax = 3; %6;
cont.mit.dfnc = 3;%5;
%     df = 5; % default df for a mit
df = 3; % default df for a mit
cont.disp = true;

options = optimset('Display','off');
% w = warning('query','last');
% id = w.identifier;
id = 'optim:fminunc:SwitchingMethod';
warning('off',id);

%         thinning = 10;
BurnIn_PCP = BurnIn/10; 
II = 10;
thinning = 1;

%% Load data 
[y, T, y_plot, data_name, time] = Load_data_empirical(data, H);  
load('IBM_new500rets.mat')
time = [time(1), 2019];
H2 = length(IBM_new);
H = H + H2;
y = [y; IBM_new];
data_name = [data_name,'_up'];

% data_stats = Plot_data(y,H, time,true,figures_path)


y_sort = sort(y(1:T));
THR_emp = y_sort(round(P_bars*T));
y_S = var(y(1:T));

hyper = 0.01; 

fprintf('Data loaded: %s\n',data_name);
 
  
        
%% Uncensored Posterior
fprintf('*** Uncensored Posterior ***\n');
%[lambda, nu, mu, omega, A, B]
theta0 = [-0.1, 6, 0, 0.7, 0.1, 0.8]
% kernel_init(theta0)
kernel_init = @(xx) - posterior_skt_gas_mex(xx, y(1:T), hyper, y_S)/T;
kernel = @(xx) posterior_skt_gas_mex(xx, y(1:T), hyper, y_S);
[mu_MLE1,~,exitflag,output,~,Sigma] = fminunc(kernel_init,mu_init,options);
[mu_MLE,~,exitflag,output,~,Sigma] = fminunc(kernel_init,mu_MLE1,options);


Sigma = inv(T*Sigma); 
r = fn_testSigma(reshape(Sigma,1,D^2)); % if r==1 then there is a problem with Sigma_C
if ~r
    Sigma_start = Sigma;   
else
    Sigma_start = nearestSPD(Sigma);      
end


prior =  @(xxx) prior_skt_gas(xxx, y, hyper, y_S);
delta = [  0.0300    0.030    0.1000    0.0200    0.1000    0.0200];
THETA = RWMH_skt_gas(kernel, prior, mu_MLE, delta, 50000, 10000, true) ;
mu_MLE = mean(THETA);
Sigma_start = cov(THETA);

%     cont.mit.CV_tol = 0.3;
cont.mit.CV_max = 1.9;

CV = cont.mit.CV_old;
while (CV(end) >= 2)
    try
%             [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
        [mit, CV] = MitISEM_new2(kernel_init, kernel, mu_init, cont, GamMat);
        [draw, accept] = IndMH_mit(mit, kernel,M,BurnIn,GamMat);
    catch
%         draw = rmvt(mu_MLE,Sigma,df,M+BurnIn);
        mit = struct('mu',mu_MLE,'Sigma',reshape(Sigma_start,1,D^2),'df', df, 'p', 1);
%                 [mit, CV] = MitISEM_new(mit, kernel, mu_init, cont, GamMat);            
        [mit, CV] = MitISEM_new2(mit, kernel, mu_init, cont, GamMat);         
        [draw, accept] = IndMH_mit(mit, kernel,M,BurnIn,GamMat);
        % mean   0.0089    6.9275    0.0938    0.0153    0.2761    0.9773
    end
end
% 
% for ii=1:D
%     subplot(2,3,ii)
%     hold on
%     histogram(draw(:,ii),50)
%     histogram(THETA(1:M,ii),50)
% end


% compute the implied volatility for the last in-sample period
fT = volatility_skt_gas(draw,y(1:T),0);  

% predictive densities
dens_post = predictive_dens_skt_gas(y(T:(T+H)), fT, draw);
 % predicitve cdfs, constant threshold for different tails
cdf_post = predictive_cdf_skt_gas(y(T:(T+H)), fT, draw, THR_emp);
C_score_post = C_ScoringRule(dens_post, cdf_post, y((T+1):(T+H)), THR_emp);

QUANT_MLE = skewtinv(P_bars,mu_MLE(2),mu_MLE(1))';
fT_MLE = volatility_skt_gas(mu_MLE,y(1:T),0);
f_MLE = volatility_path_skt_gas(mu_MLE, y);

[THR_MLE, cond_MLE] = threshold_skt_gas_varc_mle(y((T+1):(T+H)), mu_MLE, ...
    y(T), fT_MLE, QUANT_MLE);

cdf_v_post = predictive_cdf_skt_gas(y(T:(T+H)), fT, draw, THR_MLE);       
Cv_score_post = C_ScoringRule(dens_post, cdf_v_post, y((T+1):(T+H))', THR_MLE);    

results0.mit = mit;
results0.CV = CV;
results0.mu_MLE = mu_MLE;
results0.Sigma = Sigma;

results0.draw = draw;        
results0.accept = accept;
results0.dens = dens_post;
results0.cdf = cdf_post;
results0.cdf_v = cdf_v_post;
results0.C_score = C_score_post;
results0.Cv_score = Cv_score_post;
results0.fT = fT;        
  

if false
    ff = figure(123);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.44 0.4]);
    set(gcf,'defaulttextinterpreter','latex');
    xx = linspace(time(1,1),time(1,2),H);
    hold on
    plot(xx,y((T+1):(T+H)),'k')
    plot(xx,THR_MLE(1,:),'color',[0.6350    0.0780    0.1840])
    plot(xx,THR_MLE(2,:),'color',[0.9290    0.6940    0.1250])
    plot(xx,THR_MLE(3,:),'color',[  0    0.4470    0.7410])
    set(gca,'XLim',time)
    set(gca,'YLim',[-15,10])    
    plotTickLatex2D('FontSize',11);
    ll = legend('Out-of-sample data', '0.5\% MLE thr.','1\% MLE thr.','5\% MLE thr.')
    set(ll,'interpreter','latex','fontsize',11)
    set(gcf,'PaperPositionMode','auto');
    name = [figures_path,'/',model,'_',data_name,'_mle_thresholds.eps'];
    print(ff,name,'-depsc','-r0')
    name = [figures_path,'/',model,'_',data_name,'_mle_thresholds.png'];
    print(ff,name,'-dpng','-r0')    
end
    
 
name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'3.mat'];
load(name,'-regexp','results0','\w*_MLE','^THR');
draw = results0.draw;



%% CP and PCP with 5% threshold 
if arg5  
%     load(name,'results_e');
%     mit_C = results_e.mit_C; CV_C = results_e.CV_C; mu_C = results_e.mu_C; Sigma_C = results_e.Sigma_C;           
    THRES = 0.05;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S);

    results_z = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_z.THRES = THRES;                                              
end  


%% CP and PCP with 10% threshold 
if arg10  
    % mit_C = results_a.mit_C;
    THRES = 0.1;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold 10%% ***\n');
    kernel_init = @(xx) - C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S);

    results_a = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_a.THRES = THRES;        
%     results_a = results;
end    


%% CP and PCP with 10% threshold 
if arg15  
    % mit_C = results_a.mit_C;
    THRES = 0.15;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold 10%% ***\n');
    kernel_init = @(xx) - C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S);

    results_a2 = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_a2.THRES = THRES;        
%     results_a = results;
end  

%% CP and PCP with 20% threshold 
if arg20  
%     load(name,'results_b');
%     mit_C = results_b.mit_C; CV_C = results_b.CV_C; mu_C = results_b.mu_C; Sigma_C = results_b.Sigma_C;
    THRES = 0.2;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper,y_S);

    results_b = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_b.THRES = THRES;                                              
end

%% CP and PCP with 20% threshold 
if arg25  
%     load(name,'results_b');
%     mit_C = results_b.mit_C; CV_C = results_b.CV_C; mu_C = results_b.mu_C; Sigma_C = results_b.Sigma_C;
    THRES = 0.25;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper,y_S);

    results_b2 = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_b2.THRES = THRES;                                              
end

%% CP and PCP with 30% threshold 
if arg30 
%     load(name,'results_c');
%     mit_C = results_c.mit_C; CV_C = results_c.CV_C; mu_C = results_c.mu_C; Sigma_C = results_c.Sigma_C;    
    THRES = 0.3;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S);

    results_c = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_c.THRES = THRES;                                              
end

%% CP and PCP with 40% threshold 
if arg40  
%     load(name,'results_d');
%     mit_C = results_d.mit_C; CV_C = results_d.CV_C; mu_C = results_d.mu_C; Sigma_C = results_d.Sigma_C;       
    THRES = 0.4;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_gas(xx, y(1:T,1), threshold, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S);

    results_d = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_d.THRES = THRES;                                              
end    

%% CP and PCP with 50% threshold 
if arg50  
%     load(name,'results_e');
%     mit_C = results_e.mit_C; CV_C = results_e.CV_C; mu_C = results_e.mu_C; Sigma_C = results_e.Sigma_C;           
    THRES = 0.5;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_mex(xx, y(1:T,1), threshold, hyper, y_S);

    results_e = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_e.THRES = THRES;                                              
end    
    
%% CENSORED MLE PARAMETERS

prior = @(xx) prior_skt_gas(xx, y, hyper, y_S);

if arg100
    THRES = 0.1;
    quantile = skewtinv(THRES, mu_MLE(2), mu_MLE(1));
    fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',THRES);
    kernel_init = @(xx) - C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S);

    results_m_a = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_m_a.THRES = THRES;    
%     results_m_a = results;
end


if arg150
    THRES = 0.15;
    quantile = skewtinv(THRES, mu_MLE(2), mu_MLE(1));
    fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',THRES);
    kernel_init = @(xx) - C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S);

    results_m_a2 = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_m_a2.THRES = THRES;    
%     results_m_a = results;
end



if arg200
    THRES = 0.2;
    quantile = skewtinv(THRES, mu_MLE(2), mu_MLE(1));
    fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',THRES);
    kernel_init = @(xx) - C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S);

    results_m_b = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_m_b.THRES = THRES;                                              
end



if arg250
    THRES = 0.25;
    quantile = skewtinv(THRES, mu_MLE(2), mu_MLE(1));
    fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',THRES);
    kernel_init = @(xx) - C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S);

    results_m_b2 = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_m_b2.THRES = THRES;                                              
end


if arg300
    THRES = 0.3; 
    quantile = skewtinv(THRES, mu_MLE(2), mu_MLE(1));
    fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',THRES);
    kernel_init = @(xx) - C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S);
    
    results_m_c = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, 5, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_m_c.THRES = THRES;                                                                     
end
   
if arg400
%     load(name,'results_m_d');
%     mit_C = results_m_d.mit_C; CV_C = results_m_d.CV_C; mu_C = results_m_d.mu_C; Sigma_C = results_m_d.Sigma_C;           
        
    THRES = 0.4; 
    quantile = skewtinv(THRES, mu_MLE(2), mu_MLE(1));
    fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',THRES);
    kernel_init = @(xx) - C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S);

    results_m_d = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, , cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_m_d.THRES = THRES;                                                                     
end
      

if arg500
%     load(name,'results_m_e');
%     mit_C = results_m_e.mit_C; CV_C = results_m_e.CV_C; mu_C = results_m_e.mu_C; Sigma_C = results_m_e.Sigma_C;           
     
    THRES = 0.5; 
    quantile = skewtinv(THRES, mu_MLE(2), mu_MLE(1));
    fprintf('*** Censored Posterior, MLE time varying threshold, THR = %3.2f ***\n',THRES);
    kernel_init = @(xx) - C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S)/T;    
    kernel = @(xx) C_posterior_skt_gas_varc_mle_mex(xx, y(1:T,1), ...
        mu_MLE, f_MLE, quantile, hyper, y_S);

    
    results_m_e = PCP_CP_skt_gas_estimate_evaluate(kernel, kernel_init, y, T, H, draw, df, cont, GamMat, mu_init, options,...
                           M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_m_e.THRES = THRES;                                                                     
end


name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'_3.mat'];
load(name);
draw = results0.draw;
fT = results0.fT;
results_a = PCP_CP_skt_gas_evaluate(results_a, y, T, H, THR_emp, THR_MLE);
results_b = PCP_CP_skt_gas_evaluate(results_b, y, T, H, THR_emp, THR_MLE); 
results_c = PCP_CP_skt_gas_evaluate(results_c, y, T, H, THR_emp, THR_MLE);
results_d = PCP_CP_skt_gas_evaluate(results_d, y, T, H, THR_emp, THR_MLE); 
results_e = PCP_CP_skt_gas_evaluate(results_e, y, T, H, THR_emp, THR_MLE);
% results_m_d = PCP_CP_skt_gas_evaluate(results_m_d, [y;IBM_new], T, H+H2, THR_emp, THR_MLE2); 
results_m_e = PCP_CP_skt_gas_evaluate(results_m_e, y, T, H, THR_emp, THR_MLE);

        
    %% Results
    if save_on
        name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'_3.mat'];
%         name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'2.mat'];
        if exist(name,'file') %exits_on
            save(name,'-regexp','^results_','-append');%,'\w*_MLE','^THR','-append');
        else
            save(name,'-regexp','^results','\w*_MLE','^THR');
        end
    end
% end