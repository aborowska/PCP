% function results = PCP_t_garch11_empirical(sdd, data, p_bar0, p_bar1, p_bar, T, H, M, BurnIn, ...
%     mu_init, df, cont, options, partition, II, GamMat)
addpath(genpath('include/'));
% 
% clear all
% sdd = 1;
% plot_on = false;
% save_on = false; %true;

data = 103; %104; %44; % MSFT 4; % GSPC 2;
H = 1500; %1000; %2275 +248;
  
model = 'skt_agarch11';
figures_path = ['figures/',model,'/'];
parameters = {'$\\lambda$','$\\nu$','$\\mu$','$\\omega$','$\\mu_2$','$\\alpha$','$\\beta$'};
params = {'$\lambda$','$\nu$','$\mu$','$\omega$','$\mu_2$','$\alpha$','$\beta$'};

partition = 5; % INDEX OF THE BEGINNING OF THE UNCENSORED PARAMETERS SET
% quantiles of interest
p_bar05 = 0.005;
p_bar1 = 0.01;
p_bar = 0.05; 
P_bars = [p_bar05, p_bar1, p_bar];

%  theta[i] <= lambda
%  theta[i+N] <= nu
%  theta[i+2*N] <= mu
%  theta[i+3*N] <= omega
%  theta[i+4*N] <= mu2
%  theta[i+5*N] <= alpha
%  theta[i+6*N] <= beta 
mu_init = [0, 5, 0, 1, 0.1, 0.05, 0.85];
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
%     data = 22;
%     H =1523;2000;
[y, T, y_plot, data_name, time] = Load_data_empirical(data, H);  
load('IBM_new500rets.mat')
time = [time(1), 2019];
H2 = length(IBM_new);
H = H + H2;
y = [y; IBM_new];
data_name = [data_name,'_up'];

y_S = var(y(1:T));

y_sort = sort(y(1:T));
THR_emp = y_sort(round(P_bars*T));

hyper = 0.01; 

fprintf('Data loaded: %s\n',data_name);

name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'.mat'];


%% Uncensored Posterior
fprintf('*** Uncensored Posterior ***\n');    
kernel_init = @(xx) - posterior_skt_agarch11(xx, y(1:T), 0, hyper)/T;
kernel = @(xx) posterior_skt_agarch11_mex(xx, y(1:T), y_S, hyper);
% kernel(mu_init)
% kernel2 = @(xx) posterior_skt_agarch11(xx, y(1:T), y_S, hyper);
% kernel2(mu_init)

[mu_MLE1,~,exitflag,output,~,Sigma] = fminunc(kernel_init,mu_init,options);
[mu_MLE,~,exitflag,output,~,Sigma] = fminunc(kernel_init,mu_MLE1,options);
%    -0.0844    5.7338   -0.0268    0.0000    1.2304    0.0877    0.8621
% kernel(mu_MLE1)

Sigma = inv(T*Sigma);  
r = fn_testSigma(reshape(Sigma,1,D^2)); % if r==1 then there is a problem with Sigma_C
if ~r
    Sigma_start = Sigma;   
else
    Sigma_start = nearestSPD(Sigma);      
end

%         [mean_theta, median_theta, std_theta, mean_accept, Draw_MH ] = RWMH_tgarch11(kernel, ...
%             @(xx)prior_t_garch11(xx,hyper), mu_MLE, [8.5, 0.05, 10, 0.001, 0.001], M, BurnIn, 'true');


cont.mit.CV_max = 1.5;
% cont.mit.N = 2*cont.mit.N;
% CV = cont.mit.CV_old;
while (CV(end) >= 2)
    try
%             [mit, CV] = MitISEM_new(kernel_init, kernel, mu_init, cont, GamMat);
        [mit, CV] = MitISEM_new2(kernel_init, kernel, mu_init, cont, GamMat);
        [draw, accept] = IndMH_mit(mit, kernel,M,BurnIn,GamMat);
    catch
%                 draw = rmvt(mu_MLE,Sigma,df,M+BurnIn);
        mit = struct('mu',mu_MLE,'Sigma',reshape(Sigma,1,D^2),'df', df, 'p', 1);
%                 [mit, CV] = MitISEM_new(mit, kernel, mu_init, cont, GamMat);            
        [mit, CV] = MitISEM_new2(mit, kernel, mu_init, cont, GamMat);            
        [draw, accept] = IndMH_mit(mit, kernel,M,BurnIn,GamMat);
    end
end
% draw = results0.draw;
% for ii=1:7; subplot(2,4,ii); histogram(draw(:,ii)); line([mu_MLE(ii),mu_MLE(ii)],[0,1000],'color','r'); end
% mean(draw)
% median(draw)
% compute the implied volatility for the last in-sample period
hT = volatility_skt_agarch11(draw,y(1:T),y_S,0);  

% predictive densities
dens_post = predictive_dens_skt_agarch11(y(T:(T+H)), hT, draw);
 % predicitve cdfs, constant threshold for different tails
cdf_post = predictive_cdf_skt_agarch11(y(T:(T+H)), hT, draw, THR_emp);
C_score_post = C_ScoringRule(dens_post, cdf_post, y((T+1):(T+H)), THR_emp);

% QUANT_MLE = tinv(P_bars, mu_MLE(2))';
% r = sort(skewtrnd(mu_MLE(2),mu_MLE(1),[1000000,1])); 
% Q_r = r(P_bars*length(r));

mu_MLE = mean(draw);

QUANT_MLE = skewtinv(P_bars,mu_MLE(2),mu_MLE(1))';
% P_MLE = skewtcdf(QUANT_MLE,mu_MLE(2),mu_MLE(1))';
hT_MLE = volatility_skt_agarch11(mu_MLE,y(1:T),y_S,0);
[THR_MLE, cond_MLE] = threshold_skt_agarch11_varc_mle(y((T+1):(T+H)), mu_MLE, y(T), hT_MLE, QUANT_MLE);


cdf_v_post = predictive_cdf_skt_agarch11(y(T:(T+H)), hT, draw, THR_MLE);       
Cv_score_post = C_ScoringRule(dens_post, cdf_v_post, y((T+1):(T+H))', THR_MLE);
    
results0.cdf_v = cdf_v_post;
results0.Cv_score = Cv_score_post;
results0.mit = mit;
results0.CV = CV;
results0.mu_MLE = mu_MLE;
results0.Sigma = Sigma;

results0.draw = draw;        
results0.accept = accept;
results0.dens = dens_post;
results0.cdf = cdf_post;
results0.C_score = C_score_post;
results0.fT = hT;  


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
    name = [figures_path,'/',model,'_',data_name,'_mle_thresholds_new.eps'];
    print(ff,name,'-depsc','-r0')
    name = [figures_path,'/',model,'_',data_name,'_mle_thresholds_new.png'];
    print(ff,name,'-dpng','-r0')    
end



% hT = results0.fT;
% dens_post = results0.dens;
% draw = results0.draw;
%% CP and PCP with 10% threshold 
if arg10  
    THRES = 0.1;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold 10%% ***\n');
    kernel_init = @(xx) - C_posterior_skt_agarch11(xx, y(1:T,1), threshold, y_S, hyper)/T;    
    kernel = @(xx) C_posterior_skt_agarch11_mex(xx, y(1:T,1), threshold, y_S, hyper);

    results_a = PCP_CP_skt_agarch_estimate_evaluate(kernel, kernel_init, y, T, H, draw, ...
        df, cont, GamMat, mu_init, options, M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_a.THRES = THRES;
end    

%% CP and PCP with 20% threshold 
if arg20  
    THRES = 0.2;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_agarch11(xx, y(1:T,1), threshold, y_S, hyper)/T;    
    kernel = @(xx) C_posterior_skt_agarch11_mex(xx, y(1:T,1), threshold,  y_S, hyper);

    results_b = PCP_CP_skt_agarch_estimate_evaluate(kernel, kernel_init, y, T, H, draw,...
        df, cont, GamMat, mu_init, options, M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_b.THRES = THRES;

end

%% CP and PCP with 30% threshold 
if arg30  
    THRES = 0.3;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_agarch11(xx, y(1:T,1), threshold, y_S, hyper)/T;    
    kernel = @(xx) C_posterior_skt_agarch11_mex(xx, y(1:T,1), threshold,  y_S, hyper);

    results_c = PCP_CP_skt_agarch_estimate_evaluate(kernel, kernel_init, y, T, H, draw,...
        df, cont, GamMat, mu_init, options, M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_c.THRES = THRES;

end

%% CP and PCP with 40% threshold 
if arg40  
    THRES = 0.4;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_agarch11(xx, y(1:T,1), threshold, y_S, hyper)/T;    
    kernel = @(xx) C_posterior_skt_agarch11_mex(xx, y(1:T,1), threshold,  y_S, hyper);

    results_d = PCP_CP_skt_agarch_estimate_evaluate(kernel, kernel_init, y, T, H, draw,...
        df, cont, GamMat, mu_init, options, M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_d.THRES = THRES;

end

if arg50  
    THRES = 0.5;
    threshold = sort(y(1:T));
    threshold = threshold(round(THRES*T));
    fprintf('*** Censored Posterior, threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_agarch11(xx, y(1:T,1), threshold, y_S, hyper)/T;    
    kernel = @(xx) C_posterior_skt_agarch11_mex(xx, y(1:T,1), threshold,  y_S, hyper);

    results_e = PCP_CP_skt_agarch_estimate_evaluate(kernel, kernel_init, y, T, H, draw,...
        df, cont, GamMat, mu_init, options, M, BurnIn,BurnIn_PCP, THR_emp, THR_MLE);
    results_e.THRES = THRES;

end



results_a = PCP_CP_skt_agarch_evaluate(results_a, y, T, H, THR_emp, THR_MLE);
results_b = PCP_CP_skt_agarch_evaluate(results_b, y, T, H, THR_emp, THR_MLE); 
results_c = PCP_CP_skt_agarch_evaluate(results_c, y, T, H, THR_emp, THR_MLE);
results_d = PCP_CP_skt_agarch_evaluate(results_d, y, T, H, THR_emp, THR_MLE); 
results_e = PCP_CP_skt_agarch_evaluate(results_e, y, T, H, THR_emp, THR_MLE);
% results_m_d = PCP_CP_t_gas_evaluate(results_m_d, [y;IBM_new], T, H+H2, THR_emp, THR_MLE2); 
% results_m_e = PCP_CP_t_gas_evaluate(results_m_e, [y;IBM_new], T, H+H2, THR_emp, THR_MLE2);



if arg200  
    THRES = 0.2;
    threshold = skewtinv(THRES,mu_MLE(2),mu_MLE(1))';
    fprintf('*** Censored Posterior, MLE time varying threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, threshold, y_S, hyper)/T;    
    kernel = @(xx) C_posterior_skt_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, threshold,  y_S, hyper);
      
    results_m_b = PCP_CP_skt_agarch_estimate_evaluate(kernel, kernel_init, y, T, H, draw,...
        df, cont, GamMat, mu_init, options, M, BurnIn, BurnIn_PCP, THR_emp, THR_MLE);
    results_m_b.THRES = THRES;
end


if arg300  
    THRES = 0.3;
    threshold = skewtinv(THRES,mu_MLE(2),mu_MLE(1))';
    fprintf('*** Censored Posterior, MLE time varying threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, threshold, y_S, hyper)/T;    
    kernel = @(xx) C_posterior_skt_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, threshold,  y_S, hyper);
      
    results_m_c = PCP_CP_skt_agarch_estimate_evaluate(kernel, kernel_init, y, T, H, draw,...
        df, cont, GamMat, mu_init, options, M, BurnIn, BurnIn_PCP, THR_emp, THR_MLE);
    results_m_c.THRES = THRES;
end


if arg400  
    THRES = 0.4;
    threshold = skewtinv(THRES,mu_MLE(2),mu_MLE(1))';
    fprintf('*** Censored Posterior, MLE time varying threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, threshold, y_S, hyper)/T;    
    kernel = @(xx) C_posterior_skt_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, threshold,  y_S, hyper);
      
    results_m_d = PCP_CP_skt_agarch_estimate_evaluate(kernel, kernel_init, y, T, H, draw,...
        df, cont, GamMat, mu_init, options, M, BurnIn, BurnIn_PCP, THR_emp, THR_MLE);
    results_m_d.THRES = THRES;
end



if arg500  
    THRES = 0.5;
    threshold = skewtinv(THRES,mu_MLE(2),mu_MLE(1))';
    fprintf('*** Censored Posterior, MLE time varying threshold %3,2f ***\n', THRES);
    kernel_init = @(xx) - C_posterior_skt_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, threshold, y_S, hyper)/T;    
    kernel = @(xx) C_posterior_skt_agarch11_varc_mle_mex(xx, y(1:T,1), mu_MLE, threshold,  y_S, hyper);
      
    results_m_e = PCP_CP_skt_agarch_estimate_evaluate(kernel, kernel_init, y, T, H, draw,...
        df, cont, GamMat, mu_init, options, M, BurnIn, BurnIn_PCP, THR_emp, THR_MLE);
    results_m_e.THRES = THRES;
end










%% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor nu, mu and omega

fprintf('*** Partially Censored Posterior, threshold 10%%, partition = 4 ***\n');
% mit_C: joint candidate for the joint censored posterior
partition = 5;
II = 5;
draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
[draw_PC, a_PC, lnw_PC] = sim_cond_mit_MH_outloop(mit_C, draw_short,...
    partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, 10);
accept_PC = mean(a_PC); 

hT_PC = volatility_skt_agarch11(draw_PC,y(1:T),y_S,0);  


% predictive densities
dens_PC = predictive_dens_skt_agarch11(y(T:(T+H)), hT_PC, draw_PC);      
% predicitve cdfs, constant threshold for different tails
cdf_PC = predictive_cdf_skt_agarch11(y(T:(T+H)), hT_PC, draw_PC, THR_emp);
C_score_PC = C_ScoringRule(dens_PC, cdf_PC, y((T+1):(T+H)), THR_emp);

cdf_v_PC = predictive_cdf_skt_agarch11(y(T:(T+H)), hT_PC, draw_PC, THR_MLE);       
Cv_score_PC = C_ScoringRule(dens_PC, cdf_v_PC, y((T+1):(T+H))', THR_MLE);

results.draw_PC = draw_PC;
results.accept_PC = accept_PC;
results.fT_PC = hT_PC;   
results.dens_PC = dens_PC;
results.cdf_PC = cdf_PC;
results.C_score_PC = C_score_PC;
results.cdf_v_PC = cdf_v_PC;
results.Cv_score_PC = Cv_score_PC;  
    
    
        
    %% Results
    if save_on
        name = ['results/',model,'/',data_name,'/PCP_emp_',model,'_data_',data_name,'_new.mat'];
        if exist(name,'file') %exits_on
            save(name,'-regexp','^results','-append')%,'\w*_MLE','^THR','-append');
        else
            save(name,'-regexp','^results','\w*_MLE','^THR');
        end
    end
% end