load(name,'mu_MLE','mu_C','draw','draw_C',...
    'mean_draw','mean_draw_C','median_draw','median_draw_C','std_draw','std_draw_C')
%% delta = [8.5, 0.33, 80, 0.001, 0.001];
% 0.0628    0.1433    0.0817    0.5160    0.5361
% delta = [8.5, 0.33, 20, 0.003, 0.003];
% 0.0619    0.1365    0.2723    0.3028    0.3224
% delta = [5.5, 0.33, 20, 0.003, 0.003];
% 0.0955    0.1391    0.2125    0.3819    0.3967
% delta = [3.5, 0.33, 20, 0.003, 0.003];
% 0.1585    0.1375    0.2204    0.3768    0.3934
% delta = [2.0, 0.33, 15, 0.003, 0.003];
% 0.2647    0.1441    0.2339    0.3937    0.4252
% delta = [2.0, 0.2, 15, 0.003, 0.004];
% 0.2570    0.2189    0.3190    0.3040    0.2701
%%
delta = [1.5, 0.15, 15, 0.003, 0.004];
% 0.3188    0.2858    0.2918    0.3453    0.3164
kernel = @(xx) posterior_t_garch11_mex(xx, y(1:T), y_S, GamMat, hyper);      
[mean_theta, median_theta, std_theta, mean_accept, Draw_MH] = RWMH_tgarch11(kernel, ...
    @(xx)prior_t_garch11(xx,hyper), mu_MLE, delta, M, BurnIn, true);


threshold = sort(y(1:T));
threshold = threshold(round(2*p_bar*T));

%% delta_C = [5.5, 0.33, 50, 0.001, 0.001];
% 0.5353    0.3498    0.1306    0.7438    0.7522
% delta_C = [5.5, 0.33, 30, 0.001, 0.001];
% 0.5488    0.3590    0.2593    0.6701    0.6785
% delta_C = [5.5, 0.33, 20, 0.001, 0.001];
% 0.5161    0.3580    0.2588    0.7302    0.7448
% delta_C = [5.5, 0.33, 20, 0.003, 0.003];
% 0.5443    0.3477    0.3453    0.4507    0.4460
% delta_C = [5.5, 0.33, 20, 0.005, 0.005];
% 0.5311    0.3594    0.3028    0.3927    0.4095
% delta_C = [6.0, 0.33, 20, 0.005, 0.006];
% 0.5242    0.3597    0.3238    0.3632    0.3434
% delta_C = [7.0, 0.33, 20, 0.005, 0.006];
% 0.4850    0.3571    0.3176    0.3554    0.3361
% delta_C = [8.0, 0.33, 20, 0.005, 0.006];
% 0.4306    0.3296    0.3020    0.4375    0.4215
% delta_C = [9.0, 0.33, 20, 0.006, 0.007];
% 0.4405    0.3834    0.2978    0.3070    0.2972
% delta_C = [10.0, 0.33, 20, 0.006, 0.007];
% 0.4088    0.3786    0.2925    0.3088    0.2983
% delta_C = [11.0, 0.33, 20, 0.006, 0.007];
% 0.3403    0.3332    0.2728    0.4439    0.4339
%%
delta_C = [11.0, 0.33, 20, 0.007, 0.008];
%     0.3590    0.3547    0.2759    0.3461    0.3513
   
kernel = @(xx) C_posterior_t_garch11_2_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);
[mean_theta_C, median_theta_C, std_theta_C, mean_accept_C, Draw_MH_C] = RWMH_tgarch11(kernel, ...
    @(xx)prior_t_garch11(xx,hyper), mu_C, delta_C, M, BurnIn, true);


%%
kernel = @(xx) C_posterior_t_garch11_2_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);
BurnInRW = 100;
MRW = 1;
THETA_post = Draw_MH(1:MRW:end,:);

delta_pcp =  [21.0, 0.10, 10, 0.001, 0.001];

partition = 4;
prior =  @(xx)prior_t_garch11(xx,hyper);

[mean_theta_PC, median_theta_PC, std_theta_PC, mean_accept_PC, Draw_MH_PC] = ...
    Partial_RWMH_tgarch11(kernel, prior, partition, THETA_post, delta_pcp, BurnInRW, MRW, true);
  
%% Draw statistics
display(mean_draw);
display(mean_theta);

display(median_draw);
display(median_theta);

display(std_draw);
display(std_theta);



display(mean_draw_C);
display(mean_theta_C);

display(median_draw_C);
display(median_theta_C);

display(std_draw_C);
display(std_theta_C);


load(name,'draw_PC','mean_draw_PC','median_draw_PC','std_draw_PC') 

display(mean_draw_PC);
display(mean_theta_PC);

display(median_draw_PC);
display(median_theta_PC);

display(std_draw_PC);
display(std_theta_PC);

%%  Histograms
figure(244)
for ii = 1:5
    subplot(2,3,ii);
    histogram(draw(:,ii));
    title(['MitISEM Post ',params{ii}],'Interpreter','Latex');
end

figure(245)
for ii = 1:5
    subplot(2,3,ii);
    histogram(draw_C(:,ii));
    title(['MitISEM Censored ',params{ii}],'Interpreter','Latex');
end


figure(344)
for ii = 1:5
    subplot(2,3,ii);
    histogram(Draw_MH(:,ii));
    title(['MH Post ',params{ii}],'Interpreter','Latex');
end

figure(345)
for ii = 1:5
    subplot(2,3,ii);
    histogram(Draw_MH_C(:,ii));
    title(['MH Censored ',params{ii}],'Interpreter','Latex');
end

figure(246)
for ii = 1:5
    subplot(2,3,ii);
    histogram(draw_PC(:,ii));
    title(['MitISEM Partially Censored ',params{ii}],'Interpreter','Latex');
end

figure(346)
for ii = 1:5
    subplot(2,3,ii);
    histogram(Draw_MH_PC(:,ii));
    title(['MH Partially Censored ',params{ii}],'Interpreter','Latex');
end


%%  Trace Plots
figure(224)
for ii = 1:5
    subplot(2,3,ii);
    plot(draw(:,ii));
    title(['MitISEM Post ',params{ii}],'Interpreter','Latex');
end

figure(225)
for ii = 1:5
    subplot(2,3,ii);
    plot(draw_C(:,ii));
    title(['MitISEM Censored ',params{ii}],'Interpreter','Latex');
end


figure(324)
for ii = 1:5
    subplot(2,3,ii);
    plot(Draw_MH(:,ii));
    title(['MH Post ',params{ii}],'Interpreter','Latex');
end

figure(325)
for ii = 1:5
    subplot(2,3,ii);
    plot(Draw_MH_C(:,ii));
    title(['MH Censored ',params{ii}],'Interpreter','Latex');
end