clear all
close all
addpath(genpath('include/'));
load('results\ar1\ar1_1_2_1000_PCP.mat','y','mit1')

% draw = [0, 0 ,0.85];
d1 = 1;
d2 = 2;
draw1 = 0.8; % TRUE marginal draw
% draw1 = 0.85; % marginal draw 

%% Joint mixture parameters
Sigma = reshape(mit1.Sigma(1,:),d1+d2,d1+d2);
Sigma_11 = inv(Sigma(3,3));     % inv(0.003);
Sigma_12 = Sigma(3,1:2);        % [0.012,0.004];
Sigma_22 = Sigma(1:2,1:2);      % [0.1, 0.04; 0.04, 0.03];

mu_1 = mit1.mu(1,3);            % 0.8;
mu_2 = mit1.mu(1,1:2);          % [0.3, 2];
df_1 = mit1.df(1,1);            % 5;

%% Conditional parameters
df_2 = df_1 + d1;
mu_temp = (draw1 - mu_1);
mu_2_cond = mu_2 + mu_temp*Sigma_11*Sigma_12;
Sigma_22_cond =  Sigma_22 - Sigma_12'*Sigma_11*Sigma_12;
Sigma_22_cond_corr = Sigma_22_cond/df_2;
Sigma_22_cond_corr = Sigma_22_cond_corr*(df_1 + mu_temp*Sigma_11*mu_temp);

mit_cond.mu = mu_2_cond;
mit_cond.Sigma = reshape(Sigma_22_cond_corr,1,d2*d2);
mit_cond.df = df_2;
mit_cond.p = 1;

% Generation
M = 10000;
BurnIn = 1000;

draw = randn(M+BurnIn,d2);
C = chol(Sigma_22_cond_corr);
draw = draw*C;
W = df_2./chi2rnd(df_2,M+BurnIn,1);
W = sqrt(W);
draw = bsxfun(@times,draw,W);
draw = bsxfun(@plus,draw,mu_2_cond);
       
% Logkernel evaluation
T = length(y);
p_bar = 0.05;
threshold = sort(y);
threshold = threshold(2*p_bar*T);
            
% kernel = @(xx) C_posterior_ar1(xx, y, threshold);
kernel = @(xx) C_posterior_ar1_mex(xx, y, threshold);

draw_aug = [draw, ones(M+BurnIn,1)*draw1];
% tic
% lnk = kernel(draw_aug);
% toc
% tic
lnk = kernel(draw_aug);
% toc

%% For 200 draws
% Laptop -->   515.5555
% Elapsed time is 114.429607 seconds.
% Elapsed time is 0.221954 seconds.
% PC --> 8.35
% Elapsed time is 0.272951 seconds.
% Elapsed time is 0.032703 seconds.
% For 11000 draws
% PC --> 9.46
% Elapsed time is 14.827212 seconds.
% Elapsed time is 1.567201 seconds.

%% Logcandidate evaluation
x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);

% logdenisty: conditional + marginal 
lnd_c = dmvgt_mex(draw, mit_cond.mu, mit_cond.Sigma, mit_cond.df, 1, GamMat, double(1));
lnd_m =  dmvgt_mex(draw1, mu_1, inv(Sigma_11), df_1, 1, GamMat, double(1));
lnd = lnd_c + lnd_m;

% logdensity: joint
mit.mu = [mu_2,mu_1];
mit.Sigma = reshape([Sigma_22, Sigma_12';Sigma_12,inv(Sigma_11)],1,(d1+d2)^2);
mit.df = df_1;
mit.p = 1;

lnd_j = dmvgt_mex(draw_aug, mit.mu, mit.Sigma, mit.df, 1, GamMat, double(1));
% lnd should be equal to lnd_j

% Logweights and independent MH
lnw = lnk - lnd;
lnw = lnw - max(lnw);
[ind, a] = fn_MH(lnw);
accept = a/(M+BurnIn);

draw_aug = draw_aug(ind,:);
draw_aug = draw_aug(BurnIn+1:BurnIn+M,:);

subplot(2,1,1)
hist(draw_aug(:,1))
subplot(2,1,2)
hist(draw_aug(:,2))

%% MEX debug
% x = -5:0.1:5;
% cdf1 = normcdf(x,0,1);
% cdf2 = normcdf_my(x,0,1);
% cdf3 = normcdf_my2(x,0,1);
% x = x';
% cdf4 = normcdf_my_mex(x,0,1);
% cdf4 = normcdf_my_mex(0,0,1);

% theta = [0,2,0.8;0.1,1,0.9];
% threshold = 0.3;
% addpath(genpath('include/'));
% load('results\ar1\ar1_1_2_1000_PCP.mat','y')
% mex C_posterior_ar1_mex
[d_mex, T_mex] = C_posterior_ar1_mex(theta, y, threshold);
[d, T] = C_posterior_ar1(theta, y, threshold);