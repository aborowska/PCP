clear all
sigma1 = 1;
sigma2 = 2;

%% split normal, mode 0
x = linspace(-5,5,100);
fn_norm_pdf = @(xx,s) exp(-xx.^2/(2*s^2))/sqrt(2*pi*s^2);
A = sqrt(2/pi)/(sigma1+sigma2);
split = [fn_norm_pdf(x(1:50),sigma2),fn_norm_pdf(x(51:100),sigma1)];
plot(x,split)

%% iid simulation, mean 0
N = 10000;

c = 1/sqrt(2*pi);
eps =  randn(N,1);
ind = (eps>0);
eps(ind) = c + sigma1.*eps(ind);
eps(~ind) = c + sigma2.*eps(~ind);
% eps1 = c + sigma1.*eps(eps>0);
% eps2 = c + sigma2.*eps(eps<0);
% eps = [eps1;eps2];
y = eps - mean(eps);

y_sort = sort(y);
VaR_5 = y_sort(0.05*N); % -2.9320

% misspecified model: normal with unknown sigma
% Metropolis-Hastings for the parameters
M = 10000; % number of draws for preliminary and IS computations
BurnIn = 1000;

kernel_init = @(ss) -loglik_iid(ss,y)/N;
[mu,~,~,~,~,Sigma] = fminunc(kernel_init,1);
Sigma = inv(N*Sigma);
df = 5;
draw = rmvt(mu,Sigma,df,M+BurnIn);
kernel = @(ss) posterior_iid(ss,y);
lnk = kernel(draw);

x_gam = (0:0.00001:50)'+0.00001;
GamMat = gamma(x_gam);

lnd = dmvgt_mex(draw, mu, Sigma, df, 1, GamMat, double(1));
lnw = lnk - lnd;
lnw = lnw - max(lnw);
[ind, a] = fn_MH(lnw);
draw = draw(ind,:);
accept = a/(M+BurnIn);
draw = draw(BurnIn+1:BurnIn+M);    
hist(draw,20)    

y_post = sort(draw.*randn(M,1));
VaR_5_post = y_post(0.05*M); % -2.4789

y_2 = sort(2.*randn(M,1));
VaR_5_2 = y_2(0.05*M); %  -3.2886

% Misspecified model: N(0,sigma
% take values below zero (threshold = zero)
threshold = 0;
kernel_init = @(ss) - C_posterior_iid(ss, y, threshold)/N;
[mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,1);
Sigma_C = inv(N*Sigma_C);
draw_C = rmvt(mu_C,Sigma_C,df,M+BurnIn);
kernel = @(ss) C_posterior_iid(ss, y, threshold);
lnk_C = kernel(draw_C);
lnd_C = dmvgt_mex(draw_C, mu_C, Sigma_C, df, 1, GamMat, double(1));
lnw_C = lnk_C - lnd_C;
lnw_C = lnw_C - max(lnw_C);
[ind, a] = fn_MH(lnw_C);
draw_C = draw_C(ind,:);
accept_C = a/(M+BurnIn);
draw_C = draw_C(BurnIn+1:BurnIn+M);    
hist(draw_C,20)    



%% simple AR(1)
rho = 0.8;
y = zeros(N,1);

y(1,1) = eps(1,1);
for ii = 2:N
    y(ii,1) = rho*y(ii-1,1) + eps(ii);
end
