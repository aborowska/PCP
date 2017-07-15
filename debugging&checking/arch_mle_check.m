% ARCH: ML estimation MSE (unbiasedness)

clear all
close all

addpath(genpath('include/'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 


model = 'arch1'; 
partition = 3;
p_bar = 0.05;
% parameters = {'$\\mu$','$\\gamma$','$\\omega$','$\\alpha$'};
parameters = {'$\\mu$','$\\omega$','$\\mu2$','$\\alpha$'};

TRANS = 0;

for sigma2 = [1]
for T = [1000]

sigma1 = 1;
% sigma2 = 1;
c = (sigma2 - sigma1)/sqrt(2*pi); % mean of eps
kappa = 0.5*(sigma1^2 + sigma2^2 - ((sigma2-sigma1)^2)/pi); % var of eps
sigma1_k = sigma1/sqrt(kappa);
sigma2_k = sigma2/sqrt(kappa);

mu2 = 0; % gama = mu - mu2;
omega = 1;
alpha = 0.1; 

% theta  = [mu, omega, mu2, alpha]
mu_true = [0, omega, 0, alpha];
mu_init = [0, 1, 0.1, 0.05];


% T = 2500; % time series length

options = optimset('Display','off');
% w = warning('query','last');
% id = w.identifier;
id = 'optim:fminunc:SwitchingMethod';
warning('off',id);

S = 100;
MUS = zeros(S,4);
MUS_mex = zeros(S,4);

% MUS_omega = zeros(S,4);

% MUS_C_y_S = zeros(S,4);
% MUS_C_omega = zeros(S,4);

% MUS_C0_y_S = zeros(S,4);
% MUS_C0_omega = zeros(S,4);

for ss = 1:S

    % ARCH(1)
    eps = randn(T,1);
    ind = (eps>0);
    eps(ind) = c + sigma1.*eps(ind);
    eps(~ind) = c + sigma2.*eps(~ind);
    eps = eps/sqrt(kappa);  
    y = zeros(T,1);
    h_true = zeros(T,1);

    for ii = 1:T
        if (ii == 1)
            h_true(ii,1) = omega;
        else
            h_true(ii,1) = omega*(1-alpha) + alpha*(y(ii-1,1))^2;
        end
        y(ii,1) = sqrt(h_true(ii,1))*eps(ii,1);
    end

    % (Misspecified) model: GARCH(1,1) normal 
    % Uncensored Posterior
    y_S = var(y);
    if TRANS
        kernel_init = @(xx) -posterior_arch1_mex(transform_param_arch(xx, 'back'), y, y_S)/T;
        [mu,~,~,~,~,Sigma] = fminunc(kernel_init, transform_param_arch(mu_init, 'opt'),options);
        Sigma = inv(T*Sigma);
        MUS_y_S(ss,:) = transform_param_agarch(mu, 'back');
    else    
        kernel_init = @(xx) -posterior_arch1_mex(xx, y, y_S)/T;
        [mu,~,~,~,~,Sigma] = fminunc(kernel_init, mu_init,options);
        Sigma = inv(T*Sigma);
        MUS_mex(ss,:) = mu;
    end
    
    kernel_init = @(xx) -posterior_arch1(xx, y, y_S)/T;
    [mu,~,~,~,~,Sigma] = fminunc(kernel_init,mu_init,options);
    Sigma = inv(T*Sigma);
    MUS(ss,:) = mu;
     
%     % CENSORED at 10th percentile
%     threshold = sort(y(1:T));
%     threshold = threshold(2*p_bar*T);
% 
%     if TRANS
%         kernel_init = @(xx) - C_posterior_arch1_mex(transform_param_arch(xx, 'back'), y(1:T,1), threshold, y_S)/T;    
%         [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init, transform_param_arch(mu_init, 'opt'),options);
%         Sigma_C = inv(T*Sigma_C);
%         MUS_C_y_S(ss,:) = transform_param_arch(mu_C, 'back'); 
%     else    
%         kernel_init = @(xx) - C_posterior_arch1_mex(xx, y(1:T,1), threshold, y_S)/T;    
%         [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init, mu_init,options);
%         Sigma_C = inv(T*Sigma_C);
%         MUS_C_y_S(ss,:) = mu_C; 
%     end
% %     kernel_init = @(xx) - C_posterior_arch1_mex(xx, y(1:T,1), threshold, 0)/T;    
% %     [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
% %     Sigma_C = inv(T*Sigma_C);
% %     MUS_C_omega(ss,:) = mu_C; 
%     
%     % CENSORED at 0
%     threshold0 = 0;
%     if TRANS
%         kernel_init = @(xx) - C_posterior_arch1_mex(transform_param_arch(xx, 'back'), y, threshold0, y_S)/T;    
%         [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init, transform_param_arch(mu_init, 'opt'),options);
%         Sigma_C0 = inv(T*Sigma_C0);
%         MUS_C0_y_S(ss,:) = transform_param_agarch(mu_C0, 'back');
%     else 
%         kernel_init = @(xx) - C_posterior_arch1_mex(xx, y, threshold0, y_S)/T;    
%         [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init, mu_init,options);
%         Sigma_C0 = inv(T*Sigma_C0);
%         MUS_C0_y_S(ss,:) = mu_C0;     
%     end
% %     kernel_init = @(xx) - C_posterior_arch1_mex(xx, y, threshold0, 0)/T;    
% %     [mu_C0,~,~,~,~,Sigma_C0] = fminunc(kernel_init,mu_init,options);
% %     Sigma_C0 = inv(T*Sigma_C0);
% %     MUS_C0_omega(ss,:) = mu_C0;  
% end
% 
% if (sigma1 == sigma2)
%     RMSE_y_S = sqrt(sum(bsxfun(@minus,MUS_y_S,mu_true).^2)/S);
% %     RMSE_omega = sqrt(sum(bsxfun(@minus,MUS_omega,mu_true).^2)/S);
%     RMSE_C_y_S = sqrt(sum(bsxfun(@minus,MUS_C_y_S,mu_true).^2)/S);
% %     RMSE_C_omega = sqrt(sum(bsxfun(@minus,MUS_C_omega,mu_true).^2)/S);
%     RMSE_C0_y_S = sqrt(sum(bsxfun(@minus,MUS_C0_y_S,mu_true).^2)/S);
% %     RMSE_C0_omega = sqrt(sum(bsxfun(@minus,MUS_C0_omega,mu_true).^2)/S);
end

ff = figure(1);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

subplot(2,3,1)
hist([MUS_y_S(:,1),MUS_C_y_S(:,1),MUS_C0_y_S(:,1)])
% hist([MUS_y_S(:,1),MUS_omega(:,1),MUS_C_y_S(:,1),MUS_C_omega(:,1),MUS_C0_y_S(:,1),MUS_C0_omega(:,1)])
xlabel('\mu')

subplot(2,3,2)
hist([MUS_y_S(:,2),MUS_C_y_S(:,2),MUS_C0_y_S(:,2)])
% hist([MUS_y_S(:,2),MUS_omega(:,2),MUS_C_y_S(:,2),MUS_C_omega(:,2),MUS_C0_y_S(:,2),MUS_C0_omega(:,2)])
xlabel('\omega')

subplot(2,3,3)
hist([MUS_y_S(:,3),MUS_C_y_S(:,3),MUS_C0_y_S(:,3)])
% hist([MUS_y_S(:,3),MUS_omega(:,3),MUS_C_y_S(:,3),MUS_C_omega(:,3),MUS_C0_y_S(:,3),MUS_C0_omega(:,3)])
xlabel('\mu2')

subplot(2,3,4)
hist([MUS_y_S(:,4),MUS_C_y_S(:,4),MUS_C0_y_S(:,4)])
% hist([MUS_y_S(:,4),MUS_omega(:,4),MUS_C_y_S(:,4),MUS_C_omega(:,4),MUS_C0_y_S(:,4),MUS_C0_omega(:,4)])
xlabel('\alpha')

subplot(2,3,5)
hist([MUS_y_S(:,5),MUS_C_y_S(:,5),MUS_C0_y_S(:,5)])
% hist([MUS_y_S(:,5),MUS_omega(:,5),MUS_C_y_S(:,5),MUS_C_omega(:,5),MUS_C0_y_S(:,5),MUS_C0_omega(:,5)])
xlabel('\beta')


leg = legend({'Post', 'CP 10%', 'CP 0'},...
    'Position',[0.69 0.29 0.15 0.15],'units', 'normalized');
% leg = legend({'Post', 'Post omega', 'CP 10% var', 'CP 10% omega', 'CP 0 var', 'CP 0 omega' },...
%     'Position',[0.69 0.29 0.15 0.15],'units', 'normalized');
subplot(2,3,6)
axis off
if (sigma1 == sigma2)
    text(0.05,0.43,'MLE RMSEs:');
    s1 = sprintf('%6.4f, ',RMSE_y_S);
%     s2 = sprintf('%6.4f, ',RMSE_omega);

    s3 = sprintf('%6.4f, ',RMSE_C_y_S);
%     s4 = sprintf('%6.4f, ',RMSE_C_omega);
    s5 = sprintf('%6.4f, ',RMSE_C0_y_S);
%     s6 = sprintf('%6.4f, ',RMSE_C0_omega);    
    text(0.05,0.36,['y\_S: ',s1]);
%     text(0.05,0.29,['y\_omega: ',s2]);
    text(0.05,0.22,['C y\_S: ',s3]);
%     text(0.05,0.15,['C y\_omega: ',s4]);   
    text(0.05,0.08,['C0 y\_S: ',s5]);
%     text(0.05,0.01,['C0 y\_omega: ',s6]);    
end

if TRANS
    s0 = sprintf('%s, T = %i, sigma2 = %i (TRANS)',model,T,sigma2);
else
    s0 = sprintf('%s, T = %i, sigma2 = %i',model,T,sigma2);
end
suptitle(s0)

if TRANS 
    name = ['figures/',model,'/',model,'_mle_check_TRANS_T',num2str(T),'_',num2str(sigma1),'_',num2str(sigma2),'.png'];
else
    name = ['figures/',model,'/',model,'_mle_check_T',num2str(T),'_',num2str(sigma1),'_',num2str(sigma2),'.png'];
end
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r600');

end
end