function d = posterior_skt_garch11(theta, y, y_S, hyper)
%         theta = mu_init1 = theta;
%         kernel_init = @(xx) - posterior_skt_garch11(xx, y(1:T), y_S, hyper)/T;
% %         kernel = @(xx) posterior_skt_garch11_mex(xx, y(1:T), y_S, GamMat, hyper);
% kernel_init(mu_init1')
%         [mu_MLE0,~,exitflag,output,~,Sigma] = fminunc(kernel_init,mu_init0); %options);
%         [mu_MLE1,~,exitflag,output,~,Sigma] = fminunc(kernel_init,mu_init1); %options);
% lowerbounds = [-1+1e-6, 2+1e-6, -inf, 1e-8, 1e-8, 1e-8];
% upperbounds = [1-1e-6, inf, inf, 1, 1, 1];
% A = [0, 0, 1, 1, 0, 0];
% b = 1;
% [mu_MLE0, LLF1, EXITFLAG, OUTPUT, GRAD, HESSIAN] = fmincon(kernel_init, mu_init1 ,A, b, [],[],lowerbounds, upperbounds, [], options);
      
% theta = [0, 5, 0, 1, 0.1, 0.8]
% mu_MLE1 =  -0.3360    5.0052    0.0378    1.0138    0.1092    0.8908;
% mu_init = [0.5, 8, 1, 3, 0.14, 0.85]
    T = length(y);
    M = size(theta,1);
    
    lambda = theta(:,1);
    nu = theta(:,2);
    mu = theta(:,3);
    omega = theta(:,4);
    alpha = theta(:,5);
    beta = theta(:,6);

    
    prior = prior_skt_garch(M, lambda, nu, omega, alpha, beta, hyper);   

    pr = logical(prior(:,1));
%     rho = (nu-2)./nu;


%     c = gamma((nu+1)/2)/(sqrt(pi*(nu-2))*gamma(nu/2));
%     b = sqrt(1 + 3*lamda^2 - a^2);
    logc = NaN*ones(M,1);
    c = NaN*ones(M,1);
    a = NaN*ones(M,1);
    logb = NaN*ones(M,1);
    b = NaN*ones(M,1);
    tau = NaN*ones(M,1);
     
    logc(pr) = gammaln((nu(pr)+1)/2) - gammaln(nu(pr)/2) - 0.5*log(pi*(nu(pr)-2));
    c(pr) = exp(logc(pr));
    a(pr) = 4.*lambda(pr).*c(pr).*((nu(pr)-2)./(nu(pr)-1));
%     loga = log(a);
%     loga = log(4) + log(lambda) + logc + log(nu-2) - log(nu-1);
%     a = exp(loga);
    logb(pr) = 0.5.*log(1 + 3.*lambda(pr).^2 - a(pr).^2);    
    b(pr) = exp(logb(pr));
    
    tau(pr) = - a(pr)./b(pr);
    
    CONST = T*(logb + logc);
    
    if (y_S == 0) 
        h = omega;
    else        
        h = y_S*ones(M,1);
    end    
    
    d = -Inf*ones(M,1);

    for ii = 1:M        
        if prior(ii,1) % compute only for valid draws
%             scale = sqrt(rho(ii,1)*h(ii,1));
            scale = sqrt(h(ii,1));
            z = (y(1,1) - mu(ii,1))/scale;
%             d(ii,1) = log(tpdf(z,nu(ii))/scale);
            d(ii,1) = sktpdf(z, nu(ii,1), lambda(ii,1), a(ii,1), b(ii,1), tau(ii,1)) - log(scale);
            
            for jj = 2:T
                h(ii,1) = omega(ii,1)*(1-alpha(ii,1)-beta(ii,1)) + alpha(ii,1)*(y(jj-1,1)-mu(ii,1)).^2 ...
                    + beta(ii,1)*h(ii,1);
%                 scale = sqrt(rho(ii,1)*h(ii,1));
                scale = sqrt(h(ii,1));
                z = (y(jj,1)-mu(ii,1))/scale;
%                 d(ii,1) = d(ii,1) + log(tpdf(z,nu(ii,1))/scale);
                d(ii,1) = d(ii,1) + sktpdf(z, nu(ii,1), lambda(ii,1), a(ii,1), b(ii,1), tau(ii,1)) - log(scale);
                
            end  
        end
    end
    d = CONST + d + prior(:,2);
end

function lpdf = sktpdf(x, nu, lambda, a, b, tau)
% LOG  **KERNEL** OF THE SK-T DENSITY (no constant for speed)
% x = z;
 
%         indicator1 = ((data(t)-mu(t))./sqrt(h(t))<-a./b);
%         indicator2 = ((data(t)-mu(t))./sqrt(h(t))>=-a./b);
    indicator1 = (x < tau);
    indicator2 = (x >= tau);
    
%     likelihoods1 =  - ((nu+1)./2).*log(1 + (((b.*x(indicator1) + a)./(1-lambda)).^2)./(nu-2));
%     likelihoods2 =  - ((nu+1)./2).*log(1 + (((b.*x(indicator2) + a)./(1+lambda)).^2)./(nu-2)); 
% %         likelihoods1 = log(b) + log(c) - ((nu+1)./2).*log(1+1./(nu-2).*((b.*indicator1.*((data(t)-mu(t))./sqrt(h(t)))+a)./(1-lamda)).^2);
% %         likelihoods2 = log(b) + log(c) - ((nu+1)./2).*log(1+1./(nu-2).*((b.*indicator2.*((data(t)-mu(t))./sqrt(h(t)))+a)./(1+lamda)).^2);         
% %         likelihoods = - 0.5*log(h(t)) + indicator1.*likelihoods1 + indicator2.*likelihoods2; 
%     lpdf = sum(likelihoods1) + sum(likelihoods2);    

    lpdf = NaN(size(x));
    lpdf(indicator1) = - ((nu+1)./2).*log(1 + (((b.*x(indicator1) + a)./(1-lambda)).^2)./(nu-2));
    lpdf(indicator2) = - ((nu+1)./2).*log(1 + (((b.*x(indicator2) + a)./(1+lambda)).^2)./(nu-2)); 
end

function R = prior_skt_garch(M, lambda, nu, omega, alpha, beta, hyper)
    % uniform prior on alpha and beta on (0,1)
    % with restriction alpha + beta < 1
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val at the corresponding point

    c1 = ((alpha >= 0) & (alpha < 1) & (beta >= 0) & (beta < 1));
    c2 = (alpha + beta < 1);
    c3 = (omega > 0);
    c4 = (nu > 2);
    c5 = ((lambda > -1) & lambda < 1);
    
    r1 = (c1 & c2 & c3 & c4 & c5);

    r2 = -Inf*ones(M,1);
    r2(r1==true) = log(0.5) + log(hyper) - hyper*(nu(r1==true) - 2);
    
    R = [r1, r2]; %prior =  [r1, r2];
end
