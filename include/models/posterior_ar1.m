function [d, N] = posterior_ar1(theta, y)
    N = length(y);  
    M = size(theta,1);    

    mu = theta(:,1);
    sigma = theta(:,2);
    rho = theta(:,3);
  
    a1 = mu./(1-rho);
    P1 = (sigma.^2)./(1-rho.^2);
 
    d = -Inf*ones(M,1);
    r = (sigma>0) & (abs(rho)<1); % prior satisfied
    
    d(r) = log(P1(r)) + ((y(1,1) - a1(r)).^2)./P1(r); % sationary distribution for the first observation
    d(r) = d(r) + N*log(2*pi) + (N-1)*log(sigma(r).^2); % Gaussian constants
    
    for ii = 1:M
        if r(ii)
            for jj = 2:N
                d(ii) = d(ii) + ((y(jj,1) - mu(ii,1) - rho(ii,1).*y(jj-1,1)).^2)./(sigma(ii,1).^2);                      
            end
        end
    end
    d = -0.5*d;
    
%      (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);
    prior = - log(beta(20, 1.5)) + (20-1)*log((rho+1)/2) + (1.5-1)*log(1-(rho+1)/2); 
    prior = prior - log(sigma);
    d(r) = d(r) + prior(r);
end
   