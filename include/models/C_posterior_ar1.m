function [d, T] = C_posterior_ar1(theta, y, threshold)
    T = length(y);
    M = size(theta,1);
    
    mu = theta(:,1);
    sigma = theta(:,2);
    rho = theta(:,3);
    
%   Stationary distibution for the first observation
    a1 = mu./(1-rho);
    P1 = (sigma.^2)./(1-rho.^2);
    
    d = -Inf*ones(M,1);
    r = (sigma>0) & (abs(rho)<1); % prior satisfied
    
    ind = (y<threshold); % observations of interest
%     ind_T = T - sum(ind); % "censored observations"
%     P = log(1-normcdf(threshold,mu,sigma));

    if ind(1,1) 
        d(r) = - 0.5*(log(2*pi) + log(P1(r)) + ((y(1,1) - a1(r)).^2)./P1(r));
        d(r) = d(r) - 0.5*(sum(ind)-1)*(log(2*pi) + log(sigma(r).^2)); % Gaussian constant for all the remaining uncensored observations
    else
%         d(r) = log(1-normcdf(threshold,a1(r),sqrt(P1(r))));
        d(r) = log(1-normcdf_my(threshold,a1(r),sqrt(P1(r))));
%         d(r) = log(1-normcdf((threshold-a1(r))./sqrt(P1(r))));
        d(r) = d(r) - 0.5*sum(ind)*(log(2*pi) + log(sigma(r).^2)); % Gaussian constant for all the uncensored observations
    end
    
    for ii = 1:M
        if r(ii) % compute only for valid draws
            for jj = 2:T
                if ind(jj,1)
                    d(ii) = d(ii) - 0.5*(((y(jj,1) - mu(ii,1) - rho(ii,1).*y(jj-1,1)).^2)./(sigma(ii,1).^2));                      
                else  % P(y_t>=threshold|y_1,...,y_(t-1)) = 1 - P(y_t<threshold|y_1,...,y_(t-1))
%                     d(ii) = d(ii) + log(1-normcdf(threshold,mu(ii,1) + rho(ii,1).*y(jj-1,1),sigma(ii,1)));                    
                    d(ii) = d(ii) + log(1-normcdf_my(threshold,mu(ii,1) + rho(ii,1).*y(jj-1,1),sigma(ii,1)));
%                     d(ii) = d(ii) + log(1-normcdf((threshold-(mu(ii,1) + rho(ii,1).*y(jj-1,1)))/sigma(ii,1)));  
                end
            end
        end
    end
    
% %      (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);
%     prior = - log(beta(20, 1.5)) + (20-1)*log((rho+1)/2) + (1.5-1)*log(1-(rho+1)/2); 
%     prior = prior - log(sigma);
    prior = - log(sigma);
    d(r) = d(r) + prior(r);
end