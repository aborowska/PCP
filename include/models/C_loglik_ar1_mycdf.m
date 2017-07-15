function [d, T] = C_loglik_ar1_mycdf(theta, y, threshold)
    T = length(y);
    
    mu = theta(:,1);
    sigma = theta(:,2);
    rho = theta(:,3);
    
%   Stationary distibution for the first observation
    a1 = mu./(1-rho);
    P1 = (sigma.^2)./(1-rho.^2);
    
    d = -Inf;
    r = (sigma>0) & (abs(rho)<1); % prior satisfied
    
    if r
        ind = (y<threshold); % observations of interest
        if ind(1,1) 
            d = - 0.5*(log(2*pi) + log(P1) + ((y(1,1) - a1).^2)./P1);
            d = d - 0.5*(sum(ind)-1)*(log(2*pi) + log(sigma.^2)); % Gaussian constant for all the remaining uncensored observations
        else
%             d = log(1-normcdf((threshold-a1)./sqrt(P1)));
            d = log(1-normcdf_my(threshold,a1,sqrt(P1)));
            d = d - 0.5*sum(ind)*(log(2*pi) + log(sigma.^2)); % Gaussian constant for all the uncensored observations
        end

        for jj = 2:T
            if ind(jj,1)
                d = d - 0.5*(((y(jj,1) - mu - rho.*y(jj-1,1)).^2)./(sigma.^2));                      
            else  % P(y_t>=threshold|y_1,...,y_(t-1)) = 1 - P(y_t<threshold|y_1,...,y_(t-1))
%                 d = d + log(1-normcdf((threshold-(mu + rho.*y(jj-1,1)))/sigma));  
                d = d + log(1-normcdf_my(threshold,mu + rho.*y(jj-1,1),sigma));
            end
        end
    end   
end