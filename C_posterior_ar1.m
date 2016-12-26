function d = C_posterior_ar1(theta, y, threshold)
    N = length(y);
    M = size(theta,1);
    
    mu = theta(:,1);
    sigma = theta(:,2);
    rho = theta(:,3);
    
    a1 = mu./(1-rho);
    P1 = (sigma.^2)./(1-rho.^2);
    
    d = -Inf*ones(M,1);
    r = (sigma>0) & (abs(rho)<1); 
    
    ind = (y<threshold);
    ind_N = N - sum(ind);
    P = log(1-normcdf(threshold,mu,sigma)); %P(y_t>=threshold|y_1,...,y_(t-1))
    
    d(r) = - 0.5*(ind_N*log(2*pi) + ind_N*log(sigma(r).^2));
    for ii = 1:M
        if r(ii)
            if ind(1,1)
                d(ii) = d(ii) - 0.5*((y(1,1) - a1(ii,1)).^2)./P1(ii,1);
            else
                d(ii) = d(ii) + P(ii,1);
            end
            for jj = 2:N
                if ind(jj,1)
                    d(ii) = d(ii) - 0.5*((y(jj,1) - mu(ii,1) - rho(ii,1).*y(jj-1,1)).^2)./sigma(ii,1).^2;                      
                else
                    d(ii) = d(ii) + P(ii,1);                    
                end
            end
%             d(ii) = d(ii) - 0.5*(ind_N*log(2*pi) + ind_N*log(sigma(ii,1).^2));
        end
    end
    
%      (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);
    prior = - log(beta(20, 1.5)) + (20-1)*log((rho+1)/2) + (1.5-1)*log(1-(rho+1)/2); 
    prior = prior - log(sigma);
    d(r) = d(r) + prior(r);
end