function d = C_posterior_iid(theta, y, threshold)
    N = length(y);
    mu = theta(:,1);
    sigma = theta(:,2);
    M = size(theta,1);
    d = -Inf*ones(M,1);
    r = sigma>0; % prior satisfied
    
    ind = (y<threshold); % observations of interest
    ind_N = sum(ind);
    % for the % "censored observations":
    P = log(1-normcdf(threshold,mu,sigma)); %P(y_t>=threshold|y_1,...,y_(t-1))
    
    for ii = 1:M
        if r(ii)
            d(ii) = -0.5*(ind_N*log(2*pi) + ind_N*log(sigma(ii,1).^2) + sum((y(ind)-mu(ii,1)).^2)./sigma(ii,1).^2);        
            d(ii) = d(ii) + (N-ind_N)*P(ii);
        end
    end
    
    % Misspecified model: N(0,sigma)
    prior = -log(sigma);
    d(r) = d(r) + prior(r);
end