function d = C_posterior_iid(sigma, y, threshold)
    N = length(y);
    M = size(sigma,1);
    d = -Inf*ones(M,1);
    r = sigma>0; 
    
    
    ind = (y<threshold);
    ind_N = sum(ind);
    P = log(1-normcdf(threshold,0,sigma)); %P(y_t>=threshold|y_1,...,y_(t-1))

%     y_tilde = y;
%     y_tilde(ind) = threshold;
    
    d(r) = -0.5*(ind_N*log(2*pi) + ind_N*log(sigma(r).^2) + sum(y(ind).^2)./sigma(r).^2);
    d(r) = d(r) + (N-ind_N)*P(r);
    
    % Misspecified model: N(0,sigma)
    prior = -log(sigma);
    d(r) = d(r) + prior(r);
end