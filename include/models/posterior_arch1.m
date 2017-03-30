function d = posterior_arch1(theta, y, y_S)
    [M, k] = size(theta);
    T = size(y,1);
    
    mu = theta(:,1);
    omega = theta(:,2);
    if (k == 4) % different mu and mu2
        mu2 = theta(:,3);
        alpha = theta(:,4);
    else
        mu2 = mu;
        alpha = theta(:,3);
    end

    prior = double((omega > 0) & (alpha > 0) & (alpha < 1));
    
    if ((nargin == 3) && (y_S > 0))
        h = y_S.*ones(M,1);
    else
       gama = mu - mu2;
       h = omega + (alpha.*gama.^2)./(1-alpha);       
    end
    
    % stationary distribution for the first observation
    d = log(h) + ((y(1,1) - mu).^2)./h;
    % recursion - to prevent memory problems for long series
    for ii = 2:T
        h = omega.*(1-alpha) + alpha.*(y(ii-1,1) - mu2).^2;
        d = d + log(h) + ((y(ii,1) - mu).^2)./h;
    end
    
    d = d + T*log(2*pi); % add Gasian constants
    d = -0.5*d + log(prior);   
    
end