function d = posterior_iid(theta, y)
    N = length(y);  
    mu = theta(:,1);
    sigma = theta(:,2);
    M = size(theta,1);    
    d = -Inf*ones(M,1);
    r = sigma>0; 
    
    for ii = 1:M
        if r(ii)
            d(ii) = -0.5*(N*log(2*pi) + N*log(sigma(ii,1).^2) + sum((y-mu(ii,1)).^2)./sigma(ii,1).^2);        
        end
    end
    
    prior = -log(sigma);
    d(r) = d(r) + prior(r);
end