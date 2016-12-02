function d = posterior_iid(sigma, y)
    N = length(y);
    M = size(sigma,1);
    d = -Inf*ones(M,1);
    r = sigma>0; 
    d(r) = -0.5*(N*log(2*pi) + N*log(sigma(r).^2) + sum(y.^2)./sigma(r).^2);
    prior = -log(sigma);
    d(r) = d(r) + prior(r);
end