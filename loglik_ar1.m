function d = loglik_ar1(theta, y)
    N = length(y);
    M = size(theta,1);

    mu = theta(:,1);
    sigma = theta(:,2);
    rho = theta(:,3);
    
    a1 = mu./(1-rho);
    P1 = (sigma.^2)./(1-rho.^2);
    
    d = -Inf*ones(M,1);
    r = (sigma>0) & (abs(rho)<1); 
    
    for ii = 1:M
        if r(ii)
            d(ii) = ((y(1,1) - a1(ii,1)).^2)./P1(ii,1);
            for jj = 2:N
                d(ii) = d(ii) + ((y(jj,1) - mu(ii,1) - rho(ii,1).*y(jj-1,1)).^2)./sigma(ii,1).^2;                      
            end
            d(ii) = -0.5*(d(ii) + N*log(2*pi) + N*log(sigma(ii,1).^2));
        end
    end
    
    d = d/N; 
end