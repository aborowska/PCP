function d = loglik_agarch11(theta, y)
    T = length(y);
    M = size(theta,1);

    mu = theta(:,1);
    gama = theta(:,2); % "typo" on purpose: not to confuse with the gamma function
    omega = theta(:,3);
    alpha = theta(:,4);
    beta = theta(:,5);
     
    ind = 1:T-1;
    
    d = -Inf*ones(M,1);
        
    for ii = 1:M
        if mod(ii,1000) == 0
            fprintf('posterior ii = %d\n',ii);
        end
        h = zeros(T,1); 
        h(1,1) = omega(ii,1) + (alpha(ii,1)*gama(ii,1)^2)/(1-alpha(ii,1)-beta(ii,1)); % unconditional variance
        h(2:T) = omega(ii,1)*(1-alpha(ii,1)-beta(ii,1)) + alpha(ii,1)*(y(ind,1)-mu(ii,1)+gama(ii,1)).^2;
        d(ii,1) = ((y(1,1)-mu(ii,1)).^2)/h(1,1);
        
        for jj = 2:T
            h(jj,1) = h(jj,1) + beta(ii,1)*h(jj-1,1);
            d(ii,1) = d(ii,1) +  ((y(jj,1)-mu(ii,1)).^2)/h(jj,1);
        end  
        
        d(ii,1) = -0.5*(d(ii,1) + T*log(2*pi) + T*sum(log(h)));
    end

    d = d/T; 
end