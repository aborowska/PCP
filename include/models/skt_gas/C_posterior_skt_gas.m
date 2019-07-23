function d = C_posterior_skt_gas(theta, y, threshold, hyper, y_S)
    T = length(y);
    [N ,~] = size(theta);
    
    lambda = theta(:,1);    
    nu = theta(:,2);
    mu = theta(:,3);
    omega = theta(:,4);
    A = theta(:,5);
    B = theta(:,6); 
    
    prior = prior_skt_gas(N, lambda, omega, B, nu, hyper);
    pr = logical(prior(:,1));

    logc = NaN*ones(N,1);
    c = NaN*ones(N,1);
    a = NaN*ones(N,1);
    logb = NaN*ones(N,1);
    b = NaN*ones(N,1);
    tau = NaN*ones(N,1);
     
    logc(pr) = gammaln((nu(pr)+1)/2) - gammaln(nu(pr)/2) - 0.5*log(pi*(nu(pr)-2));
    c(pr) = exp(logc(pr));
    a(pr) = 4.*lambda(pr).*c(pr).*((nu(pr)-2)./(nu(pr)-1));
    logb(pr) = 0.5.*log(1 + 3.*lambda(pr).^2 - a(pr).^2);    
    b(pr) = exp(logb(pr));
    
    tau(pr) = - a(pr)./b(pr);
   
    CONST = (logb + logc);
    
    
%     f = omega./(1-B); % unconditional variance to initialize h_1
    f = log(y_S)*ones(N,1);
    
    d = -Inf*ones(N,1);
    
    sum_C = 0;
    for jj = 1:T
        if (y(jj) >= threshold) % observations of interest
            sum_C = sum_C + 1;
        end
    end
    %sum_C
    z1 = zeros(sum_C,1);
    z2 = zeros(T-sum_C,1);
    
    ind = (y<threshold); % observations of interest

    for ii = 1:N        
        if prior(ii,1) % compute only for valid draws
            h = exp(f(ii,1));
            scale = sqrt(h);
            z = (y(1,1) - mu(ii,1))/scale;  

            d(ii,1) = prior(ii,2);
            if ind(1,1) 
                tmp00 = sktpdf(z, nu(ii,1), lambda(ii,1), a(ii,1), ...
                    b(ii,1), tau(ii,1));
                d(ii,1) = d(ii,1) + tmp00 - log(scale) + CONST(ii);
            else
                z2 = (threshold - mu(ii,1))/scale;
                d(ii,1) = d(ii,1) + log(1-sktcdf(z2, nu(ii,1), lambda(ii,1), a(ii,1), b(ii,1), tau(ii,1)));
                
            end    
            ind_tau = 2*(z >= tau(ii,1)) - 1;
        
            for jj = 2:T
                nom = (nu(ii,1)+1)*b(ii,1)*z*(b(ii,1)*z+a(ii,1));
                den = (nu(ii,1)-2)*(1+ind_tau*lambda(ii,1))^2 + (b(ii,1)*z+a(ii,1))^2;
                s = 0.5*(nom/den - 1);
                f(ii,1) = omega(ii,1) + A(ii,1)*s ...
                            + B(ii,1)*f(ii,1);
                         
                h = exp(f(ii,1));               
                scale = sqrt(h);
                z = (y(jj,1)-mu(ii,1))/scale;         
                ind_tau = 2*(z >= tau(ii,1)) - 1; 
                
                if ind(jj,1)
                    tmp00 = sktpdf(z, nu(ii,1), lambda(ii,1), ...
                        a(ii,1), b(ii,1), tau(ii,1)); 
                    d(ii,1) = d(ii,1) + tmp00 - log(scale) + CONST(ii);   
                else  % P(y_t>=threshold|y_1,...,y_(t-1)) = 1 - P(y_t<threshold|y_1,...,y_(t-1))                    
                    z2 = (threshold-mu(ii,1))/scale;
                    d(ii,1) = d(ii,1) + log(1-sktcdf(z2, nu(ii,1), lambda(ii,1), a(ii,1), b(ii,1), tau(ii,1)));      %(x, nu, lambda, a, b, tau)                                  
                end
                 
            end 
        end
    end
    
%     d = d + prior(:,2);

end


function cdf = sktcdf(x, nu, lambda, a, b, tau)
% SK-T CDF 
     
    indicator1 = (x < tau);
    indicator2 = (x >= tau);
    
    y = NaN(size(x));
    cdf = NaN(size(x));
    
    temp = sqrt((nu)./(nu-2));
    y(indicator1) = (b.*x(indicator1) + a)./(1-lambda);
    y(indicator2) = (b.*x(indicator2) + a)./(1+lambda);
    y = temp.*y;
    
    tcdfy = tcdf(y,nu);
    cdf(indicator1) = (1-lambda).*tcdfy(indicator1);
    cdf(indicator2) = - lambda + (1+lambda).*tcdfy(indicator2);
 
end


function lpdf = sktpdf(x, nu, lambda, a, b, tau)
% LOG  **KERNEL** OF THE SK-T DENSITY (no constant for speed)

    indicator1 = (x < tau);
    indicator2 = (x >= tau);

    lpdf = NaN(size(x));
    lpdf(indicator1) = - ((nu+1)./2).*log(1 + (((b.*x(indicator1) + a)./(1-lambda)).^2)./(nu-2));
    lpdf(indicator2) = - ((nu+1)./2).*log(1 + (((b.*x(indicator2) + a)./(1+lambda)).^2)./(nu-2)); 
end


function R = prior_skt_gas(N, lambda, omega, B, nu, hyper)
    % uniform priors 
    
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val an the corresponding point
    
    c1 = (omega > 0);
    c3 = ((B >= 0) & (B < 1));
    c4 = (nu > 2);
    c5 = ((lambda > -1) & lambda < 1);
    
    r1 = (c1 & c3 & c4 & c5);
    r2 = -Inf*ones(N,1);
    r2(r1==true) = log(0.5) + log(hyper) - hyper*( nu(r1==true)  - 2); % exponential prior: nu~exp(1) --> p(nu)=exp(-nu) from 2 to inf  

    R = [r1, r2];
end
