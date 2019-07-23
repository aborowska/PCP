function d = C_posterior_skt_agarch11(theta, y, threshold, y_S, hyper)
    T = length(y);
    M = size(theta,1);
    
    lambda = theta(:,1);
    nu = theta(:,2);
    mu = theta(:,3);
    omega = theta(:,4);
    mu2 = theta(:,5);
    alpha = theta(:,6);
    beta = theta(:,7);

    prior = prior_skt_garch(M, lambda, nu, omega, alpha, beta, hyper);   
    pr = logical(prior(:,1));

%     c = gamma((nu+1)/2)/(sqrt(pi*(nu-2))*gamma(nu/2));
%     b = sqrt(1 + 3*lamda^2 - a^2);
    logc = NaN*ones(M,1);
    c = NaN*ones(M,1);
    a = NaN*ones(M,1);
    logb = NaN*ones(M,1);
    b = NaN*ones(M,1);
    tau = NaN*ones(M,1);
     
    logc(pr) = gammaln((nu(pr)+1)/2) - gammaln(nu(pr)/2) - 0.5*log(pi*(nu(pr)-2));
    c(pr) = exp(logc(pr));
    a(pr) = 4.*lambda(pr).*c(pr).*((nu(pr)-2)./(nu(pr)-1));
%     loga = log(a);
%     loga = log(4) + log(lambda) + logc + log(nu-2) - log(nu-1);
%     a = exp(loga);
    logb(pr) = 0.5.*log(1 + 3.*lambda(pr).^2 - a(pr).^2);    
    b(pr) = exp(logb(pr));
    
    tau(pr) = - a(pr)./b(pr);
    
    CONST = (logb + logc);
    
    
    if (y_S > 0)
        h = y_S*ones(M,1);
    else
        gama = mu - mu2;
        h = (omega + (alpha.*gama.^2))./(1-alpha-beta); % unconditional variance
    end    
    
    d = -Inf*ones(M,1);
 
    ind = (y<threshold); % observations of interest

    for ii = 1:M        
        if prior(ii,1) % compute only for valid draws
            scale = sqrt(h(ii,1));               

            d(ii,1) = prior(ii,2); 
            
            if ind(1,1) 
                z = (y(1,1) - mu(ii,1))/scale;        
                d(ii,1) = d(ii,1) + sktlogpdf(z, nu(ii,1), lambda(ii,1), a(ii,1), b(ii,1), tau(ii,1)) ...
                    - log(scale) + CONST(ii);
            else
                z = (threshold - mu(ii,1))/scale;
                d(ii,1) = d(ii,1) + log(1-sktcdf(z, nu(ii,1), lambda(ii,1), a(ii,1), b(ii,1), tau(ii,1)));
            end    
        
 
            for jj = 2:T
%                 h(ii,1) = omega(ii,1)*(1-alpha(ii,1)-beta(ii,1)) + alpha(ii,1)*(y(jj-1,1)-mu2(ii,1)).^2  + beta(ii,1)*h(ii,1);   
                h(ii,1) = omega(ii,1)*(1-alpha(ii,1)-beta(ii,1)) + ...
                    alpha(ii,1)*(y(jj-1,1)-mu2(ii,1)).^2 + beta(ii,1)*h(ii,1);
                
                scale = sqrt(h(ii,1));            
                if ind(jj,1)
                    z = (y(jj,1) - mu(ii,1))/scale;         
                    d(ii,1) = d(ii,1) + sktlogpdf(z, nu(ii,1), lambda(ii,1), a(ii,1), b(ii,1), tau(ii,1)) ...
                        - log(scale) + CONST(ii);
                else  % P(y_t>=threshold|y_1,...,y_(t-1)) = 1 - P(y_t<threshold|y_1,...,y_(t-1))
                    z = (threshold - mu(ii,1))/scale;
                    d(ii,1) = d(ii,1) + log(1-sktcdf(z, nu(ii,1), lambda(ii,1), a(ii,1), b(ii,1), tau(ii,1)));      %(x, nu, lambda, a, b, tau)              
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

function lpdf = sktlogpdf(x, nu, lambda, a, b, tau)
% LOG  **KERNEL** OF THE SK-T DENSITY (no constant for speed)

    indicator1 = (x < tau);
    indicator2 = (x >= tau);

    lpdf = NaN(size(x));
    lpdf(indicator1) = - ((nu+1)./2).*log(1 + (((b.*x(indicator1) + a)./(1-lambda)).^2)./(nu-2));
    lpdf(indicator2) = - ((nu+1)./2).*log(1 + (((b.*x(indicator2) + a)./(1+lambda)).^2)./(nu-2)); 
end

function R = prior_skt_garch(M, lambda, nu, omega, alpha, beta, hyper)
    % uniform prior on alpha and beta on (0,1)
    % with restriction alpha + beta < 1
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val at the corresponding point

    c1 = ((alpha >= 1e-6) & (alpha < 1) & (beta >= 1e-6) & (beta < 1));
    c2 = (alpha + beta < 1);
    c3 = (omega > 1e-6);
    c4 = (nu > 2);
    c5 = ((lambda > -1) & (lambda < 1));
    
    r1 = (c1 & c2 & c3 & c4 & c5);

    r2 = -Inf*ones(M,1);
    r2(r1==true) = log(0.5) + log(hyper) - hyper*(nu(r1==true) - 2);
    
    R = [r1, r2]; %prior =  [r1, r2];
end
