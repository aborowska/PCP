function d = predictive_cdf_skt_gas(y, fT, theta, threshold)
% fT is the estimated volatility for the last in-sample period T
% y is yT (the last in-sample) and H out-of-sample
%M = size(theta,1);
    
    
    H = length(y)-1;
    [P, varc] = size(threshold); % potentially P different thresholds
    if (varc == 1)
        threshold = repmat(threshold,1,H);
    end
    d = zeros(P,H);
   
    lambda = theta(:,1);
    nu = theta(:,2);
    mu = theta(:,3);
    omega = theta(:,4);
    A = theta(:,5);
    B = theta(:,6); 

    logc = gammaln((nu+1)/2) - gammaln(nu/2) - 0.5*log(pi*(nu-2));
    c = exp(logc);
    a = 4.*lambda.*c.*((nu-2)./(nu-1));
    logb = 0.5.*log(1 + 3.*lambda.^2 - a.^2);    
    b = exp(logb);
    tau = - a./b;
    
    f = fT;
    z = (y(1)-mu)./sqrt(exp(f));
    for t = 1:H
        ind_tau = 2*(z >= tau) - 1;            
        nom = (nu+1).*b.*z.*(b.*z+a);
        den = (nu-2).*(1+ind_tau.*lambda).^2 + (b.*z+a).^2;
        s = 0.5*(nom./den - 1);
        f = omega + A.*s + B.*f;  %f_{t+1}  
        scale = sqrt(exp(f));
        z = (y(t+1)-mu)./scale;

        for p = 1:P
            x = (threshold(p,t) - mu)./scale;
            dd = sktcdf(x, nu, lambda, a, b, tau); 
            d(p,t) = mean(dd);         
        end
    end
end




function cdf = sktcdf(x, nu, lambda, a, b, tau)
% SK-T CDF 
     
    ind1 = (x < tau);
    ind2 = (x >= tau);
    
    y = NaN(size(x));
    cdf = NaN(size(x));
    
    temp = sqrt((nu)./(nu-2));
    y(ind1) = (b(ind1) .*x(ind1) + a(ind1) )./(1-lambda(ind1));
    y(ind2) = (b(ind2).*x(ind2) + a(ind2))./(1+lambda(ind2));
    y = temp.*y;
    
    tcdfy = tcdf(y,nu);
    cdf(ind1) = (1-lambda(ind1) ).*tcdfy(ind1);
    cdf(ind2) = - lambda(ind2) + (1+lambda(ind2)).*tcdfy(ind2);
 
end
