function d = predictive_cdf_skt_agarch11(y, hT, theta, threshold)
% cdf_post = predictive_cdf_skt_agarch11(y(T:(T+H)), hT, draw, THR_emp);
% hT is the estimated volatility for the last in-sample period T
 
      
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
    mu2 = theta(:,5);
    alpha = theta(:,6);
    beta = theta(:,7);

    logc = gammaln((nu+1)/2) - gammaln(nu/2) - 0.5*log(pi*(nu-2));
    c = exp(logc);
    a = 4.*lambda.*c.*((nu-2)./(nu-1));
    logb = 0.5.*log(1 + 3.*lambda.^2 - a.^2);    
    b = exp(logb);
    
    tau = - a./b;
    
    h = hT;
    for t = 1:H
        h = omega.*(1-alpha-beta) + alpha.*(y(t)-mu2).^2  + beta.*h;
        for p = 1:P
            z = (threshold(p,t) - mu)./sqrt(h);
%             dd = tpdf(x,nu);
            dd = sktcdf(z, nu, lambda, a, b, tau);
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
