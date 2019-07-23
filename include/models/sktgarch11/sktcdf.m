function cdf = sktcdf(x, nu, lambda, a, b, tau)
% SK-T CDF
 
    if (nargin < 6)
        logc = gammaln((nu+1)/2) - gammaln(nu/2) - 0.5*log(pi*(nu-2));
        c = exp(logc);
        a = 4*lambda*c*((nu-2)/(nu-1));
    %     loga = log(a);
    %     loga = log(4) + log(lambda) + logc + log(nu-2) - log(nu-1);
    %     a = exp(loga);
        logb = 0.5*log(1 + 3*lambda^2 - a^2);    
        b = exp(logb);

        tau = - a./b;
    end
    
%         indicator1 = ((data(t)-mu(t))./sqrt(h(t))<-a./b);
%         indicator2 = ((data(t)-mu(t))./sqrt(h(t))>=-a./b);
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
