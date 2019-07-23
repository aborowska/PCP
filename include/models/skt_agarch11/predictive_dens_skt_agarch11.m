function d = predictive_dens_skt_agarch11(y, hT, theta)
% hT is the estimated volatility for the last in-sample period T
% y is yT (the last in-sample) and H out-of-sample
   %M = size(theta,1);
   
    H = length(y)-1;
    d = zeros(1,H);

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
    CONST = (logb + logc);
    
    tau = - a./b;

    h = hT;
    for t = 1:H
        h = omega.*(1-alpha-beta) + alpha.*(y(t)-mu2).^2  + beta.*h;
        x = (y(t+1) - mu)./sqrt(h);
        dd = sktlogpdf(x, nu, lambda, a, b, tau) - log(sqrt(h)) + CONST;
        dd = exp(dd);
        d(t) = mean(dd);         
    end
end

 

function lpdf = sktlogpdf(x, nu, lambda, a, b, tau)
% LOG  **KERNEL** OF THE SK-T DENSITY (no constant for speed)

    ind1 = (x < tau);
    ind2 = (x >= tau);

    lpdf = NaN(size(x));
    lpdf(ind1) = - ((nu(ind1)+1)./2).*log(1 + (((b(ind1).*x(ind1) + a(ind1))./(1-lambda(ind1))).^2)./(nu(ind1)-2));
    lpdf(ind2) = - ((nu(ind2)+1)./2).*log(1 + (((b(ind2).*x(ind2) + a(ind2))./(1+lambda(ind2))).^2)./(nu(ind2)-2)); 
end