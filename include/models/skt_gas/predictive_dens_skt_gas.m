function d = predictive_dens_skt_gas(y, fT, theta)
% fT is the estimated log volatility for the last in-sample period T
% y is yT (the last in-sample) and H out-of-sample
    %M = size(theta,1);
% dens_post = predictive_dens_skt_gas(y(T:(T+H)), fT, draw);

    H = length(y)-1;
    d = zeros(1,H);

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
    CONST = (logb + logc);
    tau = - a./b;

    f = fT;  
    z = (y(1,1)-mu)./sqrt(exp(f));
    for t = 1:H    
        ind_tau = 2*(z >= tau) - 1;            

        nom = (nu+1).*b.*z.*(b.*z+a);
        den = (nu-2).*(1+ind_tau.*lambda).^2 + (b.*z+a).^2;
        s = 0.5*(nom./den - 1);
        f = omega + A.*s + B.*f;    
        
        % h = omega.*(1-alpha-beta) + alpha.*(y(t)-mu).^2  + beta.*h;  
        scale = sqrt(exp(f));
        z = (y(t+1,1)-mu)./scale;
        dd = sktlogpdf(z, nu, lambda, a, b, tau) - log(scale) + CONST; 
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