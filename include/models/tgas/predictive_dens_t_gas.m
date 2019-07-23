function d = predictive_dens_t_gas(y, fT, theta)
% hT is the estimated volatility for the last in-sample period T
% y is yT (the last in-sample) and H out-of-sample
    %M = size(theta,1);

    H = length(y)-1;
    d = zeros(1,H);

    nu = theta(:,1);
    mu = theta(:,2);
    omega = theta(:,3);
    A = theta(:,4);
    B = theta(:,5);

    rho = (nu-2)./nu;
    nu_con = (nu+1)./(nu-2);
    A = A.*((nu+3)./nu);

    f = fT;
    for t = 1:H
        C = 1 + ((y(t)-mu).^2)./((nu-2).*f);              
        f = omega + A.*(nu_con.*((y(t)-mu).^2)./C - f) + B.*f;                        
        % y(jj,1) ~ mu(:,1) + sqrt(rho(:,1).*f(:,jj)).*trnd(nu);

        % h = omega.*(1-alpha-beta) + alpha.*(y(t)-mu).^2  + beta.*h;     
        x = (y(t+1) - mu)./sqrt(rho.*f);
        dd = tpdf(x,nu) - log(sqrt(rho.*f)); 
        d(t) = mean(dd);         
    end
end