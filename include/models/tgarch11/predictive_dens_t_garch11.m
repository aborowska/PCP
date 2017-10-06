function d = predictive_dens_t_garch11(y, hT, theta)
% hT is the estimated volatility for the last in-sample period T
% y is yT (the last in-sample) and H out-of-sample
   %M = size(theta,1);
   
   H = length(y)-1;
   d = zeros(1,H);
   
   nu = theta(:,1);
   mu = theta(:,2);
   omega = theta(:,3);
   alpha = theta(:,4);
   beta = theta(:,5);
    
   rho = (nu-2)./nu;
   
   h = hT;
   for t = 1:H
        h = omega.*(1-alpha-beta) + alpha.*(y(t)-mu).^2  + beta.*h;
        x = (y(t+1) - mu)./sqrt(rho.*h);
        dd = tpdf(x,nu); 
        d(t) = mean(dd);         
   end
end