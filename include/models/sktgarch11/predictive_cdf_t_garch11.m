function d = predictive_cdf_t_garch11(y, hT, theta, threshold)
% hT is the estimated volatility for the last in-sample period T
% y is yT (the last in-sample) and H out-of-sample
   %M = size(theta,1);
      
   H = length(y)-1;
   [P, varc] = size(threshold); % potentially P different thresholds
   if (varc == 1)
        threshold = repmat(threshold,1,H);
   end
   d = zeros(P,H);
   
   nu = theta(:,1);
   mu = theta(:,2);
   omega = theta(:,3);
   alpha = theta(:,4);
   beta = theta(:,5);
    
   rho = (nu-2)./nu;
   
   h = hT;
   for t = 1:H
        h = omega.*(1-alpha-beta) + alpha.*(y(t)-mu).^2  + beta.*h;
        for p = 1:P
            x = (threshold(p,t) - mu)./sqrt(rho.*h);
%             dd = tpdf(x,nu); 
            dd = tcdf(x,nu); 
            d(p,t) = mean(dd);         
        end
   end
end