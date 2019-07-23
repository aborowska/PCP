function [THR, cond] = threshold_skt_gas_varc_mle(y, theta_MLE, ...
    yT, fT_MLE, quantile)
% QUANT = tinv(P_bars, mu_MLE(1))';
% fT_MLE = volatility_skt_gas(mu_MLE,y(1:T),0);
% [THR_MLE, cond_MLE] = threshold_skt_gas_varc_mle(y((T+1):(T+H)), mu_MLE, y(T), fT_MLE, QUANT);

    P = length(quantile); % potentially P different thresholds
%     quantile = tinv(threshold_m, mu_MLE(1));

    lambda = theta_MLE(:,1);
    nu = theta_MLE(:,2);
    mu = theta_MLE(:,3);
    omega = theta_MLE(:,4);
    A = theta_MLE(:,5);
    B = theta_MLE(:,6); 

    logc = gammaln((nu+1)/2) - gammaln(nu/2) - 0.5*log(pi*(nu-2));
    c = exp(logc);
    a = 4.*lambda.*c.*((nu-2)./(nu-1));
    logb = 0.5.*log(1 + 3.*lambda.^2 - a.^2);    
    b = exp(logb);
    tau = - a./b;


    H = length(y);
    THR = zeros(P,H);
    cond = zeros(P,H);

    f_MLE = fT_MLE;
    z = (yT-mu)./sqrt(exp(f_MLE));

    for jj = 1:H
        ind_tau = 2*(z >= tau) - 1;
        nom = (nu+1).*b.*z.*(b.*z+a);
        den = (nu-2).*(1+ind_tau.*lambda).^2 + (b.*z+a).^2;
        s = 0.5*(nom./den - 1);
        f_MLE = omega + A.*s + B.*f_MLE;    
        
        scale_MLE = sqrt(exp(f_MLE));
        z = (y(jj) - mu)/scale_MLE;
       
        for p = 1:P
            if (z < quantile(p)) % observations of interest
                cond(p,jj) = 1;
            end
            THR(p,jj) = mu + scale_MLE*quantile(p);
        end
    end
    
    %if cond then dont censor
    % if cond==0 then censor
    % so censor when  (y(jj) - mu_MLE(2))/sqrt(f_MLE) >=  quantile
    % so censor when  y(jj) >= mu_MLE(2) + sqrt(f_MLE)*quantile

end