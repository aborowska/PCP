function [THR, cond] = threshold_t_gas_varc_mle(y, mu_MLE, yT, fT_MLE, quantile)
% QUANT = tinv(P_bars, mu_MLE(1))';
% fT_MLE = volatility_t_gas(mu_MLE,y(1:T),0);
% [THR_MLE, cond_MLE] = threshold_t_gas_varc_mle(y((T+1):(T+H)), mu_MLE, y(T), fT_MLE, QUANT);

    P = length(quantile); % potentially P different thresholds
%     quantile = tinv(threshold_m, mu_MLE(1));

  
    rho_MLE = (mu_MLE(1)-2)/mu_MLE(1);
    nu_con_MLE = (mu_MLE(1)+1)./(mu_MLE(1)-2);
    A_MLE = mu_MLE(4).*((mu_MLE(1)+3)./mu_MLE(1));

    H = length(y);
    THR = zeros(P,H);
    cond = zeros(P,H);
    sum_C = 0;
    
    f_MLE = fT_MLE;
    for jj = 1:H
        if (jj == 1)
            C = 1 + ((yT-mu_MLE(2)).^2)/((mu_MLE(1)-2)*f_MLE);
            f_MLE = mu_MLE(3) + A_MLE*(nu_con_MLE*((yT-mu_MLE(2)).^2)/C - f_MLE) ...
                        + mu_MLE(5)*f_MLE;
        else
            C = 1 + ((y(jj-1,1)-mu_MLE(2)).^2)/((mu_MLE(1)-2)*f_MLE);
            f_MLE = mu_MLE(3) + A_MLE*(nu_con_MLE*((y(jj-1,1)-mu_MLE(2)).^2)/C - f_MLE) ...
                        + mu_MLE(5)*f_MLE;
        end

        temp = (y(jj) - mu_MLE(2))/sqrt(rho_MLE*f_MLE);
        for p = 1:P
            if (temp < quantile(p)) % observations of interest
                cond(p,jj) = 1;
            else
                sum_C = sum_C + 1;
            end
            THR(p,jj) = mu_MLE(2) + sqrt(rho_MLE*f_MLE)*quantile(p);
        end
    end
    
    %if cond then dont censor
    % if cond==0 then censor
    % so censor when  (y(jj) - mu_MLE(2))/sqrt(f_MLE) >=  quantile
    % so censor when  y(jj) >= mu_MLE(2) + sqrt(f_MLE)*quantile

end