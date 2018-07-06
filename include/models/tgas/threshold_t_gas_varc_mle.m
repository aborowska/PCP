function [THR, cond] = threshold_t_gas_varc_mle(y, mu_MLE, fT_MLE, quantile)
            
    P = length(quantile); % potentially P different thresholds
%     quantile = tinv(threshold_m, mu_MLE(1));

%     y_S = var(y);
%y((T+1):(T+H))
    rho_MLE = (mu_MLE(1)-2)/mu_MLE(1);
    nu_con_MLE = (mu_MLE(1)+1)./(mu_MLE(1)-2);
    A_MLE = mu_MLE(4).*((mu_MLE(1)+3)./mu_MLE(1));

    H = length(y);
    THR = zeros(P,H);
    cond = zeros(P,H);
    sum_C = 0;
    for jj = 1:H
        if (jj == 1)
            f_MLE = fT_MLE;
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