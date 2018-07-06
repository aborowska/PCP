function [THR, cond] = threshold_t_garch11_varc_mle(y, mu_MLE, y_S, quantile)
            
    P = length(quantile); % potentially P different thresholds
%     quantile = tinv(threshold_m, mu_MLE(1));

%     y_S = var(y);
%y((T+1):(T+H))
    rho_MLE = (mu_MLE(1)-2)/mu_MLE(1);
    
    H = length(y);
    THR = zeros(P,H);
    cond = zeros(P,H);
    sum_C = 0;
    for jj = 1:H
        if (jj == 1)
            h_MLE = y_S;
        else
            temp = y(jj-1) - mu_MLE(2);
            h_MLE = mu_MLE(3)*(1-mu_MLE(4)-mu_MLE(5)) ...
                + mu_MLE(4)*temp^2 + mu_MLE(5)*h_MLE;
        end
        temp = (y(jj) - mu_MLE(2))/sqrt(rho_MLE*h_MLE);
        for p = 1:P
            if (temp < quantile(p)) % observations of interest
                cond(p,jj) = 1;
            else
                sum_C = sum_C + 1;
            end
            THR(p,jj) = mu_MLE(2) + sqrt(rho_MLE*h_MLE)*quantile(p);
        end
    end
    
    %if cond then dont censor
    % if cond==0 then censor
    % so censor when  (y(jj) - mu_MLE(2))/sqrt(h_MLE) >=  quantile
    % so censor when  y(jj) >= mu_MLE(2) + sqrt(h_MLE)*quantile

end