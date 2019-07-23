function [THR, cond] = threshold_skt_agarch11_varc_mle(y, mu_MLE, yT, hT_MLE, QUANT_MLE)
            
    P = length(QUANT_MLE); % potentially P different thresholds quantile = P_bars;
%     quantile = tinv(threshold_m, mu_MLE(1));

    lambda = mu_MLE(:,1);
    nu = mu_MLE(:,2);
    mu = mu_MLE(:,3);
    omega = mu_MLE(:,4);
    mu2 = mu_MLE(:,5);
    alpha = mu_MLE(:,6);
    beta = mu_MLE(:,7);
    
    
    H = length(y);
    THR = zeros(P,H);
    cond = zeros(P,H);
    sum_C = 0;
    
    h_MLE = hT_MLE;
    for jj = 1:H
        if (jj == 1)
            temp = yT - mu2;           
        else
            temp = y(jj-1) - mu2;            
        end
        h_MLE = omega*(1-alpha-beta) + alpha*temp^2 + beta*h_MLE;
        
        temp = (y(jj) - mu)/sqrt(h_MLE);
        for p = 1:P
            if (temp < QUANT_MLE(p)) % observations of interest
                cond(p,jj) = 1;
            else
                sum_C = sum_C + 1;
            end
            THR(p,jj) = mu + sqrt(h_MLE)*QUANT_MLE(p);
        end
    end
    
    %if cond then dont censor
    % if cond==0 then censor
    % so censor when  (y(jj) - mu_MLE(2))/sqrt(h_MLE) >=  quantile
    % so censor when  y(jj) >= mu_MLE(2) + sqrt(h_MLE)*quantile

end