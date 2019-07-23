function h_out = volatility_skt_agarch11(theta, y, y_S, H)
    T = length(y);
    M = size(theta,1);

%     lambda = theta(:,1);
%     nu = theta(:,2);
    mu = theta(:,3);
    omega = theta(:,4);
    mu2 = theta(:,5);
    alpha = theta(:,6);
    beta = theta(:,7);

    if (y_S > 0)
        h = y_S*ones(M,1);
    else
        gama = mu - mu2;
        h = omega + (alpha.*gama.^2)./(1-alpha-beta); % unconditional variance
    end
    
    if (nargin == 3)
        H = 1;
    end
    
    if (H > 0)
        h_out = zeros(M,H);    
    end
    
      
    for jj = 2:T
        h = omega.*(1-alpha-beta) + alpha.*(y(jj-1,1)-mu2).^2 + beta.*h;
        if (jj > T-H)
            h_out(:,jj-(T-H)) = h;
        end
    end
    
    if (H == 0)
        h_out = h;
    end
end