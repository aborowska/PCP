function h_out = volatility_garch11(theta, y, y_S, H)
    T = length(y);
    M = size(theta,1);

    mu = theta(:,1);
    omega = theta(:,2);
    alpha = theta(:,3);
    beta = theta(:,4);

    if (y_S > 0)
        h = y_S*ones(M,1);
    else
        h = omega;
    end
    
    if (nargin == 3)
        H = 1;
    end
    
    h_out = zeros(M,H);    
    for jj = 2:T
        h = omega.*(1-alpha-beta) + alpha.*(y(jj-1,1)-mu).^2 + beta.*h;
        if (jj > T-H)
            h_out(:,jj-(T-H)) = h;
        end
    end  
end