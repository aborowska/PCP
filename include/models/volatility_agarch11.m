function h_out = volatility_agarch11(theta, y, y_S, H)
    T = length(y);
    M = size(theta,1);

    mu = theta(:,1);
    gama = theta(:,2); % "typo" on purpose: not to confuse with the gamma function
    omega = theta(:,3);
    alpha = theta(:,4);
    beta = theta(:,5);

    if (y_S > 0)
        h = y_S*ones(M,1);
    else
        h = omega + (alpha.*gama.^2)./(1-alpha-beta); % unconditional variance
    end
    
    if (nargin == 3)
        H = 1;
    end
    
    h_out = zeros(M,H);    
    for jj = 2:T
        h = omega.*(1-alpha-beta) + alpha.*(y(jj-1,1)-mu+gama).^2 + beta.*h;
        if (jj > T-H)
            h_out(:,jj-(T-H)) = h;
        end
    end  
          
        
end