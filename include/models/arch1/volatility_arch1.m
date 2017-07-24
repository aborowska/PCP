function [h_out,h_all] = volatility_arch1(theta, y, y_S, H)
    T = length(y);
    M = size(theta,1);

    mu = theta(:,1);
    omega = theta(:,2); 
    mu2 = theta(:,3);
    alpha = theta(:,4);
    
    if (y_S > 0)
        h = y_S*ones(M,1);
    else
        gama = mu - mu2;
        h = omega + (alpha.*gama.^2)./(1-alpha); % unconditional variance
    end
    
    if (nargin == 3)
        H = 1;
    end
    
    h_out = zeros(M,H);  
    if (nargout == 2)
        h_all = zeros(M,T+H);
        h_all(:,1) = h;
    end
    
    for jj = 2:T
        h = omega.*(1-alpha) + alpha.*(y(jj-1,1)-mu2).^2;
        if (jj > T-H)
            h_out(:,jj-(T-H)) = h;
        end
        
        if (nargout == 2)
            h_all(:,jj) = h;
        end
    end                 
end