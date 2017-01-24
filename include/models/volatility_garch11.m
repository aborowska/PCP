function h = volatility_garch11(theta, y, y_S)
    T = length(y);
    M = size(theta,1);

    mu = theta(:,1);
    omega = theta(:,2);
    alpha = theta(:,3);
    beta = theta(:,4);

    if (nargin == 3)
        h = y_S*ones(M,1);
    else
        h = omega./(1-alpha-beta);
    end
            
    for jj = 2:T
      h = omega + alpha.*(y(jj-1,1)-mu).^2 + beta.*h;
    end  
          
        
end