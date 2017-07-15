function d = C_posterior_garch11_varc(theta, y, y_S)
    T = length(y);
    M = size(theta,1);
    
    mu = theta(:,1);
    omega = theta(:,2);
    alpha = theta(:,3);
    beta = theta(:,4);
    
    if (nargin == 4)
        h1 = y_S*ones(M,1);
    else
        h1 = omega;
    end    
    
    d = -Inf*ones(M,1);
    prior = prior_garch(M, omega, alpha, beta);   

    for ii = 1:M        
        if prior(ii,1) % compute only for valid draws
            if (y(1,1)>0)  %?????
                d(ii,1) = - 0.5*(log(2*pi) + log(h1(ii,1)) + ((y(1,1) - mu(ii,1)).^2)./h1(ii,1));
            else
                d(ii,1) = log(1-normcdf_my_mex(0,mu(ii,1),sqrt(h1(ii,1))));
            end    
        
            h = zeros(T+1,1);
            h(1,1) = h1(ii,1);
            jj = 2;
            h(jj,1) = omega(ii,1)*(1-alpha(ii,1)-beta(ii,1)) + alpha(ii,1)*(y(jj-1,1)-mu(ii,1)).^2  + beta(ii,1)*h(jj-1,1);
            
            for jj = 2:T
                % future volatility
                h(jj+1,1) = omega(ii,1)*(1-alpha(ii,1)-beta(ii,1)) + alpha(ii,1)*(y(jj,1)-mu(ii,1)).^2  + beta(ii,1)*h(jj,1);
                
                if (h(jj+1,1) <= h(jj,1))
                    d(ii,1) = d(ii,1) - 0.5*(log(2*pi) + log(h(jj,1)) + ((y(jj,1)-mu(ii,1)).^2)/h(jj,1));
                else  % P(y_t>=threshold|y_1,...,y_(t-1)) = 1 - P(y_t<threshold|y_1,...,y_(t-1))
%                     d(ii,1) = d(ii,1) + log(1-normcdf_my_mex(threshold,mu(ii,1), sqrt(h(jj,1))));  
                    CONTR = (((y(jj-1,1)-mu(ii,1)).^2) -  (beta(ii,1)./alpha(ii,1))*(h(jj,1)-h(jj-1,1)))/h(jj,1);
                    d(ii,1) = d(ii,1) + log(1-chi2cdf(CONTR,1));                    
                end
            end  
        end
    end
    d = d + prior(:,2);
end

function R = prior_garch(M, omega, alpha, beta)
    % uniform prior on alpha and beta on (0,1)
    % with restriction alpha + beta < 1
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val at the corresponding point

    c1 = ((alpha >= 0) & (alpha < 1) & (beta >= 0) & (beta < 1));
    c2 = (alpha + beta < 1);
    c3 = (omega > 0);

    r1 = (c1 & c2 & c3);

    r2 = -Inf*ones(M,1);
    r2(r1==true) = log(2);
    
    R = [r1, r2];
end
