function d = posterior_agarch11(theta, y, y_S)
    T = length(y);
    M = size(theta,1);

    mu = theta(:,1);
    gama = theta(:,2); % "typo" on purpose: not to confuse with the gamma function
    omega = theta(:,3);
    alpha = theta(:,4);
    beta = theta(:,5);
     
    prior = prior_garch(M, omega, alpha, beta);    
     
    ind = 1:T-1;
    
    d = -Inf*ones(M,1);
    
    for ii = 1:M       
        if prior(ii,1)
            h = zeros(T,1);
            if nargin == 3
                h(1,1) = y_S;
            else
                h(1,1) = omega(ii,1) + (alpha(ii,1)*gama(ii,1)^2)/(1-alpha(ii,1)-beta(ii,1)); % unconditional variance
            end
%             d(ii,1) = ((y(1,1)-mu(ii,1)).^2)/h(1,1);
            h(2:T) = omega(ii,1)*(1-alpha(ii,1)-beta(ii,1)) + alpha(ii,1)*(y(ind,1)-mu(ii,1)+gama(ii,1)).^2;
            for jj = 2:T
%                 h(jj,1) = omega(ii,1) + alpha(ii,1)*(y(jj-1,1)-mu(ii,1)+gama(ii,1)).^2 ...
%                     + beta(ii,1)*h(jj-1,1);
                h(jj,1) = h(jj,1) + beta(ii,1)*h(jj-1,1);
%                 d(ii,1) = d(ii,1) +  ((y(jj,1)-mu(ii,1)).^2)/h(jj,1);
            end  
            d(ii,1) = sum(((y-mu(ii,1)).^2)./h);
            d(ii,1) = -0.5*(d(ii,1) + T*log(2*pi) + sum(log(h)));            
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
    r2(r1==true) = log(0.5);
    
    R = [r1, r2];
end
