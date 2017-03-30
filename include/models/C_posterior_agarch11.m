function d = C_posterior_agarch11(theta, y, threshold, y_S)
    T = length(y);
    M = size(theta,1);
    
    mu = theta(:,1);
%     gama = theta(:,2); % "typo" on purpose: not to confuse with the gamma function
%     omega = theta(:,3);
    omega = theta(:,2); % "typo" on purpose: not to confuse with the gamma function
    mu2 = theta(:,3);
    alpha = theta(:,4);
    beta = theta(:,5);
     
    if (nargin == 4)
        h1 = y_S*ones(M,1);
    else
       gama = mu - mu2;
       h1 = omega + (alpha.*gama.^2)./(1-alpha-beta); % unconditional variance
    end    
    
    d = -Inf*ones(M,1);
    prior = prior_garch(M, omega, alpha, beta);   
    ind = (y<threshold); % observations of interest

    for ii = 1:M        
        if prior(ii,1) % compute only for valid draws
            if ind(1,1) 
                d(ii,1) = - 0.5*(log(2*pi) + log(h1(ii,1)) + ((y(1,1) - mu(ii,1)).^2)./h1(ii,1));
            else
                d(ii,1) = log(1-normcdf_my_mex(threshold,mu(ii,1),sqrt(h1(ii,1))));
            end    
        
            h = zeros(T,1);
            h(1,1) = h1(ii,1);
            for jj = 2:T
%                 h(jj,1) = omega(ii,1)*(1-alpha(ii,1)-beta(ii,1)) + alpha(ii,1)*(y(jj-1,1)-mu(ii,1)+gama(ii,1)).^2  + beta(ii,1)*h(jj-1,1);
                h(jj,1) = omega(ii,1)*(1-alpha(ii,1)-beta(ii,1)) + alpha(ii,1)*(y(jj-1,1)-mu2(ii,1)).^2  + beta(ii,1)*h(jj-1,1);
                if ind(jj,1)
                    d(ii,1) = d(ii,1) - 0.5*(log(2*pi) + log(h(jj,1)) + ((y(jj,1)-mu(ii,1)).^2)/h(jj,1));
                else  % P(y_t>=threshold|y_1,...,y_(t-1)) = 1 - P(y_t<threshold|y_1,...,y_(t-1))
                    d(ii,1) = d(ii,1) + log(1-normcdf_my_mex(threshold,mu(ii,1), sqrt(h(jj,1))));                    
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
    r2(r1==true) = log(0.5);
    
    R = [r1, r2];
end
