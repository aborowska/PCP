function [d, T, D, THR] = C_posterior_ar1_varc_mle(theta, y, mu_mle, threshold, fun)
    T = length(y);
    M = size(theta,1);
    
    if (nargin == 4)
        fun = @(xx) xx; % can be abs or ^2  
    end
    
    if (M < 20)
        D = zeros(M,T);
    else
        D = [];
    end
    
    mu = theta(:,1);
    sigma = theta(:,2);
    rho = theta(:,3);
    
    cond = zeros(T,1);
    THR = zeros(T,1);

    a1_mle = mu_mle(1)./(1-mu_mle(3));
    P1_mle = (mu_mle(2).^2)./(1-mu_mle(3).^2);
    for jj=1:1:T
%         // standardized variables
        if (jj==1)
            THR(jj) = (y(jj) - a1_mle)/sqrt(P1_mle);            
        else
            THR(jj) = (y(jj) - mu_mle(1) - mu_mle(3)*y(jj-1))/mu_mle(2);        
        end
        
        if (THR(jj) < threshold)        
            cond(jj) = 1;
        end
        
        if (jj==1)
            THR(jj) = a1_mle + sqrt(P1_mle)*threshold;            
        else        
            THR(jj) = mu_mle(1) + mu_mle(3)*y(jj-1) + mu_mle(2)*threshold;
        end
    end
    
%   Stationary distibution for the first observation
    a1 = mu./(1-rho);
    P1 = (sigma.^2)./(1-rho.^2);
    
    d = -Inf*ones(M,1);
    r = (sigma>0) & (abs(rho)<1); % prior satisfied

    for ii = 1:M
        if r(ii) % compute only for valid draws
%             THR = fun(threshold*a1(ii,1));
            if cond(1)
                d(ii) = - 0.5*(log(2*pi) + log(P1(ii,1)) + ((y(1,1) - a1(ii,1)).^2)./P1(ii,1));
            else
                d(ii) = log(1-normcdf_my(THR(1),a1(ii,1),sqrt(P1(ii,1))));
            end
            
            if (M < 20)
                D(ii,1) = d(ii);
            end           
            
            for jj = 2:T
                MU = mu(ii,1) + rho(ii,1).*y(jj-1,1);
%                 THR = fun(threshold*MU);
                if cond(jj)
                    CONTR = - 0.5*(log(2*pi) + log(sigma(ii,1).^2) + ((y(jj,1) - MU).^2)./(sigma(ii,1).^2));
                    d(ii) = d(ii) + CONTR;    
                else  % P(y_t>=threshold|y_1,...,y_(t-1)) = 1 - P(y_t<threshold|y_1,...,y_(t-1))
%                     d(ii) = d(ii) + log(1-normcdf_my(THR,MU,sigma(ii,1)));
                    CONTR = log(1-normcdf_my(THR(jj),MU,sigma(ii,1)));
                    d(ii) = d(ii) + CONTR;
                end
                
                if (M < 20)
%                     D(ii,jj) = d(ii) - D(ii,jj-1);
                    D(ii,jj) = CONTR;
                end                
            end
        end
    end
    
% %      (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);
%     prior = - log(beta(20, 1.5)) + (20-1)*log((rho+1)/2) + (1.5-1)*log(1-(rho+1)/2); 
%     prior = prior - log(sigma);
%     d(r) = d(r) + prior(r);
end