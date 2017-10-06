function d = C_posterior_t_garch11_2(theta, y, threshold, y_S, hyper)
    T = length(y);
    M = size(theta,1);
    
    nu = theta(:,1);
    mu = theta(:,2);
    omega = theta(:,3);
    alpha = theta(:,4);
    beta = theta(:,5);

    rho = (nu-2)./nu;

    if (y_S == 0) 
        h = omega;
    else        
        h = y_S*ones(M,1);
    end      
    
    d = -Inf*ones(M,1);
    prior = prior_t_garch(M, nu, omega, alpha, beta, hyper); 
    
    sum_C = 0;
    for jj = 1:T
        if (y(jj) >= threshold); % observations of interest
            sum_C = sum_C + 1;
        end
    end
    %sum_C
    z1 = zeros(sum_C,1);
    z2 = zeros(T-sum_C,1);
    
    ind = (y<threshold); % observations of interest

    for ii = 1:M        
        if prior(ii,1) % compute only for valid draws
            i2 = 0;
            scale = sqrt(rho(ii,1)*h(ii,1));
            d(ii,1) = 0;
            if ind(1,1) 
                z1  = (y(1,1) - mu(ii,1))/scale;
                d(ii,1) = d(ii,1) + log(tpdf(z1,nu(ii))/scale);
             else
                i2 = i2+1;
                z2(i2,1) = (threshold - mu(ii,1))/scale;
             end    
        
            for jj = 2:T
                h(ii,1) = omega(ii,1)*(1-alpha(ii,1)-beta(ii,1)) + alpha(ii,1)*(y(jj-1,1)-mu(ii,1)).^2  + beta(ii,1)*h(ii,1);
                scale = sqrt(rho(ii,1)*h(ii,1));
                if ind(jj,1)
                    z1 = (y(jj,1)-mu(ii,1))/scale;
                    d(ii,1) = d(ii,1)  + log(tpdf(z1,nu(ii))/scale);
                 else  % P(y_t>=threshold|y_1,...,y_(t-1)) = 1 - P(y_t<threshold|y_1,...,y_(t-1))
                    i2 = i2+1;
                    z2(i2,1) = (threshold-mu(ii,1))/scale;
                 end
            end
%             fprintf('d(ii,1) = %6.4f\n', d(ii,1));
%             fprintf('sum_cdf = %6.4f\n', sum(tcdf(z2,nu(ii,1))));
            d(ii,1) = d(ii,1) + sum(log(1-tcdf(z2,nu(ii,1))));
        end
    end
    d = d + prior(:,2);
end

function R = prior_t_garch(M, nu, omega, alpha, beta, hyper)
    % uniform prior on alpha and beta on (0,1)
    % with restriction alpha + beta < 1
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val at the corresponding point

    c1 = ((alpha >= 0) & (alpha < 1) & (beta >= 0) & (beta < 1));
    c2 = (alpha + beta < 1);
    c3 = (omega > 0);
    c4 = (nu > 2);
    
    r1 = (c1 & c2 & c3 & c4);

    r2 = -Inf*ones(M,1);
    r2(r1==true) = log(0.5) + log(hyper) - hyper*(nu(r1==true) - 2);
    
    R = [r1, r2];
end
