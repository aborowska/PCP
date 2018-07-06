function d = C_posterior_t_gas_varc_mle(theta, y,  mu_mle, threshold, hyper)
    T = length(y);
    M = size(theta,1);
    
    nu = theta(:,1);
    mu = theta(:,2);
    omega = theta(:,3);
    A = theta(:,4);
    B = theta(:,5);

    rho = (nu-2)./nu;
    nu_con = (nu+1)./(nu-2);
    A = A.*((nu+3)./nu);

    rho_mle = (mu_mle(1)-2)./mu_mle(1);
    nu_con_mle = (mu_mle(1)+1)./(mu_mle(1)-2);
    A_mle = mu_mle(4).*((mu_mle(1)+3)./mu_mle(1));
    
    prior = prior_t_gas(M, omega, B, nu, hyper);
    
    f = omega./(1-B); % unconditional variance to initialize h_1

    d = -Inf*ones(M,1);
    
    sum_C = 0;
    cond = zeros(T,1);
    THR = zeros(T,1);
    for jj = 1:T
        if (jj == 1)
            f_mle = mu_mle(3)/(1-mu_mle(5));
        else
            C = 1 + ((y(jj-1,1)-mu_mle(2)).^2)/((mu_mle(1)-2)*f_mle);

            f_mle = mu_mle(3) + A_mle*(nu_con_mle*((y(jj-1,1)-mu_mle(2)).^2)/C - f_mle) ...
                        + mu_mle(5)*f_mle;
        end
        THR(jj) = (y(jj) - mu_mle(2))/sqrt(rho_mle*f_mle);
        if (THR(jj) < threshold) % observations of interest
            cond(jj) = 1;
        else
            sum_C = sum_C + 1;
        end
        THR(jj) = mu_mle(2) + sqrt(rho_mle*f_mle)*threshold;
    end
    
    z1 = zeros(sum_C,1);
    z2 = zeros(T-sum_C,1);
    
    ind = (y<threshold); % observations of interest

    for ii = 1:M        
        if prior(ii,1) % compute only for valid draws
            i2 = 0;
            scale = sqrt(rho(ii,1)*f(ii,1));
            d(ii,1) = 0;
            if ind(1,1) 
                z1  = (y(1,1) - mu(ii,1))/scale;
                d(ii,1) = d(ii,1) + log(tpdf(z1,nu(ii))/scale);
             else
                i2 = i2+1;
                z2(i2,1) = (threshold - mu(ii,1))/scale;
             end    
        
            for jj = 2:T
                C = 1 + ((y(jj-1,1)-mu(ii,1)).^2)/((nu(ii,1)-2)*f(ii,1));
                
                f(ii,1) = omega(ii,1) + A(ii,1)*(nu_con(ii,1)*((y(jj-1,1)-mu(ii,1)).^2)/C - f(ii,1)) ...
                            + B(ii,1)*f(ii,1);                   
                scale = sqrt(rho(ii,1)*f(ii,1));
                if ind(jj,1)
                    z1 = (y(jj,1)-mu(ii,1))/scale;
                    d(ii,1) = d(ii,1)  + log(tpdf(z1,nu(ii))/scale);
                 else  % P(y_t>=threshold|y_1,...,y_(t-1)) = 1 - P(y_t<threshold|y_1,...,y_(t-1))
                    i2 = i2+1;
                    z2(i2,1) = (threshold-mu(ii,1))/scale;
                 end
            end 
            d(ii,1) = d(ii,1) + sum(log(1-tcdf(z2,nu(ii,1))));
        end
    end
    d = d + prior(:,2);
end

function R = prior_t_gas(M, omega, B, nu, hyper)
    % uniform priors 
    
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val an the corresponding point
    
    c1 = (omega > 0);
    c3 = ((B >= 0) & (B < 1));
    c4 = (nu > 2);
    
    r1 = (c1 & c3 & c4);
    r2 = -Inf*ones(M,1);
    r2(r1==true) = log(hyper) - hyper*( nu(r1==true)  - 2); % exponential prior: nu~exp(1) --> p(nu)=exp(-nu) from 2 to inf  

    R = [r1, r2];
end