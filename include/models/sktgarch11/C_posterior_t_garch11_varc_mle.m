function d = C_posterior_t_garch11_varc_mle(theta, y, mu_mle, threshold, ...
    y_S, hyper)
% function [d,THR,SCALE,z2] = C_posterior_t_garch11_varc_mle(theta, y, mu_mle, threshold, ...
%     y_S, hyper)
    
    T = length(y);
    M = size(theta,1);
    
    nu = theta(:,1);
    mu = theta(:,2);
    omega = theta(:,3);
    alpha = theta(:,4);
    beta = theta(:,5);

    rho = (nu-2)./nu;
    rho_mle = (mu_mle(1)-2)/mu_mle(1);
    
    if (y_S == 0) 
        h = omega;
    else        
        h = y_S*ones(M,1);
    end      
    
    d = -Inf*ones(M,1);
    prior = prior_t_garch(M, nu, omega, alpha, beta, hyper); 
    
    sum_C = 0;
    cond = zeros(T,1);
    THR = zeros(T,1);
    for jj = 1:T
        if (jj == 1)
            h_mle = y_S;
        else
            temp = y(jj-1) - mu_mle(2);
            h_mle = mu_mle(3)*(1-mu_mle(4)-mu_mle(5)) ...
                + mu_mle(4)*temp^2 + mu_mle(5)*h_mle;
        end
        THR(jj) = (y(jj) - mu_mle(2))/sqrt(rho_mle*h_mle);
        if (THR(jj) < threshold) % observations of interest
            cond(jj) = 1;
        else
            sum_C = sum_C + 1;
        end
        THR(jj) = mu_mle(2) + sqrt(rho_mle*h_mle)*threshold;
    end
    
%     z1 = zeros(sum_C,1);
    z2 = zeros(T-sum_C,1);
    

    for ii = 1:M        
        if prior(ii,1) % compute only for valid draws
            i2 = 0;
%             SCALE = zeros(T,1);
            scale = sqrt(rho(ii,1)*h(ii,1));
%             SCALE(1,1) = scale;
            d(ii,1) = 0;
%             d(ii,1) = prior(ii,2);
            if cond(1,1) 
                z1  = (y(1,1) - mu(ii,1))/scale;
                d(ii,1) = d(ii,1) + log(tpdf(z1,nu(ii))/scale);
             else
                i2 = i2+1;
                z2(i2,1) = (THR(jj) - mu(ii,1))/scale;
             end    
        
            for jj = 2:T
                h(ii,1) = omega(ii,1)*(1-alpha(ii,1)-beta(ii,1)) + alpha(ii,1)*(y(jj-1,1)-mu(ii,1)).^2  + beta(ii,1)*h(ii,1);
                scale = sqrt(rho(ii,1)*h(ii,1));
%                 SCALE(jj,1) = scale;

                if cond(jj,1)
                    z1 = (y(jj,1)-mu(ii,1))/scale;
%                     fprintf('d(%i,1) = %16.14f\n', jj-1,d(ii,1));

                    d(ii,1) = d(ii,1)  + log(tpdf(z1,nu(ii))/scale);
                 else  % P(y_t>=threshold|y_1,...,y_(t-1)) = 1 - P(y_t<threshold|y_1,...,y_(t-1))
                    i2 = i2+1;
                    z2(i2,1) = (THR(jj)-mu(ii,1))/scale;
                 end
            end
%             fprintf('d(ii,1) = %6.4f\n', d(ii,1));
%             for oo = 1:10:sum_C
%                 fprintf('z2(%i,1) = %16.14f\n',oo,z2(oo));
%                 fprintf('cdf(%i,1) = %16.14f\n',oo,tcdf(z2(oo),nu(ii,1)))
%             end
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
