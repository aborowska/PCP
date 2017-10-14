function d = posterior_t_gas(theta, y, hyper)
    % theta is Nx5, matrix of draws
    [N ,~] = size(theta);
    nu = theta(:,1);
    mu = theta(:,2);
    omega = theta(:,3);
    A = theta(:,4);
    B = theta(:,5); 
    
    rho = (nu-2)./nu;
    
    nu_con = (nu+1)./(nu-2);
    A = A.*((nu+3)./nu);

    prior = prior_t_gas(N, omega, B, nu, hyper);
    
    T = size(y,1);
                 
    f = omega./(1-B); % unconditional variance to initialize h_1

    d = -Inf*ones(N,1);
       
    for ii = 1:N               
        if (prior(ii,1)) % when all the parameter constraints are satisfied
            
%             pdf(1,1) = dmvt(y(1,1), mu(ii,1), rho(ii,1)*f(ii,1), nu(ii,1), GamMat);
%             pdf(1,1) = log(pdf(1,1));

            scale = sqrt(rho(ii,1)*f(ii,1));
            z = (y(1,1) - mu(ii,1))/scale;
            d(ii,1) = log(tpdf(z,nu(ii))/scale);
% %             d(ii,1) = log(duvt_garch(y(1,1),mu(ii,1),rho(ii,1)*f(ii,1),nu(ii,1)));             
            
            for jj = 2:T
                C = 1 + ((y(jj-1,1)-mu(ii,1)).^2)/((nu(ii,1)-2)*f(ii,1));
                
                f(ii,1) = omega(ii,1) + A(ii,1)*(nu_con(ii,1)*((y(jj-1,1)-mu(ii,1)).^2)/C - f(ii,1)) ...
                            + B(ii,1)*f(ii,1);
                        
%                 pdf(jj,1) = dmvt(y(jj,1), mu(ii,1), rho(ii,1)*f(ii,1), nu(ii,1), GamMat);
%                 pdf(jj,1) = log(pdf(jj,1));
                   
                scale = sqrt(rho(ii,1)*f(ii,1));
                z = (y(jj,1)-mu(ii,1))/scale;
                d(ii,1) = d(ii,1) + log(tpdf(z,nu(ii,1))/scale);             
            end
        end
    end
    
    d = d + prior(:,2);
end


function R = prior_t_gas(N, omega, B, nu, hyper)
    % uniform priors 
    
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val an the corresponding point
    
    c1 = (omega > 0);
    c3 = ((B >= 0) & (B < 1));
    c4 = (nu > 2);
    
    r1 = (c1 & c3 & c4);
    r2 = -Inf*ones(N,1);
    r2(r1==true) = log(hyper) - hyper*( nu(r1==true)  - 2); % exponential prior: nu~exp(1) --> p(nu)=exp(-nu) from 2 to inf  

    R = [r1, r2];
end
