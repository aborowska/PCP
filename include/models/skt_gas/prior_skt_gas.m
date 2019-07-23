function r2 = prior_skt_gas(theta, y, hyper, y_S)
    % theta is Nx5, matrix of draws
    T = length(y);
    [N ,~] = size(theta);
    
    lambda = theta(:,1);    
    nu = theta(:,2);
    mu = theta(:,3);
    omega = theta(:,4);
    A = theta(:,5);
    B = theta(:,6); 
    
    
    c1 = (omega > 0);
    c3 = ((B >= 0) & (B < 1));
    c4 = (nu > 2);
    c5 = ((lambda > -1) & lambda < 1);
    
    r1 = (c1 & c3 & c4 & c5);
    r2 = -Inf*ones(N,1);
    r2(r1==true) = log(0.5) + log(hyper) - hyper*( nu(r1==true)  - 2); % exponential prior: nu~exp(1) --> p(nu)=exp(-nu) from 2 to inf  

end
